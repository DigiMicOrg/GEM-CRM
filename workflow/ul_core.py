"""Binary structural U/L evidence. No flux estimation and no genus union."""
import math
import gzip
import re
import xml.etree.ElementTree as ET

# Only suppress these endpoints in L; retain them in reactions and U.
DEFAULT_EXCLUDED = {'cpd00001', 'cpd00067', 'cpd00002', 'cpd00008', 'cpd00009',
                    'cpd00012', 'cpd00003', 'cpd00004', 'cpd00005', 'cpd00006'}
EXCLUDED_SBO = {'SBO:0000627', 'SBO:0000628', 'SBO:0000629', 'SBO:0000630'}


def attrs(element, name, default=''):
    return next((v for k, v in element.attrib.items() if k.split('}')[-1] == name), default)


def parse_sbml(path):
    if str(path).endswith('.gz'):
        with gzip.open(path, 'rb') as handle: root = ET.parse(handle).getroot()
    else:
        root = ET.parse(path).getroot()
    model = root.find('{*}model')
    if model is None:
        raise ValueError('SBML has no model')
    params = {p.attrib['id']: float(p.attrib['value'])
              for p in model.findall('{*}listOfParameters/{*}parameter')}
    species = {}
    for s in model.findall('{*}listOfSpecies/{*}species'):
        sid = s.attrib['id']
        # ModelSEED base IDs are shared across c0/e0; do not collapse unknown IDs.
        match = re.fullmatch(r'M_(cpd\d+)_([^_]+)', sid)
        species[sid] = {'base': match[1] if match else sid,
                        'compartment': s.attrib.get('compartment', ''),
                        'name': s.attrib.get('name', sid)}
    reactions = []
    for r in model.findall('{*}listOfReactions/{*}reaction'):
        bounds = []
        for key in ('lowerFluxBound', 'upperFluxBound'):
            ref = attrs(r, key)
            if ref not in params or not math.isfinite(params[ref]):
                raise ValueError(f'{r.attrib["id"]}: missing/invalid {key}')
            bounds.append(params[ref])
        if bounds[0] > bounds[1]:
            raise ValueError('Lower bound exceeds upper bound')
        sides = []
        for side in ('Reactants', 'Products'):
            values = {}
            for s in r.findall('{*}listOf' + side + '/{*}speciesReference'):
                sid = s.attrib['species']
                if sid not in species:
                    raise ValueError('Unknown species ' + sid)
                coefficient = float(s.attrib.get('stoichiometry', 1))
                if not math.isfinite(coefficient) or coefficient <= 0:
                    raise ValueError('Invalid stoichiometry')
                values[sid] = values.get(sid, 0) + coefficient
            sides.append(values)
        seed_ids = set()
        for li in r.findall('.//{*}li'):
            m = re.search(r'seed\.reaction[/:](rxn\d+)', attrs(li, 'resource'))
            if m:
                seed_ids.add(m[1])
        gpr = r.find('{*}geneProductAssociation')
        genes = sorted({attrs(g, 'geneProduct') for g in r.findall('.//{*}geneProductRef')})
        # An association is model evidence, not proof of direct homology.
        origin = ''
        for element in r.iter():
            if element.tag.split('}')[-1] in ('origin', 'gs.origin'):
                origin = ''.join(element.itertext()).strip()
            for key, value in element.attrib.items():
                if key.split('}')[-1] in ('origin', 'gs.origin'):
                    origin = value
        reactions.append(dict(id=r.attrib['id'], lb=bounds[0], ub=bounds[1],
                              reactants=sides[0], products=sides[1],
                              sbo=r.attrib.get('sboTerm', ''), seed=seed_ids,
                              genes=';'.join(genes), gpr_xml=ET.tostring(gpr, encoding='unicode') if gpr is not None else '',
                              origin=origin or 'unknown'))
    return species, reactions


def directed(reactions):
    for reaction in reactions:
        if reaction['ub'] > 0:
            yield reaction, 'forward', reaction['reactants'], reaction['products']
        if reaction['lb'] < 0:
            yield reaction, 'reverse', reaction['products'], reaction['reactants']


def extract(species, reactions, transporter_rows, excluded=None):
    excluded = DEFAULT_EXCLUDED if excluded is None else set(excluded)
    links = {}
    for t in transporter_rows:
        base = re.sub(r'_e0$', '', re.sub(r'^EX_', '', t.get('exid', '').strip()))
        for rid in re.split(r'\s*,\s*', t.get('rea', '').strip()):
            if base and rid:
                links.setdefault(base, set()).add(rid)
    imports, exports, boundary = {}, {}, {}
    consumption, intracellular = {}, []
    trace = []
    def record(reaction, direction):
        return dict(reaction=reaction['id'], direction=direction,
                    genes=reaction.get('genes', ''), origin=reaction.get('origin', 'unknown'))
    for r, direction, left, right in directed(reactions):
        # Boundary output is detected from oriented stoichiometry, not a
        # presumed sign convention for EX_ reaction IDs.
        if (r['sbo'] == 'SBO:0000627' or r['id'].startswith('R_EX_')) and not right and len(left) == 1:
            sid = next(iter(left))
            if species[sid]['compartment'] == 'e0':
                boundary.setdefault(species[sid]['base'], []).append(record(r, direction))
        if r['sbo'] in EXCLUDED_SBO or not left or not right:
            continue
        lc = {(species[s]['base'], species[s]['compartment']): v for s, v in left.items()}
        rc = {(species[s]['base'], species[s]['compartment']): v for s, v in right.items()}
        for base, seeds in links.items():
            if not (r['seed'] & seeds):
                continue
            if lc.get((base, 'e0'), 0) > rc.get((base, 'e0'), 0) and rc.get((base, 'c0'), 0) > lc.get((base, 'c0'), 0):
                imports.setdefault(base, []).append(record(r, direction))
            if lc.get((base, 'c0'), 0) > rc.get((base, 'c0'), 0) and rc.get((base, 'e0'), 0) > lc.get((base, 'e0'), 0):
                exports.setdefault(base, []).append(record(r, direction))
        if not all(species[s]['compartment'] == 'c0' for s in set(left) | set(right)):
            continue
        # Only net-consumed/net-produced metabolites count; catalysts do not.
        consumed = {species[s]['base'] for s in left if left[s] > right.get(s, 0)}
        produced = {species[s]['base'] for s in right if right[s] > left.get(s, 0)}
        for base in consumed:
            consumption.setdefault(base, []).append(record(r, direction))
        intracellular.append((r, direction, consumed, produced))
    uptake = set(imports) & set(consumption)
    exportable = set(exports) & set(boundary)
    pairs = set()
    for base in sorted(uptake):
        trace.append(dict(kind='U', substrate=base, product='',
                          import_evidence=imports[base], reaction_evidence=consumption[base],
                          export_evidence=[], boundary_evidence=[]))
    for r, direction, consumed, produced in intracellular:
        for a in sorted((consumed & uptake) - excluded):
            for b in sorted((produced & exportable) - excluded - {a}):
                pairs.add((a, b))
                trace.append(dict(kind='L', substrate=a, product=b,
                                  import_evidence=imports[a], reaction_evidence=[record(r, direction)],
                                  export_evidence=exports[b], boundary_evidence=boundary[b]))
    return dict(uptake=sorted(uptake), pairs=sorted(pairs), evidence=trace,
                metabolites=sorted({s['base'] for s in species.values() if s['compartment'] in ('c0', 'e0')}))


def select_hits(rows, identity=98.0, coverage=95.0, score_delta=0.0):
    """Best eligible HSP per genome, then retain co-best genomes per ASV.

    qcovhsp is coverage of this alignment; fragmented HSPs are not summed.
    These are configurable screening thresholds, not species boundaries.
    """
    best = {}
    for row in rows:
        if float(row['pident']) < identity or float(row['qcovhsp']) < coverage:
            continue
        key = (row['qseqid'], row['accession'])
        rank = tuple(float(row[x]) for x in ('bitscore', 'qcovhsp', 'pident'))
        if key not in best or rank > best[key][0]:
            best[key] = (rank, row)
    top = {}
    for (asv, _), (rank, _) in best.items():
        top[asv] = max(top.get(asv, -math.inf), rank[0])
    return [row for (asv, _), (rank, row) in sorted(best.items()) if rank[0] >= top[asv] - score_delta]
