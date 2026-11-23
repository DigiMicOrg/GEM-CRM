import math

SUBSTRATES = {'D-glucose': 'cpd00027', 'L-glutamine': 'cpd00053', 'L-glutamate': 'cpd00023', 'L-asparagine': 'cpd00132', 'L-leucine': 'cpd00107', 'L-phenylalanine': 'cpd00066'}
REF_BIOMASS_C = 41.5
REF_PROTEIN_FRACTION = 0.55

def number(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return float('nan')

def gmean(xs):
    xs = [x for x in xs if math.isfinite(x) and x > 0]
    return math.exp(sum((math.log(x) for x in xs)) / len(xs)) if xs else None

def record_values(e, biomass_C, protein_fraction):
    """Convert one audited whole-cell uptake record to carbon-normalised units."""
    carbon = number(e['carbon_atoms'])
    K = number(e['K_value']) * {'uM': 0.001, 'mM': 1}.get(e['K_unit'], float('nan')) * carbon
    rate = number(e['Vmax_value']) * 0.06
    conversion = 'dry_mass'
    if e['Vmax_unit'] == 'nmol/min/mgProtein':
        if protein_fraction is None:
            return None
        rate *= protein_fraction
        conversion = 'model_biomass_protein_fraction'
    elif e['Vmax_unit'] != 'nmol/min/mgDW':
        return None
    vmax = rate * carbon / biomass_C
    if not (math.isfinite(K) and math.isfinite(vmax) and (K > 0) and (vmax > 0)):
        return None
    return (K, vmax, vmax / K, conversion)

def usable_record(e):
    return e['measurement_type'] == 'whole_cell_uptake' and e['decision'].startswith('candidate')

def resource_prior(resource, evidence):
    """Substrate-level kinetic scale from all audited records for this substrate."""
    pool = []
    for e in evidence:
        if SUBSTRATES.get(e['substrate']) != resource or not usable_record(e):
            continue
        x = record_values(e, REF_BIOMASS_C, REF_PROTEIN_FRACTION)
        if x:
            pool.append(x)
    if not pool:
        return None
    ks = [x[0] for x in pool if math.isfinite(x[0]) and x[0] > 0]
    return {'K': gmean([x[0] for x in pool]), 'K_min': min(ks) if ks else None, 'K_max': max(ks) if ks else None, 'u': gmean([x[2] for x in pool]), 'n': len(pool), 'evidence_ids': ';'.join((e['record_id'] for e in evidence if SUBSTRATES.get(e['substrate']) == resource and usable_record(e))), 'source': ';'.join((e['source_url'] for e in evidence if SUBSTRATES.get(e['substrate']) == resource and usable_record(e)))}