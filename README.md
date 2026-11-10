# MiCRM 16S workflow

Go from raw 16S FASTQ files to a "MiCRM simulated abundance vs 16S observed abundance" comparison bar chart in one command.
No data, caches or pre-built GEMs are included. All required reference databases, genomes and GEMs are downloaded or built automatically on the first run and cached, so they can be reused directly when you switch to new data.

## Directory layout

```
micrm_workflow/
├── main.py              # One-command entry point; runs the 6 steps in order
├── config.toml          # All tunable parameters (paths, threads, thresholds, MiCRM settings)
├── environment.yml      # conda environment (R/DADA2, gapseq, barrnap, mash, NCBI datasets)
├── workflow/            # Function files for each step
│   ├── denoise.py / denoise.R    # Step 1: DADA2 denoising + SILVA classification to genus
│   ├── genomes.py                # Step 3: NCBI reference genomes + 16S Mash selection + similar-genus/family fallback
│   ├── gems.py                   # Step 4: build GEMs with gapseq (resumable)
│   ├── parameterize.py           # Step 5: binary uptake capability, resource panel, K
│   ├── maintenance.py            # Step 5: cell size -> m
│   ├── kinetic_helpers.py + resources/kinetics_evidence.csv   # Literature evidence for K
│   ├── micrm.py                  # Step 6: U/L estimation, abundance constraints, BDF simulation, plotting
│   ├── estimate.py, network.py, entropy_solver.py, bdf_solver.py   # Original micrm_shared core (unchanged)
│   ├── ul_core.py                # SBML parsing (originally vendored, unchanged)
│   ├── databases.py, common.py   # Tool checks, reference database download, common functions
├── 16data/              # <- Put the raw FASTQ files here
├── databases/           # Auto-downloaded: SILVA training set, Madin trait table
├── genomes/             # Auto-cached: NCBI taxonomy queries, genomes, 16S, selection records
├── gems/                # Auto-built: gapseq models (one folder per accession) + U/L computation cache
└── results/             # Per-community results, figures/, summary.csv
```

## Installation

```bash
conda env create -f environment.yml
conda activate micrm16s
python main.py --check      # Check tools and pre-download SILVA (about 140 MB) and the Madin table
```

An internet connection is required: NCBI (taxonomy + datasets), Zenodo (SILVA), GitHub (Madin table), and the sequence libraries that gapseq downloads by itself on first run. If NCBI requests are rate-limited, you can request an API key and put it in `ncbi_api_key` in `config.toml`, or set the environment variable `NCBI_API_KEY`.

## Running

```bash
# Put the FASTQ files into 16data/
python main.py
```

- **One FASTQ = one community**, and the file name is the community name. `.fastq(.gz)` and `.fq(.gz)` are supported, and files may also be placed in subfolders.
- Paired-end data (`xxx_R1/_R2`, `xxx_1/_2`, `xxx_R1_001`) are detected automatically and only R1 is used, consistent with the original workflow.
- After an interruption, just run `python main.py` again; completed parts are skipped: already denoised communities, already selected genomes, already built GEMs, and simulations whose inputs have not changed.
- Common options:
  - `--steps params,simulate`: rerun only parameterization and simulation (e.g. after changing `[micrm]` in `config.toml`).
  - `--rerun`: force the simulation to rerun.
  - `--communities A,B`: run only the specified communities.
- gapseq modelling of several dozen genera can take from a few hours to a day, so it is recommended to run on a server with `nohup python main.py > run.log 2>&1 &`. The number of threads is set in `[compute]`.

## What each step does

1. **Denoising (DADA2, forward reads only)**: truncLen = min(modal read length, 230). If the quality scores are nearly constant, `noqualErrfun` is used automatically. After a first failure, it retries once with relaxed parameters (no truncation, maxEE 5, noqual). Finally, ASVs are classified to genus with SILVA 138.2. Only ASVs with Kingdom = Bacteria and a genus assignment go into modelling.
2. **Genus summary**: collect the ASV sequences of each genus across communities, together with the family assigned by SILVA.
3. **Reference genomes**: tried in the following order, stopping at the first hit:
   1. A genus of the same name in NCBI (restricted to Bacteria)
   2. Aliases for combined labels (Escherichia-Shigella -> Escherichia + Shigella)
   3. Similar genus: extract the genus name from the label, e.g. `[Eubacterium] xx group` -> Eubacterium, `Clostridium sensu stricto 1` -> Clostridium, and hyphenated combinations are split into individual genera
   4. **Family fallback**: reference genomes of other genera in the same family, one per genus
   5. If none is found, skip

   Within each tier, the search order is: NCBI reference -> type material -> any annotated assembly, taking up to `max_candidates` candidates. The 16S sequences of the candidate genomes are extracted with barrnap and compared with the ASVs of that genus using Mash, and the closest candidate is selected. The origin is recorded in the `source` column of `results/genus_references.csv`.
4. **GEM**: gapseq find -> find-transport -> draft -> medium -> fill. Each stage has a checkpoint, stored in `gems/<accession>/`, and GEMs that have already been built are reused directly.
5. **Parameter preparation** (one set per community, stored in `results/<community>/02_inputs/`):
   - **Binary uptake matrix** `binary_uptake.csv`: a substance is marked 1 if the GEM has an extracellular exchange and a transport reaction for it, it contains carbon, and it is not CO2/HCO3-.
   - **Resource panel**: the union of the resources that can be taken up by all modelled genera in the community; all of them are supplied.
   - **K**: from literature records of whole-cell uptake, converted to mM-C, taking the geometric mean for the same substrate. Substrates without evidence fall back to the global geometric mean (0.01262 mM-C). The same resource uses the same K across all genera.
   - **m**: cell volume is calculated from the cell diameter and length in Madin et al. 2020, and m = 0.01449·V^-0.09. Genera without trait data are assigned a random value within the known volume range (seed 42). This is a per-genus table lookup, **not inferred from genes**.
6. **Simulation** (`results/<community>/03_micrm/`):
   - **U**: for each genus × resource pair, supply 1 unit of carbon, find the first feasible retention fraction in a prior order, then use LP to find the minimum carbon flux burden B, with u = 10/(1+B).
   - **L**: take the maximum output entropy over the same feasible flux polytope as the leakage allocation. `binary_leakage.csv` lists every substrate -> product combination.
   - **Abundance constraints**: each row of U is multiplied by abundance / mean abundance; L is reallocated as 90% abundance-weighted demand + 10% uniform, keeping the total carbon unchanged. In some samples, missing information during genome reconstruction may lead to an incomplete representation of the metabolic potential of dominant consumers, causing them to become disadvantaged in simulations, in direct contrast to observations. The underlying assumption is that more abundant taxa are likely to consume more resources.
   - **MiCRM**: C0 = 1, total R0 = 1e-4, total ρ = 5 (both distributed uniformly over the resource panel), turnover 0.1, integrated with BDF to steady state.
   - **Plotting**: draw a bar chart comparing the relative abundances of the modelled genera in the community, and collect them into `results/figures/<community>.png`.

## Main outputs

| File | Content |
|---|---|
| `results/summary.csv` | Per-community status, number of genera, number of resources, coverage, JSD |
| `results/figures/<community>.png` | MiCRM vs 16S comparison bar chart |
| `results/genus_references.csv` | Genome selected for each genus, source (genus / alias / related_genus / family_surrogate), Mash distance |
| `results/gem_catalog.csv` | genus -> GEM path and status |
| `results/maintenance.csv` | Volume, m and source for each genus |
| `results/<community>/01_dada2/` | asvs.csv, counts.csv, qc.csv, quality.png |
| `results/<community>/02_inputs/` | genera / resources / observed / kinetics / maintenance / binary_uptake |
| `results/<community>/03_micrm/` | parameters.npz (U, L, K, m, R0, ρ), parameter_audit.json, simulation.json, abundance_comparison.csv, comparison.png, metadata.json |

## Notes

- `kinetics_evidence.csv` is the literature evidence for K (part of the method), so it is shipped together with the code.
- A family-fallback result is a "representative close-relative genome" and is less reliable than a same-genus match. You can filter out those rows in `genus_references.csv`, or set `family_fallback = false` in the config.
- Communities with fewer than 2 usable GEMs are marked as skipped, with the reason written to summary.csv.
