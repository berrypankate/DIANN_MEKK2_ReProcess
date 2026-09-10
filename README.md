# DIA-NN MAP3K2 / MAP3K3 Analysis

This project tests how different DIA-NN protein-inference settings affect peptide mapping and protein quantification for the closely related proteins MAP3K2 (MEKK2) and MAP3K3 (MEKK3).

## Project Goal

MAP3K2 and MAP3K3 have similar protein sequences, so some peptides can map to both proteins.

The goal is to determine which DIA-NN settings best distinguish MAP3K2 from MAP3K3 while keeping useful high-confidence peptides for protein quantification.

MAP3K2 depletion provides a biological check: MAP3K2-specific peptide evidence should disappear after depletion, while MAP3K3 evidence can remain.

## DIA-NN Conditions

Three conditions have been analyzed:

- **Baseline** — standard DIA-NN protein inference
- **PF (`--proteoforms`)** — enables proteoform confidence scoring
- **NPI (`--no-prot-inf`)** — disables DIA-NN protein inference and uses protein groups from the spectral library

A fourth condition combines both parameters:

- **PF + NPI** — `--proteoforms --no-prot-inf`

## Current Findings

The main biological result is consistent across Baseline, PF, and NPI:

- MAP3K2-specific peptides are present in control samples.
- MAP3K2-specific peptides disappear after MAP3K2 depletion.
- MAP3K3-specific peptides remain.
- One peptide, `QVQFDPDSPETSK`, maps to both MAP3K2 and MAP3K3.

This shared peptide has `Proteotypic = 0`, meaning it cannot uniquely identify either protein.

Baseline and PF assign the shared peptide to the MAP3K3 protein group. NPI keeps it in a combined MAP3K2/MAP3K3 protein group.

This shows that changing protein-inference settings can change protein grouping and protein-level quantification even when the detected peptide sequences remain the same.

## Repository Structure

- `analysis/` — Python scripts for peptide and protein-level analysis
- `bin/` — DIA-NN batch processing and finalization scripts
- `config/` — experiment configuration files
- `config/grid/` — configs for PF, NPI, and combined parameter testing

Large RAW files, DIA-NN quant files, spectral libraries, and parquet outputs are not stored in GitHub.

## Main Analysis Scripts

`analysis/compare_inference_modes.py` compares Baseline, PF, and NPI results.

`analysis/inspect_inference_rows.py` examines individual peptide rows, protein-group assignments, proteotypic status, and protein-level abundance differences.

`analysis/run_diann_inference_grid.sh` runs the DIA-NN inference conditions.

## Next Steps

1. Complete the PF + NPI condition.
2. Compare all four DIA-NN modes.
3. Prioritize `Proteotypic = 1` peptides for confident paralog-specific evidence.
4. Compare protein grouping and abundance across conditions.
5. Identify the DIA-NN setting that gives reliable MAP3K2/MAP3K3 mapping while retaining useful quantitative peptide evidence.