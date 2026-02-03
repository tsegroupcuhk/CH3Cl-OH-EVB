# forcefield

This directory contains the EVB and classical force-field parameters used in the simulations.

## Files

- `para.json`  
  JSON file with the EVB parameters (state energies, couplings, switching functions, etc.) for the CH<sub>3</sub>Cl + OH system.

- `off.xml`  
  Force-field definition in CHARMM/OpenMM-compatible XML format.

- `toppar.xml`, `toppar.str`  
  CHARMM-style topology and parameter files (XML and STR formats).

- `computed_energy_with_this_evb_model.txt`  
  Sanity check/validation file listing reference configurations and the corresponding energies computed with this EVB model.

## Notes

- Copies of these files also appear in various `inits/` directories under `md_simulations/` so that each example can be run independently.
- If you modify any parameters, ensure consistency between `para.json` and the associated topology/parameter files.
