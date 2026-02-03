# evb

This directory contains the standalone implementation of the empirical valence bond (EVB) model used in the paper.

## Files

- `evb.py`  
  Core EVB Hamiltonian and energy/force evaluation routines.

- `parser.py`  
  Utilities for reading parameter files (e.g. `para.json`, force-field XML/STR) and constructing EVB systems.

- `check_forces.py`  
  Consistency checks for energies and forces (e.g. finite-difference tests).

- `utils.py`  
  Helper functions used across the EVB code (I/O, unit conversions, etc.).

- `__init__.py`  
  Allows importing the module via `import evb`.

Compiled `.pyc` files in `__pycache__/` are automatically generated and can be ignored.

## Usage

The same EVB code is copied into several `md_simulations/*/eqm/evb/` and `md_simulations/*/pull-1/evb/` subdirectories 
to ensure self-contained examples. In most cases, you do not need to edit this code; you only call it via the provided MD scripts 
(e.g., `runSMD.py`, `eqm.py`).

To use the EVB module directly in Python:

```python
import evb
# See evb.evb and evb.parser for available functions and classes.
```

For details, see comments inside `evb.py` and the simulation-specific READMEs in `md_simulations/`.
