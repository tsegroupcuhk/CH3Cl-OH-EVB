# md_simulations/reactive_slab

This directory contains EVB-based simulations of the reaction at the air–water interface and at selected depths within a water slab.

There are two types of folders:

```text
reactive_slab/
├── _z_*          # Aggregated folders for each z-position (with README only)
├── _ecn_z_*      # Aggregated folders for ECN-resolved data at selected z (with README only)
├── ecn_z_10/     # Full simulation setup with ECN-resolved windows (example, subfolder ecn_2/)
└── z_-6/         # Full simulation setup at a given z (example, similar to bulk layout)
```

The folders with a leading underscore (`_z_*`, `_ecn_z_*`) primarily contain README files that describe how the data in `data/reactive_pmfs/` map to the various interfacial positions and effective coordination numbers. 
The **full simulation setups** are contained in:

- `z_-6/`
- `ecn_z_10/ecn_2` 

## Typical directory structure (e.g. `z_-6/`)
```text
z_-6/
├── config_bank/        # Restart configurations for multiple SMD trajectories
├── CV/                 # Resulting CV time series files
├── eqm/                # Equilibration scripts and EVB code
├── pull-1/             # Example SMD at this z
├── template/           # Template for generating more SMD folders
└── generate_running_folder_1-150.sh
```
The internal layout and workflow are analogous to the bulk case (`md_simulations/reactive_bulk/`):

1. Equilibrate the system at the desired z-constraint (`eqm/`).
2. Use `generate_running_folder_1-150.sh` to create multiple SMD replicas.
3. Run SMD (`pull-*` folders) and collect CV and work data.
4. Analyze work distributions and reconstruct PMFs.
5. Compare to the processed PMFs in `data/reactive_pmfs/pmf_slab_z_*.dat`.

## ECN-resolved simulations (e.g. `ecn_z_10/`)
Within `ecn_z_10/`, subfolders such as `_ecn_0`, `_ecn_1`, `_ecn_3`, `_ecn_4`, and `ecn_2/` correspond to windows at different effective coordination numbers (ECN):
- Subdirectories starting with `_ecn_*` contain the documentation about the ECN selection and mapping to processed PMFs.
- The `ecn_2/` folder contains a full simulation setup analogous to `z_-6/`, but with additional ECN constraints applied.

The resulting ECN-resolved PMFs are stored in:
- `data/reactive_pmfs/ecn_z_10/pmf_slab_ecn_*.dat`, etc.

Please refer to the main text for the definition of ECN and its physical interpretation.

