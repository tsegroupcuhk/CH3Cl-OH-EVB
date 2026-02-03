# qm

This directory contains the quantum chemistry calculations used to:

1. Determine stationary points (reactant, product, transition state) for the CH<sub>3</sub>Cl + OH reaction.
2. Generate a reference potential energy surface (PES) along the intrinsic reaction coordinate (IRC).
3. Provide input for parameterizing and validating the EVB model.

```text
qm/
├── cluster_irc/             # IRC and single-point energies along the reaction path
└── stationary_points_coord/ # Optimized geometries and logs for stationary points
```

### `cluster_irc/`
- `irc.gjf`, irc.log`\
Gaussian input and output for the IRC calculation.

- `irc_points_coord/`\
Cartesian coordinates (`*.xyz`) of selected points along the IRC.

- `single_point_along_irc_dlpno/`\
  - `single_point_orca.inp`: ORCA input file used to compute TightPNO DLPNO-CCSD(T)/aug-cc-pVQZ single points on top of the IRC geometries.
  - `pes_irc_dlpno.dat`: Resulting PES (energies vs IRC point index), also copied to data/qm_pes/.

### `stationary_points_coord/`
See the local `README` for details. In brief, this folder provides:

- Gaussian inputs and logs for:
  - Reactant optimization (`gaussian_reactant_opt/`)
  - Product optimization (`gaussian_product_opt/`)
  - Transition state optimization (`gaussian_TS_opt/`)
- Convenience `.xyz` files for reactant, product, and TS geometries.

These data were used to define the EVB states, validate activation barriers, and compare solution-phase energetics with the EVB model.

