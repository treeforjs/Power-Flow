# Cross-Section Assets

This directory stores reusable collision cross-section provenance and tables.

- `prab_26_040401_reactions.yaml` is extracted from `PhysRevAccelBeams.26.040401.pdf` and records the refs. 75-80 reaction mapping.
- `cross_section_manifest.yaml` is the loader manifest used by `mhdlab.cross_sections.CrossSectionLibrary`.
- `nist/*.csv` contains NIST BEB/BEQ starter tables converted from Angstrom^2 to m^2.

Older measured ion-impact and charge-exchange tables from PRAB refs. 77-80 are listed in the manifest with `file: null` until digitized/tabulated data are added.
