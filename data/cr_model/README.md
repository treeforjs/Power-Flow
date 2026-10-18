# Collisional-Radiative Model Data

This directory keeps CR-model provenance separate from the run configuration.

- `references.yaml` stores the hydrogen/H2 CR references supplied from the Yacora help-page bibliography.
- `yacora_hydrogen_reactions.yaml` captures the Yacora hydrogen reaction overview in a code-readable form.
- `molecular_state_requirements.yaml` records the publication requirement for explicit rovibrational populations in every retained molecular neutral, cation and anion, including polyatomic state labels and missing collision/formation data. It is a requirements inventory, not an executable or complete population model.
- `rovibrational_manifest.yaml` contains the legacy H2/H2+ data-slot subset and links to the network-wide requirements. Those legacy rows do not establish molecular coverage for the publication model.
- `molat_h2/` supplies the retained H2 B/C-to-X bound-bound radiative transitions; it contains no solved populations or state-to-state collision rates.

Rows with `file: null` are intentional placeholders for curated or digitized tables. Tables marked `can_feed_bolsig: true` are candidates for BOLSIG/LXCat-style electron collision decks once state-resolved data are added.
