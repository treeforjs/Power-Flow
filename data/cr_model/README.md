# Collisional-Radiative Model Data

This directory keeps CR-model provenance separate from the run configuration.

- `references.yaml` stores the hydrogen/H2 CR references supplied from the Yacora help-page bibliography.
- `yacora_hydrogen_reactions.yaml` captures the Yacora hydrogen reaction overview in a code-readable form.
- `rovibrational_manifest.yaml` defines where state-resolved H2/H2+ rovibrational cross sections, rates, and transition probabilities should be placed.

Rows with `file: null` are intentional placeholders for curated or digitized tables. Tables marked `can_feed_bolsig: true` are candidates for BOLSIG/LXCat-style electron collision decks once state-resolved data are added.
