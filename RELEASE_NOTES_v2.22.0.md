# Reticula 2.22.0

Feature release: TEM particle sizing (edge particles, rods, agglomerate tagger).

- **Edge particles** cut by the image edge are not measured; the number excluded is shown.
- **Agglomerate tagger**: click an object on a micrograph, tag it agglomerate or rod. The model is refitted on all tags and
  recounts the picture, with five-fold cross-validation shown. Used from 20 tags per kind; below that a fixed rule applies.
- **Rods**: length and width per particle; film rim excluded; large compact clumps flagged as agglomerates.
- **Precise mode, sub-20 nm** (experimental, opt-in): model fit per particle. Not validated on real images yet.
- **Annotation filter** for scale bars and labels (opt-in).
- New info document "TEM particle sizing (v2.21)" with five figures.

Known limits:
- Precise mode and rod measurement have no manual reference yet.
- The agglomerate tagger is not yet placed on the TEM page.
- On the 100 and 200 nm sphere series, particles are clusters; use the sphere mode there.

Full list: app/CHANGELOG.md.
