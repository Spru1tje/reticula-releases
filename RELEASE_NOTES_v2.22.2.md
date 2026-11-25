# Reticula 2.22.2 — hotfix: TEM scale bar, circle exclusion, per-sample tags

## Fixes
- **Scale bar on TEM micrographs**: the contrast stretch of the reader erased thin scale bars, so e.g. `x.tif` showed
  "Scale bar not found". Calibration now reads the unstretched grey values; the bar (372 px for 100 nm on `x.tif`) and
  its label are found again.
- **Automatic threshold on textured films**: capped at the 99.5th percentile of the contrast; textured film images no
  longer give zero particles.

## New
- **Exclusion circles** on the rods measurement: click the centre of an area (for example empty carbon film), pick a
  radius and press "Add circle". Particles whose centre lies inside are not counted. "Clear all circles" removes them.
- **Rods**: particles cut by the picture edge and particles in the film rim are not counted; a minimum area can be set.
- **Agglomerate tags per sample**: a new image resets the page; the tags of a sample are kept and only deleted after an
  explicit button press.

## Known limits
- Precise sphere mode is still **experimental**: the radius is not validated against manual sizing on real micrographs.
  Do not publish results from it.
- The OCR of the scale-bar label can misread digits (e.g. 100 as 400); the metadata calibration is used when they
  disagree, and a warning is shown.
