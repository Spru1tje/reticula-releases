# Reticon v2.20: professional icon overhaul (plan)

Goal: scientific, publication-grade vector icons in the style of BioRender: restrained, precise, detailed, consistent. Not playful.

## Style rules (to be enforced by the generator, not by hand)
- One stroke weight per icon family (1.5 px on the 64 grid for outlines; no decorative doodles).
- Two-tone palette only: neutral slate for structure, one accent (steel blue) for the sample or the signal. No rainbow.
- Real proportions: vessels, instruments and crystals drawn from measured geometry, not caricature.
- Every icon is pure vector, no masks, images or text inside the artwork. Label text is added by the app.
- Legible at 24 px: detail is removed, not shrunk. A QA step rejects icons whose smallest feature is below 1.5 px.

## Generator
- Parametric families (silhouette x content x detail), so 5000 icons come from a small set of audited building blocks.
- Each family has a reference spec and a golden SVG; a regression test compares new output to the golden file.
- Variant rules are kept from the current library (6 shapes per icon: base + 5 variants).

## Phases
1. Style sheet and generator core, with 20 reference icons reviewed by eye at 24, 48 and 96 px.
2. Families: glassware, instruments, crystals and materials, spectra, process flow, biology, energy devices (about 40 families).
3. Scale to 5000 icons in batches of 500, each batch QA-checked (geometry, stroke, contrast, variants).
4. Replace the current library in the app; update help documents and tests in the same release.
