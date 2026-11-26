# Icon style guide (Reticula, 2.20 proposal, branch reticon/biorender-style)

Scientific, restrained, for figures meant for publication. Not playful. Modelled on
the clean, flat look of scientific illustration tools (BioRender-like), built only from
our own drawings. No BioRender or other commercial icon sets are copied.

## Grid and stroke

- Grid: 64 x 64 units for all Reticon icons (decided by the coordinator; supersedes the 24 and 60 grids in 2.19.0).
- Stroke: one line weight per icon, 1.5 units. Round caps and joins.
- Minimum detail: no feature smaller than 1.5 units. Detail is removed for small sizes, not shrunk.

## Palette (fixed tokens, same on light and dark)

- ink `#1F2D3D`: outlines and text-like marks;
- accent `#2E6F95` (steel blue): the one element the icon is about;
- muted `#8FA3B5`: secondary parts (axes, supports, grid lines, shadow);
- paper `#FFFFFF`: fill for solids. Nothing else.

## Shadow

- Flat, soft offset shadow only: a copy of the icon's main shapes, moved by (2, 2),
  in muted at opacity 0.3: filled where the shape is a filled solid, stroked where it is a line.
  Drawn first (behind the artwork).
- No gradients, no blur filters, no SVG filter elements (they break in PNG export).

## Not allowed

- No gradients, glossy highlights, sunbursts, emoji-like forms.
- No text or labels inside the artwork. No embedded images or masks.
- More than one accent element per icon.

## Files and naming

- Source: one SVG per icon, `app/assets/icons_proof/<domain>_<object>.svg` for the proof set.
- Names: lowercase, underscores, domain first (`tem_lattice`, `pxrd_pattern`, `mof_pore_frame`).
- Each icon's source URL and license are recorded in the index (to be added with the library).

## Quality check (per icon)

1. viewBox is 0 0 64 64; the file is valid XML.
2. One stroke weight (1.5); round caps and joins.
3. Only the four palette tokens; no gradients, filters or text.
4. Smallest feature at least 1.5 units; still readable at 24 px.
5. Contrast of ink on paper at least 3:1 (ink on white is about 14:1).
6. Before/after contact sheet at 24, 48, 96 and 128 px.
7. Tests pass with the CI command before any commit.

## Export

- PNG export at 24, 48, 96 and 192 px from the SVG source (renderer to be fixed in the pipeline).
- Proof renders are made with headless Edge; the contact sheet is `icons_proof/contact_sheet.png`.
- Every icon renders to PNG and passes the checks above before it is added to the library.
