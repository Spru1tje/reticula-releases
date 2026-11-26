# Reticon proof set: source index

Branch `reticon/biorender-style`. All icons below are drawn for Reticula from scratch.
No external artwork, no BioRender or other commercial sets, no CC0 sources used.

| Icon | Domain | Source | License | Grid |
|---|---|---|---|---|
| pxrd_pattern | PXRD | own drawing | own work (Reticula) | 64 |
| tem_instrument | TEM | own drawing | own work (Reticula) | 64 |
| tem_lattice | TEM | own drawing | own work (Reticula) | 64 |
| sem_particle | SEM | own drawing | own work (Reticula) | 64 |
| mof_pore_frame | MOF / porous materials | own drawing | own work (Reticula) | 64 |
| sorption_isotherm | Sorption | own drawing | own work (Reticula) | 64 |
| nmr_spectrum | NMR | own drawing | own work (Reticula) | 64 |
| ftir_spectrum | FTIR | own drawing | own work (Reticula) | 64 |
| zeta_potential | Zetasizer | own drawing | own work (Reticula) | 64 |
| structure_cell | Structure (unit cell) | own drawing | own work (Reticula) | 64 |

Exports: `app/assets/icons_proof/png/<name>_<size>.png` at 24, 48, 96 and 192 px,
made with `app/tools/export_icons.py` (headless Edge).
Contact sheet: `app/assets/icons_proof/contact_sheet.png`.

Migration of the 24-grid navigation icons (`app/assets/nav_icons`, still used by the app and left unchanged):
- nmr -> nmr_spectrum; ftir -> ftir_spectrum; zeta -> zeta_potential; structure_cif -> structure_cell.
- sorption -> sorption_isotherm (already in this set; no new drawing).
