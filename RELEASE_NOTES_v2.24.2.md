# Reticula 2.24.2

## Correcties
- Versie 2.24.2. De tag v2.24.1 is mislukt op de versie-controle (app/mfs/__init__.py stond op 2.24.0); er is niets gepubliceerd.
- Reticon-index-test werkt nu op Windows: de index staat vast op LF via .gitattributes.

## Nieuw
- Reticon: de CIF-symmetrieparser gebruikt geen eval meer (safe affine parser).
- PXRD: Vegard-paneel met vier bronnen voor de eindleden.

## Bekende beperkingen
- Zie RELEASE_NOTES_v2.24.0.md voor de Vegard-beperkingen (√3-supercel strict xfail; upload-route en Estimate x niet getest).
