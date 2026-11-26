# Reticula 2.23.0

## Nieuw
- **Info-PDF's voor alle 24 modules** (bouwscript `tools/build_info_pdf.py`, reproduceerbaar: twee bouwslagen geven identieke bytes).
- **TEM: waarschuwing bij metadata-afwijking.** Als het OCR-label van de schaalbalk meer dan 2 % afwijkt van de metadata-pixelgrootte, toont de app de metadata-waarde en de gemeten balklengte. Labels en waarden veranderen niet.
- **TEM: scale-check en extra regressietests** (OCR-fix, cirkels, randuitsluiting, drempelcap, per-sample tags).
- **Reticon: stijlgids en icoonexport** (64-raster, platte schaduw, 24/48/96/192 px, fallback-renderer).
- **Structuur: gemmi als reserve-lezer voor CIF-fasen** (pymatgen eerst, gemmi bij fout).
- **Zetasizer 8.02: test voor zeta-only `.dts`-bestanden** (alleen met lokaal bestand).

## Correcties
- Helpdoc PXRD: grenswaarde kristalgrootte ("treat the size as a lower limit").
- Helpdoc: de standaardzoekopdracht is 500 punten en 2000 stappen; de figuur met 150/300/5 is de studie-instelling.
- Helpdoc: onbewezen claim over de HKUST-1-keten verwijderd.
- Grenstekst kristalgrootte in de UI.

## Bekende beperkingen
- Precieze deeltjesmodus blijft **experimenteel**: sub-20 nm-resultaten niet publiceren; geen validatie tegen handmatige metingen.
- Rods-overlap: aanrakende staven worden als één object gemeten (bekende beperking, vastgelegd in een test).
- Structuur: MIL-53 O–C-validatie en HKUST-1-keten zijn niet gevalideerd.
- Rigaku `.HO`-reader: energie-as niet gekalibreerd; geen elementtoewijzing.
- XRF-piekzoeker: experimenteel; alleen synthetisch getest.

## Zetasizer
- De test met het lokale ZrO2-bestand (privé, niet in de repo) wordt overgeslagen als het bestand ontbreekt; de verwachte waarden horen bij een ander monster en zijn niet bevestigd.

## Rigaku .HO-reader
- Experimenteel: de energie-as is niet gekalibreerd; geen elementtoewijzing. Beide spectrum-kopieën worden gelabeld; referentie is een gebruikerskeuze.
