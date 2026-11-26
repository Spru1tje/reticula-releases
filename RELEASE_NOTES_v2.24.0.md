# Reticula 2.24.0

## Nieuw
- **PXRD: Vegard-schatting uit meerdere pieken**; tetragonale en hexagonale eindleden (lineair).
- **Structuur: indexering met cel-controle**, niet-I-afwezigheden, en een controle van cel-only versus de structuurroute.
- **PL-vensters: minimum aantal tellingen per venster (100, synthetisch)**, gedeelde ratio-helper en signaalvlaggen; reproduceerbaar controlescript voor de vensterfit.
- **Reticon: navigatienaam-test**, Edge-versie in het contactvel en de log, palet-tokens bij export.
- **TEM: regressietests** (testbestand met LF-regeleinden), polydispersiteitstest en edge-profieltest voor de precieze bolmodus (synthetisch). De testomgeving achter de CI-totalen is gedocumenteerd.
- **Info-documenten**: de figuurbesturing heet overal "🎛️ Customise this figure"; een test controleert dat getallen in de info-documenten een bron hebben (plafond 3 bronloze getallen).

## Correcties
- Helpdoc: verouderde naam "Figure controls" vervangen (Getting started, PXRD).
- Helpdoc: onbronnige getallen uit de tekst verwijderd of gecorrigeerd (PL-rss, UV/Vis-tolerantie).

## Bekende beperkingen
- De precieze deeltjesmodus blijft **experimenteel**; sub-20 nm-resultaten niet publiceren.
- TEM-helptekst en TEM-PDF zijn nog concept; drie TEM-getallen zonder bron wachten op een besluit.
- Playwright-tests (browser-e2e) draaien alleen met `playwright` geïnstalleerd.
