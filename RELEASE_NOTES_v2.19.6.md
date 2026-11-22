# Reticula 2.19.6

## Sliders with a typed number
- Every slider in every module now has a number box next to it. The slider stays; typing a number in the box moves the slider, and dragging the slider updates the box. A range gets two boxes.
- Typed values outside the range are set to the nearest end of the range; a range is kept in order.
- The spectrum selector of the PL thermometry page takes a typed index as well.
- Help text: "Ranges and numbers" in the Getting started document.

## Checks
- Known-answer tests: `tests/test_v2196_numentry.py` (typed number to slider, slider to box, clamping, ordered ranges, index box, values set by other code).
- Unit suite on the release commit: 788 passed, 50 skipped, 7 excluded OCR tests (the SEM scale-bar check is not part of the gate).
