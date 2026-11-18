## Reticula 2.16.0

**Downloads**
- Windows 10/11: `Reticula_Setup_2.16.0.exe` (not code-signed - SmartScreen: *More info* → *Run anyway*)
- macOS, Apple Silicon: `Reticula_2.16.0_macOS_arm64.dmg` (not notarised - first start: right-click → *Open*; see `README_macOS.txt`); the Mac build kit is attached as well.
- `SHA256SUMS.txt` to verify the downloads. Everything else: run from source (see the README).

**New**
- **PXRD, almost automatic**: one click identifies the phases (library of 18 reference structures) and refines cell, peak shape, atom positions and displacement parameters, with bond restraints that keep the geometry sensible. Validated on two real laboratory patterns.
- **Results you can publish**: bond lengths / angles with esd, coordination, bond-valence sums, a built-in **3D viewer**, a refined **CIF with esd**, a structure figure and a Methods paragraph.
- **Unknown phase**: indexing finds the unit cell, Le Bail checks it.
- **Home page redesigned**, a **problem report** button, **sample summary with the source of every number**, and **cross-checks between techniques**.

**Limits, stated plainly**: the structure refinement needs a starting structure (library, your CIF, or an analogue) – it refines, it does not solve unknown structures from scratch. Monoclinic / triclinic indexing and ab-initio structure solution are not included. Atom positions from laboratory data are meaningful for heavy-atom frameworks and robust fragments; light atoms (H, C, O next to heavy metals) are held by the restraints.

Full list: `app/CHANGELOG.md`.
