"""v2.18.3: the info documents of the Help and info page (one per module). English. Every abbreviation is written out where it is first used.
Sources: the user guide (app/README.md), the program's own code and tests. Advice that is general practice (not something the program checks)
is marked 'General practice'. A document = {"title", "subtitle", "sections": [(kind, text), ...]}, kinds: h p b n f k (see helppdf.py)."""
from __future__ import annotations
import os

GP = "General practice (not checked by the program): "
ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets", "docs")


def _doc(title, sub, *sections):
    return {"title": title, "subtitle": sub, "sections": [s for s in sections]}


def H(t): return ("h", t)
def P(t): return ("p", t)
def B(t): return ("b", t)
def N(t): return ("n", t)
def F(t): return ("f", t)
def K(t): return ("k", t)
def I(file, caption, height=0.3): return ("i", {"path": os.path.join(ASSETS, file), "caption": caption, "height": height})


def terms(*items):
    return [H("Terms used in this document")] + [B(i) for i in items]


DOCS: dict[str, dict] = {}

DOCS["Getting started"] = _doc(
    "Getting started with Reticula", "Read this first: it applies to every module.",
    H("What Reticula does"),
    P("Reticula turns raw instrument files into journal-style figures and does the usual analysis on the way. It runs on your computer; your data are never uploaded. "
      "Every technique has its own module in the sidebar. They all work the same way."),
    *terms("Supporting information (SI): the extra material that goes with a paper.",
           "Crystallographic information file (CIF): the standard text file with a crystal structure.",
           "Arbitrary units (a.u.): intensities that have no absolute scale.",
           "Dots per inch (dpi): the resolution of a picture."),
    H("Tutorial: from file to figure"),
    N("Pick a technique in the sidebar (grouped by Structure, Porosity, Spectroscopy, Thermal, Microscopy, Kinetics / Bio, Tools)."),
    N("Upload one or more files, or tick 'Use the bundled EXAMPLE data' to try the module first. The example files are synthetic."),
    N("Rename the samples and put them in the order you want in the sample list. Colours follow that order (colour-blind-safe Okabe-Ito palette by default). "
      "A colour rule in the sidebar keeps a colour for a sample in every figure."),
    N("Choose the processing (baseline, smoothing, normalisation, offset) and the figure width (ACS single column 3.25 in, double column 7 in, Wiley 8.5 or 17.5 cm)."),
    N("Open 'Figure controls' under the figure to change title, axes, units, ticks, legend and colours. Presets can be saved as JSON; the sidebar holds a global house style."),
    N("Download the figure as SVG, PDF, TIFF (600 dpi), PNG (600 dpi), the processed data as CSV, or everything as one ZIP. "
      "If your browser blocks downloads, use 'Save to folder'."),
    N("Use 'Add to combined figure' to put panels of different techniques into one figure, and tick 'Collect figures for the SI bundle' to get the SI export."),
    H("Ranges and numbers"),
    P("Every slider (for example a wavelength range, a fit window, a shift tolerance or a histogram bin count) has a number box next to it. "
      "Drag the slider, or type the value in the box and press Enter; the slider and the box always show the same value. "
      "A typed value outside the slider range is set to the nearest end of the range, and a range given as two boxes is kept in order."),
    H("Detail levels"),
    P("The sidebar switch 'Detail level' has three settings. Standard shows only what a publishable figure needs. Advanced adds fit options, corrections, ticks and presets. "
      "Developer adds diagnostics. You can set it for all modules or for one module. "
      "Not every page follows the level yet: the Structure (CIF) tables and two Microscopy sections do; the other technique pages are being audited."),
    H("Simple mode"),
    P("Under a figure, Figure controls → Style & colours, 'Simple mode' makes the figure black only and tells the samples apart by shape: one black symbol per sample for data "
      "points, one line style per sample for continuous curves, grey bars with a hatch per series. Legends follow. Switch it off to get the colours back."),
    H("Style rules"),
    B("Arial or Helvetica 7–9 pt, inward ticks on four sides, no grid, no title by default (each can be changed in Figure controls)."),
    B("Units in parentheses, for example 2θ (°) or Wavenumber (cm⁻¹). Unit conversion only relabels the ticks; the CSV keeps the original units."),
    H("Limits"),
    B("The file readers were built from published format descriptions and tested on synthetic and some real files. If an import fails, export plain text as described in the Export guide; two-column text always works."),
    B("If the Arial font is missing, Liberation Sans or DejaVu Sans is used: check fonts in the SVG/PDF before submission."),
)

DOCS["PXRD"] = _doc(
    "PXRD – powder X-ray diffraction", "Patterns, simulated references, crystallite size, Rietveld and structure refinement.",
    H("What it does"),
    P("Plots one or more patterns (stacked or overlaid), simulates a reference pattern from a crystal structure, and can fit the pattern: crystallite size from peak broadening, "
      "Rietveld / Le Bail refinement, a one-click identification and structure refinement, and (new in 2.18.4) indexing and structure solution for an unknown compound."),
    *terms("Powder X-ray diffraction (PXRD): intensity versus the diffraction angle 2θ of a powder.",
           "Full width at half maximum (FWHM): the width of a peak at half its height.",
           "Crystallographic information file (CIF): the text file that holds a crystal structure.",
           "Rietveld refinement: fitting the whole pattern with a crystal structure. Le Bail fit: the same without a structure; the intensities are free, only cell and peak shape are fitted.",
           "Residual weighted profile (Rwp) and goodness of fit (GoF): numbers that tell how well the calculated pattern matches the measured one.",
           "Estimated standard deviation (esd): the uncertainty of a refined value.",
           "Instrument profile: the peak width the diffractometer itself adds. Caglioti parameters U, V, W (Gaussian part) and X, Y (Lorentzian part) describe it as a function of angle."),
    H("Tutorial: a figure"),
    N("Upload the patterns (STOE WinXPOW .xy/.dat/.txt/.csv; Rigaku .rasx/.ras). Order and rename them in the sample list."),
    N("Choose the processing: background removal (asymmetric least squares, rolling ball or polynomial), smoothing, normalisation, offset."),
    N("Optional: add a simulated pattern from a CIF (Cu Kα1, pseudo-Voigt peaks) as a grey trace or sticks, and highlight bands."),
    N("Adjust the figure in Figure controls and download."),
    H("Tutorial: crystallite size and microstrain (opt-in)"),
    N("Switch on 'Crystallite size and microstrain'. Peaks are fitted automatically."),
    N("Choose the instrumental broadening: typical laboratory value (0.05°), a standard you measured (LaB₆ or Si), or a value you know."),
    N("Read the Scherrer, Williamson–Hall and Halder–Wagner results with their uncertainties."),
    F(r"$D = \dfrac{K\,\lambda}{\beta\,\cos\theta}$"),
    P("D is the crystallite size, K the shape factor (about 0.9), λ the wavelength, β the sample broadening in radians and θ the Bragg angle. "
      "The sample broadening is the measured width with the instrument's contribution removed (Gaussian: β² = β_obs² − β_inst²; Lorentzian: β = β_obs − β_inst)."),
    K("A too small instrumental width gives a too large size. The program says so when the width was only assumed (0.05°) and when the sample width is as narrow as the instrument (size above about 100 nm: treat the size as a lower limit)."),
    H("Tutorial: Rietveld / Le Bail (opt-in)"),
    N("Switch on 'Rietveld / Le Bail refinement'. Choose the phases from the library or upload CIF files. Or press 'Identify the phases'; or the 'Do everything automatically' button."),
    N("Choose the method: Rietveld (structure known), Rietveld + refine the atoms, or Le Bail. Choose the 2θ range."),
    N("Press Refine. A pattern takes about a second. Read the figure, Rwp, GoF, cell ± esd and weight fractions."),
    F(r"$R_{wp} = \sqrt{\dfrac{\sum_i w_i\,(y_i^{obs} - y_i^{calc})^2}{\sum_i w_i\,(y_i^{obs})^2}}$"),
    P("Advanced level adds: background polynomial order, monochromator angle, instrument offset (sample displacement, zero shift or both), preferred orientation, and bond-restraint strength."),
    H("Background, glass holder and standard (Advanced)"),
    B("The background is a polynomial (default order 6, range 2–14). On top of it, broad humps between 15° and 35° 2θ model the amorphous signal of a glass sample holder; this is automatic and can be switched off."),
    B("You can load a measured background (an empty holder pattern); it is scaled and subtracted in the refinement."),
    B("Instrument profile: refined from the pattern by default. Or fixed from a Si standard (SRM 640) measured on the same instrument: the standard is fitted with Le Bail at its certificate lattice parameter (5.4312 Å, not refined) and U, V, W, X, Y are then kept fixed for the sample."),
    H("Solving the structure of an unknown compound (new in 2.18.4)"),
    P("For a compound that is in no library and has no CIF, the expander 'Solve the crystal structure of an unknown compound' follows the route of the field and shows every decision. "
      "You give the chemical formula (the electron-density peaks cannot name the elements by themselves) and, if you have one, a silicon reference pattern."),
    *terms("Systematic absences: reflections that a space group forbids; their absence in the pattern points to the space group.",
           "Charge flipping: a method that finds the missing phases of the diffraction amplitudes: random phases are improved in a loop that makes the electron density positive and atom-like.",
           "Formula units per cell (Z): how many times the chemical formula fits in one unit cell.",
           "R(F): the agreement between the amplitudes of the model and those taken from the Le Bail fit; smaller is better.",
           "Instrument zero shift: the angle error of the diffractometer, measured on a standard such as silicon."),
    N("Data check: number of peaks, angle range and step. Resolution to d ≈ 1.2 Å (Cu Kα: 2θ ≈ 80°) is needed to place atoms; below that only the heavy atoms appear."),
    N("Unit cell: the indexing proposes several cells; each one is fitted with Le Bail on the real pattern. The lowest Rwp wins; a cell that fits no better than a smaller cell of the same lattice is not wanted."),
    N("Space group: the isolated reflections that a group forbids must be absent. Groups that agree are ranked (most confirmed absences first; centrosymmetric before non-centrosymmetric). The best few are tried until the model fits well."),
    N("Structure solution: charge flipping with the space group imposed, on 12 (default) random starts in parallel threads, with a progress bar. The amplitudes of overlapping reflections are shared and re-partitioned from the density itself."),
    N("Atoms: the density peaks become symmetry-independent positions; your formula gives the elements. The program chooses the number of formula units, the elements and the positions that agree best with the intensities (R(F)), adds missing atoms from a difference Fourier map and refines the positions."),
    N("Refinement: Rietveld refinement of the model with soft bond restraints (the same module as the manual refinement). You get the 3D structure, the CIF, the fit, and a ZIP with the pattern files, reflection lists and a report of every decision."),
    F(r"$R(F) = \dfrac{\sum_k \left|\sqrt{G_k^{obs}} - K\sqrt{G_k^{calc}}\,e^{-B s_k^2}\right|}{\sum_k \sqrt{G_k^{obs}}}$"),
    P("G is the total intensity of one group of reflections with the same d-spacing after the Lorentz-polarisation correction, K a scale factor, B an overall displacement parameter and s = sin θ / λ."),
    B("What is checked: the zero shift from the silicon reference is applied; the cell is chosen on the real pattern, not on the peak list alone; space groups are tested with the matched-filter intensities of isolated reflections (a bump below 0.3 % of the strongest line is not counted as a reflection); the model must agree with the composition."),
    B("Measured on synthetic patterns of known structures (public-domain COD files, Cu Kα1, 2θ 8–90°): the atoms were recovered for SrTiO₃, rutile, ZnO, anatase and tetragonal ZrO₂. Harder cases (several light atoms next to a heavy atom, low symmetry, large cells) are less reliable: see the limits."),
    B("Limits: a model from a powder pattern is a starting model. Light atoms and hydrogen are often missing; guest molecules and solvent in the pores are not found; monoclinic cells need the option 'Also try monoclinic cells'; triclinic cells are not searched; non-standard cell choices of monoclinic groups are not generated; "
      "the number of formula units is estimated from the cell volume when you do not give it. " + GP + "check bond lengths and the coordination of the metal, and compare the density with the measured one before using the structure."),
    H("What you get after a refinement (new in 2.18.4)"),
    P("Every Rietveld or Le Bail result, also one started by hand, ends with the same outputs under the figure:"),
    B("A 3D viewer of each phase in the refined cell (rotate, show or hide the cell box and bonds)."),
    B("A CIF per phase with the refined cell. With plain Rietveld or Le Bail the atoms are the starting model and were NOT refined; the file says so in a comment. With 'Rietveld + refine the atoms' the refined positions are written."),
    B("Pattern files as plain text for other programs: observed, calculated, difference and background as .xy, one table (CSV) with all of them, and the reflection positions (ticks)."),
    B("The reflection list (h k l, d-spacing, 2θ, intensity) and a .hkl file for programs that read it."),
    B("One ZIP with all of the above."),
    H("Indexing: finding the cell of an unknown phase (improved in 2.18.4)"),
    P("Indexing finds the unit cell from the peak positions alone. The quality number is the de Wolff figure of merit M20: the inverse of the average difference between the observed and calculated line positions, "
      "scaled by the number of possible lines (above 10 is reliable, above 20 good)."),
    F(r"$M_{20} = \dfrac{Q_{20}}{2\,\langle\varepsilon\rangle\,N_{20}}$"),
    P("Q20 is the Q-value (1/d²) of the twentieth line, ⟨ε⟩ the mean absolute difference between observed and calculated Q over the first 20 lines, and N20 the number of calculated lines up to Q20. "
      "The calculation lets the program fit the zero shift of the pattern together with the cell, tests the centring types (P, I, F, C, A, B and R-type rules), and counts M20 over the lines that are indexed, "
      "so that one impurity line does not destroy the figure. Cubic, tetragonal, hexagonal / trigonal, orthorhombic and (option 'Also try monoclinic cells', slower) monoclinic cells are searched, "
      "in stages from small to larger Miller indices, with a progress bar."),
    B("The program calls a cell 'reliable' when M20 is at least 20 and at most a quarter of the lines are left unindexed. " + GP + "also check that the Le Bail fit with this cell is good; a peak list with many impurity lines gives wrong cells."),
    B(GP + "put a Si or LaB₆ standard in the sample, or measure it right before: a zero shift or sample displacement of 0.05° changes the cell enough to spoil the indexing."),
    H("What is checked"),
    B("The sample displacement and the zero shift are not refined together by default (they are strongly correlated and bias the cell); 'both' needs a standard."),
    B("The program warns when the instrument broadening was not removed from the apparent size, and when the model misses something (large difference curve)."),
    B("Checked on synthetic data: the fixed Si profile gives the right cell for a second material. The fitted profile is less reliable above about 80° 2θ: check it against your own standard."),
    H("Limits"),
    B("Binary STOE .raw files need the optional xylib package; export .xy instead."),
    B("Indexing covers cubic, tetragonal, hexagonal and orthorhombic cells, and monoclinic cells on request; triclinic cells are not searched."),
    B("Only Si is built in as a standard profile; LaB₆ is not yet."),
)

DOCS["Structure (CIF)"] = _doc(
    "Structure (CIF) – crystal structure drawings and tables", "From a CIF to figures, SI tables and a simulated pattern.",
    H("What it does"),
    P("Reads a CIF and draws the structure, makes the tables for the SI, simulates a powder pattern and gives some pore and surface-contact quantities."),
    *terms("Crystallographic information file (CIF): the text file that holds a crystal structure (from SHELXL, Olex2, COD or CSD).",
           "Oak Ridge thermal-ellipsoid plot (ORTEP): atoms drawn as ellipsoids of their thermal motion (30, 50 or 70 % probability).",
           "Asymmetric unit (ASU): the smallest part of the cell from which symmetry builds the rest.",
           "Estimated standard deviation (esd): the uncertainty of a bond length or angle.",
           "Hirshfeld surface: a surface around a molecule that shows how it touches its neighbours."),
    H("Tutorial"),
    N("Upload a CIF."),
    N("Choose the drawing style (ORTEP, ball-and-stick, stick, space-filling), the view (along an axis, a [uvw] direction, a (hkl) plane or best view), and what to show (ASU, molecules, packing, supercell, hydrogen bonds)."),
    N("Download the figure, or the SI tables (crystal data, bond lengths, angles, torsions with esd) as DOCX or XLSX."),
    N("Optional: simulated powder pattern with hkl labels and overlay on your measured pattern; MOF pore quantities (Di, Df, void fraction, accessible surface area); Hirshfeld fingerprint (experimental)."),
    H("Limits"),
    B("The drawings here have no colour controls, because the colours carry information; Simple mode does not apply."),
    B("A NU-1000 CIF is not bundled; load your own."),
)

DOCS["Structure solution (unknown compound)"] = _doc(
    "Structure solution from powder data",
    "Charge flipping, metal test, net choice, rigid-body linker search and carbon recovery. Validated on simulated UiO-66 data only (see the limits).",
    H("What it does"),
    P("For a compound that is in no library and has no CIF, the solver tries to find where the metal atoms and the linker sit in the cell. "
      "It starts from indexed powder intensities (|F|², the squared structure-factor amplitudes). The inputs are the measured reflections, the cell and the space group."),
    K("Current scope: frameworks with Zr₆ clusters (fcu net) and a benzene-dicarboxylate linker are validated. The metal test also routes paddlewheel frameworks to tbo, but their linker search is not validated."),
    I("struct_fig_pipeline.png", "Figure 1. The route. Steps 1–3 and 6 use only the data. Step 5 needs a linker shape and a topology from the library.", 0.2),
    H("Step by step"),
    N("Charge flipping (Oszlanyi–Suto): the observed |F|² are turned into a density map by iteration. The strongest density peaks are the heavy atoms."),
    N("Metal test: the 48 strongest peaks are checked for pairs. Paddlewheel metals (Cu₂) give nearest-neighbour distances of 2.3–2.9 Å; a single-metal framework does not (UiO-66 Zr: about 3.5 Å)."),
    N("Net choice: a paired metal set selects tbo, the only paired net in the library. Otherwise fcu, pcu or dia are ranked by fit."),
    N("Rigid-body search: the linker (ideal shape, centre, three rotation angles, and one radius per Zr face for the μ₃-O atoms) is moved by simulated annealing to minimise the R factor."),
    N("Carbon recovery: the linker images at the best pose give the carbon positions. These are compared with the reference to report the matched fraction."),
    I("struct_fig_metal_pairs.png", "Figure 2. Nearest-neighbour distances of the heavy-metal peaks. HKUST-1 (teal) has a clear group at 2.5 Å, the Cu–Cu pair of a paddlewheel. The UiO-66 Zr control (grey) lies at about 3.5 Å. The shaded band is the paired window.", 0.26),
    H("Choosing the net from intensities"),
    P("A powder pattern shows only the reflections that are present. Mixed-parity reflections (for example 100 or 110 in a face-centred cell) are absent for fcu but present for pcu."
      "Without these absences, a list of observed peaks cannot tell the two nets apart. The matcher therefore includes all reflections up to the resolution, with zero where nothing is observed."),
    I("struct_fig_absences.png", "Figure 3. Node-model intensities of fcu and pcu for the lowest reflections. The absent mixed-parity reflections are what separates the two nets.", 0.25),
    F(r"$R_{\mathrm{net}} = \sqrt{\sum_h (|F_o| - s\,|F_c|)^2 \,/\, \sum_h |F_o|^2}$"),
    P("This is the R factor on amplitudes, with a fitted scale s. The node model places one scatterer per node and has no atomic detail, so these values only rank nets; they are not a fit quality."),
    I("struct_fig_net_ranking.png", "Figure 4. R values of the candidate nets. For UiO-66, fcu beats dia by 0.09 (confident). For HKUST-1 the metal test selects tbo.", 0.27),
    H("Rigid-body search"),
    P("The linker is treated as a rigid molecule. The cost is the scaled R factor on squared intensities, with the symmetry images of the linker. "
      "The search has two stages: a grid of orientations (evaluated in parallel), then simulated annealing from the best three to five orientations, with translation, rotation and face radii free."),
    F(r"$R_{F^2} = \sum_h |F_o^2 - k\,F_c^2| \,/\, \sum_h F_o^2$"),
    P("Over three seeds, a 150-point orientation grid with 300 steps and five refined starts gives R 0.196 (spread 0.001). The program's default search is a 500-point grid with 2000 steps."),
    I("struct_fig_speed.png", "Figure 5. Left: the same lead-sulfate refinement with float64 and float32 arithmetic (the program uses float32 in the peak profile): 89 s versus 48 s, Rwp 0.086712 versus 0.086710. Right: R over three seeds for four search settings. The grid of 150 points, 300 steps and five starts is the study setting, not the program's default (500 points, 2000 steps).", 0.27),
    H("What is checked"),
    K("Charge flipping recovers the Zr and O framework from simulated UiO-66 powder data with 5 % intensity noise (Zr at least 90 %, O at least 80 % matched within 0.5–0.7 Å). The two regression tests take about 8 s together."),
    B("The rigid-body fit on simulated UiO-66 data reaches R about 0.20 (the true model gives 0.19). A wrong topology gives 0.6 or more."),
    B("Carbon recovery: 100 % of the 192 carbons are within 0.7 Å (median 0.05 Å) for the UiO-66 pipeline, with Zr from the charge-flip peaks and not from the CIF."),
    I("struct_fig_validation.png", "Figure 6. UiO-66 carbon recovery. The pattern is simulated from the reference CIF, so this checks the method, not a measurement.", 0.2),
    B("Net matcher: UiO-66 selects fcu (R 0.690, against dia 0.783 and pcu 0.932). HKUST-1 selects tbo through the metal test. Five regression tests cover the metal test, the complete reflection list and the routing."),
    B("The refinement that follows the solver (Rietveld and structure model) is checked against real laboratory patterns: fluorapatite agrees with the GSAS-II reference within 0.002 in fractional coordinates, and lead sulfate passes its tests."),
    H("Limits (read before publishing)"),
    K("The solver has not been validated on a measured pattern. All its results come from simulated data generated from the same reference CIF. Measured data of a Zr₆ MOF and an independent test structure are still needed."),
    B("Sensitivity: a Zr position error of 0.07 Å raises R from 0.23 to 0.41; 0.04 Å (the measured charge-flip error) gives 0.23. The 0.35 pass threshold was set after seeing the results."),
    B("The tbo net nodes were taken from the HKUST-1 CIF, so the HKUST-1 net choice is not an independent test."),
    B("HKUST-1 carbon recovery and MIL-53 have not been run with this pipeline. Only UiO-66 has the full chain."),
    B("The chance level of the search has not been measured against random restarts. The 5–10 % figure in an earlier note is an estimate."),
    B("The symmetry shortcut in the cost (signed permutations) is valid only for cubic symmetry. The solver does not check this, so use it only for cubic cells."),
    P(GP + "Before publication: report Rwp and the R factors per run, use at least three seeds, compare with an established method (for example TOPAS or GSAS-II) on the same data, and check the final CIF with checkCIF."),
    H("Terms"),
    B("|F|² and |F|: squared and absolute structure-factor amplitude of one reflection (hkl); Fo observed, Fc calculated."),
    B("CIF: Crystallographic Information File, the standard text format of a crystal structure."),
    B("R factor: the disagreement between observed and calculated intensities; a lower value means a better fit. Rwp: the weighted profile R factor of a Rietveld fit."),
    B("Net (topology): the way the nodes of a framework connect. fcu: face-centred cubic; pcu: primitive cubic; dia: diamond; tbo: the net of HKUST-1."),
    B("Zr₆ cluster: six zirconium atoms forming an octahedron. Cu₂ paddlewheel: two copper atoms bridged by four carboxylate groups. μ₃-O: an oxygen bridging three metals (the face of a Zr₆ octahedron)."),
    B("Charge flipping: an iterative density method (Oszlanyi–Suto) that finds heavy atoms from |F|² without a starting model."),
    B("GSAS-II and TOPAS: established Rietveld programs used as reference methods. checkCIF: the IUCr validation service for a CIF."),
)

DOCS["TEM particle sizing (v2.21)"] = _doc(
    "TEM particle sizing: modes, edge particles, annotations, rods and agglomerates",
    "How the particle measurement works, what it has been checked against, and where it is not yet validated.",
    H("What it does"),
    P("Measures the size of particles on a TEM micrograph after the scale is set. Round particles are measured with a circle model, "
      "elongated particles with their length and width, and clumps (agglomerates) are recognised and can be left out. "
      "Every number is reported together with its method, so the origin of each result can be stated in a paper."),
    K("Status: the default circle model agrees with a manual sizing on the test images (209 nm manual; 209.4 nm measured on 52 particles). "
      "The other modes are not yet validated on real images, see the limits at the end."),
    *terms("TEM: transmission electron microscopy. Micrograph: the image. Scale bar: the bar with its length in nm, used for calibration.",
           "Agglomerate: a clump of particles that touches or fuses, measured as one object.",
           "Solidity: the area of an object divided by the area of its convex hull; 1 means convex, lower values mean irregular or clumped.",
           "Aspect ratio: length divided by width. Rod: an elongated particle, length well above its width.",
           "Equivalent diameter: the diameter of a circle with the same area.",
           "Voronoi cell: the area of the image that lies closer to one particle centre than to any other."),
    H("The modes"),
    B("Automatic (default, 'spheres'): round particles are fitted with a circle on their real edge. Touching and overlapping spheres are split. "
      "On the test images this gives 209.4 nm against 209 nm manual."),
    B("Outline only: the thresholded outline of every object. Fast (0.3 s on a 370 × 282 px image). On the same test images it gives 218.5 nm (+4.5 %). "
      "On a dense synthetic network of known spheres it underestimates by about 38 %, so use it only for non-round particles."),
    B("Precise, sub-20 nm (experimental): each particle is fitted with a blurred absorbing-sphere model on its own Voronoi cell. "
      "On synthetic spheres of known radius the error is about 0 %, but on the real 20k image it gives 256.4 nm (+23 %). "
      "Not for publication until calibrated on a manual sizing of small particles. Time: about 6 s on the 370 × 282 px test image."),
    I("tem_fig_modes.png", "Figure 1. Left: real image, mean diameter per mode against the manual 209 nm. Right: synthetic spheres of known radius, "
      "our own model (not an independent check). The precise mode is not yet validated on real images.", 0.3),
    H("Edge particles"),
    P("A particle cut by the edge of the image is not measured, because its size is unknown. The circle model used to judge a particle "
      "fitted only the visible arc, so a cut particle could look complete. The rule therefore also checks the rays from the centre: when a ray "
      "leaves the image before it finds the particle edge, the particle is cut and is excluded. The number of excluded particles is shown in the module."),
    H("Annotations (scale bar, labels)"),
    P("The scale bar, the panel letter and the text next to the bar are dark objects that are not particles. The annotation filter finds "
      "the scale bar (a long, almost black bar), removes the text in the band just above it, and removes dark panel letters in the outer 12 % of the image. "
      "It is opt-in: it removed real dark particles on other images, so it is not on by default."),
    I("tem_fig_annotations.png", "Figure 2. The publication micrograph (scale bar 100 nm = 108 px). Left: without the filter the letter, the scale bar and "
      "the text are counted. Right: with the annotation filter. Some fragments at the ends of the bar remain.", 0.3),
    H("Elongated particles (rods): length and width"),
    P("For rod or worm-like particles each object is measured by its second moments: the length is the full major axis and the width the full minor axis "
      "of the ellipse with the same moments. The image is first normalised by the local film level, so the grey film texture does not give objects. "
      "Objects in the dark rim of the carbon-film holes are flagged 'film edge' and not counted. Objects that are large and compact, or very irregular, are flagged as "
      "agglomerates and not counted (large = more than max(150 nm², 6 × the median area of the image))."),
    I("tem_fig_rods_overlay.png", "Figure 3. Overlay on the first rod-type micrograph. Red: counted rods. Green: agglomerates, not counted. "
      "The rule removes the large clumps; chains of rods that are still connected can remain.", 0.3),
    I("tem_fig_rod_table.png", "Figure 4. Median length and width for the nine micrographs of the archive, each with its own scale. "
      "On the 100 and 200 nm series the particles are clusters of spheres; those values are not single-particle sizes.", 0.3),
    H("Agglomerate tagger (learned)"),
    P("You can tag objects yourself: click an object on the overlay and mark it 'agglomerate' or 'rod / particle'. The model is refitted on all tags "
      "and then recounts the whole picture, so agglomerates that resemble your tags are also left out. "
      "The tags are stored per object (its measured features), so tags from several micrographs are combined."),
    I("tem_fig_tagger.png", "Figure 5. The tagging loop. The model is a logistic regression on six measured features "
      "(area, length, width, aspect, solidity, relative area) with an L2 penalty.", 0.22),
    B("The learned model is used only with at least 20 tags of each kind. Below that, the fixed rule above is used."),
    B("Its quality is shown as five-fold cross-validated accuracy and balanced accuracy on your tags, with the coefficients per feature."),
    B("A refit takes a few milliseconds; the segmentation is cached, so a click takes about 0.03 s on the test image."),
    K("Not yet on the TEM page: the tagger is a module and is not yet placed on the TEM page. Its first real use, on the rod images, is still to be checked."),
    H("Limits and what is needed for publication"),
    B("The circle model is checked against one manual sizing of 209 nm on the test images. Sub-20 nm particles have no manual reference yet."),
    B("The precise mode and the rod measurement are not validated against manual sizes. They are experimental."),
    B("The learned model can only use the six measured features. It cannot see that a clump is next to other clumps."),
    B("Images 5 to 9 of the archive (100 and 200 nm series) are clusters of spheres; they need the particle-level model, not the rod measurement."),
    P(GP + "For a publication: measure 20 to 30 particles by hand on each micrograph type, report the method per number, state the number of excluded edge and "
         "agglomerate objects, and keep the scale calibration in the same units. Validate the learned model on tags from images that were not used for training."),
)

DOCS["N2 sorption (BET/PSD)"] = _doc(
    "N₂ sorption – isotherms, BET and pore size distribution", "From the instrument export to surface area, micropore volume and pore sizes.",
    H("What it does"),
    P("Plots nitrogen isotherms and calculates the BET surface area with the Rouquerol criteria, the t-plot micropore volume and external area, and pore size distributions."),
    *terms("Brunauer–Emmett–Teller (BET) theory: the monolayer capacity from the adsorption isotherm gives the surface area.",
           "Relative pressure (p/p₀): the pressure divided by the saturation pressure.",
           "Pore size distribution (PSD): how the pore volume is spread over pore widths.",
           "Density functional theory (DFT), Horváth–Kawazoe (HK), Saito–Foley (SF), Barrett–Joyner–Halenda (BJH): methods to calculate a PSD from the isotherm.",
           "t-plot (with the Harkins–Jura thickness): separates micropore volume from external area.",
           "Standard temperature and pressure (STP): the reference for gas volumes (cm³ g⁻¹ STP)."),
    H("Tutorial"),
    N("Upload the MicroActive Excel export (.xls/.xlsx, full report or the Isotherm Tabular table) or a CSV. The 'File structure' expander shows which tables were found."),
    N("Plot the isotherms (adsorption filled, desorption open; linear and/or logarithmic p/p₀)."),
    N("Read the BET table: the app's Rouquerol BET next to the value in the file, and the app's value on the file's p/p₀ range. Check the Rouquerol criteria panel."),
    N("Choose a PSD method (DFT kernel fit, HK, SF or BJH) and read the t-plot values."),
    F(r"$\dfrac{1}{v\,[(p_0/p)-1]} = \dfrac{c-1}{v_m\,c}\,\dfrac{p}{p_0} + \dfrac{1}{v_m\,c}$"),
    P("v is the adsorbed quantity, v_m the monolayer capacity and c a constant; the straight-line fit over a p/p₀ window gives v_m, from which the area follows."),
    H("What is checked"),
    P("The Rouquerol criteria (ISO 9277) choose the fitting window: (1) the constant c must be positive; (2) the quantity v·(1 − p/p₀) must keep increasing with p/p₀ in the window; "
      "(3) the relative pressure of the monolayer must lie inside the window; (4) it must agree with the value predicted by the BET equation. The app also shows the range of areas over all windows."),
    B("Ultramicroporous MOFs (MOF-801, UiO-66) often cannot meet criterion 3 or 4. The app flags this; then report a 'BET-equivalent area' and the p/p₀ range used."),
    B("The instrument usually fits a fixed window (0.05–0.30). For microporous MOFs the Rouquerol window lies at much lower p/p₀, so the two values differ by design."),
    H("Limits"),
    B("The DFT kernel is a carbon-slit N₂ 77 K kernel; it is not an oxide NLDFT/QSDFT kernel, so PSDs of oxides are approximate."),
    B("MicroActive .SMP/.QPS files cannot be read; use the Excel export."),
)

DOCS["FTIR / Raman"] = _doc(
    "FTIR and Raman – vibrational spectra", "Stacked spectra, band labels, band fitting and difference spectra.",
    H("What it does"),
    P("Plots infrared and Raman spectra, processes them, labels bands and fits overlapping bands."),
    *terms("Fourier-transform infrared spectroscopy (FTIR), here with attenuated total reflection (ATR): infrared absorption of a solid pressed on a crystal.",
           "Fourier-transform Raman (FT-Raman): scattered light with a 1064 nm laser.",
           "Transmittance (T) and absorbance (A): A = −log₁₀ T.",
           "Wavenumber (cm⁻¹): the horizontal axis, the number of waves per centimetre.",
           "Asymmetric least squares (ALS): a baseline-removal method."),
    H("Tutorial"),
    N("Choose the spectrum type: ATR-FTIR (transmittance), ATR-FTIR (absorbance) or FT-Raman."),
    N("Upload the files: Bruker OPUS binaries (.0, .1, …), OMNIC .spa (one spectrum per file), .dpt, or CSV/TXT. Header-less OMNIC CSV exports are read as one spectrum."),
    N("If a text file holds several spectra (a wavenumber column plus one column per spectrum), tick 'The file(s) contain several spectra'. Each column then becomes a spectrum, named by its header."),
    N("Set the range, the baseline (none for FTIR, ALS for Raman), normalisation, and highlight bands ('lo-hi:label')."),
    N("Advanced: normalise to a reference band, automatic peak labels, a zoomed inset, difference spectra, and a band-fit (assignment table and simultaneous multi-band fit)."),
    F(r"$A = -\log_{10} T$"),
    H("What is checked"),
    B("OMNIC .spa files were checked against their own CSV export: same points and intensities (wavenumber within 0.001 cm⁻¹); the .spa holds the full measured range."),
    B("An OMNIC .spg (a series) is not read."),
)

DOCS["UV/Vis"] = _doc(
    "UV/Vis – absorbance, assays, band gaps", "Spectra, Beer–Lambert calibration, assay presets, kinetics and Tauc band gaps.",
    H("What it does"),
    P("Plots ultraviolet–visible spectra, and calculates concentrations, assay results, single-wavelength kinetics and optical band gaps."),
    *terms("Ultraviolet–visible spectroscopy (UV/Vis): light absorption between about 200 and 800 nm.",
           "Absorbance (A), transmittance (%T) and reflectance (%R).",
           "Limit of detection (LOD) and limit of quantification (LOQ): the smallest amounts that can be detected or measured reliably.",
           "Kubelka–Munk function: turns diffuse reflectance into a quantity proportional to absorption.",
           "Tauc plot: a plot to read the optical band gap E_g from the absorption edge.",
           "Bicinchoninic acid (BCA), o-phthalaldehyde (OPA): reagents of protein and amine assays."),
    H("Tutorial"),
    N("Upload spectra (Implen NanoPhotometer Excel, Cary or generic CSV, JCAMP-DX)."),
    N("Choose the quantity (A, %T or %R), blank, smoothing and the maximum wavelength (λmax)."),
    N("Beer–Lambert: enter standards, fit the calibration (with LOD and LOQ) and read unknowns. Pick an assay preset (BCA, Bradford, Lowry, ninhydrin, OPA)."),
    N("Band gap: use Kubelka–Munk and a Tauc plot on reflectance data."),
    F(r"$A = \varepsilon\, c\, l$"),
    F(r"$F(R) = \dfrac{(1-R)^2}{2R}\qquad (F(R)\,h\nu)^{n}\ \mathrm{versus}\ h\nu$"),
    P("ε is the molar absorption coefficient, c the concentration, l the path length; R the reflectance, hν the photon energy. The exponent n depends on the type of transition (2 for direct, 1/2 for indirect)."),
    B(GP + "read the band gap from the linear part of the Tauc plot and state which n you used."),
    H("Automatic direct or indirect band gap (v2.18.11)"),
    P("The library can choose between the direct (n = 2) and the indirect (n = 1/2) Tauc form by the better automatic linear edge (higher R²; direct on a tie). Checked on a synthetic direct edge at 3.0 eV, which returns direct and 3.00 eV."),
    B("General practice: the rule is a heuristic. Look at both Tauc plots and state which transition you assume."),

)

DOCS["Photoluminescence"] = _doc(
    "Info implementatie ratiometric thermometry in Reticula", "Photoluminescence module: spectra, thermal quenching and ratiometric thermometry.",
    H("1. The idea in plain words"),
    P("Some luminescent materials (for example Er³⁺-doped upconversion particles or phosphors) emit light in two bands whose energy levels lie close together. "
      "When the material is heated, the upper level is populated more strongly, so the band at the shorter wavelength grows relative to the other one. "
      "We therefore do not use the absolute brightness (it depends on laser power, sample position and concentration) but the RATIO of the two bands. "
      "That ratio is the thermometer: ratiometric luminescence thermometry."),
    H("2. Terms, defined before they are abbreviated"),
    B("Photoluminescence (PL): light emitted by a material after it was excited with light."),
    B("Intensity of band 1 and band 2 (I1, I2): the area under the shorter- and the longer-wavelength emission band."),
    B("Ratio (Δ): I1 divided by I2."),
    B("Energy gap (ΔE): the energy difference between the two emitting levels, in cm⁻¹ or eV."),
    B("Boltzmann constant (kB): 0.695 cm⁻¹ per kelvin; the temperature T is in kelvin (K)."),
    B("Relative sensitivity (Sr): how many percent the ratio changes per kelvin."),
    B("Temperature uncertainty (σT): the smallest temperature change that can be told from the noise, in K."),
    B("Full width at half maximum (FWHM): the width of a band at half of its height."),
    B("Bayesian information criterion (BIC): a score that compares two models and punishes extra parameters; the lower score wins."),
    B("Coefficient of determination (R²) and reduced chi-squared (χ²red): numbers that tell how well a fit describes the data."),
    B("Standard deviation (SD), arbitrary units (a.u.), supporting information (SI): used in the usual way."),
    H("3. The formulas"),
    F(r"$\Delta = \dfrac{I_1}{I_2} = \alpha\,\exp\!\left(-\dfrac{\Delta E}{k_B T}\right)\ (+\,C)$"),
    B("α is a scale factor of the material; C is an optional constant (tick 'Add a constant C')."),
    P("Relative sensitivity (percent per kelvin):"),
    F(r"$S_r = 100\%\cdot\left|\dfrac{1}{\Delta}\dfrac{\partial\Delta}{\partial T}\right| = 100\%\cdot\dfrac{\Delta E}{k_B T^{2}}$"),
    P("Temperature uncertainty with photon counts:"),
    F(r"$\sigma_T = \dfrac{1}{S_r}\sqrt{\dfrac{1}{I_1}+\dfrac{1}{I_2}}$"),
    P("If the error of the ratio itself is known (from the global fit), the general form is σT = (σΔ / Δ) / Sr. The absolute sensitivity Sa = |dΔ/dT| = (Δ − C)·ΔE / (kB T²) is calculated as well."),
    B("Consequence: Sr falls with 1/T². A thermometer is more sensitive at low temperature, and a larger ΔE gives a higher Sr."),
    H("4. The steps (tab 'Thermometry (band ratio)')"),
    N("Load the spectra: one table (wavelength plus one column per temperature) or separate files. The temperature is read from the column header or file name ('293.15 K', '25C', 'T=300K'); numbers above 150 are read as K, smaller ones as °C. At least 3 temperatures are needed."),
    N("Find the bands: automatically (the two strongest separate bands, also hidden shoulders) or choose two integration windows yourself."),
    N("Determine the ratio at every temperature (section 6)."),
    N("Fit the Boltzmann law to Δ versus T. Result: ΔE (cm⁻¹ and eV, with error), α, R²."),
    N("Calculate Sr and σT per temperature and draw the figures: emission map, Δ versus T with fit, Sr(T), σT(T). They can go into the combined figure."),
    H("4b. Which files are read (new in 2.18.4)"),
    B("Edinburgh Instruments `.FS` files (written by the Fluoracle software): a single scan, or a 'Multiple Scans' file with one scan per temperature. Every scan becomes one spectrum; the sample name, "
      "the detector (for example near-infrared photomultiplier tube, NIR PMT) and the excitation source are read from the file, and the temperature is taken from the sample name ('..._25C...')."),
    B("Text and Excel files from other instruments (Horiba, Ocean Optics, Hitachi …): two columns per file, or one table with one column per temperature."),
    B("Wavelength ranges: visible (VIS, below about 750 nm), visible plus near-infrared (VIS-NIR) and near-infrared (NIR, above about 850 nm) are recognised from the data and shown in the background block."),
    H("5a. Background subtraction before the fit (new in 2.18.4)"),
    P("Before any ratio is formed, the background that the experiment adds to the emission is removed. It has several sources: the detector pedestal (dark counts; in the NIR several thousand counts under a "
      "lanthanide spectrum is normal), stray light and laser leakage, second-order light of the grating, and broad emission of the solvent, cuvette or host. A pedestal that is left in changes the ratio, "
      "and the ratio changes with temperature in a way that looks like a physical effect, so the removal must be standard, visible and the same for every spectrum of the series."),
    P("Setting 'Background' under the spectra (default: Automatic). What the automatic mode does, step by step:"),
    N("Find the signal-free windows of the whole series: a region is free when no spectrum of the series rises above its own lower envelope (statistics-sensitive non-linear iterative peak clipping, SNIP) "
      "by more than a few noise levels. Short runs and the edges of each run are dropped. For a lanthanide series measured from 900 to 1250 nm these are the stretches between and beside the emission bands."),
    N("Fit a polynomial in these windows only (Chebyshev basis, Tukey-bisquare re-weighting so that stray peaks do not pull the fit). Order 0 is a pure pedestal, order 1 adds a slope, order 2 a curve. "
      "The order is chosen with the Bayesian information criterion (BIC) added up over all spectra; a higher order must lower the BIC by more than 10 to be used, and the order is capped at 2 "
      "because a higher order bends in the gap between windows (inside an emission band, where nothing checks it)."),
    N("Subtract the SAME plan (windows and order) from every spectrum; each spectrum gets its own coefficients. The figure shows the raw spectra, the background (dashed) and the windows (green)."),
    N("Check: the mean left in the signal-free windows should be zero within one noise level (σ). The block shows this number; a warning appears if it is not."),
    P("Other settings, for when you know better: Manual windows (you type the windows, and choose polynomial, rubber band, ModPoly, arPLS or SNIP), Dark / background file (a measured blank with the same "
      "detector and settings, scaled in the signal-free windows and subtracted), Constant: minimum, Constant: mean of a region, None. 'Why this method?' lists all candidates with bias and scatter in the windows."),
    B("When no signal-free window exists (emission over the whole range) the rubber band is used: the straight-segment floor under the data (lower convex hull). It never rises above the data, so it cannot swallow a broad band; "
      "it removes only the floor."),
    B("Spikes (cosmic rays, hot pixels of array detectors): a Hampel filter replaces points that differ from the median of their neighbours by more than 8 noise levels. Off by default; not needed for a PMT scan."),
    B("The global fit of section 6 still fits its own small polynomial under the bands (section 5). The two layers do not fight: the first removes the instrument background over the whole range, the second "
      "absorbs what is left under the bands."),
    P("What influences ratiometric thermometry, and what Reticula does about it:"),
    B("Detector pedestal and slope (corrected, above)."),
    B("Spectral response of the spectrometer (grating and detector efficiency): corrected when you load a correction file or a calibration lamp (measured and certified spectrum). A response that does not change during the series multiplies every ratio by the same factor: it changes α and the absolute ratio, NOT ΔE or Sᵣ. It matters when the calibration is moved to another instrument, and for colour coordinates."),
    B("Water absorption near 975 and 1195 nm in the NIR changes band ratios in aqueous samples: a note appears for NIR spectra. " + GP + "compare with a dry sample or D₂O."),
    B("Self-absorption and reabsorption at high concentration, laser heating, drift of laser power and detector with time: not corrected. " + GP + "check that the ratio does not change with power or with repeated scans."),
    B("Photon-counting noise: enters σT (section 3). The size of the background raises the noise: a high pedestal lowers the precision even when it is subtracted correctly."),
    H("5. The background polynomial"),
    P("Under the emission bands there is always a background: detector offset, stray light, or the tail of a neighbouring emission. It is described by a polynomial in wavelength x:"),
    F(r"$\mathrm{background}(x) = b_0 + b_1 x + b_2 x^2 + b_3 x^3$"),
    P("The order is the highest power used. Order 0 is a flat offset, order 1 a straight sloping line, order 2 a curve, order 3 an S-shaped curve."),
    B("Default: order 1 (a straight line). In the Advanced level the slider 'Background polynomial order' allows 0 to 3."),
    B("Every spectrum gets its own background over the fit window (the bands plus 3 FWHM on each side), built from Chebyshev polynomials (a numerically stable form of the same curve)."),
    B("The coefficients are solved exactly together with the band areas in every step of the fit (the problem is linear in them)."),
    B("The integration-window method uses no polynomial (see section 6 for its background)."),
    P(GP + "keep order 1. Use 0 for a flat baseline. Go to 2 only if the residual shows a smooth curve. A high order can absorb band tails and shift the ratio: check that the ratio does not change when the order changes by one."),
    H("6. Two ways to determine the ratio"),
    P("A. Global fit (recommended). All spectra of the series are described together by the SAME bands (pseudo-Voigt shapes with unit area). Position and width may drift linearly with temperature; areas and background of each spectrum are solved exactly. "
      "An extra 'nuisance' band (a shoulder or third peak) can be added so that it does not contaminate the ratio. The simple model (fixed shape) and the drifting model are compared with the BIC; the better one is used."),
    P("B. Integration windows. The classic method: the area in two windows above a straight line through the window ends. Drawback: overlapping bands are cut at an arbitrary valley and the result depends on your windows."),
    P("'Compare both' shows the two next to each other. If they differ strongly, the bands overlap and the global fit is preferred."),
    P("C. Manifold fit (new in 2.18.4), for crowded bands such as the Stark-split lines of a lanthanide. Every line between the limits of the fit range is found (maxima and shoulders from the second derivative) and gets its own pseudo-Voigt shape, shared by the series. "
      "You drag two ranges: the lines whose centre lies in range 1 are summed into I₁, those in range 2 into I₂, the others are nuisance lines. Lines left in the residual are added while the Bayesian information criterion improves. "
      "When the two-band fit flags more than 30 % of the spectra, the page proposes this method."),
    F(r"$I_1 = \sum_{k \in R_1} A_k,\quad I_2 = \sum_{k \in R_2} A_k,\quad \sigma^2_{I_1} = \sum_{k,l \in R_1} \mathrm{cov}(A_k, A_l)$"),
    P("The error of the ratio uses the full covariance of the line areas (overlapping lines are anti-correlated, which the sum partly cancels)."),
    B("The integration windows integrate above ZERO when the standard background subtraction is on (a straight line through the window ends would cut into neighbouring lines), and above the straight line through the window ends when it is off. Their counting error is calculated from the raw counts."),
    B("The ratio is put the right way round automatically: when I₁/I₂ falls with temperature, Δ = I₂/I₁ is used (band 2 is then the thermally populated level) and the page says so."),
    H("6b. Quality report and confidence (new in 2.18.4)"),
    P("Above the results the page shows a verdict – HIGH, MODERATE or LOW confidence – collected from every check, each with the reason and what to do. LOW means: do not report these numbers before the problem is solved; they are still shown, marked."),
    B("Data: detector saturation (a maximum at the converter limit of an array detector, or three or more points with exactly the maximum; photon-counting detectors have no converter limit), settings that differ inside the series."),
    B("Information: the change of the ratio over the temperature range compared with its error (slope of ln Δ against 1/T divided by its error, t). Below 3 there is no reliable thermometer; 3–10 is weakly determined. The ratio must be monotonic (rank correlation)."),
    B("Fit: the share of spectra that the band model does not describe; the Boltzmann residuals (reduced χ²) and the runs test of their signs (long runs of one sign = a curvature that the law does not have)."),
    B("Physics: the fitted ΔE against the spacing of the two bands (centroids in cm⁻¹; for two thermally coupled levels they should agree within about 30 %), and against the usual range of thermally coupled levels (General practice: 200–2000 cm⁻¹; below 200 cm⁻¹ is Stark-level thermometry with a low sensitivity)."),
    B("An independent method: the integration windows give their own ΔE, which must agree with the fit within twice the combined uncertainty."),
    B("Experiment: heating against cooling branch (from the measurement order), reproducibility between runs (one file = one run), laser heating, the simulation self-check (consistency only)."),
    H("6c. The uncertainty budget (new in 2.18.4)"),
    P("The ± of ΔE is the square root of the sum of squares of every component:"),
    F(r"$u(\Delta E)^2 = u_\mathrm{stat}^2 + u_{T,\mathrm{rand}}^2 + u_{T,\mathrm{cal}}^2 + \sum_j (\Delta E_j - \Delta E)^2 + u_\mathrm{runs}^2 \,(+\, u_\mathrm{Birge}^2)$"),
    B("Statistical: Monte Carlo – the ratios are redrawn from their errors 2500 times and the Boltzmann law (a weighted straight line of ln Δ against 1/T) is refitted every time; u is half the 15.9–84.1 % range of the results. The ratio errors (version 2.18.5) include the uncertainty of the non-linear line shapes (centres, widths, η): their covariance from the global fit is propagated through the sensitivity of each ratio to them. Without this term the ratio errors were 10–20 % too small (checked against 60 repeated simulated series). At low counts (a few hundred counts per point) the reported ΔE error is about 15–25 % larger than the scatter of repeated simulations; this is the safe direction, and the cause (the first-order propagation of a ratio of noisy areas) is not removed by changing the fit loss or the weights."),
    B("Temperature sensor: the same Monte Carlo with a random error on every temperature, and with one common offset (sensor calibration, sample not at the sensor temperature). You type both values; zero means not counted."),
    B("Processing: the Boltzmann law is refitted on the ratios obtained with every reasonable alternative – other fit range or margin, the other peak model, the other background subtractions, wider and narrower windows. The difference in ΔE is what the CHOICE costs."),
    B("Runs: the standard deviation of ΔE between independent runs divided by the square root of their number. Birge ratio (version 2.18.5): the errors are scaled by √χ² only when the excess scatter is significant (χ² above 1 by more than two standard deviations of χ²); scaling on every χ² above 1 biased the errors upwards."),
    B("Ratio noise (version 2.18.6): when the reduced chi-squared (the fit misfit divided by the stated error) is above 10, the counting error is too small to be used as a weight. The fit is then repeated with a constant relative error equal to the measured scatter of the ratio, and the report shows a 'Ratio noise' warning. General practice, not checked by the program: look at the integration time, the spectral-response correction and whether the spectra were averaged or normalised before the ratio was taken."),
    B("Poisson maximum likelihood (version 2.18.6): on counts the ratio is fitted to the raw counts of the two integration windows, with the background counts from the background plan. Each point is weighted by its Poisson error (the variance of a count equals the count), no subtraction is needed before the ratio, and no point is dropped. The uncertainty of the energy gap comes from the observed information (the curvature of the likelihood). The result is shown next to the subtraction result as a check. General practice, not checked by the program: the windows should contain the whole band."),
    B("Emission-masked baseline (version 2.18.6): when no signal-free window exists, the baseline is a straight line (order 1 at most) fitted only outside the emission. The emission is found from the data by iteration. This replaces the rubber band as the default, because SNIP and arPLS remove the tails of a broad band (measured area errors +13 % and −26 %). Order 2 is not used: it extrapolates into the band (measured −30 %). The rubber band stays as an alternative. A safety check falls back to the rubber band, and says so in the notes, when the emission covers the whole range and the masked baseline would rise above the spectrum by more than three noise units. General practice, not checked by the program: a blank or solvent spectrum is the reliable background for broad emission."),
    B("Publication checks (version 2.20.4): the report now also tests (1) whether the simple Boltzmann law is enough or a constant offset is needed, by AIC (a difference above 2 prefers the offset, above 10 strongly); (2) whether a spectral response correction was applied; and (3) the accuracy against an independent reference thermometer, if you type its temperatures in the validation box (bias, RMS and largest deviation; 1 K ok, 3 K warn, above that fail, as a convention). Without a reference the report says so: the absolute accuracy is then not shown. The PDF has a Methods section (what was done, in plain steps, with a figure of the analysis chain), a model figure (both laws and their AIC) and, when references are given, a validation figure with the 1 K and 3 K bands."),
    B("Analysis report (PDF, version 2.18.6): next to the CSV and the JSON record, the thermometry tab offers a PDF with the verdict, the settings, the equations, the calibration figure, the table, the uncertainty budget, every check with its advice, and the limitations. The report states what the program computed; advice the program does not check is marked as general practice."),
    B("Two temperature uncertainties are given: the resolution σT,res = (σΔ/Δ)/Sᵣ, the precision of ONE reading; and the calibration uncertainty σT,cal, the accuracy of the curve T(Δ) at that temperature (Monte Carlo over ratios and sensor). Report both."),
    F(r"$T(\Delta) = \dfrac{\Delta E}{k_B \ln\!\left(\alpha / (\Delta - C)\right)},\qquad \sigma_{T,\mathrm{res}} = \dfrac{\sigma_\Delta / \Delta}{S_r}$"),
    H("6d. Experiment checks and the self-check (new in 2.18.4)"),
    B("Heating against cooling: when the series was measured up and back down, the cooling ratios are compared with the heating branch (interpolated in ln Δ against 1/T) and turned into an apparent temperature difference ΔT = Δ(ln Δ)·k_BT²/ΔE. More than twice the resolution (at least 1 K) is flagged."),
    B("Laser heating: spectra at one set temperature and several powers (power in the file name) are analysed with exactly the same method; the temperature read from the calibration is fitted linearly against the power. The slope is the heating per watt."),
    B("Simulation self-check (version 2.18.5): six series like yours are simulated from the same band model with your fitted ΔE and analysed again. This tests the statistics and the code path (bias, and whether the error bars contain the known value in about 95 % of the runs). It cannot detect a wrong model: if the band model is wrong for your sample, the check still passes."),
    B("Analysis record: a JSON file with the program version, the input files and their MD5 fingerprints, every setting, the results, the budget and the checks – for the supporting information and to reproduce the result later."),
    H("7. What is checked"),
    B("Error of the ratio = statistical error (including the covariance of the two areas) + systematic error. The systematic part is the spread of refits with narrower or wider windows, with the other model, and (new in 2.18.4) with the other reasonable background subtractions "
      "(the other polynomial orders and shortened windows, or the other flexible methods); the part due to the background choice is shown on its own column. Total = √(statistical² + systematic²)."),
    B("The fit works on the background-subtracted spectra; the statistical weights still use the raw counts (signal plus background), because the Poisson noise of the background does not disappear when it is subtracted."),
    B("The integration windows integrate above zero after the standard background subtraction, and above a straight line through the window ends without it (section 6)."),
    B("Quality flags per spectrum: R², χ²red and the autocorrelation of the residuals. Badly fitting spectra are marked: do not trust those points."),
    B("The program looks for structure left in the residuals and proposes an extra band only if the BIC improves (at most 3 extra bands)."),
    B("The ratio must be positive for every spectrum; otherwise a warning appears and the Boltzmann fit stops."),
    B("σT follows from the statistical error of every ratio: from the raw counts for photon counts, from the point-to-point noise of the spectrum otherwise (leave 'counts' unticked for intensities in a.u.)."),
    B("Saturation: the extras panel measures the power dependence (slope n < 1 at high power means saturation)."),
    H("8. The other tabs of the module"),
    B("Spectra: coloured by temperature, waterfall, emission map; nm, eV or cm⁻¹ axes; peak position and FWHM."),
    B("Temperature dependence: thermal quenching (activation energy E_a ± its error, T₅₀ ± its error by Monte Carlo over the correlated parameters, Arrhenius plot). A warning appears when fewer than three temperatures lie on the drop of the intensity, and when T₅₀ is an extrapolation."),
    B("Spectra (new in 2.18.4): a data-check panel (settings read from the files, saturation), optional labels of typical lanthanide / Cr³⁺ transition ranges (a guide, never an assignment), normalisation per second when the dwell time is in the files, and the removal of the excitation line and its second order (scattered light) before anything else."),
    B("Colour (Advanced): CIE 1931 x, y, the colour temperature and (new in 2.18.4) the dominant wavelength and colour purity against a chosen white point; purple colours give the complementary wavelength. Colour needs response-corrected spectra; the page says so."),
    H("8b. Decay curves and lifetimes (new in 2.18.4)"),
    P("Choose 'Decay curves (lifetimes)' at the top. Files: time in the first column (unit from the header), one column per decay. Three models: a tail fit after the peak (the excitation is short compared with the lifetime), rise and decay (the emitting level is fed with a delay, as in upconversion), and reconvolution with a measured instrument response function (IRF) for lifetimes close to the instrument response."),
    F(r"$I(t) = \sum_i a_i\, e^{-(t-t_0)/\tau_i} + B,\qquad \bar\tau_\mathrm{amp} = \dfrac{\sum a_i \tau_i}{\sum a_i},\qquad \bar\tau_\mathrm{int} = \dfrac{\sum a_i \tau_i^2}{\sum a_i \tau_i}$"),
    B("Weights: photon counts → 1/counts (Neyman χ², slightly biased below about 10 counts per channel). The amplitudes and the background are solved exactly (non-negative least squares) for every trial set of lifetimes."),
    B("The number of exponentials: 1 to 3 are fitted, and a higher number is kept only if the information criterion drops by more than 10, the lifetimes are separated (ratio above 1.5) and every amplitude is above 2 %. Quality flags from χ² and the Durbin–Watson statistic of the residuals."),
    B("Lifetime against temperature: τ(T) = τ₀/(1 + C·exp(−Eₐ/k_BT)) (Mott–Seitz) with Eₐ ± error and the relative sensitivity Sᵣ = |dτ/dT|/τ."),
    B(GP + "state which average lifetime you report; keep the count rate below about 5 % of the repetition rate (pile-up); measure the IRF with a scatterer at the emission wavelength settings."),
    H("9. Practical tips"),
    B("A crowded band system (for example the overlapping Stark lines of Yb³⁺ between 940 and 1050 nm) is not two bands plus a few extras: the two-band fit flags most spectra (wavy residual). Use the manifold fit and drag the two ranges; then read the quality report."),
    B(GP + "use as many temperatures as possible over a wide range; a Boltzmann fit through 3 points is weak."),
    B(GP + "correct the spectrometer response (detector and grating) when the calibration must work on another instrument or when you report colour; it does not change ΔE or Sᵣ."),
    B(GP + "measure the series up and back down, repeat it in two or three runs, and fill in the temperature-sensor uncertainty: these are the checks that make a calibration trustworthy."),
    B(GP + "laser heating falsifies the temperature: check that the ratio does not change with power."),
    B(GP + "report ΔE, Sr at a reference temperature and σT, not only the curve."),
    H("10. What changed in version 2.18.6, and the gain"),
    P("Three changes make the ratio and the background more reliable. The gains below are measured on synthetic spectra with a known answer (see the figures). They are not claims about one particular material."),
    I("pl_fig_ratio.png", "Figure 1. Ratio fit on 300 synthetic series with a true energy gap of 700 cm⁻¹ (80–400 K, 50 background counts per window). Left: share of the results inside the 1σ and 2σ intervals (ideal: about 68 % and 95 %). Right: bias of the energy gap.", 0.3),
    B("New: a Poisson maximum-likelihood fit on the raw counts of the two windows. Before: the background was subtracted, points with zero or negative signal were dropped, and a logarithm was fitted with weights. Gain: the bias fell from about −2.7 % to about 0 %; the share inside 1σ rose from about 29 % to about 69 %, and inside 2σ from about 57 % to about 95 %. The error bars now have the right size."),
    B("New: the maximum-likelihood result is shown next to the subtraction result as a check. If the two disagree by more than twice their combined uncertainty, the program warns."),
    B(GP + "the maximum-likelihood fit needs counts (not normalised intensities), and the windows should contain the whole band."),
    I("pl_fig_background.png", "Figure 2. Area error of one band with two background methods (40 synthetic spectra each, the true area is known). Shaded: ±2 %. Case D (a concave background) is not corrected by either method.", 0.3),
    B("New: when no signal-free window exists, the baseline is a straight line (order 1 at most), fitted only outside the emission. The emission is found from the data by iteration. Gain: on a broad band over the whole range the area error fell from about +6 % (rubber band) to about −2 %; on sharp peaks on a slope from about +78 % to about +1 %; on noisy spectra from about +89 % to about 0 %."),
    B("Limit: a concave background (a dip in the middle) is corrected by neither method (about +715 % and +485 % in Figure 2). The program does not yet detect this case. For such a background, use a blank or solvent spectrum."),
    B(GP + "a blank or solvent spectrum, measured under the same settings, is the reliable background for broad emission."),
    K("Summary: where a signal-free window exists, nothing changes. Without one, the emission-masked polynomial replaces the rubber band as the default. The energy gap is reported with the maximum-likelihood uncertainty, and every analysis can be saved as a PDF report."),
    H("11. Version 2.19.5: a faster and more exact manifold fit"),
    P("The manifold fit (every Stark line of the two manifolds fitted in all spectra at once) took 11 s on the 7-spectrum test series and 74 s on a crowded "
      "series of 13 spectra with 10 lines. It is now 12–40 times faster, finds the right minimum where the old fit got stuck, and reports error bars "
      "that match the real scatter. Every number below was measured on one computer with the real code of the older versions; the figures are made by "
      "tools/make_pl_speed_figures.py from tools/pl_speed_2195_data.json."),
    *terms("Jacobian: the matrix of the derivatives of every residual with respect to every shape parameter (centre, width, Lorentzian share η, drifts). The optimiser needs it in every step.",
           "Variable projection: the areas and the background are solved exactly for given shapes, so the optimiser only sees the shape parameters (Golub and Pereyra, 1973).",
           "Non-negative least squares (NNLS): a least-squares solution in which no area may be negative.",
           "Trust region reflective (TRF): the bounded least-squares optimiser of SciPy that the fit uses.",
           "Basic Linear Algebra Subprograms (BLAS): the library that does the matrix arithmetic.",
           "Karush–Kuhn–Tucker (KKT) conditions: the conditions that prove a solution with bounds is the optimum.",
           "Bayesian information criterion (BIC): the number that decides between the shared and the drifting shape, and whether a line is added.",
           "Residual sum of squares (rss): how far the fit is from the data; lower is better."),
    I("pl_fig_speed_versions.png", "Figure 3. Time of the whole manifold analysis of the thermometry tab (line search, both shape models, robustness refits with other windows, the other model and two other backgrounds) for five data sets, same computer, same settings.", 0.33),
    I("pl_fig_speed_table.png", "Table 1. The same times with the result of each version: the energy gap ΔE of the Boltzmann fit and the chosen shape model. The test series was made without drift: the old code chose the drifting shape only because its shared fit got stuck.", 0.2),
    H("12. Why it was slow"),
    P("A profile of the crowded series showed 77 % of the time in the Jacobian. It was made by finite differences: for every one of the 30–50 shape parameters "
      "the whole series was evaluated once more, and every evaluation solved an NNLS problem for every spectrum. One analysis needed 169 Jacobians, "
      "4 829 series evaluations and 36 456 NNLS solutions. The optimisations below remove this work, one after the other (Figure 4)."),
    I("pl_fig_speed_steps.png", "Figure 4. Time after each optimisation step (cumulative, median of three runs; the start point is v2.19.4). The numbers 1–8 are the steps in the list below.", 0.3),
    N("Exact Jacobian. The derivative of the projected residual is written out (variable projection) with the analytic derivatives of the pseudo-Voigt: one evaluation per step "
      "instead of one per parameter, and exact instead of approximate. The areas are solved by a QR decomposition; NNLS is used only when an area would come out negative "
      "(when every free area is ≥ 0, the free solution IS the NNLS solution). All spectra on a common wavelength grid are handled in one batched array operation. "
      "BLAS is limited to one thread during the fit: the matrices are small (a few hundred × ten) and several threads made each operation up to 5 times slower. "
      "The optimiser moves offsets from the start values, so its first step is one scale unit (the old first step could widen a line to 7 nm and switch another off, see Figure 5). "
      "The robustness refits run at the same time on several cores, and the other shape model of the automatic choice is reused instead of fitted again."),
    F(r"$J_p = -\left(P^{\perp}\,\dfrac{\partial M}{\partial \theta_p}\,a + (M^{+})^{T}\left(\dfrac{\partial M}{\partial \theta_p}\right)^{T} r\right)$"),
    P("M is the design matrix (one column per line plus the background), a the areas, r the residual, P⊥ the projection on the part of the data that M cannot describe, M⁺ the pseudo-inverse and θₚ one shape parameter."),
    N("Warm starts: a trial with one line more and the robustness refits start from the main solution instead of from the first guess."),
    N("Lines switched off by NNLS in some spectra: the derivatives of all spectra are still computed in one batch; only those spectra are projected one by one."),
    N("Start values that lie on a bound (η = 0 for a Gaussian line) are moved 0.1 % inside, because the optimiser scales its steps by the distance to the bound and crawled along it "
      "(80 evaluations without converging, now 8). A warm-started fit that does not converge is repeated from the ordinary start and the lower cost wins, so a warm start can only save time. "
      "Parameters held on a bound are left out of the error propagation (section 13)."),
    N("The drift columns of the Jacobian follow exactly from the others (∂M/∂s = t·∂M/∂c and ∂M/∂g = t·∂M/∂w, t is a constant per spectrum), so only three columns per line are computed. "
      "The pseudo-Voigt is computed in place; a difference column for a switched-off line recomputes only that line; the robustness refits skip the error terms they do not use."),
    N("Spectra with a negative free area are solved in one batch without that line, and the solution is accepted only if it passes the KKT test (then it is the unique NNLS solution; "
      "checked against NNLS on 1 400 spectra: equal to 5·10⁻¹⁵). A cold fit starts from a fit of the mean spectrum (one spectrum: a fraction of the cost). "
      "The convergence tolerance on the cost is 10⁻⁸ instead of 10⁻⁹ (8–14 % fewer evaluations; the ratios change by less than 3·10⁻⁹ of their error)."),
    N("The projection and the derivatives are written without temporary arrays (the same numbers to 10⁻¹⁵)."),
    N("Both shape models start from the same warm start and run at the same time. The drifting model contains the shared one, so it may never end above it; "
      "if it does, it is repeated from the shared solution. On the 41-spectrum series a trial with a spurious line no longer crawls (74 evaluations, now 29)."),
    B("On the page: with a fixed shape model the fit of the automatic choice is reused (it had already fitted both models); the laser-heating check refit starts warm."),
    H("13. Two errors of the old fit that were found and fixed"),
    P("Stuck fits. With the old start the first optimiser step could be huge: on the test series with the window 500–575 nm one line widened to 7.0 nm (true 3.8 nm) and its neighbour "
      "was switched off (area 0); the fit stayed there (rss 131 916 instead of 2 603; the ratio at 293 K was 0.2051, true 0.2012; up to 2.1 % off). "
      "Over 60 noisy copies of the test series the old fit got stuck in 8, which made the true scatter of its ratio 4.5 times larger (Figure 6A; worst ratio 4.2 % off, now 0.66 %)."),
    I("pl_fig_speed_stuck.png", "Figure 5. Left: one spectrum of the test series with the old (stuck) and the new fit, and their residuals. Right: the ratio of every spectrum against the true ratio.", 0.3),
    P("Error bars. The old finite-difference step was relative (SciPy's diff_step): for a parameter at zero (η of a Gaussian line on its bound, a drift near zero) the step became about 10⁻²¹ "
      "and the derivative was rounding noise (up to 10⁷ instead of about 10) or exactly zero. The covariance of the shape parameters, and with it the ratio errors, was then erratic: "
      "0.65–0.74 times the true scatter of the old results, and 4.1 times at one temperature (Figure 6B). With the exact Jacobian the errors were 1.33–1.43 times the true scatter, "
      "because a parameter held on a bound was treated as free. It does not scatter (in repeated noisy series it stays on the bound), so it is now left out: 1.05–1.16 times the true scatter. "
      "For lines with η inside its range nothing changes (0.89–1.02 times, Figure 6C)."),
    I("pl_fig_speed_errors.png", "Figure 6. Monte Carlo over noisy copies of the test series (known answer: the real scatter). A: scatter of the ratio. B, C: reported error divided by the real scatter; ideal 1, shaded: ±2 standard errors of the scatter itself.", 0.28),
    H("14. What stayed the same, and how it is checked"),
    B("The same lines in every data set, the same ratios within their errors, the same ΔE (Table 1). The test series now gets the shared shape and the systematic errors no longer come from stuck refits "
      "(test series 4.7·10⁻³ → 1.5·10⁻⁴; 41 spectra 5.9·10⁻³ → 3.3·10⁻⁵)."),
    B("Warm against cold starts over 38 noisy series (30 test series, 8 crowded series): the same model and lines in every case, the BIC never worse, ratios equal to 1.1·10⁻⁴ of their error; 19.5 s against 52.9 s."),
    B("Tests with known answers (tests/test_v2195_plfit.py): the exact Jacobian against finite differences (common grid, different grids, lines switched off, drift), the ratio gradient, "
      "noiseless data recovered exactly, the stuck minimum, the Monte Carlo error bars, the NNLS shortcut against NNLS, warm against cold, serial against parallel."),
    B("Switches for comparisons: MFS_PL_FD_JAC=1 runs the old finite-difference path, MFS_PL_WARM=0 switches the warm starts off, MFS_NO_PARALLEL=1 uses one core, MFS_PL_BLAS_THREADS=keep leaves BLAS alone."),
    H("15. Tried and rejected"),
    B("A compact singular value decomposition inside SciPy's optimiser (about 12 % faster): it needs a change inside SciPy; rejected for robustness. The decomposition is now the largest single cost (about 25 %)."),
    B("Two BLAS threads in the sequential part, an iterative trust-region solver (LSMR), difference rows only for the spectra where a line is off: all measured slower (LSMR also gave a worse BIC)."),
    B("Normal equations (Cholesky) instead of QR, or single precision: they square the condition number of overlapping lines; rejected for precision."),
    B("Stopping a trial fit early by its BIC: during a slow convergence the BIC still changed by about 2, which can change which lines are accepted; rejected."),
    B("Allowing negative areas instead of NNLS: it changes the physics; rejected."),
    B(GP + "compare times only on the same computer with the same settings; on one core the parallel parts run one after the other (the results are identical)."),
    K("Summary: the manifold fit is 12–40 times faster (crowded series 74 s → 1.9 s, 41 spectra 95 s → 3.1 s), it no longer gets stuck in a wrong minimum, and its ratio errors now match the real scatter."),
)

DOCS["1H NMR (digest)"] = _doc(
    "¹H NMR – spectra, assignment and digestion analysis", "From a Bruker folder to a labelled figure, a peak list, and linker : modulator ratios.",
    H("What it does"),
    P("Processes a one-dimensional proton NMR spectrum, finds the multiplets, assigns peaks to atoms of a structure you draw, and (in the digestion mode) turns integrals into ratios and concentrations."),
    *terms("Nuclear magnetic resonance (NMR): the radio-frequency response of nuclei in a magnetic field; ¹H is the proton.",
           "Chemical shift (δ, in ppm): the position of a signal relative to a reference.",
           "Coupling constant (J, in Hz): the splitting of a signal by neighbouring protons.",
           "Quantitative NMR (qNMR): concentrations from integrals compared with an internal standard.",
           "Trimethylsilylpropanoic acid (TMSP): a common internal standard (here TMSP-d4).",
           "Metal–organic framework (MOF): a porous crystalline material that is digested (dissolved) before measurement."),
    H("Tutorial: characterisation"),
    N("Upload a Bruker experiment folder as .zip (fid or processed pdata/1/1r), a JCAMP-DX file, a TopSpin 'totxt' export, or x,y text."),
    N("Phase, baseline and referencing are automatic; adjust if needed. Multiplets and J couplings are found."),
    N("Draw the structure in the built-in editor (or paste a molfile or SMILES). The program predicts shifts (nmrshiftdb2 table) and assigns peaks to atoms; drag to correct."),
    N("Check impurities against the Fulmer 2010 and Gottlieb solvent tables. Make the labelled figure and the ACS-style peak list for the SI."),
    P("'Simple mode (black)' draws the whole figure in black."),
    H("Tutorial: digestion and qNMR"),
    N("Switch to the digestion mode. Choose the integration windows for linker, modulator and the internal standard."),
    N("Read linker : modulator ratios and concentrations versus the standard. They feed the defect calculator."),
    F(r"$n_i \propto \dfrac{A_i}{N_{H,i}}\qquad \dfrac{n_i}{n_{std}} = \dfrac{A_i\,/\,N_{H,i}}{A_{std}\,/\,N_{H,std}}$"),
    P("A is the integral of a signal and N_H the number of protons it represents."),
    H("What is checked"),
    B("Phasing, baseline and referencing are automatic; the shift prediction is a database lookup, not a quantum calculation: check assignments by eye."),
    B(GP + "use a long enough relaxation delay (D1) for quantitative integrals."),
)

DOCS["XPS"] = _doc(
    "XPS – X-ray photoelectron spectroscopy", "Region spectra, background, doublet fitting and atomic percentages.",
    H("What it does"),
    P("Plots XPS regions (binding energy axis reversed), subtracts a background, fits doublets and gives atomic percentages."),
    *terms("X-ray photoelectron spectroscopy (XPS): the energy of electrons emitted by X-rays tells the elements and their chemical state.",
           "Binding energy (BE): the energy axis, in eV.",
           "Shirley and Tougaard backgrounds: two common models of the step background under a peak.",
           "Spin–orbit doublet: a pair of peaks (for example Zr 3d, Hf 4f, S 2p) with a fixed splitting and area ratio.",
           "Relative sensitivity factor (RSF): the correction that turns peak areas into amounts (Scofield values by default)."),
    H("Tutorial"),
    N("Upload CSV/TXT region files (CasaXPS, Avantage) or simple VAMAS (.vms)."),
    N("Choose the region, the background (Shirley or Tougaard) and the components. For a doublet the splitting and the area ratio are fixed."),
    N("Read peak positions and areas; give or edit RSFs to get atomic percent."),
    F(r"$x_i = \dfrac{I_i / \mathrm{RSF}_i}{\sum_j I_j / \mathrm{RSF}_j}\times 100\%$"),
    B("VAMAS support is minimal (single-block regions). Use the multi-region EXAMPLE to try the module."),
)

DOCS["ICP-OES / XRF"] = _doc(
    "ICP-OES and XRF – elemental analysis", "From concentrations to weight percent and molar ratios.",
    H("What it does"),
    P("Reads elemental analysis tables and calculates weight percent, mean ± standard deviation bars and molar ratios (also Zr : Hf : S with propagated uncertainty)."),
    *terms("Inductively coupled plasma optical emission spectrometry (ICP-OES): elements in a digested sample, as mg/L.",
           "X-ray fluorescence (XRF): elements in a solid from their characteristic X-rays.",
           "Standard deviation (SD): the spread over replicates.",
           "Metal node: the metal cluster of a MOF (for example M₆)."),
    H("Tutorial"),
    N("Upload a Syngistix CSV (long or wide) or any table with sample, element and value."),
    N("Enter the digested mass and the volume to convert mg/L to weight percent."),
    N("Choose the bars (one panel per element or grouped) and the molar ratios (pick the reference metal, per n atoms)."),
    N("The Zr : Hf : S case gives Hf/(Zr + Hf) and S per M₆ node with propagated uncertainty. The values go to the sample store and the defect calculator."),
    F(r"$w\,[\mathrm{wt\%}] = \dfrac{c\,[\mathrm{mg/L}]\cdot V\,[\mathrm{L}]}{m\,[\mathrm{mg}]}\times 100\%$"),
)

DOCS["TGA"] = _doc(
    "TGA – thermogravimetric analysis", "Mass loss curves, derivative, steps and the oxide residue.",
    H("What it does"),
    P("Plots mass versus temperature with the derivative, labels steps, normalises to the oxide residue and passes values to the defect calculator."),
    *terms("Thermogravimetric analysis (TGA): the mass of a sample while it is heated.",
           "Derivative thermogravimetry (DTG): the rate of mass change; its minima mark the steps.",
           "Oxide residue: the metal oxide (MO₂) that remains after combustion in air.",
           "Quadrupole mass spectrometry (QMS): a mass spectrometer coupled to the TGA to see what leaves."),
    H("Tutorial"),
    N("Upload Netzsch Proteus ASCII, TA Universal Analysis / TRIOS text, Mettler STARe export, or generic CSV (mg or %)."),
    N("Choose the DTG display (second panel or twin axis) and optional renormalisation to the oxide residue = 100 %."),
    N("Read the automatic step labels; set the plateau if needed. The values W(200 °C), plateau and end go to the defect calculator. A TGA–MS overlay is available."),
    B("The plateau is picked automatically at the DTG minimum between 250 and 480 °C. Check it, and override it in the defect tab if dehydroxylation and combustion overlap."),
    B(GP + "measure in air so that the end product is the oxide MO₂."),
)

DOCS["TGA"]["sections"] += [
    H("Uncertainty of a mass-loss step (v2.19.2)"),
    P("A step is the difference of two balance readings, so its mass-loss uncertainty from the balance resolution is √2 × resolution. With a 10 mg sample and 0.001 mg resolution the uncertainty is 0.0141 % of the initial mass; the mass lost (in mg) is the step % times the initial mass."),
    B("General practice: the balance resolution is a lower limit; buoyancy, purge-gas flow and sample weighing add to it."),
    B("Checked on synthetic data with a known answer (tests/test_v2192_items.py)."),
]
DOCS["DSC"] = _doc(
    "DSC – differential scanning calorimetry", "Peaks, glass transitions, cycles and heat capacity.",
    H("What it does"),
    P("Plots heat flow versus temperature or time and calculates transition temperatures, enthalpies, glass transitions, cycle comparisons and heat capacity."),
    *terms("Differential scanning calorimetry (DSC): the heat flow into a sample while it is heated or cooled.",
           "Onset, peak and endset: the start, top and end of a thermal event.",
           "Enthalpy (ΔH): the heat of a transition, from the area of the peak (J/g).",
           "Glass transition (Tg) and heat-capacity step (Δcp).",
           "Heat capacity (cp): heat needed per gram and kelvin; measured against a sapphire reference (ASTM E1269).",
           "Simultaneous thermal analysis (STA): TGA and DSC together."),
    H("Tutorial"),
    N("Upload Mettler STARe, TA TRIOS / Universal Analysis, Netzsch Proteus, PerkinElmer Pyris text exports, or generic CSV/XLSX."),
    N("Choose the units, the exo direction (up or down) and the segment colours."),
    N("Pick a peak and a baseline (linear or sigmoidal) to get onset, peak, endset and ΔH; or Tg (midpoint, inflection, Δcp)."),
    N("Full cycle: heating and cooling segments are found from dT/dt; compare Heat 1 with Heat 2; see the supercooling ΔT = Tm − Tc."),
    N("Heat capacity: load the blank, the sapphire and the sample; the ratio method gives cp."),
    F(r"$\Delta H = \dfrac{1}{\beta\, m}\int (\Phi - \Phi_{base})\,dT$"),
    P("Φ is the heat flow, β the heating rate and m the mass."),
)

DOCS["SEM / TEM / EDX"] = _doc(
    "SEM, TEM and EDX – images, scale bars and particle sizes", "Calibrated images, scale bars, particle statistics and element maps.",
    H("What it does"),
    P("Shows electron-microscope images with a correct scale bar, measures particle sizes and draws element maps."),
    *terms("Scanning electron microscopy (SEM) and transmission electron microscopy (TEM).",
           "High-resolution TEM (HRTEM) and bright-field TEM: image modes of a TEM.",
           "Energy-dispersive X-ray spectroscopy (EDX): element maps and spectra from X-rays.",
           "Pixel size (nm/px): the calibration that turns pixels into nanometres.",
           "Instruments named in the documentation: JEOL JEM-1010 (bright-field TEM, 100 kV, tungsten filament, CCD camera, 0.4 nm resolution, up to 500 000×) and JEOL JEM-2200FS (HRTEM, 200 kV field-emission gun, image-lens spherical-aberration corrector, electron energy loss spectrometer and EDX)."),
    H("Tutorial"),
    N("Upload the ORIGINAL image: TIFF (with the databar), TEM .dm3/.dm4/.ser/.emd/.mrc, or PNG/JPG with a manual pixel size. EDX maps are greyscale images, one per element."),
    N("Check the calibration panel: it compares the file metadata, the scale bar burned into the info strip (measured, with text recognition) and the magnification. If all fail, click the two bar ends and type its value."),
    N("Set the scale bar (value or length, position, colour, thickness, font, box), build the (a)–(d) panel grid, and use matched bars or one calibration for a batch."),
    N("Particle sizing: automatic recognition with background flattening; click to exclude or to measure; log-normal histogram. EDX maps can be overlaid on the SE image; line scans are available."),
    B("JEM-1010: save the camera TIFF at full resolution with the scale bar on (Olympus SIS / analySIS / MegaView III files are read). JEM-2200FS: save the raw Digital Micrograph or Velox file, not a screenshot."),
    B("The 127 mm magnification estimate is only valid for Philips/FEI XL SEMs."),
)

DOCS["N2 sorption (BET/PSD)"]["sections"] += [
    H("Gurvich total pore volume (v2.19.2)"),
    P("Gurvich rule (Gurvich 1915): the nitrogen uptake at a relative pressure near 0.95 is taken to fill all pores; the STP gas volume (cm³ per gram at 0 °C and 1 atm) is converted to liquid nitrogen with molar volume 34.67 cm³/mol, using 22414 cm³ STP per mol."),
    B("General practice: the Gurvich volume depends on the chosen relative pressure (0.95 is conventional); compare it with the BET or density-functional values before you report a pore volume."),
    B("Checked on synthetic data with a known answer (tests/test_v2192_items.py: 22.414 cm³ STP/g gives 0.03467 cm³/g)."),
]
DOCS["Kinetics & adsorption"] = _doc(
    "Kinetics and adsorption", "Rate constants, activation parameters, turnover frequencies and isotherms.",
    H("What it does"),
    P("Fits concentration-versus-time data and adsorption isotherms."),
    *terms("Observed rate constant (k_obs) and half-life (t½).",
           "Activation energy (E_a), activation enthalpy (ΔH‡) and entropy (ΔS‡).",
           "Turnover frequency (TOF): conversions per open site per hour.",
           "Confidence interval (CI): the range that holds the true value with 95 % probability.",
           "Langmuir and Freundlich: two adsorption isotherm models."),
    H("Tutorial"),
    N("Upload a table: sample, replicate, time, concentration; or raw NMR integrals with an internal standard."),
    N("Read the first-order fits (k_obs ± SD, t½, R²), the C/C₀ plots and the k_obs bars. Blank subtraction is available."),
    N("With several temperatures: Arrhenius and Eyring analysis with confidence intervals; initial rates; TOF per open site (open sites from the defect calculator); pD–rate profile (pKa)."),
    N("Adsorption table (sample, Ce, qe): Langmuir and Freundlich fits."),
    F(r"$C = C_0\,e^{-k t}\qquad t_{1/2} = \dfrac{\ln 2}{k}$"),
    F(r"$k = A\,e^{-E_a/RT}\qquad k = \dfrac{k_B T}{h}\,e^{\Delta S^\ddagger/R}\,e^{-\Delta H^\ddagger/RT}$"),
    F(r"$q = \dfrac{q_{max} K C_e}{1 + K C_e}\ \ (\mathrm{Langmuir})\qquad q = K_F C_e^{1/n}\ \ (\mathrm{Freundlich})$"),
    H("Avrami, isotherm model comparison and heating-rate activation energy (v2.18.12)"),
    P("Avrami (Johnson and Mehl 1939; Avrami 1939, 1940): on the page as an expander under the crystallisation heading. The result is n with its standard error (1 SE), k, the half time t½ and R². Check the plot: the Avrami form assumes random nucleation and growth."),
    P("Model comparison (under the Langmuir figure): Langmuir, Sips (Sips 1948) and Toth (Toth 1971) are ranked by the corrected Akaike criterion AICc (Hurvich and Tsai 1989); the lowest AICc is the best model. Differences below 2 are not a real preference."),
    P("Activation energy from heating rates (under the Arrhenius heading): Kissinger (1957) uses the peak temperatures Tp at several heating rates β; Ozawa (1965) and Flynn and Wall (1966) use the temperatures at one fixed conversion. Both need temperatures in kelvin; a value below 273.15 gives a warning."),
    B("General practice: a single Kissinger or Ozawa–Flynn–Wall energy assumes one simple reaction. For several steps, use the energy per conversion, not one number."),
    B("Depends on: the example tables are synthetic; replace them with your data. The energies are apparent values, not proof of a mechanism."),
    H("Avrami and isotherm model checks (library, v2.18.11)"),
    P("Avrami (JMAK) crystallisation kinetics: X(t) = 1 − exp(−k tⁿ) (Avrami 1939, 1940), fitted in the linear form ln(−ln(1 − X)) = ln k + n ln t between 5 % and 95 % transformed. The fit returns k, n, the half time and R²."),
    P("Isotherms: Langmuir, Sips and Toth (Toth 1971) are fitted and ranked by the corrected Akaike criterion AICc (Hurvich and Tsai 1989); the lowest AICc is the best model. Differences below 2 are not a real preference."),
    B("Checked on synthetic data with known parameters (tests/test_v21811_known_answers.py). General practice: n of Avrami is a guide to the mechanism, not proof of it."),
    B("Checked on synthetic data with known answers (tests/test_v21811_known_answers.py, tests/test_v21812_known_answers.py)."),
    H("Additional methods and file checks (v2.19.1)"),
    P("Kissinger–Akahira–Sunose (KAS, Akahira and Sunose 1971): an isoconversional activation energy at one conversion from ln(β / T²) against 1 / T. It uses no fixed constant, so it is a second check next to Ozawa–Flynn–Wall (OFW) at the same conversion. Both need the temperatures in kelvin."),
    P("Vegard's law (Vegard 1921): the composition x of a solid solution from its lattice parameter, a = x a_A + (1 − x) a_B. With a bowing parameter b the bowed form a = x a_A + (1 − x) a_B − b x (1 − x) is used (Denton and Ashcroft 1991). A value outside the two end members is reported, not guessed."),
    P("Sauerbrey relation (Sauerbrey 1959): the mass per area of a thin, rigid film on a quartz crystal microbalance from the frequency shift, Δm/A = −Δf √(μq ρq) / (2 f0²). For 5 MHz quartz this is about 17.7 ng per hertz and square centimetre."),
    P("Birch–Murnaghan third-order equation of state (Birch 1947; Murnaghan 1944): V0, B0 (bulk modulus) and B0′ (its pressure derivative) fitted to pressure–volume data, with standard errors."),
    P("Randles–Sevcik (Randles 1948; Ševčík 1948): the cyclic-voltammetry peak current ip = 0.4463 n F A C √(n F v D / RT) of a reversible, diffusion-controlled couple; the inverse gives the diffusion coefficient D. Specific capacity (mAh/g) = I t / m and coulombic efficiency = discharge / charge capacity are also in the library."),
    P("File reading: a heat-flow column in mW is converted to W, a wavelength column in µm, Å or m to nm, and a tab-separated export with decimal commas (1,5) is read as 1.5. Detection names a JCAMP-DX file before parsing; JCAMP-DX is detected but not read yet."),
    B("General practice: KAS and OFW assume one kinetic model per conversion; Vegard is an approximation (check the bowing parameter); Sauerbrey holds only for thin, rigid films, and a warning is given above 2 % of the crystal frequency; Randles–Sevcik assumes a reversible, planar, diffusion-controlled process."),
    B("Depends on: the constants used (F, R, quartz modulus and density) and the temperature (298.15 K for Randles–Sevcik unless you give another). Check the file's units column before you trust a converted value."),
    B("Checked on synthetic data with known answers (tests/test_v2191_known_answers.py)."),
    H("Friedman, heating-rate checks, binding and enzyme curves (v2.19.2)"),
    P("Friedman (1964) differential isoconversional method: at one conversion α the slope of ln(β · dα/dT) against 1 / T gives the activation energy Ea, with no model assumed. The curve is taken as α against temperature for each heating rate β (kelvin per second). It is sensitive to noise, so smooth the curve first."),
    P("Heating-rate check (before Kissinger or KAS): the rates should span at least a factor of 2 and neighbouring rates should differ by more than 10 %; the peak temperature should rise with the rate. A warning is shown otherwise."),
    P("Kissinger–Friedman cross-check: the two activation energies must agree within 15 % (relative difference |a − b| / mean). A larger difference points to wrong peak picks or too much smoothing."),
    P("Hill plot (Hill 1910): the fractional occupancy θ = xⁿ / (Kⁿ + xⁿ) gives K (the concentration at half occupancy) and the Hill coefficient n (n above 1 means cooperative binding). Michaelis–Menten (1913) and Lineweaver–Burk (1934): v = Vmax S / (Km + S), fitted by non-linear least squares; the double-reciprocal line is shown for comparison."),
    B("General practice: a Hill coefficient from fewer than about 6 concentrations, or from data that do not reach both plateaus, is not reliable; check the residuals."),
    B("Checked on synthetic data with known answers (tests/test_v2192_items.py: K = 10, n = 2.5; Vmax = 2, Km = 5; Friedman Ea = 120 kJ/mol)."),

)

DOCS["SDS-PAGE (Image Lab)"] = _doc(
    "SDS-PAGE – protein cleavage from gel bands", "Band percentages versus time.",
    H("What it does"),
    P("Reads the band table of a gel and plots how much of each band (intact protein and fragments) there is over time."),
    *terms("Sodium dodecyl sulfate polyacrylamide gel electrophoresis (SDS-PAGE): proteins separated by size in a gel.",
           "Image Lab: the Bio-Rad software that quantifies the bands (Gel Doc EZ)."),
    H("Tutorial"),
    N("Upload the Image Lab 'Lane and Bands' table (CSV or XLSX)."),
    N("Map each lane to a time point."),
    N("Read the band percentage versus time for each band."),
)

DOCS["Zeta / DLS"] = _doc(
    "Zeta potential and dynamic light scattering", "Size, polydispersity and surface charge from Malvern Zetasizer data.",
    H("What it does"),
    P("Reads Zetasizer data (CSV, Excel or the .dts file), and plots zeta potentials, size distributions and the extra quantities of the instrument report."),
    *terms("Dynamic light scattering (DLS): particle size from the flicker of scattered light.",
           "Z-average: the intensity-weighted mean hydrodynamic diameter (nm), from a cumulant fit of the correlation function.",
           "Polydispersity index (PDI): the width of the size distribution (a small number means uniform).",
           "Zeta potential (ζ, in mV): the electrical potential at the particle's slipping plane, from the electrophoretic mobility.",
           "Zeta deviation: the width (standard deviation) of the zeta distribution of one record.",
           "Conductivity (mS/cm) and count rate (kcps, thousands of counts per second)."),
    H("Tutorial"),
    N("Upload the export or the .dts file. For .dts files of software 8.x every record is read; size records give names, the intensity distribution, Z-average and PDI as stored by the software; zeta records give ζ, deviation, conductivity, count rate, the three zeta peaks (mean, area %, SD) and the zeta distribution."),
    N("Untick runs (records) you want to leave out of the averages and graphs."),
    N("Read the table (also saved for the SI), choose what to show in the figures (Z-average, PDI, ζ, ± SD), and download the figures."),
    N("The Methods text is generated from the data; fill in the bracketed items yourself."),
    B("Checked on real files against the Zetasizer reports: Z-average, PDI, ζ, deviation, conductivity and the peaks are equal to the reported values."),
    B("The zeta distribution is shown only for records whose stored block reproduces the software's ζ."),
    B("Software 7.11 files: Z-average and PDI are calculated from the stored correlation function; ζ needs the Excel/CSV export."),
    F(r"$d_H = \dfrac{k_B T}{3\pi\,\eta\,D}$"),
    P("The Stokes–Einstein relation behind the Z-average: d_H is the hydrodynamic diameter, η the viscosity and D the diffusion coefficient."),
    H("Mobility to zeta potential (library, v2.18.11)"),
    P("Converts an electrophoretic mobility (µm·cm/(V·s)) into ζ with Henry's function (Henry 1931) and the approximation of Ohshima (1994, J. Colloid Interface Sci. 168, 269). The Debye length follows from the ionic strength of a 1:1 electrolyte; the parameter ka (particle radius over Debye length) decides the regime."),
    B("Smoluchowski form (ζ = η μ / ε): valid for ka above about 300, where Henry's function is within about 1 % of 1.5. Smoluchowski gives 12.8 mV per 1 µm·cm/(V·s) in water at 25 °C."),
    B("Hückel form (f = 1) for ka below about 1 (small particles, thick double layer)."),
    B("Checked: the Debye length of a 1 mM 1:1 electrolyte in water at 25 °C is 9.6 nm (tests in tests/test_v21811_known_answers.py)."),
    B("v2.18.12: on the Zeta page, the expander '⚡ Zeta potential from the electrophoretic mobility' gives ζ from the mobility, ionic strength, radius and temperature. Depends on: the radius (for ka), the ionic strength (for the Debye length) and the temperature (for viscosity and permittivity of water)."),
    B("General practice: report the ionic strength and the medium with the ζ value; ζ is not the surface potential."),

)

DOCS["Defect calculator"] = _doc(
    "Defect calculator – missing linkers in MOFs", "Every method the data allows, side by side, with a consensus and the experimental-ratio method of Lázaro 2020.",
    H("What it does"),
    P("Takes the results from TGA, NMR, ICP, N₂ sorption and titration (or typed values) and calculates the number of linkers per metal node by every possible method. "
      "The experimental-ratio method (i) additionally handles multivariate MOFs, partly decomposing components and composites."),
    *terms("Metal–organic framework (MOF) and its metal node (for example M₆O₄(OH)₄ with M = Zr or Hf).",
           "Linkers per node (L), missing linkers (x = ideal − L), modulator or cap per node (m): a modulator is a monocarboxylic acid (formate, acetate, benzoate) that sits where a linker is missing.",
           "Open sites (OS): metal positions free for catalysis.",
           "Thermogravimetric analysis (TGA): the mass of a sample while it is heated in air; the residue at the end is the metal oxide (for example ZrO₂).",
           "Dehydroxylated (DH) MOF: the framework just before its linker burns off, after the solvent, the modulators and the water of the node (M₆O₄(OH)₄ → M₆O₆) have left.",
           "Nuclear magnetic resonance (NMR): the molar ratios of linker and modulators in the digested MOF (¹H NMR; high-performance liquid chromatography (HPLC) or ICP can also give them).",
           "Inductively coupled plasma (ICP) and energy-dispersive X-ray spectroscopy (EDX): metal ratios and chloride.",
           "Building unit (SBU, secondary building unit): the metal–oxo cluster of the node, here counted per metal atom.",
           "Experimental ratio (R_exp): the mass of the structure divided by the mass of its residue, read from the TGA curve.",
           "Multivariate (MTV) MOF: a MOF with several linkers; MTVM: with several modulators.",
           "Composite: a MOF with a second material (an oxide, a polymer, a loaded molecule).",
           "Consensus: the median of the quantitative methods, with the minimum, maximum and standard deviation."),
    H("Tutorial"),
    N("Fill the data, or pick the sample in the sample store: TGA (W at 200 °C, plateau, residue), NMR integrals, ICP and nitrogen results of that sample are pulled in automatically; you can still edit them. Choose the framework and node type."),
    N("In 'Species' list the main linker, the modulators (role cap) and any guest. Give the NMR integrals. Mark a species 'at plateau' if it survives to the dehydroxylated plateau. "
      "For a component that only partly burns away, set its 'residue fraction' (for example 0.5)."),
    N("Read the table of methods, the consensus and the minimal formula. Open each method for its equations, your own numbers step by step, pitfalls and the DOI."),
    N("Method (i) runs by itself when the TGA and NMR data are there. Only for a composite, or to choose the species that balances the charge, open '5b · Composite and charge balance'."),
    N("Download the PDF report, the CSV table or the JSON input."),
    H("The methods (a)–(h)"),
    B("a: TGA plateau (ideal-normalised). b: M-200 (mass at 200 °C + NMR ratios + Cl with charge balance) and b′ the missing-cluster model. c: TGA plateau corrected with the NMR modulator ratio. "
      "d: experimental minimal formula (oxide residue, quantitative NMR, Cl, water by difference). e: ICP metal : linker ratio and mass balance. f: porosity versus the ideal pore volume (semi-quantitative). "
      "g: potentiometric titration. h: open-site rule."),
    F(r"$\mathrm{OS} = (12 - CN) + t\,x - \sum_i z_i\,n_i - n_{Cl}$"),
    P("CN is the ideal connectivity, t the number of carboxylates per linker, z_i the charge and n_i the amount of cap i. For 12-connected nodes this reduces to OS = 2x − m − Cl."),
    H("Method (i): the experimental-ratio method (Lázaro 2020, new in 2.18.4)"),
    P("Idea in plain words. Instead of reading step heights from the TGA curve and deciding which step belongs to which species, the method uses ONE number: how many times the mass of the structure "
      "is larger than the mass of its metal-oxide residue. Combined with the molar ratios from NMR and the charge balance of the node, that single ratio fixes the composition exactly. "
      "The oxygen that the linker's carboxylate gives to the residue is counted correctly, which older plateau methods miss."),
    F(r"$R_{exp} = \dfrac{M_w(\mathrm{structure})}{M_w(\mathrm{residue})}$"),
    F(r"$x = \dfrac{R\,M_{res} - M_{SBU} - \sum n_j M_j - (q - \sum n_j z_j)\,M_{bal}/z_{bal}}{\sum_i r_i M_i - \left(\sum_i r_i z_i\right) M_{bal}/z_{bal} - R\sum_i r_i f_i M_i}$"),
    P("All amounts are per metal atom and are multiplied by the metals per node at the end. x is the main linker per metal, M_res the residue per metal (ZrO₂), M_SBU the building unit per metal, "
      "r_i the molar amount of component i relative to the main linker (from NMR), z_i its charge, f_i the fraction of it that stays in the residue (0 = burns away), n_j a fixed counter-ion measured by EDX or ICP "
      "(chloride), q the charge available per metal, and M_bal, z_bal the mass and charge of the one species that closes the charge balance."),
    P("Two stages are solved, each with its own equation:"),
    B("(i1) the desolvated structure at about 200 °C: building unit still hydroxylated, modulators present. The balancing species is hydroxide by default; you can choose oxide, chloride or none "
      "(then the charge that is left per node is reported)."),
    B("(i2) the dehydroxylated structure just before the linker burns off: modulators gone (except species marked 'at plateau'), oxide (O²⁻) compensates the missing linkers: ZrO(L)ₓO₁₋ₓ. "
      "Option 'no compensation': ZrO(L)ₓ."),
    H("Composites and special cases (section 5b of the page)"),
    B("Inert second material (an oxide or metal that does not decompose, for example TiO₂@MOF or ZrO₂@MOF): give the share of the residue that comes from the MOF (from ICP or EDX); "
      "the second material is removed from both the structure mass and the residue."),
    B("A second material that burns away completely before the linker (polymer, loaded molecule): the DH stage (i2) gives the MOF composition; the start ratio then gives the weight percent of the second material."),
    F(r"$w_{comp} = 1 - \dfrac{M_{MOF}}{R_{exp}\,M_{res}}$"),
    B("A second material that burns away together with the linker: give its weight percent (from ICP, UV-Vis or HPLC); the MOF mass is the measured mass times (1 − w)."),
    B("Several linkers (multivariate MOF) or several modulators: list them all with their NMR integrals. A linker that only partly decomposes gets a residue fraction (for example 0.5)."),
    B("Other nodes: every node of the node library works (for example Zn₄O of MOF-5, Ti₈ of MIL-125, M₃O of MIL-100, paddlewheel of HKUST-1, ZIF-8); the building unit and the charge come from the node."),
    H("What is checked"),
    B("The implementation reproduces the Supporting Information of the paper: Tables S9 and S10 (pristine and defective UiO-66 with formate, both with and without oxygen compensation), S14/S15 (multivariate MOF with a partly decomposing linker), S23 (inert TiO₂ composite) and S24/S25 (polymer that burns before the linker). Tests on those values are part of the program."),
    B("It also reproduces the error the paper describes for the older plateau method: 1.57 linkers instead of 2 for the most defective structure without oxygen compensation."),
    B("Every result shows the check R_calc = R_exp (the solution must give back the measured ratio) and the charge left per node."),
    B("A warning appears when the balancing species comes out negative (the NMR table or the stage does not fit) and when the two stages disagree by more than 0.15 linkers."),
    B("Method (i) is shown next to the methods (a)–(h) and is not part of the consensus, so earlier results stay the same."),
    H("What it needs, and what it does not do"),
    B("Needs: the TGA residue and the mass at 200 °C (i1) and/or at the plateau (i2), measured in air; the NMR ratios if there are modulators or more linkers."),
    B("The nature of the balancing species (hydroxide, chloride, oxide) must come from another technique (ICP, EDX, FT-IR). The program only reports what is left."),
    B("Not implemented: MOFs grown on MOFs (sections S8 of the paper), a second material that burns away after the DH MOF or that decomposes only partly (S11), and amorphous materials."),
    H("Other checks of the calculator"),
    B("The code reproduces published values (Valenzano 2011, Shearer 2014/2016, De Roo 2025) and recovers the compositions used to build the synthetic examples."),
    B("Methods f and b′ are semi-quantitative and are excluded from the consensus. The Klet titration is only checked for its conventions, not against raw data."),
    B("Extra-framework ZrO₂ in composites is removed from the TGA masses before the methods a–d and (i)."),
    B(GP + "measure the TGA in air at a slow heating rate, weigh the sample after the same drying as the NMR sample, and check the plateau at the minimum of the derivative curve."),
)
DOCS["Figure composer"] = _doc(
    "Combined figure – several panels in one figure", "Panels from different techniques with journal templates.",
    H("What it does"),
    P("Puts graphs of different techniques into one multi-panel figure with panel letters."),
    *terms("Panel: one graph of the combined figure.", "Journal template: preset column width, font size and label style (ACS, RSC, Wiley)."),
    H("Tutorial"),
    N("Under any graph click 'Add to combined figure'."),
    N("Open 'Combined figure' and arrange the panels (1–4 columns). Choose the panel letters (A/B/C, a/b/c or a)) and their place."),
    N("Choose a journal template; the program checks the width and the minimum text size."),
    N("Export as vector graphics."),
    B("Panels use the default style of each figure; the composer has its own controls."),
)

DOCS["Reticon"] = _doc(
    "Reticon – schemes with lab and materials icons", "A drawing editor for graphical abstracts and method schemes.",
    H("What it does"),
    P("A scheme editor with about 800 icons (CC0, bundled) for lab equipment, materials, instruments and symbols. Everything is offline."),
    *terms("Scalable vector graphics (SVG), portable document format (PDF): vector export formats.",
           "Portable network graphics (PNG), tagged image file format (TIFF): picture export formats with a chosen resolution (dpi)."),
    H("Tutorial"),
    N("The editor opens with the module. Press the prominent 'Full screen' button (or Ctrl+Shift+M) to work in a large window."),
    N("Drag icons from the library onto the canvas. Widen the library by dragging its edge (double-click resets it)."),
    N("Add arrows (connect icons), text and boxes; align, group and distribute items; undo with Ctrl+Z."),
    N("Start from a template or from the content tab (methods schematic, your data, molecule, structure, combined figure)."),
    N("Export to SVG, PDF, PNG or TIFF at the width and dpi of your journal. Schemes are saved with the project file."),
)

DOCS["Batch processing"] = _doc(
    "Batch processing – one preset for many files", "A whole folder to a set of figures.",
    H("What it does"),
    P("Applies one set of figure settings to every file of a folder."),
    *terms("Preset: a saved set of figure settings (JSON) or the house style."),
    H("Tutorial"),
    N("Choose the technique (PXRD, FTIR, Raman, XPS, TGA, DSC, N₂ sorption or SEM/TEM particles) and a folder (sub-folders included), an upload or the example folder."),
    N("Choose the preset, the journal width and the layout (one figure per file, all overlaid or stacked)."),
    N("Run. Download the ZIP (PNG, TIFF, PDF, SVG and CSV) or write to an output folder. A log lists every file."),
)

DOCS["SI export bundle"] = _doc(
    "SI export bundle – everything for the supporting information", "Figures, data, parameters and tables in one ZIP.",
    H("What it does"),
    P("Collects every figure shown in the session into one ZIP."),
    *terms("Supporting information (SI): the extra material that goes with a paper."),
    H("Tutorial"),
    N("Keep 'Collect figures for the SI bundle' ticked in the sidebar while you work."),
    N("Open the SI export page and choose what to include."),
    N("Download the ZIP: figures (PDF, PNG, SVG), data CSVs, SI_parameters.xlsx (figure list, analysis and figure settings, results per sample) and the defect SI table (DOCX/XLSX)."),
    N("'Dossier per material' makes one ZIP per material (figures, data, methods, results). The cross-checks between techniques (for example crystallite size versus particle size) are listed here too."),
    B("Every number in the sample summary carries its page, data files, time and program version."),
)

DOCS["Sample database"] = _doc(
    "Sample database – your samples and their results", "A local list of sample IDs linked to protocols, results and data files.",
    H("What it does"),
    P("A local database (SQLite, in the Reticula data folder) of your samples."),
    *terms("Sample ID: the code of a sample; data files are linked when it appears in their names (<ID>-B<batch>_…).",
           "Digital object identifier (DOI): the permanent link of a source paper."),
    H("Tutorial"),
    N("Load your Excel sample list (a 'Sample ID' column; optional Name, Paper/Project, Protocol family, Parent / starting material, Protocol basis, Main source (DOI), protocol files)."),
    N("Point it at your data folders. Browse the samples with their protocol files, family, parent and children, source DOI, results from the sample store and linked data files; add notes."),
    B("A synthetic example list is bundled."),
)

DOCS["Export guide"] = _doc(
    "Export guide – which file to export from each instrument", "Open the 'Export guide' module for the full table.",
    H("What it does"),
    P("Lists for every instrument which file to export, which formats Reticula reads, and the menu path in the instrument software. It includes the JEOL JEM-1010 and JEM-2200FS electron microscopes and the Zetasizer."),
    *terms("Comma-separated values (CSV): a plain text table; the safest export format.",
           "Transmission electron microscope (TEM) and its high-resolution mode (HRTEM).",
           "Dynamic light scattering (DLS): the Zetasizer size measurement."),
    H("Tip"),
    B("If an import fails, export plain text as described there; two-column text always works."),
)
