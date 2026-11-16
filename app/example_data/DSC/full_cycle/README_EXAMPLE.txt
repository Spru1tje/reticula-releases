Full heating/cooling cycle EXAMPLE built from the user's REAL Mettler STARe heating curve (SG_DSC_1.1, 1 K/min,
25 -> 70 °C) + a SIMULATED 5 min isothermal at 70 °C + a SIMULATED cooling 70 -> 25 °C at 1 K/min with a
crystallisation exotherm (onset 40.8 °C, ΔHc = +98 J/g). The cooling part is synthetic and labelled as EXAMPLE.

Writes example_data/DSC/full_cycle/:
  SG_DSC_full-cycle_Mettler_EXAMPLE.txt       one curve block, continuous time (case a: one file, all segments)
  SG_DSC_heating_Mettler_REAL.txt             the real heating export (case b, part 1)
  SG_DSC_cooling_Mettler_EXAMPLE.txt          the simulated cooling alone, t restarting at 0, later time stamp (case b, part 2)
  SG_DSC_3-curves_Mettler_EXAMPLE.txt         three 'Curve Name:' blocks heat / iso / cool, t restarting (case c)
Truth of the simulated crystallisation: onset 40.8 °C, peak 39.6 °C (programme T), ΔHc = +98.0 J/g.
