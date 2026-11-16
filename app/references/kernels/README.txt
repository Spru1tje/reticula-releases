DFT-N2-77K-carbon-slit.csv
  N2 at 77 K, slit-shaped carbon pores, 0.4-9.2 nm. Kernel isotherms in mmol per g per cm3/g of pore volume
  versus relative pressure. Taken unchanged from pyGAPS 4.6.1 (pygaps/data/kernels), MIT licence
  (LICENSE_pyGAPS_MIT.txt; Iacomi & Llewellyn, Adsorption 2019, doi:10.1007/s10450-019-00168-5).
  pyGAPS does not document how this kernel was generated; Reticula labels it "DFT kernel (pyGAPS)" and does
  NOT call it NLDFT or QSDFT.

No openly licensed NLDFT/QSDFT kernel for N2 on oxides (cylindrical pores) was found. For oxide / cylindrical
pores Reticula therefore uses documented classical methods instead:
  - Horvath-Kawazoe (slit) and Saito-Foley (cylinder) micropore analysis, with the Cheng-Yang correction,
    ported from pyGAPS (MIT) - adsorbent parameters "Carbon (HK)" or "AlSi oxide ion" (Saito & Foley 1991);
  - BJH (Barrett-Joyner-Halenda, cylinder, Kelvin + Harkins-Jura thickness) for mesopores.
Commercial NLDFT/QSDFT kernels (Micromeritics, Anton Paar/Quantachrome) are proprietary; their PSDs can be
imported from the instrument export instead.
