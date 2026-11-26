"""Zetasizer .dts with software 8.02 that holds ZETA records only (no size-intensity distribution).

The file is private measurement data and is NOT in the repository. Set MFS_TEST_DTS to its path to run this test;
otherwise it is skipped. This test checks the structure of the reader output only; the measured values are not in the repository.

Known and unknown (checked against the bytes of the file, values not recorded here):
  * Zeta, zeta deviation and count rate are read from the zeta records; zeta deviation of one record and count rate of
    two records could not be located as bytes by a direct search. The origin of those values is unknown.
  * The file contains no intensity size distribution (no 70-point size axis 0.4-10 000 nm in any record), so `dist` is None.
    The zeta distribution is read separately (zdist) and is checked below.
"""
import os

import numpy as np
import pytest


def test_zeta_only_8_02_file_reads_zeta_and_no_size_distribution():
    """Zeta-only .dts (software 8.02): three records with names and zeta values; dist is None, zdist is present."""
    path = os.environ.get("MFS_TEST_DTS")
    if not path or not os.path.exists(path):
        pytest.skip("no zeta-only 8.02 .dts file (MFS_TEST_DTS); private measurement data is not part of the repository")
    from mfs import readers as R

    with open(path, "rb") as fh:
        z = R.read_zetasizer(fh.read(), "x.dts")
    sm = z["summary"]
    assert len(sm) == 3
    assert all(isinstance(s, str) and s.strip() for s in sm["sample"])
    assert np.isfinite(sm["zeta"].astype(float)).all()
    assert np.isfinite(sm["zeta_sd"].astype(float)).all()
    assert z["dist"] is None                      # no intensity size distribution in this file (see module docstring)
    assert z.get("zdist") is not None             # the zeta distribution is present and read
