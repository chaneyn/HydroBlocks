import sys
from pathlib import Path
from unittest.mock import MagicMock
import numpy as np
import pytest
sys.path.append('../')
from pyRichards import richards


def _net_mass_m3s_from_hdiv_mm_s(hdiv_mm_s, area_m2):
    # Convert [mm/s] divergence back to volumetric [m3/s].
    print("[_net_mass_m3s_from_hdiv_mm_s] hdiv_mm_s=", hdiv_mm_s, flush=True)
    print("[_net_mass_m3s_from_hdiv_mm_s] area_m2=", area_m2, flush=True)
    return np.sum(hdiv_mm_s * area_m2 / 1000.0)


def _configure_richards_case(area):
    print("[_configure_richards_case] area=", np.asarray(area, dtype=float), flush=True)
    obj = richards.richards(nhru=2, nsoil=1, vsp_flag=True)

    obj.dem1 = np.array([11.0, 10.0], dtype=float)
    obj.w = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)
    obj.dx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)
    obj.area = np.asarray(area, dtype=float)
    obj.af = 1.0
    obj.flag_sat = True

    # Keep both units saturated to isolate numerical conservation behavior.
    obj.theta[:] = np.array([[0.45], [0.45]], dtype=float)
    obj.thetar[:] = np.array([[0.05], [0.05]], dtype=float)
    obj.thetas[:] = np.array([[0.45], [0.45]], dtype=float)
    obj.b[:] = np.array([[4.0], [4.0]], dtype=float)
    obj.satpsi[:] = np.array([[-0.3], [-0.3]], dtype=float)
    obj.ksat[:] = np.array([[1.0e-5], [1.0e-5]], dtype=float)
    obj.m[:] = np.array([2.0, 2.0], dtype=float)
    obj.dz[:] = np.array([[1.0], [1.0]], dtype=float)

    print("[_configure_richards_case] theta=", obj.theta[:, 0], flush=True)
    print("[_configure_richards_case] dem1=", obj.dem1, flush=True)
    print("[_configure_richards_case] w=", obj.w, flush=True)
    print("[_configure_richards_case] dx=", obj.dx, flush=True)

    return obj


def _configure_richards_hbands_case(area):
    print("[_configure_richards_hbands_case] area=", np.asarray(area, dtype=float), flush=True)
    obj = richards.richards_hbands(nhru=2, nhband=2, nsoil=1, vsp_flag=True)

    obj.ncsbasins = 1
    obj.dem1hband = np.array([11.0, 10.0], dtype=float)
    obj.w = {"Basin1": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)}
    obj.dx = {"Basin1": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)}
    obj.area = np.asarray(area, dtype=float)
    obj.af = 1.0
    obj.flag_sat = True

    # Keep both units saturated to isolate numerical conservation behavior.
    obj.theta[:] = np.array([[0.45], [0.45]], dtype=float)
    obj.thetar[:] = np.array([[0.05], [0.05]], dtype=float)
    obj.thetas[:] = np.array([[0.45], [0.45]], dtype=float)
    obj.b[:] = np.array([[4.0], [4.0]], dtype=float)
    obj.satpsi[:] = np.array([[-0.3], [-0.3]], dtype=float)
    obj.ksat[:] = np.array([[1.0e-5], [1.0e-5]], dtype=float)
    obj.m[:] = np.array([2.0, 2.0], dtype=float)
    obj.dz[:] = np.array([[1.0], [1.0]], dtype=float)

    print("[_configure_richards_hbands_case] theta=", obj.theta[:, 0], flush=True)
    print("[_configure_richards_hbands_case] dem1hband=", obj.dem1hband, flush=True)
    print("[_configure_richards_hbands_case] w=", obj.w["Basin1"], flush=True)
    print("[_configure_richards_hbands_case] dx=", obj.dx["Basin1"], flush=True)

    return obj


def test_update_numba_richards_conserves_mass_equal_area():
    print("\n[test_update_numba_richards_conserves_mass_equal_area] start", flush=True)
    obj = _configure_richards_case(area=[100.0, 100.0])

    print("[test_update_numba_richards_conserves_mass_equal_area] calling update_numba(flag=True)", flush=True)
    obj.update_numba(vsp_flag=True)
    print("[test_update_numba_richards_conserves_mass_equal_area] hdiv=", obj.hdiv[:, 0], flush=True)
    net_mass = _net_mass_m3s_from_hdiv_mm_s(obj.hdiv[:, 0], obj.area)
    print("[test_update_numba_richards_conserves_mass_equal_area] net_mass=", net_mass, flush=True)

    assert np.isfinite(net_mass)
    assert net_mass == pytest.approx(0.0, abs=1e-12)
    print("[test_update_numba_richards_conserves_mass_equal_area] PASS", flush=True)


def test_update_numba_richards_conserves_mass_unequal_area():
    print("\n[test_update_numba_richards_conserves_mass_unequal_area] start", flush=True)
    obj = _configure_richards_case(area=[100.0, 200.0])

    print("[test_update_numba_richards_conserves_mass_unequal_area] calling update_numba(flag=True)", flush=True)
    obj.update_numba(vsp_flag=True)
    print("[test_update_numba_richards_conserves_mass_unequal_area] hdiv=", obj.hdiv[:, 0], flush=True)
    net_mass = _net_mass_m3s_from_hdiv_mm_s(obj.hdiv[:, 0], obj.area)
    print("[test_update_numba_richards_conserves_mass_unequal_area] net_mass=", net_mass, flush=True)

    # Physical expectation: global net should remain zero regardless of area heterogeneity.
    assert np.isfinite(net_mass)
    print("[test_update_numba_richards_conserves_mass_unequal_area] asserting mass closure against zero", flush=True)
    assert net_mass == pytest.approx(0.0, abs=1e-12)
    print("[test_update_numba_richards_conserves_mass_unequal_area] PASS", flush=True)


def test_update_numba_richards_hbands_conserves_mass_equal_area():
    print("\n[test_update_numba_richards_hbands_conserves_mass_equal_area] start", flush=True)
    obj = _configure_richards_hbands_case(area=[100.0, 100.0])

    print("[test_update_numba_richards_hbands_conserves_mass_equal_area] calling update_numba(flag=True)", flush=True)
    obj.update_numba(vsp_flag=True)
    print("[test_update_numba_richards_hbands_conserves_mass_equal_area] hdiv=", obj.hdiv[:, 0], flush=True)
    net_mass = _net_mass_m3s_from_hdiv_mm_s(obj.hdiv[:, 0], obj.area)
    print("[test_update_numba_richards_hbands_conserves_mass_equal_area] net_mass=", net_mass, flush=True)

    assert np.isfinite(net_mass)
    assert net_mass == pytest.approx(0.0, abs=1e-12)
    print("[test_update_numba_richards_hbands_conserves_mass_equal_area] PASS", flush=True)


def test_update_numba_richards_hbands_conserves_mass_unequal_area():
    print("\n[test_update_numba_richards_hbands_conserves_mass_unequal_area] start", flush=True)
    obj = _configure_richards_hbands_case(area=[100.0, 200.0])

    print("[test_update_numba_richards_hbands_conserves_mass_unequal_area] calling update_numba(flag=True)", flush=True)
    obj.update_numba(vsp_flag=True)
    print("[test_update_numba_richards_hbands_conserves_mass_unequal_area] hdiv=", obj.hdiv[:, 0], flush=True)
    net_mass = _net_mass_m3s_from_hdiv_mm_s(obj.hdiv[:, 0], obj.area)
    print("[test_update_numba_richards_hbands_conserves_mass_unequal_area] net_mass=", net_mass, flush=True)

    # Physical expectation: global net should remain zero regardless of area heterogeneity.
    assert np.isfinite(net_mass)
    print("[test_update_numba_richards_hbands_conserves_mass_unequal_area] asserting mass closure against zero", flush=True)
    assert net_mass == pytest.approx(0.0, abs=1e-12)
    print("[test_update_numba_richards_hbands_conserves_mass_unequal_area] PASS", flush=True)

def _configure_richards_case_no_vsp(area):
    print("[_configure_richards_case_no_vsp] area=", np.asarray(area, dtype=float), flush=True)
    obj = richards.richards(nhru=2, nsoil=1, vsp_flag=False)

    obj.dem1 = np.array([11.0, 10.0], dtype=float)
    obj.w = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)
    obj.dx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)
    obj.area = np.asarray(area, dtype=float)
    obj.af = 1.0

    obj.theta[:] = np.array([[0.45], [0.45]], dtype=float)
    obj.thetar[:] = np.array([0.05, 0.05], dtype=float)
    obj.thetas[:] = np.array([0.45, 0.45], dtype=float)
    obj.b[:] = np.array([4.0, 4.0], dtype=float)
    obj.satpsi[:] = np.array([-0.3, -0.3], dtype=float)
    obj.ksat[:] = np.array([1.0e-5, 1.0e-5], dtype=float)
    obj.m[:] = np.array([2.0, 2.0], dtype=float)
    obj.dz[:] = np.array([[1.0], [1.0]], dtype=float)

    print("[_configure_richards_case_no_vsp] theta=", obj.theta[:, 0], flush=True)
    print("[_configure_richards_case_no_vsp] dem1=", obj.dem1, flush=True)
    print("[_configure_richards_case_no_vsp] w=", obj.w, flush=True)
    print("[_configure_richards_case_no_vsp] dx=", obj.dx, flush=True)

    return obj


def _configure_richards_hbands_case_no_vsp(area):
    print("[_configure_richards_hbands_case_no_vsp] area=", np.asarray(area, dtype=float), flush=True)
    obj = richards.richards_hbands(nhru=2, nhband=2, nsoil=1, vsp_flag=False)

    obj.ncsbasins = 1
    obj.dem1hband = np.array([11.0, 10.0], dtype=float)
    obj.w = {"Basin1": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)}
    obj.dx = {"Basin1": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)}
    obj.area = np.asarray(area, dtype=float)
    obj.af = 1.0

    obj.theta[:] = np.array([[0.45], [0.45]], dtype=float)
    obj.thetar[:] = np.array([0.05, 0.05], dtype=float)
    obj.thetas[:] = np.array([0.45, 0.45], dtype=float)
    obj.b[:] = np.array([4.0, 4.0], dtype=float)
    obj.satpsi[:] = np.array([-0.3, -0.3], dtype=float)
    obj.ksat[:] = np.array([1.0e-5, 1.0e-5], dtype=float)
    obj.m[:] = np.array([2.0, 2.0], dtype=float)
    obj.dz[:] = np.array([[1.0], [1.0]], dtype=float)

    print("[_configure_richards_hbands_case_no_vsp] theta=", obj.theta[:, 0], flush=True)
    print("[_configure_richards_hbands_case_no_vsp] dem1hband=", obj.dem1hband, flush=True)
    print("[_configure_richards_hbands_case_no_vsp] w=", obj.w["Basin1"], flush=True)
    print("[_configure_richards_hbands_case_no_vsp] dx=", obj.dx["Basin1"], flush=True)

    return obj


def test_update_numba_richards_nonvsp_conserves_mass_equal_area():
    print("\n[test_update_numba_richards_nonvsp_conserves_mass_equal_area] start", flush=True)
    obj = _configure_richards_case_no_vsp(area=[100.0, 100.0])

    print("[test_update_numba_richards_nonvsp_conserves_mass_equal_area] calling update_numba(flag=False)", flush=True)
    obj.update_numba(vsp_flag=False)
    print("[test_update_numba_richards_nonvsp_conserves_mass_equal_area] hdiv=", obj.hdiv[:, 0], flush=True)

    net_mass = _net_mass_m3s_from_hdiv_mm_s(obj.hdiv[:, 0], obj.area)
    print("[test_update_numba_richards_nonvsp_conserves_mass_equal_area] net_mass=", net_mass, flush=True)

    assert np.isfinite(net_mass)
    assert net_mass == pytest.approx(0.0, abs=1e-12)
    print("[test_update_numba_richards_nonvsp_conserves_mass_equal_area] PASS", flush=True)


def test_update_numba_richards_nonvsp_conserves_mass_unequal_area():
    print("\n[test_update_numba_richards_nonvsp_conserves_mass_unequal_area] start", flush=True)
    obj = _configure_richards_case_no_vsp(area=[100.0, 200.0])

    print("[test_update_numba_richards_nonvsp_conserves_mass_unequal_area] calling update_numba(flag=False)", flush=True)
    obj.update_numba(vsp_flag=False)
    print("[test_update_numba_richards_nonvsp_conserves_mass_unequal_area] hdiv=", obj.hdiv[:, 0], flush=True)

    net_mass = _net_mass_m3s_from_hdiv_mm_s(obj.hdiv[:, 0], obj.area)
    print("[test_update_numba_richards_nonvsp_conserves_mass_unequal_area] net_mass=", net_mass, flush=True)

    assert np.isfinite(net_mass)
    print("[test_update_numba_richards_nonvsp_conserves_mass_unequal_area] asserting mass closure against zero", flush=True)
    assert net_mass == pytest.approx(0.0, abs=1e-12)
    print("[test_update_numba_richards_nonvsp_conserves_mass_unequal_area] PASS", flush=True)


def test_update_numba_richards_hbands_nonvsp_conserves_mass_equal_area():
    print("\n[test_update_numba_richards_hbands_nonvsp_conserves_mass_equal_area] start", flush=True)
    obj = _configure_richards_hbands_case_no_vsp(area=[100.0, 100.0])

    print("[test_update_numba_richards_hbands_nonvsp_conserves_mass_equal_area] calling update_numba(flag=False)", flush=True)
    obj.update_numba(vsp_flag=False)
    print("[test_update_numba_richards_hbands_nonvsp_conserves_mass_equal_area] hdiv=", obj.hdiv[:, 0], flush=True)

    net_mass = _net_mass_m3s_from_hdiv_mm_s(obj.hdiv[:, 0], obj.area)
    print("[test_update_numba_richards_hbands_nonvsp_conserves_mass_equal_area] net_mass=", net_mass, flush=True)

    assert np.isfinite(net_mass)
    assert net_mass == pytest.approx(0.0, abs=1e-12)
    print("[test_update_numba_richards_hbands_nonvsp_conserves_mass_equal_area] PASS", flush=True)


def test_update_numba_richards_hbands_nonvsp_conserves_mass_unequal_area():
    print("\n[test_update_numba_richards_hbands_nonvsp_conserves_mass_unequal_area] start", flush=True)
    obj = _configure_richards_hbands_case_no_vsp(area=[100.0, 200.0])

    print("[test_update_numba_richards_hbands_nonvsp_conserves_mass_unequal_area] calling update_numba(flag=False)", flush=True)
    obj.update_numba(vsp_flag=False)
    print("[test_update_numba_richards_hbands_nonvsp_conserves_mass_unequal_area] hdiv=", obj.hdiv[:, 0], flush=True)

    net_mass = _net_mass_m3s_from_hdiv_mm_s(obj.hdiv[:, 0], obj.area)
    print("[test_update_numba_richards_hbands_nonvsp_conserves_mass_unequal_area] net_mass=", net_mass, flush=True)

    assert np.isfinite(net_mass)
    print("[test_update_numba_richards_hbands_nonvsp_conserves_mass_unequal_area] asserting mass closure against zero", flush=True)
    assert net_mass == pytest.approx(0.0, abs=1e-12)
    print("[test_update_numba_richards_hbands_nonvsp_conserves_mass_unequal_area] PASS", flush=True)