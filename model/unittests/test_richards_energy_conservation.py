import sys
from pathlib import Path
from unittest.mock import MagicMock
import numpy as np
import pytest
sys.path.append('../')
from pyRichards import richards
from pyRichards import mssubsurface


def _load_module_under_test():
    # mssubsurface imports mpi4py at module import time. Provide a lightweight mock
    # when mpi4py is unavailable in the active test environment.
    if "mpi4py" not in sys.modules:
        mpi4py_mock = MagicMock()
        mpi4py_mock.MPI = MagicMock()
        sys.modules["mpi4py"] = mpi4py_mock
        sys.modules["mpi4py.MPI"] = mpi4py_mock.MPI


MODULE = _load_module_under_test()

RHO_W = 1000.0
C_W = 4186.0


def _integrated_heat_power_j_s(hdiv_heat_j_m2_s, area_m2):
    # Convert heat divergence [J/m2/s] to integrated power [J/s].
    # Sign convention follows the model output; net over all HRUs should close to ~0.
    power = np.sum(hdiv_heat_j_m2_s * area_m2 )
    print("[_integrated_heat_power_j_s] hdiv_heat=", hdiv_heat_j_m2_s, flush=True)
    print("[_integrated_heat_power_j_s] area=", area_m2, flush=True)
    print(f"[_integrated_heat_power_j_s] net_power_J_s={power:+.6e}", flush=True)
    return power


def _configure_richards_for_heat_test(area):
    obj = richards.richards(nhru=2, nsoil=1, vsp_flag=True)

    obj.dem1 = np.array([11.0, 10.0], dtype=float)
    obj.w = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)
    obj.dx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)
    obj.area = np.asarray(area, dtype=float)
    obj.af = 1.0
    obj.flag_sat = True

    # Saturated state so lateral flow and advective heat transport are active.
    obj.theta[:] = np.array([[0.45], [0.45]], dtype=float)
    obj.thetar[:] = np.array([[0.05], [0.05]], dtype=float)
    obj.thetas[:] = np.array([[0.45], [0.45]], dtype=float)
    obj.b[:] = np.array([[4.0], [4.0]], dtype=float)
    obj.satpsi[:] = np.array([[-0.3], [-0.3]], dtype=float)
    obj.ksat[:] = np.array([[1.0e-5], [1.0e-5]], dtype=float)
    obj.m[:] = np.array([2.0, 2.0], dtype=float)
    obj.dz[:] = np.array([[1.0], [1.0]], dtype=float)

    temperature = np.array([[300.0], [280.0]], dtype=float)

    print("[_configure_richards_for_heat_test] area=", obj.area, flush=True)
    print("[_configure_richards_for_heat_test] dem1=", obj.dem1, flush=True)
    print("[_configure_richards_for_heat_test] theta=", obj.theta[:, 0], flush=True)
    print("[_configure_richards_for_heat_test] temperature=", temperature[:, 0], flush=True)

    return obj, temperature


def _configure_richards_hbands_for_heat_test(area):
    obj = richards.richards_hbands(nhru=2, nhband=2, nsoil=1, vsp_flag=True)

    obj.ncsbasins = 1
    obj.dem1hband = np.array([11.0, 10.0], dtype=float)
    obj.w = {"Basin1": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)}
    obj.dx = {"Basin1": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)}
    obj.area = np.asarray(area, dtype=float)
    obj.af = 1.0
    obj.flag_sat = True
    # Saturated state so lateral flow and advective heat transport are active.
    obj.theta[:] = np.array([[0.45], [0.45]], dtype=float)
    obj.thetar[:] = np.array([[0.05], [0.05]], dtype=float)
    obj.thetas[:] = np.array([[0.45], [0.45]], dtype=float)
    obj.b[:] = np.array([[4.0], [4.0]], dtype=float)
    obj.satpsi[:] = np.array([[-0.3], [-0.3]], dtype=float)
    obj.ksat[:] = np.array([[1.0e-5], [1.0e-5]], dtype=float)
    obj.m[:] = np.array([2.0, 2.0], dtype=float)
    obj.dz[:] = np.array([[1.0], [1.0]], dtype=float)

    temperature = np.array([[300.0], [280.0]], dtype=float)

    print("[_configure_richards_hbands_for_heat_test] area=", obj.area, flush=True)
    print("[_configure_richards_hbands_for_heat_test] dem1hband=", obj.dem1hband, flush=True)
    print("[_configure_richards_hbands_for_heat_test] theta=", obj.theta[:, 0], flush=True)
    print("[_configure_richards_hbands_for_heat_test] temperature=", temperature[:, 0], flush=True)

    return obj, temperature


def test_calculate_advective_heat_divergence_from_q_conserves_energy_directly():
    print("\n[test_calculate_advective_heat_divergence_from_q_conserves_energy_directly] start", flush=True)

    # Two-node antisymmetric lateral flow field [mm/s].
    # q[i, j] is flux contribution from i toward j in model convention.
    q = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=float)
    temperature = np.array([300.0, 280.0], dtype=float)

    print("[direct] q=", q, flush=True)
    print("[direct] temperature=", temperature, flush=True)

    hdiv_heat = richards.calculate_advective_heat_divergence_from_q(q, temperature, RHO_W, C_W)
    print("[direct] hdiv_heat [J/m2/s]=", hdiv_heat, flush=True)

    # For equal volumes, global energy closure is simple sum ≈ 0.
    net = np.sum(hdiv_heat)
    print(f"[direct] net_heat_divergence={net:+.6e}", flush=True)

    assert np.all(np.isfinite(hdiv_heat))
    assert net == pytest.approx(0.0, abs=1e-12)
    print("[test_calculate_advective_heat_divergence_from_q_conserves_energy_directly] PASS", flush=True)


def test_update_numba_richards_conserves_energy():
    print("\n[test_update_numba_richards_conserves_energy] start", flush=True)
    obj, temperature = _configure_richards_for_heat_test(area=[100.0, 200.0])

    print("[richards] calling update_numba(vsp_flag=True, temperature=..., rho_w=..., c_w=...)", flush=True)
    obj.update_numba(vsp_flag=True, temperature=temperature, rho_w=RHO_W, c_w=C_W, hdiv_heat = np.zeros((2,1)))

    print("[richards] hdiv [mm/s]=", obj.hdiv[:, 0], flush=True)
    print("[richards] hdiv_heat [J/m2/s]=", obj.hdiv_heat[:, 0], flush=True)

    net_heat_power = _integrated_heat_power_j_s(
        hdiv_heat_j_m2_s=obj.hdiv_heat[:, 0],
        area_m2=obj.area,
    )

    # Non-trivial transport should exist with temperature gradient and active flow.
    assert np.max(np.abs(obj.hdiv_heat[:, 0])) > 0.0
    # Domain-integrated heat source/sink should close to zero.
    assert net_heat_power == pytest.approx(0.0, abs=1e-8)
    print("[test_update_numba_richards_conserves_energy] PASS", flush=True)


def test_update_numba_richards_hbands_conserves_energy():
    print("\n[test_update_numba_richards_hbands_conserves_energy] start", flush=True)
    obj, temperature = _configure_richards_hbands_for_heat_test(area=[100.0, 200.0])

    print("[richards_hbands] calling update_numba(vsp_flag=True, temperature=..., rho_w=..., c_w=...)", flush=True)
    obj.update_numba(vsp_flag=True, temperature=temperature, rho_w=RHO_W, c_w=C_W, hdiv_heat = np.zeros((2,1)))

    print("[richards_hbands] hdiv [mm/s]=", obj.hdiv[:, 0], flush=True)
    print("[richards_hbands] hdiv_heat [J/m2/s]=", obj.hdiv_heat[:, 0], flush=True)

    net_heat_power = _integrated_heat_power_j_s(
        hdiv_heat_j_m2_s=obj.hdiv_heat[:, 0],
        area_m2=obj.area,
    )

    assert np.max(np.abs(obj.hdiv_heat[:, 0])) > 0.0
    assert net_heat_power == pytest.approx(0.0, abs=1e-8)
    print("[test_update_numba_richards_hbands_conserves_energy] PASS", flush=True)


def test_update_numba_richards_unsaturated_has_near_zero_heat_transport():
    print("\n[test_update_numba_richards_unsaturated_has_near_zero_heat_transport] start", flush=True)
    obj, temperature = _configure_richards_for_heat_test(area=[100.0, 200.0])

    # Force clearly unsaturated conditions below cutoff: theta <= (1-eps)*theta_s.
    obj.theta[:] = np.array([[0.30], [0.30]], dtype=float)
    print("[richards-unsat] adjusted theta=", obj.theta[:, 0], flush=True)
    print("[richards-unsat] saturation cutoff=", (1.0 - 0.01) * obj.thetas[:, 0], flush=True)
    print("[richards-unsat] temperature=", temperature[:, 0], flush=True)

    print("[richards-unsat] calling update_numba(vsp_flag=True, temperature=..., rho_w=..., c_w=...)", flush=True)
    obj.update_numba(vsp_flag=True, temperature=temperature, rho_w=RHO_W, c_w=C_W)

    print("[richards-unsat] hdiv [mm/s]=", obj.hdiv[:, 0], flush=True)
    print("[richards-unsat] hdiv_heat [J/m2/s]=", obj.hdiv_heat[:, 0], flush=True)

    net_heat_power = _integrated_heat_power_j_s(
        hdiv_heat_j_m2_s=obj.hdiv_heat[:, 0],
        area_m2=obj.area,
    )

    max_abs_hdiv_heat = np.max(np.abs(obj.hdiv_heat[:, 0]))
    print(f"[richards-unsat] max_abs_hdiv_heat={max_abs_hdiv_heat:.6e} J/m2/s", flush=True)
    print(f"[richards-unsat] net_heat_power={net_heat_power:+.6e} J/s", flush=True)

    # In unsaturated layers, lateral flow is suppressed, so heat advection should be negligible.
    assert max_abs_hdiv_heat == pytest.approx(0.0, abs=1e-12)
    assert net_heat_power == pytest.approx(0.0, abs=1e-8)
    print("[test_update_numba_richards_unsaturated_has_near_zero_heat_transport] PASS", flush=True)