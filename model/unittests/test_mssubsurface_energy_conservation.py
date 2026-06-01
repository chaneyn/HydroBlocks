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
    power = np.sum(hdiv_heat_j_m2_s * area_m2)
    print("[_integrated_heat_power_j_s] hdiv_heat=", hdiv_heat_j_m2_s, flush=True)
    print("[_integrated_heat_power_j_s] area=", area_m2, flush=True)
    print(f"[_integrated_heat_power_j_s] net_power_J_s={power:+.6e}", flush=True)
    return power


def _reconstruct_net_heat_cluster(flows_m3_s, temperatures_k):
    """Rebuild the per-cluster heat power balance implied by the signed flow tensor."""
    nclusters, _, nsoil = flows_m3_s.shape
    net_heat_cluster = np.zeros((nclusters, nsoil), dtype=float)

    for layer in range(nsoil):
        for source in range(nclusters):
            for target in range(source+1, nclusters):
                q_link = flows_m3_s[source, target, layer]
                if np.abs(q_link) == 0.0:
                    continue

                if q_link > 0.0:
                    upwind = source
                    downwind = target
                else:
                    upwind = target
                    downwind = source

                heat_power = np.abs(q_link) * RHO_W * C_W * temperatures_k[upwind, layer]
                net_heat_cluster[upwind, layer] += heat_power
                net_heat_cluster[downwind, layer] -= heat_power

    print("[_reconstruct_net_heat_cluster] net_heat_cluster=", net_heat_cluster, flush=True)
    return net_heat_cluster


def test_compute_intermediate_hdiv_heat_conserves_energy_inside_same_cid():
    print("\n[test_compute_intermediate_hdiv_heat_conserves_energy_inside_same_cid] start", flush=True)

    cls = mssubsurface.mssubsurface
    ss = cls.__new__(cls)

    # 2 intermediate units mapped 1:1 to 2 HRUs in the same CID.
    ss.farea_gw = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
        ],
        dtype=float,
    )

    # Antisymmetric volumetric exchange [m3/s] for one layer.
    # Positive link from unit 0 to unit 1 and equal reverse sign.
    ss.inter_unit_flow_m3s = np.zeros((2, 2, 1), dtype=float)
    ss.inter_unit_flow_m3s[0, 1, 0] = 2.0e-6
    ss.inter_unit_flow_m3s[1, 0, 0] = -2.0e-6

    ss.temp_gw = np.array(
        [
            [300.0],
            [280.0],
        ],
        dtype=float,
    )
    area_hrus = np.array([100.0, 200.0], dtype=float)
    dz_hrus = np.array([[1.0], [1.0]], dtype=float)

    print("[intermediate] farea_gw=", ss.farea_gw, flush=True)
    print("[intermediate] flows layer0=", ss.inter_unit_flow_m3s[:, :, 0], flush=True)
    print("[intermediate] temperatures=", ss.temp_gw[:, 0], flush=True)

    hdiv_heat = ss.compute_intermediate_hdiv_heat(
        area_hrus=area_hrus,
        dz_hrus=dz_hrus,
    )

    print("[intermediate] hdiv_heat layer0 [J/s/m2]=", hdiv_heat[:, 0], flush=True)

    net_power = _integrated_heat_power_j_s(hdiv_heat[:, 0], area_hrus)

    # Active transport should be present with nonzero flow and temperature gradient.
    assert np.max(np.abs(hdiv_heat[:, 0])) > 0.0
    # Energy conservation inside same CID: no net source/sink.
    assert net_power == pytest.approx(0.0, abs=1e-16)

    print("[test_compute_intermediate_hdiv_heat_conserves_energy_inside_same_cid] PASS", flush=True)


def test_compute_regional_hdiv_heat_conserves_energy_across_multiple_cids():
    print("\n[test_compute_regional_hdiv_heat_conserves_energy_across_multiple_cids] start", flush=True)

    cls = mssubsurface.mssubsurface

    # Shared regional exchange tensor across 2 regional units (each representing one CID).
    flows = np.zeros((2, 2, 1), dtype=float)
    flows[0, 1, 0] = 3.0e-6
    flows[1, 0, 0] = -3.0e-6

    reg_temperature = np.array(
        [
            [302.0],  # CID 1 unit
            [285.0],  # CID 2 unit
        ],
        dtype=float,
    )

    # Local object for CID 1
    ss_cid1 = cls.__new__(cls)
    ss_cid1.regional_inter_unit_flow_m3s_cross = flows
    ss_cid1.reg_temperature_gw = reg_temperature
    ss_cid1.this_cid = np.array([0], dtype=np.int64)
    ss_cid1.farea_gw = np.array([[1.0]], dtype=float)

    # Local object for CID 2
    ss_cid2 = cls.__new__(cls)
    ss_cid2.regional_inter_unit_flow_m3s_cross = flows
    ss_cid2.reg_temperature_gw = reg_temperature
    ss_cid2.this_cid = np.array([1], dtype=np.int64)
    ss_cid2.farea_gw = np.array([[1.0]], dtype=float)

    area_cid1 = np.array([120.0], dtype=float)
    dz_cid1 = np.array([[1.0]], dtype=float)
    area_cid2 = np.array([180.0], dtype=float)
    dz_cid2 = np.array([[1.0]], dtype=float)

    print("[regional] flows layer0=", flows[:, :, 0], flush=True)
    print("[regional] reg_temperature=", reg_temperature[:, 0], flush=True)

    hdiv_heat_cid1 = ss_cid1.compute_regional_hdiv_heat(area_hrus=area_cid1, dz_hrus=dz_cid1)
    hdiv_heat_cid2 = ss_cid2.compute_regional_hdiv_heat(area_hrus=area_cid2, dz_hrus=dz_cid2)

    print("[regional] hdiv_heat_cid1 [J/s/m2]=", hdiv_heat_cid1[:, 0], flush=True)
    print("[regional] hdiv_heat_cid2 [J/s/m2]=", hdiv_heat_cid2[:, 0], flush=True)

    power_cid1 = _integrated_heat_power_j_s(hdiv_heat_cid1[:, 0], area_cid1)
    power_cid2 = _integrated_heat_power_j_s(hdiv_heat_cid2[:, 0], area_cid2)
    net_power_global = power_cid1 + power_cid2

    print(f"[regional] power_cid1={power_cid1:+.6e} J/s", flush=True)
    print(f"[regional] power_cid2={power_cid2:+.6e} J/s", flush=True)
    print(f"[regional] net_power_global={net_power_global:+.6e} J/s", flush=True)

    # Each CID should receive/shed heat, but total across CIDs should close.
    assert np.max(np.abs(hdiv_heat_cid1[:, 0])) > 0.0
    assert np.max(np.abs(hdiv_heat_cid2[:, 0])) > 0.0
    assert net_power_global == pytest.approx(0.0, abs=1e-8)

    print("[test_compute_regional_hdiv_heat_conserves_energy_across_multiple_cids] PASS", flush=True)


def test_compute_enthalpy_flux_conserves_energy_through_cluster_to_hru_mapping():
    print("\n[test_compute_enthalpy_flux_conserves_energy_through_cluster_to_hru_mapping] start", flush=True)

    clusters = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.75, 0.25, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )

    temperatures = np.array(
        [
            [300.0, 299.0],
            [290.0, 288.0],
            [280.0, 281.0],
        ],
        dtype=float,
    )

    flows = np.zeros((3, 3, 2), dtype=float)

    # Layer 0: 0 <-> 1 and 1 <-> 2 exchanges, antisymmetric by construction.
    flows[0, 1, 0] = 2.0e-6
    flows[1, 0, 0] = -2.0e-6
    flows[1, 2, 0] = 1.5e-6
    flows[2, 1, 0] = -1.5e-6

    # Layer 1: different pattern to ensure the conservation check is layer-wise.
    flows[0, 2, 1] = 1.0e-6
    flows[2, 0, 1] = -1.0e-6
    flows[1, 0, 1] = 0.5e-6
    flows[0, 1, 1] = -0.5e-6

    area_hrus = np.array([120.0, 180.0, 150.0, 200.0], dtype=float)
    dz_hrus = np.array(
        [
            [1.0, 0.8],
            [1.0, 0.8],
            [1.0, 0.8],
            [1.0, 0.8],
        ],
        dtype=float,
    )

    print("[entropy] clusters=", clusters, flush=True)
    print("[entropy] temperatures=", temperatures, flush=True)
    print("[entropy] flows layer0=", flows[:, :, 0], flush=True)
    print("[entropy] flows layer1=", flows[:, :, 1], flush=True)
    print("[entropy] area_hrus=", area_hrus, flush=True)
    print("[entropy] dz_hrus=", dz_hrus, flush=True)

    expected_net_heat_cluster = _reconstruct_net_heat_cluster(flows, temperatures)

    for layer in range(flows.shape[2]):
        layer_net = np.sum(expected_net_heat_cluster[:, layer])
        print(f"[entropy] layer={layer} expected_net_heat_cluster_sum={layer_net:+.6e} J/s", flush=True)
        assert layer_net == pytest.approx(0.0, abs=1e-12)

    hdiv_heat = mssubsurface.compute_enthalpy_flux(clusters, flows, temperatures, area_hrus, dz_hrus)
    print("[entropy] hdiv_heat=", hdiv_heat, flush=True)

    for layer in range(flows.shape[2]):
        integrated_power = _integrated_heat_power_j_s(hdiv_heat[:, layer], area_hrus)
        print(f"[entropy] layer={layer} integrated_power={integrated_power:+.6e} J/s", flush=True)
        assert integrated_power == pytest.approx(0.0, abs=1e-12)

    total_power = np.sum(hdiv_heat * area_hrus[:, np.newaxis] * dz_hrus)
    print(f"[entropy] total_power={total_power:+.6e} J/s", flush=True)
    assert total_power == pytest.approx(0.0, abs=1e-12)

    # Conservation at the cluster level must survive the HRU projection.
    assert np.allclose(expected_net_heat_cluster.sum(axis=0), 0.0, atol=1e-12)
    assert np.all(np.isfinite(hdiv_heat))

    print("[test_compute_enthalpy_flux_conserves_energy_through_cluster_to_hru_mapping] PASS", flush=True)
