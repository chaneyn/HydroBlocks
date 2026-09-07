import sys
import numpy as np
import numpy.testing as npt

sys.path.append('../')

from pyRichards.advectivetransport import AdvectiveTransport


def test_advect_tracer_step_conserves_mass_and_moves_downstream():
    """
    A single-layer, two-node pair should conserve tracer mass and move the
    mass from the donor node to the receiver node for a closed exchange.
    """
    print("\n[test_advect_tracer_step_conserves_mass_and_moves_downstream] start", flush=True)
    adv = AdvectiveTransport(
        nhrus=2,
        nrisfus=2,
        nsoil=1,
        reg_ids={0: np.arange(2, dtype=int)},
        comm=None,
    )

    area = np.array([100.0, 100.0])
    theta = np.array([[0.40], [0.40]], dtype=float)
    dz = np.array([[1.0], [1.0]], dtype=float)
    tracer_mass0 = np.array([[40.0], [120.0]], dtype=float)

    q = np.zeros((2, 2, 1), dtype=float)
    q[0, 1, 0] = 0.20  # 0 -> 1

    print("[tracer_step] tracer_mass0=", tracer_mass0, flush=True)
    print("[tracer_step] q=", q, flush=True)

    tracer_mass1 = adv.advect_tracer_step(tracer_mass0, q, area, theta, dz, dt=1.0)

    Vw_before = area[:, None] * theta * dz
    Vw_after = Vw_before.copy()
    Vw_after[0, 0] -= 0.20
    Vw_after[1, 0] += 0.20

    mass0 = tracer_mass0
    mass1 = tracer_mass1

    print("[tracer_step] tracer_mass1=", tracer_mass1, flush=True)
    print("[tracer_step] mass0=", np.sum(mass0), flush=True)
    print("[tracer_step] mass1=", np.sum(mass1), flush=True)

    npt.assert_allclose(np.sum(mass1), np.sum(mass0), rtol=1e-12, atol=1e-12)
    assert mass1[0, 0] < mass0[0, 0]
    assert mass1[1, 0] > mass0[1, 0]
    print("[test_advect_tracer_step_conserves_mass_and_moves_downstream] PASS", flush=True)


def test_aggregate_and_redistribute_mass_round_trip_preserves_mass():
    """
    Tracer mass should round-trip between RISFU and HRU space without loss.
    """
    print("\n[test_aggregate_and_redistribute_mass_round_trip_preserves_mass] start", flush=True)
    adv = AdvectiveTransport(
        nhrus=3,
        nrisfus=3,
        nsoil=1,
        reg_ids={0: np.arange(3, dtype=int)},
        comm=None,
    )

    farea = np.eye(3, dtype=float)

    tracer_mass_risfu = np.array([[12.0], [24.0], [36.0]], dtype=float)

    print("[redistribute] tracer_mass_risfu=", tracer_mass_risfu, flush=True)
    print("[redistribute] farea=", farea, flush=True)

    tracer_mass_hrus = adv.redistribute_mass_hrus(tracer_mass_risfu, farea)
    tracer_mass_risfu_back = adv.aggregate_mass_risfu(tracer_mass_hrus, farea)

    mass_risfu = np.sum(tracer_mass_risfu)
    mass_hrus = np.sum(tracer_mass_hrus)

    print("[redistribute] tracer_mass_hrus=", tracer_mass_hrus, flush=True)
    print("[redistribute] tracer_mass_risfu_back=", tracer_mass_risfu_back, flush=True)
    print("[redistribute] mass_risfu=", mass_risfu, flush=True)
    print("[redistribute] mass_hrus=", mass_hrus, flush=True)

    npt.assert_allclose(mass_hrus, mass_risfu, rtol=1e-12, atol=1e-12)
    npt.assert_allclose(tracer_mass_risfu_back, tracer_mass_risfu, rtol=1e-12, atol=1e-12)
    print("[test_aggregate_and_redistribute_mass_round_trip_preserves_mass] PASS", flush=True)


def test_compute_reg_tracer_returns_local_rows_only():
    """
    compute_reg_tracer should return the updated tracer mass only for the
    requested local regional-unit indices.
    """
    print("\n[test_compute_reg_tracer_returns_local_rows_only] start", flush=True)
    adv = AdvectiveTransport(
        nhrus=3,
        nrisfus=3,
        nsoil=1,
        reg_ids={0: np.arange(3, dtype=int)},
        comm=None,
    )

    area = np.array([100.0, 100.0, 100.0], dtype=float)
    theta = np.array([[0.45], [0.45], [0.45]], dtype=float)
    dz = np.array([[1.0], [1.0], [1.0]], dtype=float)
    c0 = np.array([[1.0], [2.0], [3.0]], dtype=float)

    q = np.zeros((3, 3, 1), dtype=float)
    q[0, 1, 0] = 0.10
    q[1, 2, 0] = 0.08

    print("[regional] c0=", c0, flush=True)
    print("[regional] q=", q, flush=True)

    c_expected = adv.advect_tracer_step(c0, q, area, theta, dz, dt=1.0)
    c_local = adv.compute_reg_tracer(c0, q, np.array([1, 2], dtype=int), area, theta, dz, dt=1.0)

    print("[regional] c_expected=", c_expected, flush=True)
    print("[regional] c_local=", c_local, flush=True)

    npt.assert_allclose(c_local, c_expected[[1, 2], :], rtol=1e-12, atol=1e-12)
    print("[test_compute_reg_tracer_returns_local_rows_only] PASS", flush=True)
