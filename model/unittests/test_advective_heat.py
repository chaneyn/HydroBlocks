import sys
import numpy as np
import numpy.testing as npt
sys.path.append('../')
from pyRichards import richards
from pyRichards import mssubsurface

def test_advective_heat_matches_compute_enthalpy_flux():
    # Synthetic domain: 3 HRUs, single soil layer
    n = 3
    area = np.array([1000.0, 1000.0, 1000.0])  # m2
    dz = np.array([0.5, 0.5, 0.5])  # m (layer thickness)
    temps = np.array([280.0, 285.0, 290.0])  # K

    # Define signed volumetric flows [m3/s] (flows[i,j]: from i to j)
    flows = np.array([
        [0.0, 0.001, -0.001],
        [-0.001, 0.0, 0.001],
        [0.001, -0.001, 0.0]
    ])

    # Convert flows to q [mm/s] using the same sign convention inverted by calculate_advective_heat_divergence_from_q
    # q[i,j] = -flows[i,j] * 1000 / area[i]
    q_mm_s = -flows * 1000.0 / area[:, np.newaxis]

    print("\n[test_advective_heat_matches_compute_enthalpy_flux] start", flush=True)
    print("[advective] flows (m3/s)=", flows, flush=True)
    print("[advective] q (mm/s)=", q_mm_s, flush=True)

    rho_w = 1000.0
    c_w = 4186.0

    # Call refactored function (layer vector)
    rhs_ref = richards.calculate_advective_heat_divergence_from_q(q_mm_s, temps, rho_w, c_w)
    print("[advective] rhs_ref (J/s/m2)=", rhs_ref, flush=True)

    # Build inputs for compute_enthalpy_flux (clusters identity -> one cluster per HRU)
    clusters = np.eye(n)
    flows_3d = np.zeros((n, n, 1))
    flows_3d[:, :, 0] = flows
    temps_2d = temps.reshape((n, 1))
    dz_hrus = dz.reshape((n, 1))

    hdiv_heat = mssubsurface.compute_enthalpy_flux(clusters, flows_3d, temps_2d, area, dz_hrus)
    rhs_expected = hdiv_heat[:, 0]
    print("[advective] rhs_expected (J/s/m2)=", rhs_expected, flush=True)

    # compute_enthalpy_flux counts ordered pairs (k,m) and (m,k) separately
    # in the flows tensor; in our q-based reconstruction each link appears
    npt.assert_allclose(rhs_ref, rhs_expected, rtol=1e-8, atol=1e-10)
    print("[test_advective_heat_matches_compute_enthalpy_flux] PASS", flush=True)

    # Energy conservation: net volumetric heat divergence over closed domain ~ 0
    print("\n[test_mass_conversion_and_conservation] start", flush=True)
    npt.assert_allclose(np.sum(rhs_ref), 0.0, atol=1e-16)
    npt.assert_allclose(np.sum(rhs_expected), 0.0, atol=1e-16)
    print("[test_mass_conservation] PASS", flush=True)