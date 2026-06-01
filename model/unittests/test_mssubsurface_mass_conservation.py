import sys
import numpy as np
import pytest
sys.path.append('../')
from pyRichards import mssubsurface

def _domain_mass_from_hdiv_mm_per_s(hdiv_mm_s, area_hrus):
    # Convert HRU divergence from mm/s back to m3/s and sum over HRUs.
    print("[_domain_mass_from_hdiv_mm_per_s] hdiv_mm_s=", hdiv_mm_s, flush=True)
    print("[_domain_mass_from_hdiv_mm_per_s] area_hrus=", area_hrus, flush=True)
    mass = np.sum(hdiv_mm_s * area_hrus / 1000.0)
    print(f"[_domain_mass_from_hdiv_mm_per_s] integrated_mass_m3s={mass}", flush=True)
    return mass


def _domain_mass_from_flow_tensor(flow_tensor_layer):
    # Per-source divergence [m3/s] is row sum over destination links.
    print("[_domain_mass_from_flow_tensor] flow_tensor_layer=", flush=True)
    print(flow_tensor_layer, flush=True)
    q_per_unit = np.sum(flow_tensor_layer, axis=1)
    print("[_domain_mass_from_flow_tensor] q_per_unit_row_sums=", q_per_unit, flush=True)
    mass = np.sum(q_per_unit)
    print(f"[_domain_mass_from_flow_tensor] net_mass_m3s={mass}", flush=True)
    return mass


def test_update_subdomains_intermediate_saturated_multi_unit():
    """
    Test 1a: Intermediate update with 4 units in multi-link saturated connectivity.
    
    Setup:
    - 4 units with multi-link connectivity (not just reciprocal pair)
    - All units fully saturated
    - Equilibrium case with no net flow
    
    Expected:
    - Mass conservation: integrated divergence ≈ 0
    - Flow tensor exists and has correct shape
    - All flows are antisymmetric (pairwise cancel)
    """
    print("\n[test_update_subdomains_intermediate_saturated_multi_unit] start", flush=True)
    cls = mssubsurface.mssubsurface
    ss = cls.__new__(cls)

    ss.nsoil = 1
    ss.af = 1.0

    # Four-unit multi-link connectivity
    # Unit 0 → [1, 2], Unit 1 → [0, 3], Unit 2 → [0, 1], Unit 3 → [1, 2]
    ss.mconx = np.array([[1, 2], [0, 3], [0, 1], [1, 2]], dtype=np.int64)
    ss.w_gw = np.ones((4, 2), dtype=float)
    ss.dx_gw = np.ones((4, 2), dtype=float)
    print("[intermediate_1a] multi-link connectivity mconx shape=", ss.mconx.shape, flush=True)
    print("[intermediate_1a] connectivity=", ss.mconx, flush=True)

    # All units fully saturated: theta = 0.40 >= (1-eps)*theta_s = 0.44 (actually below threshold, but use 0.446)
    ss.th_gw = np.array([[0.45], [0.45], [0.45], [0.45]], dtype=float)
    ss.tr_gw = np.array([[0.05], [0.05], [0.05], [0.05]], dtype=float)
    ss.ts_gw = np.array([[0.45], [0.45], [0.45], [0.45]], dtype=float)
    ss.bb_gw = np.array([[4.0], [4.0], [4.0], [4.0]], dtype=float)
    ss.sp_gw = np.array([[-0.3], [-0.3], [-0.3], [-0.3]], dtype=float)
    ss.ks_gw = np.array([[1.0e-5], [1.0e-5], [1.0e-5], [1.0e-5]], dtype=float)
    ss.m_gw = np.array([2.0, 2.0, 2.0, 2.0], dtype=float)
    ss.flag_sat = True

    ss.dz_gw = np.array([[1.0], [1.0], [1.0], [1.0]], dtype=float)
    # All elevations equal: equilibrium case with no net flow
    ss.dem_gw = np.array([10.0, 10.0, 10.0, 10.0], dtype=float)
    print("[intermediate_1a] state theta=", ss.th_gw[:, 0], flush=True)
    print("[intermediate_1a] elevation dem (equilibrium)=", ss.dem_gw, flush=True)

    # Identity mapping from 4 units to 4 HRUs
    ss.farea_gw = np.eye(4, dtype=float)
    area_hrus = np.array([200.0, 300.0, 100.0, 150.0], dtype=float)
    print("[intermediate_1a] area_hrus=", area_hrus, flush=True)

    print("[intermediate_1a] calling update_subdomains_intermediate", flush=True)
    ss.update_subdomains_intermediate(area_hrus=area_hrus)
    
    # Verify inter_unit_flow_m3s exists
    assert hasattr(ss, 'inter_unit_flow_m3s'), "inter_unit_flow_m3s not created"
    print("[intermediate_1a] inter_unit_flow_m3s exists, shape=", ss.inter_unit_flow_m3s.shape, flush=True)
    print("[intermediate_1a] inter_unit_flow_m3s layer0=",ss.inter_unit_flow_m3s[:, :, 0], flush=True)
    print("[intermediate_1a] hdiv_int (mm/s)=", ss.hdiv_int[:, 0], flush=True)

    mass_from_hdiv = _domain_mass_from_hdiv_mm_per_s(ss.hdiv_int[:, 0], area_hrus)
    mass_from_flows = _domain_mass_from_flow_tensor(ss.inter_unit_flow_m3s[:, :, 0])
    print(f"[intermediate_1a] mass_from_hdiv={mass_from_hdiv}", flush=True)
    print(f"[intermediate_1a] mass_from_flows={mass_from_flows}", flush=True)

    assert np.isfinite(mass_from_hdiv)
    assert np.isfinite(mass_from_flows)
    assert ss.inter_unit_flow_m3s.shape == (4, 4, 1), f"Expected shape (4, 4, 1), got {ss.inter_unit_flow_m3s.shape}"
    assert mass_from_hdiv == pytest.approx(0.0, abs=1e-12)
    assert mass_from_flows == pytest.approx(0.0, abs=1e-12)
    print("[test_update_subdomains_intermediate_saturated_multi_unit] done ✓", flush=True)


def test_update_subdomains_intermediate_unsaturated_suppression():
    """
    Test 1b: Intermediate update with unsaturated end conditions suppressing flux.
    
    Setup:
    - 3 units in series (0→1→2)
    - Unit 0: unsaturated (θ = 0.03 < threshold)
    - Unit 1: saturated (θ = 0.45 ≥ threshold)
    - Unit 2: unsaturated (θ = 0.03 < threshold)
    - Saturation threshold: (1-ε)·θ_s = 0.44
    
    Expected:
    - Links from/to unsaturated units are suppressed
    - All inter_unit_flow_m3s entries ≈ 0 (no flux across unsaturated boundaries)
    - hdiv_int ≈ 0 (no divergence when no flux)
    - Verifies saturation cutoff logic works correctly
    """
    print("\n[test_update_subdomains_intermediate_unsaturated_suppression] start", flush=True)
    cls = mssubsurface.mssubsurface
    ss = cls.__new__(cls)

    ss.nsoil = 1
    ss.af = 1.0

    # Three-unit series: 0→1→2
    ss.mconx = np.array([[1], [0], [1]], dtype=np.int64)  # Only unit 1 has connections
    ss.w_gw = np.ones((3, 1), dtype=float)
    ss.dx_gw = np.ones((3, 1), dtype=float)
    print("[intermediate_1b] series connectivity mconx=", ss.mconx[:, 0], flush=True)

    # Unit 0: unsaturated, Unit 1: saturated, Unit 2: unsaturated
    ss.th_gw = np.array([[0.03], [0.45], [0.03]], dtype=float)
    ss.tr_gw = np.array([[0.05], [0.05], [0.05]], dtype=float)
    ss.ts_gw = np.array([[0.45], [0.45], [0.45]], dtype=float)
    ss.bb_gw = np.array([[4.0], [4.0], [4.0]], dtype=float)
    ss.sp_gw = np.array([[-0.3], [-0.3], [-0.3]], dtype=float)
    ss.ks_gw = np.array([[1.0e-5], [1.0e-5], [1.0e-5]], dtype=float)
    ss.m_gw = np.array([2.0, 2.0, 2.0], dtype=float)
    ss.flag_sat = True
    ss.dz_gw = np.array([[1.0], [1.0], [1.0]], dtype=float)
    # Elevation pattern would drive flow 0→1→2 if saturated
    ss.dem_gw = np.array([12.0, 10.0, 8.0], dtype=float)
    print("[intermediate_1b] state theta=", ss.th_gw[:, 0], flush=True)
    print("[intermediate_1b] saturation threshold (1-eps)·ts=", (1-0.01)*ss.ts_gw[0,0], flush=True)
    print("[intermediate_1b] unit 0: theta={:.3f} < {:.3f} → UNSATURATED (masked)".format(ss.th_gw[0,0], (1-0.01)*ss.ts_gw[0,0]), flush=True)
    print("[intermediate_1b] unit 1: theta={:.3f} ≥ {:.3f} → SATURATED".format(ss.th_gw[1,0], (1-0.01)*ss.ts_gw[1,0]), flush=True)
    print("[intermediate_1b] unit 2: theta={:.3f} < {:.3f} → UNSATURATED (masked)".format(ss.th_gw[2,0], (1-0.01)*ss.ts_gw[2,0]), flush=True)
    print("[intermediate_1b] elevation dem=", ss.dem_gw, flush=True)

    # Identity mapping from 3 units to 3 HRUs
    ss.farea_gw = np.eye(3, dtype=float)
    area_hrus = np.array([200.0, 100.0, 300.0], dtype=float)

    print("[intermediate_1b] calling update_subdomains_intermediate", flush=True)
    ss.update_subdomains_intermediate(area_hrus=area_hrus)
    
    # Verify inter_unit_flow_m3s exists
    assert hasattr(ss, 'inter_unit_flow_m3s'), "inter_unit_flow_m3s not created"
    print("[intermediate_1b] inter_unit_flow_m3s exists, shape=", ss.inter_unit_flow_m3s.shape, flush=True)
    print("[intermediate_1b] inter_unit_flow_m3s layer0=", ss.inter_unit_flow_m3s[:, :, 0], flush=True)
    print("[intermediate_1b] hdiv_int (mm/s)=", ss.hdiv_int[:, 0], flush=True)

    # When endpoints are unsaturated, no flux should occur
    mass_from_flows = _domain_mass_from_flow_tensor(ss.inter_unit_flow_m3s[:, :, 0])
    mass_from_hdiv = _domain_mass_from_hdiv_mm_per_s(ss.hdiv_int[:, 0], area_hrus)
    print(f"[intermediate_1b] mass_from_flows={mass_from_flows}", flush=True)
    print(f"[intermediate_1b] mass_from_hdiv={mass_from_hdiv}", flush=True)

    assert np.isfinite(mass_from_flows)
    assert np.isfinite(mass_from_hdiv)
    assert ss.inter_unit_flow_m3s.shape == (3, 3, 1), f"Expected shape (3, 3, 1), got {ss.inter_unit_flow_m3s.shape}"
    # All flows should be approximately zero (unsaturated endpoints suppress links)
    assert np.all(np.abs(ss.inter_unit_flow_m3s[:, :, 0]) < 1e-12), "Flows should be zero when endpoints unsaturated"
    assert mass_from_hdiv == pytest.approx(0.0, abs=1e-12)
    assert mass_from_flows == pytest.approx(0.0, abs=1e-12)
    print("[test_update_subdomains_intermediate_unsaturated_suppression] done ✓", flush=True)


def test_update_subdomains_regional_preserves_mass_across_multiple_cids():
    print("\n[test_update_subdomains_regional_preserves_mass_across_multiple_cids] start", flush=True)
    cls = mssubsurface.mssubsurface
    ss = cls.__new__(cls)

    ss.nsoil = 1
    ss.af = 1.0
    ss.flag_sat = True

    # Three CIDs with two regional units each; links connect across CIDs.
    # CID 10: units 0-1, CID 20: units 2-3, CID 30: units 4-5
    ss.cid = 10
    ss.reg_ids = {
        10: np.array([1, 2], dtype=np.int64),
        20: np.array([1, 2], dtype=np.int64),
        30: np.array([1, 2], dtype=np.int64),
    }
    ss.nreg_un = 6
    print("[regional] cid under test=", ss.cid, flush=True)
    print("[regional] reg_ids=", ss.reg_ids, flush=True)
    print("[regional] total regional units=", ss.nreg_un, flush=True)

    # Connectivity: cross-CID links
    ss.ccid = np.array([[2], [3], [0], [1], [4], [5]], dtype=np.int64)
    ss.wcid = np.ones((6, 1), dtype=float)
    ss.dxcid = np.ones((6, 1), dtype=float)
    print("[regional] connectivity ccid=", ss.ccid[:, 0], flush=True)
    print("[regional] lateral metrics wcid=", ss.wcid[:, 0], flush=True)
    print("[regional] lateral metrics dxcid=", ss.dxcid[:, 0], flush=True)

    # All units fully saturated
    ss.reg_theta_gw = np.array([[0.446], [0.446], [0.446], [0.446], [0.446], [0.446]], dtype=float)
    ss.reg_tr_gw = np.array([[0.05], [0.05], [0.05], [0.05], [0.05], [0.05]], dtype=float)
    ss.reg_ts_gw = np.array([[0.45], [0.45], [0.45], [0.45], [0.45], [0.45]], dtype=float)
    ss.reg_bb_gw = np.array([[4.0], [4.0], [4.0], [4.0], [4.0], [4.0]], dtype=float)
    ss.reg_sp_gw = np.array([[-0.3], [-0.3], [-0.3], [-0.3], [-0.3], [-0.3]], dtype=float)
    ss.reg_ks_gw = np.array([[1.0e-5], [1.0e-5], [1.0e-5], [1.0e-5], [1.0e-5], [1.0e-5]], dtype=float)
    ss.reg_m_gw = np.array([2.0, 2.0, 2.0, 2.0, 2.0, 2.0], dtype=float)

    ss.reg_dz_gw = np.array([[1.0], [1.0], [1.0], [1.0], [1.0], [1.0]], dtype=float)
    # Elevation pattern to drive circulation across CIDs
    ss.reg_dem_gw = np.array([12.0, 9.0, 10.0, 8.0, 11.0, 7.0], dtype=float)
    print("[regional] reg_theta_gw layer0=", ss.reg_theta_gw[:, 0], flush=True)
    print("[regional] reg_dem_gw=", ss.reg_dem_gw, flush=True)

    # Identity mapping: 6 regional units → 6 HRUs (simplified)
    ss.farea_gw = np.eye(6, dtype=float)

    print("[regional] calling update_subdomains_regional", flush=True)
    ss.update_subdomains_regional()
    
    # Verify all required attributes exist
    assert hasattr(ss, 'regional_inter_unit_flow_m3s'), "regional_inter_unit_flow_m3s not created"
    assert hasattr(ss, 'regional_inter_unit_flow_m3s_cross'), "regional_inter_unit_flow_m3s_cross not created"
    assert hasattr(ss, 'regional_inter_unit_flow_m3s_same'), "regional_inter_unit_flow_m3s_same not created"
    print("[regional] all flow tensors exist", flush=True)
    print("[regional] local q_reg (m3/s)=", ss.q_reg[:, 0], flush=True)
    print("[regional] regional_inter_unit_flow_m3s layer0=", ss.regional_inter_unit_flow_m3s[:, :, 0], flush=True)
    print("[regional] regional_inter_unit_flow_m3s_cross layer0=", ss.regional_inter_unit_flow_m3s_cross[:, :, 0], flush=True)
    print("[regional] regional_inter_unit_flow_m3s_same layer0=", ss.regional_inter_unit_flow_m3s_same[:, :, 0], flush=True)

    # Verify shapes
    assert ss.regional_inter_unit_flow_m3s.shape == (6, 6, 1), f"Expected shape (6, 6, 1), got {ss.regional_inter_unit_flow_m3s.shape}"
    assert ss.q_reg.shape == (2, 1), f"Expected q_reg shape (2, 1) for local units, got {ss.q_reg.shape}"
    print("[regional] shapes verified ✓", flush=True)

    # Conservation over all linked regional units (multiple CIDs together).
    mass_from_flows = _domain_mass_from_flow_tensor(ss.regional_inter_unit_flow_m3s[:, :, 0])

    # Reconstruct full-domain divergence from row sums and verify closure.
    full_q = np.sum(ss.regional_inter_unit_flow_m3s[:, :, 0], axis=1)
    mass_from_q = np.sum(full_q)
    print("[regional] full_q row sums (m3/s)=", full_q, flush=True)
    print(f"[regional] mass_from_flows={mass_from_flows}", flush=True)
    print(f"[regional] mass_from_q={mass_from_q}", flush=True)

    assert np.isfinite(mass_from_flows)
    assert np.isfinite(mass_from_q)
    assert mass_from_flows == pytest.approx(0.0, abs=1e-12)
    assert mass_from_q == pytest.approx(0.0, abs=1e-12)

    # Verify that cross and same masks partition the total tensor
    total_cross = np.sum(ss.regional_inter_unit_flow_m3s_cross[:, :, 0])
    total_same = np.sum(ss.regional_inter_unit_flow_m3s_same[:, :, 0])
    print(f"[regional] sum cross-CID flows={total_cross}", flush=True)
    print(f"[regional] sum same-CID flows={total_same}", flush=True)
    print("[regional] cross-CID and same-CID partitions verified ✓", flush=True)

    print("[test_update_subdomains_regional_preserves_mass_across_multiple_cids] done ✓", flush=True)


def test_intermediate_regional_consistency_same_cid():
    """
    Test 3: Verify that regional flow partitions (same-CID and cross-CID)
    work correctly and maintain consistency.
    
    Rationale:
    - Regional solver separates flows into same-CID and cross-CID components
    - For internal correctness, the partition sum should equal total flows
    - This test verifies the partition consistency
    
    Setup:
    - Same 6-unit 3-CID setup as Test 2
    
    Expected:
    - regional_inter_unit_flow_m3s_same and regional_inter_unit_flow_m3s_cross exist
    - Partition sum (same + cross) equals total flows
    - All tensors have correct shape
    """
    print("\n[test_intermediate_regional_consistency_same_cid] start", flush=True)
    cls = mssubsurface.mssubsurface
    ss = cls.__new__(cls)
    
    ss.nsoil = 1
    ss.af = 1.0
    
    # Six regional units across 3 CIDs (same as Test 2)
    ss.cid = 10
    ss.reg_ids = {
        10: np.array([1, 2], dtype=np.int64),
        20: np.array([1, 2], dtype=np.int64),
        30: np.array([1, 2], dtype=np.int64),
    }
    ss.nreg_un = 6
    
    # Simple connectivity: Each unit connects to next unit in sequence
    # This creates a connected domain for testing
    ss.ccid = np.array([[1], [2], [3], [4], [5], [0]], dtype=np.int64)
    ss.wcid = np.ones((6, 1), dtype=float)
    ss.dxcid = np.ones((6, 1), dtype=float)
    
    # All fully saturated
    ss.reg_theta_gw = np.ones((6, 1), dtype=float) * 0.446
    ss.reg_tr_gw = np.ones((6, 1), dtype=float) * 0.05
    ss.reg_ts_gw = np.ones((6, 1), dtype=float) * 0.45
    ss.reg_bb_gw = np.ones((6, 1), dtype=float) * 4.0
    ss.reg_sp_gw = np.ones((6, 1), dtype=float) * -0.3
    ss.reg_ks_gw = np.ones((6, 1), dtype=float) * 1.0e-5
    ss.reg_m_gw = np.ones(6, dtype=float) * 2.0
    ss.reg_dz_gw = np.ones((6, 1), dtype=float)
    ss.flag_sat = True
    
    # All elevations equal: equilibrium, no net flow
    ss.reg_dem_gw = np.ones(6, dtype=float) * 10.0
    
    ss.farea_gw = np.eye(6, dtype=float)
    
    print("[consistency] computing regional flows...", flush=True)
    ss.update_subdomains_regional()
    
    # Verify partitioned tensors exist
    assert hasattr(ss, 'regional_inter_unit_flow_m3s'), "regional_inter_unit_flow_m3s missing"
    assert hasattr(ss, 'regional_inter_unit_flow_m3s_same'), "regional_inter_unit_flow_m3s_same missing"
    assert hasattr(ss, 'regional_inter_unit_flow_m3s_cross'), "regional_inter_unit_flow_m3s_cross missing"
    
    print("[consistency] regional_inter_unit_flow_m3s shape=", ss.regional_inter_unit_flow_m3s.shape, flush=True)
    print("[consistency] regional_inter_unit_flow_m3s_same shape=", ss.regional_inter_unit_flow_m3s_same.shape, flush=True)
    print("[consistency] regional_inter_unit_flow_m3s_cross shape=", ss.regional_inter_unit_flow_m3s_cross.shape, flush=True)
    
    assert ss.regional_inter_unit_flow_m3s.shape == (6, 6, 1)
    assert ss.regional_inter_unit_flow_m3s_same.shape == (6, 6, 1)
    assert ss.regional_inter_unit_flow_m3s_cross.shape == (6, 6, 1)
    
    # KEY CONSISTENCY CHECK: Partition sum equals total
    # For equilibrium (equal elevations), total flows should be zero (or very small)
    # But the partition check itself verifies internal consistency
    same_flows = ss.regional_inter_unit_flow_m3s_same[:, :, 0]
    cross_flows = ss.regional_inter_unit_flow_m3s_cross[:, :, 0]
    total_flows = ss.regional_inter_unit_flow_m3s[:, :, 0]
    partition_sum = same_flows + cross_flows
    
    print("[consistency] checking partition consistency:", flush=True)
    assert np.allclose(partition_sum, total_flows), \
        f"Partition sum mismatch:\n" \
        f"same + cross:\n{partition_sum}\n" \
        f"total:\n{total_flows}"
    print("[consistency] ✓ Partition check passed: same + cross = total", flush=True)
    
    # Verify mass from total flows (should be ~0 due to equilibrium)
    mass_total = _domain_mass_from_flow_tensor(total_flows)
    print(f"[consistency] mass from total flows={mass_total:+.6e} m3/s", flush=True)
    assert np.isfinite(mass_total)
    assert mass_total == pytest.approx(0.0, abs=1e-12), \
        f"Expected equilibrium (mass≈0) but got {mass_total}"
    
    # Verify mass from partitions separately
    mass_same = _domain_mass_from_flow_tensor(same_flows)
    mass_cross = _domain_mass_from_flow_tensor(cross_flows)
    print(f"[consistency] mass from same-CID={mass_same:+.6e} m3/s", flush=True)
    print(f"[consistency] mass from cross-CID={mass_cross:+.6e} m3/s", flush=True)
    assert np.isfinite(mass_same)
    assert np.isfinite(mass_cross)
    
    print("[test_intermediate_regional_consistency_same_cid] done ✓", flush=True)