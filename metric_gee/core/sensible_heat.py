"""
Module 9: Sezuvchan issiqlik oqimi H — to'liq iterativ Monin-Obukhov.

TO'LIQ ITERATSIYA MANTIQI (har iteratsiyada):
  1. rah yangilash  (stability correction)
  2. Anchor nuqtada rah, rho sampling  → dt_hot, dt_cold yangilash
  3. a, b qayta hisoblash
  4. dT xaritasi yangilash              ← asosiy yangilik
  5. Konvergensiya tekshiruvi

Formulalar (Allen 2007):
  F.28: H = ρ·Cp·dT / rah
  F.29: dT = a·TsDEM + b
  F.38: u* = k·u200 / (ln(200/zom) - ψm200)
  F.39: rah = (ln(z2/z1) - ψh2 + ψh01) / (u*·k)
  F.40: L = -ρ·Cp·u*³·Ts / (k·g·H)
  F.41-45: ψm, ψh (Paulson 1970 unstable | Webb 1970 stable)
"""
import ee
import math
from metric_gee.config.constants import (
    VON_KARMAN, Z_BLEND, Z1, Z2, CP, GRAVITY,
)

MIN_ITER = 3    # kamida shu marta iteratsiya
MAX_ITER = 15   # maksimum


class SensibleHeatFlux:
    """H to'liq iterativ hisoblash — anchor qayta kalibrlash bilan."""

    # ================================================================
    # ASOSIY METOD
    # ================================================================

    def calc_h_iterative(self, anchor_data, ts_dem, ts,
                         aero, rn, g, geometry, verbose=True):
        """
        To'liq iterativ H hisoblash.

        Har iteratsiyada rah yangilanadi → anchor da sample →
        dt_hot, dt_cold → a, b → dT xarita → H → L → ψ → rah → ...

        Args:
            anchor_data : dict (_run_cimec dan)
                'h_cold', 'h_hot'         — anchor H [W/m²]  (ee.Number)
                'cold_mask', 'hot_mask'   — piksel masklari  (ee.Image)
                'cold_ts', 'hot_ts'       — anchor Ts [K]    (ee.Number)
                'cold_rah', 'hot_rah'     — boshlang'ich rah (ee.Number)
                'cold_rho', 'hot_rho'     — boshlang'ich rho (ee.Number)
            ts_dem   : ee.Image — lapse corrected Ts (CIMEC uchun)
            ts       : ee.Image — actual Ts (L uchun)
            aero     : dict — u200, zom, u_star, rah, rho_air
            rn, g    : ee.Image — energiya balansi
            geometry : ee.Geometry — anchor sampling uchun

        Returns:
            H_final  : ee.Image  [W/m²]
        """
        # ── Anchor qiymatlari ──────────────────────────────────────
        h_cold    = anchor_data['h_cold']
        h_hot     = anchor_data['h_hot']
        # Mask Landsat UTM 300m da bo'lishi kerak (_run_cimec da reproject qilingan)
        cold_mask = anchor_data['cold_mask']
        hot_mask  = anchor_data['hot_mask']
        cold_ts   = anchor_data['cold_ts']
        hot_ts    = anchor_data['hot_ts']

        # Ts farqi — nolga bo'lishdan saqlash
        ts_diff     = hot_ts.subtract(cold_ts)
        ts_diff_safe = ee.Number(
            ee.Algorithms.If(ts_diff.abs().lt(0.5), ee.Number(0.5), ts_diff)
        )

        # ── Boshlang'ich aerodynamika ──────────────────────────────
        u200    = aero['u200']
        zom     = aero['zom']
        u_star  = aero['u_star']
        rah     = aero['rah']
        rho_air = aero['rho_air']

        # ── Boshlang'ich a, b (neutral rah bilan) ─────────────────
        cold_rah = anchor_data['cold_rah']
        hot_rah  = anchor_data['hot_rah']
        cold_rho = anchor_data['cold_rho']
        hot_rho  = anchor_data['hot_rho']

        dt_cold_ee = h_cold.multiply(cold_rah).divide(cold_rho.multiply(CP))
        dt_hot_ee  = h_hot.multiply(hot_rah).divide(hot_rho.multiply(CP))

        a = dt_hot_ee.subtract(dt_cold_ee).divide(ts_diff_safe)
        b = dt_hot_ee.subtract(a.multiply(hot_ts))

        # Boshlang'ich dT xaritasi
        dT = ts_dem.multiply(a).add(b).rename('dT')

        # ── Iteratsiya ─────────────────────────────────────────────
        prev_dt_hot  = None
        prev_rah_hot = None

        for i in range(MAX_ITER):

            # F.28: H = ρ·Cp·dT / rah
            H_safe = rho_air.multiply(CP).multiply(dT).divide(rah)
            H_safe = H_safe.where(H_safe.abs().lt(1.0), ee.Image.constant(1.0))

            # F.40: L = -ρ·Cp·u*³·Ts / (k·g·H)
            L = (rho_air.multiply(CP)
                 .multiply(u_star.pow(3))
                 .multiply(ts)
                 .divide(ee.Image.constant(VON_KARMAN * GRAVITY).multiply(H_safe))
                 .multiply(-1))

            # L singularity: |L| ≥ 2 m
            is_stable = L.gt(0)
            L = L.abs().max(2.0).where(is_stable,
                L.abs().max(2.0)).where(is_stable.Not(),
                L.abs().max(2.0).multiply(-1))
            L = L.clamp(-2000, 2000)

            # F.41-45: Stability corrections (Paulson 1970 / Webb 1970)
            psi_m200, psi_h2, psi_h01 = self._calc_stability(L)

            # F.38: u* yangilash
            denom_m = (ee.Image.constant(Z_BLEND).divide(zom)
                       .log().subtract(psi_m200).max(0.5))
            u_star = u200.multiply(VON_KARMAN).divide(denom_m).max(0.01)

            # F.39: rah yangilash — fizik chegaralar [5, 500] s/m
            rah_new = (ee.Image.constant(math.log(Z2 / Z1))
                       .subtract(psi_h2).add(psi_h01)
                       .divide(u_star.multiply(VON_KARMAN)))
            rah_new = rah_new.where(rah_new.lt(0), ee.Image.constant(5.0))
            rah     = rah_new.clamp(5.0, 500.0).rename('rah')

            # ── Anchor rah sampling → a, b yangilash ──────────────
            def _masked_mean(img, mask):
                """Maskadagi o'rtacha — null-safe. Mask Landsat UTM 300m da."""
                for scale in [300, 1000]:
                    res = img.updateMask(mask).reduceRegion(
                        ee.Reducer.mean(), geometry, scale, maxPixels=1e9
                    ).getInfo()
                    if res:
                        val = next((v for v in res.values() if v is not None), None)
                        if val is not None:
                            return ee.Number(val)
                # Fallback: maskasiz butun qism
                res = img.reduceRegion(
                    ee.Reducer.mean(), geometry, 300, maxPixels=1e9
                ).getInfo()
                val = next((v for v in res.values() if v is not None), 0.0) if res else 0.0
                return ee.Number(val)

            cold_rah_new = _masked_mean(rah,     cold_mask)
            hot_rah_new  = _masked_mean(rah,     hot_mask)
            cold_rho_new = _masked_mean(rho_air, cold_mask)
            hot_rho_new  = _masked_mean(rho_air, hot_mask)

            # dt_hot, dt_cold yangilash (yangi rah bilan)
            dt_cold_ee = h_cold.multiply(cold_rah_new).divide(cold_rho_new.multiply(CP))
            dt_hot_ee  = h_hot.multiply(hot_rah_new).divide(hot_rho_new.multiply(CP))

            # a, b yangilash
            a = dt_hot_ee.subtract(dt_cold_ee).divide(ts_diff_safe)
            b = dt_hot_ee.subtract(a.multiply(hot_ts))

            # dT xaritasi yangilash ← ASOSIY YANGILIK
            dT = ts_dem.multiply(a).add(b).rename('dT')

            # dT fizik chegarasi (±20% anchor oralig'idan)
            dt_min  = dt_cold_ee.min(dt_hot_ee)
            dt_max  = dt_cold_ee.max(dt_hot_ee)
            margin  = dt_max.subtract(dt_min).multiply(0.2)
            dT = dT.clamp(dt_min.subtract(margin), dt_max.add(margin))

            # ── Konvergensiya tekshiruvi (bitta server call) ───────
            conv = ee.Dictionary({
                'dt_hot':  dt_hot_ee,
                'rah_hot': hot_rah_new,
                'a_val':   a,
            }).getInfo()

            dt_hot_v  = float(conv.get('dt_hot',  0.0))
            rah_hot_v = float(conv.get('rah_hot', 0.0))
            a_val     = float(conv.get('a_val',   0.0))

            if verbose:
                b_val = b.getInfo()
                print(f"    Iter {i+1:2d}: "
                      f"dT_hot={dt_hot_v:.4f} K  "
                      f"rah_hot={rah_hot_v:.3f} s/m  "
                      f"a={a_val:.5f}  b={b_val:.3f}")

            if prev_dt_hot is not None and (i + 1) >= MIN_ITER:
                if (abs(dt_hot_v  - prev_dt_hot)  < 0.01 and
                        abs(rah_hot_v - prev_rah_hot) < 0.1):
                    if verbose:
                        print(f"    ✅ Converged at iter {i + 1}")
                    break

            prev_dt_hot  = dt_hot_v
            prev_rah_hot = rah_hot_v

        # ── Oxirgi H + fizik limit ─────────────────────────────────
        H_final = rho_air.multiply(CP).multiply(dT).divide(rah).rename('H')

        rn_g = rn.subtract(g)
        H_final = H_final.min(rn_g)                  # H ≤ Rn-G  (LE ≥ 0)
        # ✅ KERAK (Allen 2007: "ET can exceed Rn in arid locations"):
        H_final = H_final.max(rn_g.multiply(-1.0))   # LE ≤ 2.0*(Rn-G)

        return H_final

    # ================================================================
    # STABILITY CORRECTIONS
    # ================================================================

    def _calc_stability(self, L):
        """
        Monin-Obukhov barqarorlik tuzatmalari.

        Unstable (L < 0): Paulson (1970) — F.41-43
        Stable   (L > 0): Webb (1970)    — F.44-45

        Returns:
            psi_m200 : ee.Image — u* uchun (200m)
            psi_h2   : ee.Image — rah uchun (2m)
            psi_h01  : ee.Image — rah uchun (0.1m)
        """
        # x parametrlari (unstable, L < 0).
        # Stable (L > 0) piksellar uchun L_neg → katta manfiy → x≈1 (neytral).
        # Bu psi_unst ni masked qilmaydi; .where() keyinroq psi_st bilan almashadi.
        L_neg = L.where(L.gte(0), ee.Image.constant(-1e6))
        x200 = ee.Image.constant(1).subtract(
            ee.Image.constant(16 * Z_BLEND).divide(L_neg)).pow(0.25)
        x2   = ee.Image.constant(1).subtract(
            ee.Image.constant(16 * Z2).divide(L_neg)).pow(0.25)
        x01  = ee.Image.constant(1).subtract(
            ee.Image.constant(16 * Z1).divide(L_neg)).pow(0.25)

        # F.41: ψm(200m) unstable
        psi_m_unst = (x200.add(1).divide(2).log().multiply(2)
                      .add(x200.pow(2).add(1).divide(2).log())
                      .subtract(x200.atan().multiply(2))
                      .add(math.pi / 2))

        # F.42a: ψh(2m) unstable
        psi_h2_unst = x2.pow(2).add(1).divide(2).log().multiply(2)

        # F.42b: ψh(0.1m) unstable
        psi_h01_unst = x01.pow(2).add(1).divide(2).log().multiply(2)

        # F.44: ψm(200m) stable — -5·z/L, z=200m  ← to'g'ri!
        psi_m_st   = ee.Image.constant(-5 * Z2).divide(L)
        # F.45: ψh(2m) stable
        psi_h2_st  = ee.Image.constant(-5 * Z2).divide(L)
        # F.45: ψh(0.1m) stable
        psi_h01_st = ee.Image.constant(-5 * Z1).divide(L)

        is_unstable = L.lt(0)

        psi_m200 = psi_m_unst.where(is_unstable.Not(),   psi_m_st)
        psi_h2   = psi_h2_unst.where(is_unstable.Not(),  psi_h2_st)
        psi_h01  = psi_h01_unst.where(is_unstable.Not(), psi_h01_st)

        # Neytral (|L| < 0.01) → 0
        is_neutral = L.abs().lt(0.01)
        psi_m200 = psi_m200.where(is_neutral, 0)
        psi_h2   = psi_h2.where(is_neutral,   0)
        psi_h01  = psi_h01.where(is_neutral,  0)

        # Psi juda katta bo'lmasin
        psi_m200 = psi_m200.clamp(-10, 10)
        psi_h2   = psi_h2.clamp(-10, 10)
        psi_h01  = psi_h01.clamp(-10, 10)

        return psi_m200, psi_h2, psi_h01
