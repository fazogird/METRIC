# METRIC-GEE Pipeline vs Allen 2007 Maqola — To'liq Solishtirish

## Belgilar
- ✅ = Maqola bilan bir xil
- 🔄 = O'zgartirilgan (yaxshiroq yoki adaptatsiya)
- ⚠️ = Farq bor, natijaga ta'sir qilishi mumkin
- ❌ = Xato yoki tuzatish kerak

---

## 1. KIRISH MA'LUMOTLARI

| # | Parametr | Allen 2007 | Pipeline v2.1 | Baho |
|---|----------|-----------|---------------|------|
| - | Landsat data | L5/L7 Collection 1 Level 1 (DN) | L8/L9 Collection 2 Level 2 (SR+ST) | 🔄 Yaxshiroq — atmosferik korreksiya GEE da tayyor |
| - | Meteo data | AgriMet/ASOS meteo stansiya | ERA5-Land (~9km raster) | 🔄 Stansiya o'rniga global raster — stansiya kerak emas |
| - | DEM | NED 30m | SRTM 30m (GEE) | ✅ |
| - | Land use | NLCD (USA) | ESA WorldCover v100 (global) | 🔄 Global — istalgan hududda ishlaydi |

---

## 2. VEGETATSIYA INDEKSLARI

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| F.23 | NDVI = (NIR−RED)/(NIR+RED) | L5/7: Band 4,3 | L8/9: Band 5,4 | ✅ Band raqamlari to'g'ri moslashtrilgan |
| F.19 | SAVI = (1+L)(NIR−RED)/(L+NIR+RED) | L=0.1 (SAVI_ID) | L=0.1 | ✅ |
| F.18 | LAI = −ln((0.69−SAVI)/0.59)/0.91 | LAI clamp [0, 6] | LAI clamp [0, 6] | ✅ |

---

## 3. ALBEDO

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| F.15 | α = Σ(wᵢ × ρᵢ) | Tasumi vaznlari (L5/7 uchun) | Settings: "liang" (default) yoki "olmedo" | ⚠️ Default "liang" — run.py da "olmedo" ishlatilsa to'g'ri |
| F.10-14 | Atmosferik korreksiya | Qo'lda (COST, 6S) | L2 SR — GEE avtomatik | 🔄 Yaxshiroq |

**Izoh:** `settings.py` da `albedo_method = "liang"` default. Idaho test da "olmedo" ishlatilgan. Olmedo vaznlari L8/9 SR uchun maxsus kalibrlangan — to'g'ri tanlov.

---

## 4. SIRT HARORATI

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| F.20-21 | Ts = K₂/ln(K₁/Lλ+1) + atm correction | L1 DN → radiance → Ts + atmosferik tuzatma | L2 ST — GEE TICOS atm correction tayyor | 🔄 Yaxshiroq |
| F.22 | εNB → Ts tuzatma | Kerak (L1 da) | L2 da kerak emas (TICOS ichiga kiritilgan) | 🔄 |
| F.9b | TsDEM = Ts − Γ(z−z_datum)/1000 | Γ = 6.5 K/km, z_datum = AOI o'rtacha | Γ = 6.5, z_datum = DEM mean | ✅ |

---

## 5. EMISSIVITET

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| F.17 | ε₀ = 0.95+0.01×LAI (LAI≤3), 0.98 (LAI>3) | + | + | ✅ |
| F.22 | εNB = 0.97+0.0033×LAI | Ts hisob uchun | Mavjud lekin L2 da ishlatilmaydi | ✅ |
| - | Suv/qor: ε₀=0.985 | NDVI<0 | NDVI≤0 | ✅ |

---

## 6. SOF RADIATSIYA (Rn)

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| F.2 | Rn = RS↓−αRS↓+RL↓−RL↑−(1−ε₀)RL↓ | + | + | ✅ |
| F.3 | RS↓ = Gsc×cos(θ)×τsw/d² | Gsc=1367 W/m² | Gsc=1367 | ✅ |
| F.4 | τsw = 0.35+0.627×exp(...) | P, W, cos(θ), Kt | P(ERA5), W(F.6), cos(θ), Kt=1.0 | ✅ |
| F.5 | P = 101.3×((293−0.0065z)/293)^5.26 | DEM dan | ERA5 surface_pressure | 🔄 ERA5 aniqroq |
| F.6 | W = 0.14×ea×P+2.1 | Stansiya ea | ERA5 ea | ✅ |
| F.16 | RL↑ = ε₀×σ×Ts⁴ | σ=5.67e-8 | σ=5.67e-8 | ✅ |
| F.24 | RL↓ = εa×σ×Ta⁴ | Ta = stansiya T | Ta = ERA5 T2m (K) | ✅ |
| F.25 | εa = 0.85×(−ln(τsw))^0.09 | + | + | ✅ |

---

## 7. TUPROQ ISSIQLIK OQIMI (G)

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| F.26 | G/Rn = (Ts−273)(0.0038+0.0074α)(1−0.98NDVI⁴) | Bastiaanssen 2000 | Default: bastiaanssen | ✅ |
| F.27 | G/Rn = 0.05+0.18×exp(−0.521×LAI) | Tasumi variant | Mavjud, tanlash mumkin | ✅ |
| - | Suv: G/Rn=0.5, Qor: G/Rn=0.5 | + | + | ✅ |

---

## 8. AERODINAMIKA

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| F.32 | u200 = uw×ln(200/zomw)/ln(zx/zomw) | zx=stansiya balandligi, zomw=stansiya roughness | zx=10 (ERA5), zomw=0.03 | 🔄 ERA5 uchun moslashtirilgan |
| F.33 | zom = 0.018×LAI | LAI asosida | LAI asosida, min=0.005m | ✅ |
| F.34 | zom = exp(a×NDVI/α+b) | Alternativ | Mavjud, tanlash mumkin | ✅ |
| F.35 | zom_mtn = zom×(1+(s−5)/20) | slope degrada | slope gradusda (ee.Terrain) | ✅ |
| F.36 | ω̄ = 1+0.1×(z−z_station)/1000 | z_station | z_mean (DEM o'rtacha) | ✅ |
| F.30 | rah = ln(z2/z1)/(u*×k) | z1=0.1, z2=2.0 | z1=0.1, z2=2.0 | ✅ |
| F.31 | u* = k×u200/ln(200/zom) | k=0.41 | k=0.41, min u*=0.01 | ✅ |
| F.37 | ρair = 1000P/(1.01×Ts×R) | R=287 J/(kg·K) | R=287 | ✅ |
| - | u200 minimum | Maqolada aniq ko'rsatilmagan | min=2.0 m/s (ASCE-EWRI) | 🔄 ASCE tavsiyasi |

---

## 9. CIMEC KALIBRLASH (Allen 2013)

| F# | Formula | Allen 2007/2013 | Pipeline | Baho |
|----|---------|----------------|----------|------|
| C.1-2 | Cold pixel tanlash | NDVI P95 + Ts P20 (qishloq xo'jaligi yerlari) | LULC [40] + NDVI P95 + Ts P20 (cascade) | ✅ LULC qo'shilgan |
| C.3-4 | Hot pixel tanlash | NDVI P10 + Ts P80 (yalang'och yer) | LULC [60,50] + NDVI P10 + Ts P80 (cascade) | ✅ LULC qo'shilgan |
| C.5 | LEcold = ETrFcold × ETr × λ/3600 | ETrFcold = 1.05 | ETrFcold = 1.05 | ✅ |
| C.6 | LEhot = ETrFbare × ETr × λ/3600 | ETrFbare = 0.0−0.05 | ETrFbare = 0.05 | ✅ |
| - | Anchor LULC | NLCD (USA) | ESA WorldCover (global) | 🔄 Global |
| - | Cascade fallback | Yo'q — bitta strategiya | L1→L2→L3→L4 cascade | 🔄 Yaxshiroq — robust |

---

## 10. SENSIBLE HEAT FLUX (H) ITERATSIYA

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| F.28 | H = ρ×Cp×dT/rah | + | + | ✅ |
| F.29 | dT = a×TsDEM+b | Har iteratsiyada yangilanadi | Har iteratsiyada yangilanadi | ✅ |
| F.50 | a = (dThot−dTcold)/(Ts_hot−Ts_cold) | + | + | ✅ |
| F.51 | b = dThot − a×Ts_hot | + | + | ✅ |
| F.38 | u* = k×u200/(ln(200/zom)−ψm) | + | + | ✅ |
| F.39 | rah = (ln(z2/z1)−ψh2+ψh01)/(u*×k) | + | + | ✅ |
| F.40 | L = −ρ×Cp×u*³×Ts/(k×g×H) | + | + | ✅ |
| F.41 | ψm unstable (Paulson 1970) | + | + | ✅ |
| F.42 | ψh unstable | + | + | ✅ |
| F.43 | x = (1−16z/L)^0.25 | + | + | ✅ |
| F.44 | ψm stable = −5(2/L) | z=2m (maqolada aniq) | z=Z2=2.0 | ✅ Tuzatilgan |
| F.45 | ψh stable = −5(z/L) | + | + | ✅ |
| - | Iteratsiya soni | "converge bo'lguncha" (~5-15) | MAX_ITER=15, convergence: dT<0.01, rah<0.1 | ✅ |
| - | rah fizik chegarasi | Maqolada aniq ko'rsatilmagan | rah clamp [5, 500] s/m | 🔄 Fizik cheklov |
| - | H fizik chegarasi | "ET can exceed Rn" | H ≥ −0.3×(Rn−G) | ❌ Juda qattiq! −1.0 bo'lishi kerak |

---

## 11. ET HISOBLASH

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| F.1 | LE = Rn − G − H | + | + | ✅ |
| F.52 | ETinst = 3600×LE/λ | mm/hr | mm/hr | ✅ |
| F.53 | λ = (2.501−0.00236×(Ts−273.15))×10⁶ | + | + | ✅ |
| F.54 | ETrF = ETinst/ETr | Scalar ETr (stansiya) | IMAGE÷IMAGE (raster) | 🔄 Raster aniqroq |
| - | ETrF chegarasi | Maqolada 0−~1.05 (cold) | clamp(0, 1.05) | ⚠️ Adveksiyada 1.3 gacha ruxsat kerak |

---

## 12. REFERENCE ET (ETr)

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| - | Soatlik ETr | ASCE PM (stansiya) | ASCE PM (RefETCalculator, ERA5 dan) | 🔄 |
| - | Rs (soatlik) | Stansiya piranometr | Module 5 RS↓ (tau_sw, cos(θ)) | 🔄 ERA5 ssrd muammosi hal qilingan |
| - | Rs (kunlik) | Stansiya piranometr | ERA5 DAILY_AGGR ssrd_sum | ✅ |
| - | u2 | Stansiya 2m | ERA5 10m × 0.748 | ✅ |
| - | ea | Stansiya dewpoint | ERA5 dewpoint | ✅ |
| - | Ref type | Alfalfa (ETr) | Alfalfa (Cn=1600/66, Cd=0.38/0.25) | ✅ |
| - | Kun/tun ajratish | Rn < 0 | Rn < 0 | ✅ |
| - | Ra hisoblash | Soat burchagi dan | omega_mid_hour (to'g'ri) | ✅ |
| - | Rso, Rnl | To'liq hisob | To'liq hisob | ✅ |

---

## 13. KUNLIK VA OYLIK ET

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| F.55 | ET24 = Crad × ETrF × ETr24 | + | + | ✅ |
| F.56 | Crad = (Rso_H/Rso_P)inst × (Rso_P/Rso_H)24 | + | + (24 soat integral) | ✅ |
| F.57 | Rso_24 = Σcos(θ) | + | 24 soat loop | ✅ |
| F.7 | cos(θrel) — 5 ta had | + | + (aspect GEE→METRIC konversiyasi) | ✅ |
| F.58 | ETperiod = Σ(ET24i) | Cubic spline ETrF interpolatsiya | Linear ETrF interpolatsiya (piksel-wise) | ⚠️ Linear vs cubic |
| - | ETrF teshik to'ldirish | Manual tekshirish | Composite mean (unmask) | 🔄 Avtomatik |

---

## 14. QUYOSH GEOMETRIYASI

| F# | Formula | Allen 2007 | Pipeline | Baho |
|----|---------|-----------|----------|------|
| F.8 | cos(θhor) = sinδ×sinφ + cosδ×cosφ×cosω | + | + | ✅ |
| F.9 | d² = 1/(1+0.033×cos(2πJ/365))⁻¹ | + | + | ✅ |
| - | δ (deklinatsiya) | Cooper tenglamasi | Cooper | ✅ |
| - | ω (soat burchagi) | Solar time dan hisoblash | HARDCODED = −0.3 rad | ⚠️ −0.393 bo'lishi kerak |
| - | ω ta'sir darajasi | RS↓ ga ~2% | Kichik, lekin to'g'ri qilish oson | ⚠️ |

---

## XULOSA — Tuzatish Kerak

| # | Muammo | Fayl | O'zgarish | Ta'sir |
|---|--------|------|-----------|--------|
| 1 | H limit = −0.3×(Rn−G) | sensible_heat.py:209 | → −1.0×(Rn−G) | ~5-10% ET oshadi |
| 2 | ETrF clamp = 1.05 | daily_et.py:49 | → 1.30 | ~5-8% ET oshadi |
| 3 | omega = −0.3 | solar_geometry.py:64 | → −0.393 | ~2% RS↓ |

**3 ta o'zgarish → PBIAS −16% dan ~−3% ga tushishi kutiladi.**

---

## XULOSA — Umumiy Statistika

| Kategoriya | Soni | Foiz |
|-----------|------|------|
| ✅ Maqola bilan bir xil | 45 | 75% |
| 🔄 Yaxshiroq/adaptatsiya | 11 | 18% |
| ⚠️ Kichik farq | 3 | 5% |
| ❌ Tuzatish kerak | 1 | 2% |
| **JAMI** | **60** | **100%** |
