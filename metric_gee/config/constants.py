"""
Fizik konstantalar — METRIC energiya balansi uchun.
Barcha qiymatlar SI birlikda.
"""

# Stefan-Boltzmann constantasi (W m⁻² K⁻⁴)
SIGMA = 5.67e-8

# Von Karman constantasi (dimensionless)
VON_KARMAN = 0.41

# Quyosh constantasi (W m⁻²) — F.3
G_SC = 1367

# Havo issiqlik sig'imi (J kg⁻¹ K⁻¹)
CP = 1004

# Gravitatsiya tezlanishi (m s⁻²) — F.40
GRAVITY = 9.807

# Suv zichligi (kg m⁻³) — F.52
RHO_WATER = 1000

# Gaz constantasi (J kg⁻¹ K⁻¹) — F.37
R_AIR = 287

# Blending height (m) — METRIC 200m
Z_BLEND = 200

# dT hisoblash balandliklari (m) — F.30
Z1 = 0.1   # pastki
Z2 = 2.0   # yuqori

# Yalang'och tuproq minimal zom (m) — F.33
ZOM_MIN = 0.005

# Standart atmosfera temperatura (K) — F.5
T_STD = 293

# Standart bosim (kPa) — F.5
P_STD = 101.3

# Lapse rate default (K/km)
LAPSE_RATE_DEFAULT = 6.5

# Landsat 8/9 Collection 2 Level 2 — scale factors
L8_SR_SCALE = 0.0000275
L8_SR_OFFSET = -0.2
L8_ST_SCALE = 0.00341802
L8_ST_OFFSET = 149.0

# CIMEC thresholds — Allen 2013
CIMEC_NDVI_COLD_PCT = 5      # top 5% NDVI
CIMEC_TS_COLD_PCT = 20       # coldest 20% TsDEM
CIMEC_NDVI_HOT_PCT = 10      # lowest 10% NDVI
CIMEC_TS_HOT_PCT = 20        # hottest 20% TsDEM
CIMEC_TS_TOLERANCE = 0.2     # K
CIMEC_ALBEDO_TOLERANCE = 0.02
CIMEC_ETRF_COLD_DEFAULT = 1.05
CIMEC_NDVI_BARE = 0.15
CIMEC_NDVI_FULL = 0.80
CIMEC_CV_THRESHOLD = 0.15    # NDVI homogeneity filter

# Iteratsiya limiti — Mod.9
MAX_ITERATIONS = 30
H_CONVERGENCE = 0.01         # W/m²
