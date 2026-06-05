"""
METRIC-GEE Test — Module 0-3
Pipeline modullardan import qilib test.

Ishlatish (Colab):
  !tar xzf metric_gee_v1.1.tar.gz
  %cd /content   (yoki metric_gee ning parent papkasi)
  %run metric_gee/test_mod0_to_3.py
"""
import ee
import sys
import os

# Parent papkani sys.path ga qo'shish (import ishlashi uchun)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ee.Initialize()
print("✅ GEE initialized")

# =============================================
# Pipeline modullarni import
# =============================================
from metric_gee.config.settings import Settings
from metric_gee.config.constants import LAPSE_RATE_DEFAULT
from metric_gee.inputs.landsat import LandsatLoader
from metric_gee.inputs.era5 import ERA5Loader
from metric_gee.inputs.dem import DEMLoader
from metric_gee.inputs.landuse import LandUseLoader
from metric_gee.core.vegetation import VegetationIndices
from metric_gee.core.albedo import AlbedoCalculator
from metric_gee.core.surface_temp import SurfaceTemperature
from metric_gee.core.emissivity import EmissivityCalculator

print("✅ Barcha modullar import qilindi")

# =============================================
# CONFIG
# =============================================
cfg = Settings(
    cloud_cover_max=20,
    albedo_method="olmedo",
)

GEOMETRY = ee.Geometry.Rectangle([-114.45, 42.35, -114.15, 42.55])
DATE = "2024-07-15"

print(f"📍 Idaho Twin Falls | 📅 {DATE}")

# =============================================
# MOD 0: INPUTS
# =============================================
print("\n" + "="*60)
print("MOD 0: Inputs — Landsat, DEM")
print("="*60)

landsat = LandsatLoader(cfg)
dem_loader = DEMLoader(cfg)

image = landsat.get_image(GEOMETRY, DATE)
dem = dem_loader.get_terrain(GEOMETRY)
z_datum = dem_loader.get_mean_elevation(GEOMETRY)

# Landsat metadata
props = image.getInfo()['properties']
print(f"  Satellite: {props.get('SPACECRAFT_ID', '?')}")
print(f"  Date: {props.get('DATE_ACQUIRED', '?')}")
print(f"  Cloud: {props.get('CLOUD_COVER', 0):.1f}%")
print(f"  Sun elev: {props.get('SUN_ELEVATION', 0):.1f}°")

# SR band check
sr_stats = image.select(['B_BLUE','B_GREEN','B_RED','B_NIR','B_SWIR1','B_SWIR2']).reduceRegion(
    reducer=ee.Reducer.mean(), geometry=GEOMETRY, scale=30, maxPixels=1e9).getInfo()

print(f"  SR means: ", end="")
for k, v in sr_stats.items():
    print(f"{k}={v:.4f}  ", end="")
print()

# Ts check
ts_check = image.select('B_THERMAL').reduceRegion(
    reducer=ee.Reducer.mean().combine(ee.Reducer.minMax(), sharedInputs=True),
    geometry=GEOMETRY, scale=30, maxPixels=1e9).getInfo()
print(f"  Ts: {ts_check['B_THERMAL_min']:.1f} / {ts_check['B_THERMAL_mean']:.1f} / {ts_check['B_THERMAL_max']:.1f} K")
print(f"  z_datum: {z_datum.getInfo():.0f} m")
print("✅ Mod 0 OK")

# =============================================
# MOD 1: VEGETATION
# =============================================
print("\n" + "="*60)
print("MOD 1: Vegetation — NDVI, SAVI, LAI")
print("="*60)

veg = VegetationIndices(cfg)
ndvi, savi, lai = veg.calc_all(image)

veg_stats = ndvi.addBands(savi).addBands(lai).reduceRegion(
    reducer=ee.Reducer.mean().combine(
        ee.Reducer.percentile([5, 95]), sharedInputs=True),
    geometry=GEOMETRY, scale=30, maxPixels=1e9).getInfo()

print(f"  NDVI: p5={veg_stats['NDVI_p5']:.3f}  mean={veg_stats['NDVI_mean']:.3f}  p95={veg_stats['NDVI_p95']:.3f}")
print(f"  SAVI: mean={veg_stats['SAVI_mean']:.3f}")
print(f"  LAI:  mean={veg_stats['LAI_mean']:.2f}  p95={veg_stats['LAI_p95']:.2f}")
print("✅ Mod 1 OK")

# =============================================
# MOD 2: ALBEDO
# =============================================
print("\n" + "="*60)
print("MOD 2: Albedo — Olmedo vs Liang")
print("="*60)

alb = AlbedoCalculator(cfg)

albedo_olmedo = alb.calc_albedo(image, method='olmedo')
albedo_liang = alb.calc_albedo(image, method='liang')
albedo_diff = albedo_olmedo.subtract(albedo_liang).rename('diff')

alb_stats = albedo_olmedo.addBands(albedo_liang).addBands(albedo_diff).reduceRegion(
    reducer=ee.Reducer.mean().combine(ee.Reducer.minMax(), sharedInputs=True),
    geometry=GEOMETRY, scale=30, maxPixels=1e9).getInfo()

print(f"  Olmedo: min={alb_stats['albedo_min']:.4f}  mean={alb_stats['albedo_mean']:.4f}  max={alb_stats['albedo_max']:.4f}")
print(f"  Liang:  min={alb_stats['albedo_min_1']:.4f}  mean={alb_stats['albedo_mean_1']:.4f}  max={alb_stats['albedo_max_1']:.4f}")
print(f"  Farq:   mean={alb_stats['diff_mean']:.4f} (Olmedo - Liang)")
print("✅ Mod 2 OK")

# =============================================
# MOD 3: SURFACE TEMPERATURE
# =============================================
print("\n" + "="*60)
print("MOD 3: Ts → TsDEM (lapse correction)")
print("="*60)

st = SurfaceTemperature(cfg)

ts = st.calc_ts(image)
ts_dem = st.calc_ts_dem(ts, dem, z_datum)

ts_stats = ts.addBands(ts_dem).reduceRegion(
    reducer=ee.Reducer.mean().combine(ee.Reducer.minMax(), sharedInputs=True),
    geometry=GEOMETRY, scale=30, maxPixels=1e9).getInfo()

print(f"  Ts:    min={ts_stats['Ts_min']:.1f}  mean={ts_stats['Ts_mean']:.1f}  max={ts_stats['Ts_max']:.1f} K")
print(f"  TsDEM: min={ts_stats['TsDEM_min']:.1f}  mean={ts_stats['TsDEM_mean']:.1f}  max={ts_stats['TsDEM_max']:.1f} K")

ts_range = ts_stats['TsDEM_max'] - ts_stats['TsDEM_min']
print(f"  TsDEM range: {ts_range:.1f} K (CIMEC uchun >15K yaxshi)")

# =============================================
# MOD 4: EMISSIVITY
# =============================================
print("\n" + "="*60)
print("MOD 4: Emissivity — ε₀, εNB")
print("="*60)

emiss = EmissivityCalculator()
e0, enb = emiss.calc_emissivity(lai, ndvi)

em_stats = e0.addBands(enb).reduceRegion(
    reducer=ee.Reducer.mean().combine(ee.Reducer.minMax(), sharedInputs=True),
    geometry=GEOMETRY, scale=30, maxPixels=1e9).getInfo()

print(f"  ε₀:  min={em_stats['e0_min']:.4f}  mean={em_stats['e0_mean']:.4f}  max={em_stats['e0_max']:.4f}")
print(f"  εNB: min={em_stats['enb_min']:.4f}  mean={em_stats['enb_mean']:.4f}  max={em_stats['enb_max']:.4f}")
print("✅ Mod 4 OK")

# =============================================
# XULOSA
# =============================================
print("\n" + "="*60)
print("XULOSA — Mod 0-4 pipeline test")
print("="*60)
print(f"""
  Satellite:  {props.get('SPACECRAFT_ID', '?')} | {props.get('DATE_ACQUIRED', '?')}
  NDVI mean:  {veg_stats['NDVI_mean']:.3f}
  Albedo:     {alb_stats['albedo_mean']:.4f} (Olmedo)
  Ts mean:    {ts_stats['Ts_mean']:.1f} K ({ts_stats['Ts_mean']-273.15:.1f} °C)
  TsDEM range:{ts_range:.1f} K
  ε₀ mean:    {em_stats['e0_mean']:.4f}
  
  ✅ Barcha modullar pipeline dan import qilib ishladi!
  📌 Keyingi: Mod 5-6 (Rn, G) → ERA5 integratsiya test
""")
