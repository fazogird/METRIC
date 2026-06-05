"""
test_idaho.py — Idaho validatsiya + OpenET solishtirish.

WRS-2 Path=40, Row=30 — Boise/Magic Valley, Idaho hududi.
Geometry WRS-2 tile dan avtomatik olinadi.
"""
import ee
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ee.Initialize()
print("✅ GEE initialized\n")

from pipeline import METRICPipeline
from metric_gee.config.settings import Settings

# ===== WRS-2 TILE GEOMETRY =====
WRS_PATH = 40
WRS_ROW  = 30

print(f"📍 WRS-2 Path={WRS_PATH} Row={WRS_ROW} — tile geometry olinmoqda...")

print(f"📍 WRS-2 Path={WRS_PATH} Row={WRS_ROW} — tile geometry olinmoqda...")

# Landsat 8 to'plamidan kerakli Path va Row ga tegishli 1-rasmni topamiz
dummy_image = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2") \
    .filter(ee.Filter.eq('WRS_PATH', WRS_PATH)) \
    .filter(ee.Filter.eq('WRS_ROW', WRS_ROW)) \
    .first()

# O'sha rasmning tayyor chegarasini olamiz
GEOMETRY = dummy_image.geometry()

# Tile bounds ni ko'rsatish
bounds = GEOMETRY.bounds().getInfo()['coordinates'][0]
lons = [p[0] for p in bounds]
lats = [p[1] for p in bounds]
print(f"  Tile bbox: lon=[{min(lons):.2f}, {max(lons):.2f}]  lat=[{min(lats):.2f}, {max(lats):.2f}]")

START = "2024-07-01"
END   = "2024-07-31"

# ===== CONFIG =====
cfg = Settings(
    cloud_cover_max=15,
    albedo_method="olmedo",
    g_method="bastiaanssen",
    zom_method="lai",
    lapse_rate=6.5,
    etrf_bare=0.05,
    timezone_lon_deg=-105,   # Idaho UTC-7 → -7×15 = -105
    wrs_path=WRS_PATH,
    wrs_row=WRS_ROW,
)

# ===== METRIC PIPELINE =====
print(f"\n{'='*50}")
print(f"METRIC Pipeline  |  {START} → {END}")
print("="*50)

pipe = METRICPipeline(settings=cfg)
results = pipe.run_monthly(GEOMETRY, START, END)

if 'error' in results:
    print(f"❌ {results['error']}")
    sys.exit(1)

monthly = results['monthly_et']

print(f"\n{'='*50}")
print("NATIJA  (NDVI≥0.4 qishloq xo'jaligi piksellar)")
print("="*50)
print(f"  Oylik ET:   {monthly['et_monthly']:.1f} mm")
print(f"  Oylik ETr:  {monthly['etr_monthly']:.1f} mm")
print(f"  ETrF (w.):  {monthly['etrf_weighted']:.3f}")

# ETrF tasvir qiymatlarini ko'rsatish
print("\n  Tasvirlar bo'yicha ETrF:")
for r in results.get('etrf_results', []):
    print(f"    {r['date']}: "
          f"ETrF_agri={r['etrf_mean']:.3f}  "
          f"P75={r.get('etrf_agri_p75',0):.3f}  "
          f"tile_P90={r.get('etrf_tile_p90',0):.3f}  "
          f"ETr={r.get('etr_hourly_val',0):.3f} mm/hr  "
          f"LE={r['stats']['LE']:.0f} W/m²")

metric_et = monthly['et_monthly']

# ===== OpenET VALIDATSIYA =====
print(f"\n{'='*50}")
print("OpenET Validatsiya  (CONUS Ensemble v2)")
print("="*50)
try:
    openet_col = ee.ImageCollection("OpenET/ENSEMBLE/CONUS/GRIDMET/MONTHLY/v2_0")
    openet_img = openet_col.filterBounds(GEOMETRY).filterDate(START, END).first()

    # Band nomi: et_ensemble_mad (mm/oy)
    ET_BAND = 'et_ensemble_mad'
    openet_stats = openet_img.select(ET_BAND).reduceRegion(
        ee.Reducer.mean()
                  .combine(ee.Reducer.percentile([25, 50, 75]), sharedInputs=True),
        GEOMETRY, 30, maxPixels=1e9
    ).getInfo()

    openet_mean = openet_stats.get(f'{ET_BAND}_mean') or openet_stats.get(ET_BAND) or 0
    p25  = openet_stats.get(f'{ET_BAND}_p25',  0)
    p50  = openet_stats.get(f'{ET_BAND}_p50',  0)
    p75  = openet_stats.get(f'{ET_BAND}_p75',  0)

    diff     = metric_et - openet_mean
    diff_pct = diff / max(openet_mean, 1) * 100

    print(f"  METRIC:  {metric_et:.1f} mm/oy")
    print(f"  OpenET:  {openet_mean:.1f} mm/oy  (P25={p25:.0f} P50={p50:.0f} P75={p75:.0f})")
    print(f"  Farq:    {diff:+.1f} mm  ({diff_pct:+.1f}%)")

    if abs(diff_pct) < 10:
        print("  ✅ Yaxshi  (<10%)")
    elif abs(diff_pct) < 20:
        print("  ⚠️ O'rtacha (10-20%)")
    else:
        print("  ❌ Katta farq (>20%) — kalibrlash kerak")

except Exception as e:
    print(f"  ⚠️ OpenET xato: {e}")
