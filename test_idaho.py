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
    timezone_lon_deg=114.0,   # Idaho UTC-7 → -7×15 = -105
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
print("NATIJA  (butun tile)")
print("="*50)
print(f"  ET mean:  {monthly['et_monthly_mean']:.1f} mm")

# ETrF tasvir qiymatlarini ko'rsatish
print("\n  Tasvirlar bo'yicha:")
for r in results.get('etrf_results', []):
    print(f"    {r['date']}: "
          f"ETrF={r['stats']['ETrF']:.3f}  "
          f"ETinst={r['stats']['ETinst']:.3f} mm/hr  "
          f"LE={r['stats']['LE']:.0f} W/m²")

metric_et = monthly['et_monthly_mean']

# ===== OpenET Validatsiya =====
print(f"\n{'='*50}")
print("OpenET Validatsiya")
print("="*50)
try:
    openet = (ee.ImageCollection("OpenET/ENSEMBLE/CONUS/GRIDMET/MONTHLY/v2_0")
        .filterBounds(GEOMETRY).filterDate(START, END).first())
    openet_stats = openet.select('et_ensemble_mad').reduceRegion(
        ee.Reducer.mean().combine(
            ee.Reducer.percentile([25, 50, 75, 90]), sharedInputs=True),
        GEOMETRY, 300, maxPixels=1e9).getInfo()

    openet_mean = openet_stats.get('et_mean', 0)
    diff = metric_et - openet_mean
    diff_pct = diff / max(openet_mean, 1) * 100

    print(f"  METRIC:  {metric_et:.1f} mm/oy")
    print(f"  OpenET:  {openet_mean:.1f} mm/oy")
    print(f"  Farq:    {diff:+.1f} mm ({diff_pct:+.1f}%)")

    if abs(diff_pct) < 15:
        print("  ✅ Yaxshi (<15%)")
    elif abs(diff_pct) < 25:
        print("  ⚠️ O'rtacha (15-25%)")
    else:
        print("  ❌ Katta farq (>25%)")
except Exception as e:
    print(f"  ⚠️ OpenET: {e}")