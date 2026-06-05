"""
run.py — METRIC-GEE ishga tushirish.
Faqat sozlamalar va geometry. Barcha mantiq pipeline.py da.
"""
import ee
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ee.Initialize()
print("✅ GEE initialized\n")

from pipeline import METRICPipeline
from metric_gee.config.settings import Settings

# ===== SOZLAMALAR =====
cfg = Settings(
    cloud_cover_max=15,
    albedo_method="olmedo",
    g_method="bastiaanssen",
    zom_method="lai",
    lapse_rate=6.5,
    etrf_bare=0.05,
    timezone_lon_deg=-105,   # Idaho UTC-7 → -7×15
)

# ===== STUDY AREA =====
GEOMETRY = ee.Geometry.Rectangle([-114.0, 42.35, -113.6, 42.55])
START = "2024-07-01"
END = "2024-07-31"

# ===== PIPELINE =====
pipe = METRICPipeline(settings=cfg)
results = pipe.run_monthly(GEOMETRY, START, END)

# ===== NATIJA =====
if 'error' in results:
    print(f"\n❌ {results['error']}")
else:
    m = results['monthly_et']
    print(f"\n{'='*50}")
    print(f"NATIJA:")
    print(f"  Oylik ET:  {m['et_monthly']:.1f} mm")
    print(f"  Oylik ETr: {m['etr_monthly']:.1f} mm")
    print(f"  ETrF:      {m['etrf_weighted']:.3f}")
