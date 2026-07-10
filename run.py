"""
run.py — METRIC-GEE Idaho pipeline.

Natijalar saqlash:
  outputs/{REGION}_{YYYY-MM}/
    daily_et.csv        — kunlik ETrF, ETr, ET
    etrf_images.csv     — har Landsat tasvir uchun anchor + ETrF
    summary.json        — oylik statistika
    pipeline.log        — to'liq konsol chiqishi

GEE Drive export:
  Export toggol qiling: EXPORT_TO_DRIVE = True
  Papka: EXPORT_FOLDER  (Drive ichida)
"""
import ee
import sys
import os
import json
import csv
import io
from datetime import datetime

# ─────────────────────────────────────────────────────────────
# SOZLAMALAR — shu yerni o'zgartiring
# ─────────────────────────────────────────────────────────────

WRS_PATH = 40
WRS_ROW  = 30

START = "2024-07-01"
END   = "2024-07-31"

CLOUD_MAX      = 50     # % — Landsat metadata (butun scene)
ROI_CLOUD_MAX  = 20     # % — cropland ustida piksel-darajasida bulut
LAPSE_RATE     = 6.5    # K/km
ETRF_BARE      = 0.05
TIMEZONE_LON   = -114.0  # Idaho UTC-7 → -7 × 15 = -105

# Albedo usuli — qo'lda tanlang:
#   'liang'   — Liang (2001), barcha vaznlar musbat, universal tavsiya
#   'olmedo'  — Olmedo (2012), L8/9 OLI, manfiy green vazn bor
#   'ke'      — Ke (2016), L8, B_UB (coastal/ultra-blue) kanal bilan
#   'tasumi'  — Tasumi (2008), METRIC Idaho, F.15 ruhida
#   'avg3'    — olmedo + liang + ke o'rtachasi, eng barqaror
ALBEDO_METHOD  = "olmedo"

EXPORT_TO_DRIVE = True   # ET rasterini Google Drive ga export qilish
EXPORT_SCALE    = 30     # m  (katta tile uchun 100 yoki 300 ishlating)
EXPORT_FOLDER   = "METRIC_output_2026"  # Drive papkasi

# Chiqish papkasi — skript bilan bir joyda
BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
OUT_ROOT    = os.path.join(BASE_DIR, "outputs")

# ─────────────────────────────────────────────────────────────
# LOG TIZIMI — konsol + fayl bir vaqtda
# ─────────────────────────────────────────────────────────────

class Tee:
    """stdout ni ekranga ham, faylga ham yozadi."""
    def __init__(self, *streams):
        self.streams = streams
    def write(self, data):
        for s in self.streams:
            s.write(data)
    def flush(self):
        for s in self.streams:
            s.flush()

# ─────────────────────────────────────────────────────────────
# ASOSIY FUNKSIYA
# ─────────────────────────────────────────────────────────────

def main():
    # 1. Output papkasi yaratish
    month_tag = datetime.strptime(START, "%Y-%m-%d").strftime("%Y-%m")
    region_tag = f"P{WRS_PATH:03d}R{WRS_ROW:03d}"
    out_dir = os.path.join(OUT_ROOT, f"{region_tag}_{month_tag}")
    os.makedirs(out_dir, exist_ok=True)

    log_path = os.path.join(out_dir, "pipeline.log")
    log_file = open(log_path, "w", encoding="utf-8")
    sys.stdout = Tee(sys.__stdout__, log_file)

    print(f"{'='*60}")
    print(f"METRIC-GEE  |  Idaho WRS Path={WRS_PATH} Row={WRS_ROW}")
    print(f"Sana:  {START}  →  {END}")
    print(f"Chiqish: {out_dir}")
    print(f"{'='*60}\n")

    # 2. GEE ishga tushirish
    ee.Initialize()
    print("GEE initialized\n")

    from pipeline import METRICPipeline
    from metric_gee.config.settings import Settings
    from metric_gee.utils.export import to_drive

    # 3. WRS tile geometry
    print(f"WRS-2 tile geometriyasi olinmoqda...")
    dummy = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
             .filter(ee.Filter.eq('WRS_PATH', WRS_PATH))
             .filter(ee.Filter.eq('WRS_ROW', WRS_ROW))
             .first())
    geometry = dummy.geometry()

    bounds = geometry.bounds().getInfo()['coordinates'][0]
    lons = [p[0] for p in bounds]
    lats = [p[1] for p in bounds]
    print(f"  Tile bbox: lon=[{min(lons):.2f}, {max(lons):.2f}]  "
          f"lat=[{min(lats):.2f}, {max(lats):.2f}]")

    # 4. Sozlamalar
    cfg = Settings(
        cloud_cover_max  = CLOUD_MAX,
        roi_cloud_max    = ROI_CLOUD_MAX,
        albedo_method    = ALBEDO_METHOD,
        g_method         = "bastiaanssen",
        zom_method       = "lai",
        lapse_rate       = LAPSE_RATE,
        etrf_bare        = ETRF_BARE,
        timezone_lon_deg = TIMEZONE_LON,
        wrs_path         = WRS_PATH,
        wrs_row          = WRS_ROW,
        etrf_interpolation="article",   # ← 
        rl_down_temp="ts",       # pyMETRIC usuli (default)
        # rl_down_temp="t_cold", # cold piksel Ts
        # rl_down_temp="era5",   # ERA5 Ta (eski usul)
    )

    # 5. Pipeline
    print(f"\n{'='*60}")
    print(f"Pipeline ishga tushdi...")
    print(f"{'='*60}\n")
    t0 = datetime.now()

    pipe = METRICPipeline(settings=cfg)
    results = pipe.run_monthly(geometry, START, END)

    elapsed = (datetime.now() - t0).total_seconds()

    if 'error' in results:
        print(f"\nXATO: {results['error']}")
        sys.stdout = sys.__stdout__
        log_file.close()
        sys.exit(1)

    # ─── 6. Natijalar yig'ish ──────────────────────────────
    monthly     = results['monthly_et']
    daily_res   = results.get('daily_results', [])
    etrf_res    = results.get('etrf_results', [])
    et_mean     = monthly.get('et_monthly_mean', 0) or 0

    # ─── 7. Kunlik CSV ────────────────────────────────────
    daily_csv = os.path.join(out_dir, "daily_et.csv")
    with open(daily_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["date", "etrf", "etr_mm", "et_mm", "is_image"])
        w.writeheader()
        img_dates = set(r['date'] for r in etrf_res)
        for row in daily_res:
            w.writerow({
                "date":     row['date'],
                "etrf":     f"{row['etrf']:.4f}",
                "etr_mm":   f"{row['etr']:.3f}",
                "et_mm":    f"{row['et']:.3f}",
                "is_image": "1" if row['date'] in img_dates else "0",
            })
    print(f"\n  Kunlik CSV saqlandi: {daily_csv}")

    # ─── 8. Tasvir ETrF CSV ───────────────────────────────
    img_csv = os.path.join(out_dir, "etrf_images.csv")
    with open(img_csv, "w", newline="", encoding="utf-8") as f:
        fields = ["date", "etrf_tile", "etr_mm_hr", "le_wm2",
                  "rn_wm2", "g_wm2", "h_wm2"]
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for r in etrf_res:
            s = r.get('stats', {})
            w.writerow({
                "date":      r['date'],
                "etrf_tile": f"{r.get('etrf_mean', 0):.4f}",
                "etr_mm_hr": f"{r.get('etr_hourly_val', 0):.4f}",
                "le_wm2":    f"{s.get('LE', 0):.2f}",
                "rn_wm2":    f"{s.get('Rn', 0):.2f}",
                "g_wm2":     f"{s.get('G', 0):.2f}",
                "h_wm2":     f"{s.get('H', 0):.2f}",
            })
    print(f"  Tasvir ETrF CSV saqlandi: {img_csv}")

    # ─── 9. Summary JSON ──────────────────────────────────
    summary = {
        "region":        f"WRS Path={WRS_PATH} Row={WRS_ROW}",
        "start":         START,
        "end":           END,
        "run_time_sec":  round(elapsed, 1),
        "n_images":      len(etrf_res),
        "image_dates":   [r['date'] for r in etrf_res],
        "et_monthly_mm": round(et_mean, 2),
        "etrf_images": [
            {
                "date":       r['date'],
                "etrf_tile":  round(r.get('etrf_mean', 0), 4),
                "etr_mm_hr":  round(r.get('etr_hourly_val', 0), 4),
                "le_wm2":     round(r.get('stats', {}).get('LE', 0), 1),
                "rn_wm2":     round(r.get('stats', {}).get('Rn', 0), 1),
                "g_wm2":      round(r.get('stats', {}).get('G', 0), 1),
                "h_wm2":      round(r.get('stats', {}).get('H', 0), 1),
            }
            for r in etrf_res
        ],
    }

    summary_path = os.path.join(out_dir, "summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"  Summary JSON saqlandi: {summary_path}")

    # ─── 10. Konsol natija ────────────────────────────────
    print(f"\n{'='*60}")
    print(f"OYLIK NATIJA  ({region_tag}  {month_tag})")
    print(f"{'='*60}")
    print(f"  Tasvir soni:  {len(etrf_res)} ta")
    print(f"  ET mean:      {et_mean:.1f} mm/oy")
    print(f"  Ish vaqti:    {elapsed/60:.1f} daq")

    print("\n  Tasvirlar:")
    for r in etrf_res:
        s = r.get('stats', {})
        print(f"    {r['date']}:  "
              f"ETrF={r.get('etrf_mean', 0):.3f}  "
              f"Rn={s.get('Rn', 0):.0f}  G={s.get('G', 0):.0f}  "
              f"H={s.get('H', 0):.0f}  LE={s.get('LE', 0):.0f} W/m²  "
              f"ETr={r.get('etr_hourly_val', 0):.3f} mm/hr")

    # ─── 11. Google Drive Export ──────────────────────────
    if EXPORT_TO_DRIVE:
        et_img = monthly.get('et_monthly_image')
        if et_img is not None:
            desc = f"METRIC_ET_{region_tag}_{month_tag}"
            try:
                task = to_drive(
                    image       = et_img.rename('ET_monthly_mm'),
                    description = desc,
                    folder      = EXPORT_FOLDER,
                    geometry    = geometry,
                    scale       = EXPORT_SCALE,
                )
                print(f"\n  GEE Drive export boshlandi:")
                print(f"    Task: {desc}")
                print(f"    Papka: {EXPORT_FOLDER}/")
                print(f"    Scale: {EXPORT_SCALE} m")
                print(f"    Status: {task.status()['state']}")

                # Task ID ni saqlash
                task_info = {
                    "task_id":   task.id,
                    "task_name": desc,
                    "scale_m":   EXPORT_SCALE,
                    "folder":    EXPORT_FOLDER,
                    "status":    task.status()['state'],
                }
                task_path = os.path.join(out_dir, "export_task.json")
                with open(task_path, "w", encoding="utf-8") as f:
                    json.dump(task_info, f, indent=2)
                print(f"    Task ID saqlandi: {task_path}")
            except Exception as e:
                print(f"\n  Export xato: {e}")
        else:
            print("\n  ESLATMA: et_monthly_image yo'q, export o'tkazildi.")

    print(f"\n{'='*60}")
    print(f"Barcha fayllar: {out_dir}")
    print(f"{'='*60}\n")

    # stdout ni qaytarish
    sys.stdout = sys.__stdout__
    log_file.close()


if __name__ == "__main__":
    main()
