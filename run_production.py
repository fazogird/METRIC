"""
run_production.py — METRIC-GEE Ko'p tile × Ko'p oy batch runner.

Drive naming pattern:
    METRIC_monthly_ET_{YYYY}-{MM}_P{path:03d}_R{row:03d}
    Misol: METRIC_monthly_ET_2025-04_P040_R030

Lokal output:
    outputs/{YYYY}-{MM}_P{path:03d}_R{row:03d}/
      daily_et.csv      — kunlik ETrF, ETr, ET
      summary.json      — oylik statistika + export task ID
      export_task.json  — Drive task ID
      pipeline.log      — to'liq konsol

SKIP_EXISTING = True → allaqachon bajarilgan kombinatsiyalar o'tkaziladi.
"""
import ee
import os
import sys
import json
import csv
import calendar
from datetime import datetime


# ================================================================
# SOZLAMALAR
# ================================================================

# WRS-2 tile (Path, Row) ro'yxati
TILES = [
    # (156, 33),     # Qashqadaryo 156_033
    # (155, 33),      # Qashqadaryo 155_033 
    # (156, 32),   # Buxoro 156_032
    # (155, 32),   # Samarqand 155_032
    # (155, 33)
    # (42, 30),
    # (42, 29),
    (152, 32),
    (153, 32),
]

# ════════════════════════════════════════
#Samarqand 
# ════════════════════════════════════════
# tile_groups = {
#     'west_156': ['156_032', '156_033'],   # ← G'arb (Bukhara/Navoi tomon)
#     'east_155': ['155_032', '155_033'],   # ← Sharq (Samarqand tomon)
# }
        
        
# Mavsumiy oylar: [(year, month), ...]
SEASONAL_MONTHS = [
    (2026, 3),    # Mart
    (2026, 4),    # Aprel
    (2026, 5),    # May
#     (2026, 6),    # Iyun
#     (2026, 7),    # Iyul
#     (2026, 8),    # Avgust
#     (2026, 9),    # Sentyabr
#     (2026, 10),   # Oktyabr
]

# ── Pipeline sozlamalari ──────────────────────────────────────────
CLOUD_MAX      = 70     # % — Landsat metadata (butun scene)
ROI_CLOUD_MAX  = 20     # % — cropland ustida piksel-darajasida bulut
LAPSE_RATE     = 6.5    # K/km
ETRF_BARE      = 0.05
TIMEZONE_LON   = 66.0  # -105.0  # Idaho UTC-7 → -7 × 15 = -105.0
                          # O'zbekiston UTC+5 → +5 × 15 = +75.0  Qashqadaryo 66.0

# ── Albedo usuli ─────────────────────────────────────────────────
#   'liang'   — Liang (2001), barcha vaznlar musbat, universal
#   'olmedo'  — Olmedo (2012), L8/9 OLI
#   'ke'      — Ke (2016), L8, B_UB kanal
#   'tasumi'  — Tasumi (2008), METRIC Idaho
#   'avg3'    — olmedo + liang + ke o'rtachasi
ALBEDO_METHOD  = "avg3"

# ── Export ───────────────────────────────────────────────────────
EXPORT_TO_DRIVE = True
EXPORT_SCALE    = 30        # m  (katta tile uchun 100 yoki 300)
EXPORT_FOLDER   = "METRIC_Fargona_2026" # Google Drive folder nomi

# ── Xizmat ───────────────────────────────────────────────────────
SKIP_EXISTING   = True      # True → summary.json mavjud bo'lsa o'tkazib yuborish

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_ROOT = os.path.join(BASE_DIR, "outputs")


# ================================================================
# YORDAMCHI FUNKSIYALAR
# ================================================================

def make_tag(path: int, row: int, year: int, month: int) -> str:
    """METRIC_monthly_ET_2025-04_P040_R030"""
    return f"METRIC_monthly_ET_{year}-{month:02d}_P{path:03d}_R{row:03d}"


def month_dates(year: int, month: int):
    """(start_str, end_str) — oyning birinchi va oxirgi kuni."""
    last = calendar.monthrange(year, month)[1]
    return f"{year}-{month:02d}-01", f"{year}-{month:02d}-{last:02d}"


def get_tile_geometry(path: int, row: int):
    """WRS-2 tile geometriyasini Landsat katalogidan olish."""
    dummy = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
             .filter(ee.Filter.eq('WRS_PATH', path))
             .filter(ee.Filter.eq('WRS_ROW',  row))
             .first())
    return dummy.geometry()


def is_done(out_dir: str) -> bool:
    """summary.json mavjudligi → allaqachon bajarilgan."""
    return os.path.exists(os.path.join(out_dir, "summary.json"))


# ================================================================
# BITTA TILE × OY — TO'LIQ PIPELINE
# ================================================================

def run_tile_month(path: int, row: int, year: int, month: int,
                   log_fh=None) -> dict | None:
    """
    Bitta (path, row) × (year, month) kombinatsiyasi uchun pipeline.

    Returns:
        dict — natija (et_monthly_mm, n_images, ...) yoki None (xatolik)
    """
    from pipeline import METRICPipeline
    from metric_gee.config.settings import Settings
    from metric_gee.utils.export import to_drive

    tag        = make_tag(path, row, year, month)
    start, end = month_dates(year, month)
    out_dir    = os.path.join(OUT_ROOT, f"{year}-{month:02d}_P{path:03d}_R{row:03d}")
    os.makedirs(out_dir, exist_ok=True)

    # Mahalliy log fayl + konsol
    tile_log_path = os.path.join(out_dir, "pipeline.log")
    tile_log_fh   = open(tile_log_path, "w", encoding="utf-8")

    class Tee:
        def __init__(self, *streams):
            self.streams = streams
        def write(self, data):
            for s in self.streams:
                s.write(data)
        def flush(self):
            for s in self.streams:
                s.flush()

    tee = Tee(sys.__stdout__, tile_log_fh)
    if log_fh:
        tee = Tee(sys.__stdout__, tile_log_fh, log_fh)
    sys.stdout = tee

    def log(msg: str):
        print(msg)

    try:
        log(f"\n{'─'*60}")
        log(f"  {tag}")
        log(f"  Sana: {start}  →  {end}")
        log(f"  Albedo: {ALBEDO_METHOD}  |  Scale: {EXPORT_SCALE}m")
        log(f"{'─'*60}")

        t0 = datetime.now()

        # Tile geometriyasi
        geometry = get_tile_geometry(path, row)
        bounds   = geometry.bounds().getInfo()['coordinates'][0]
        lons     = [p[0] for p in bounds]
        lats     = [p[1] for p in bounds]
        log(f"  Bbox: lon=[{min(lons):.2f}, {max(lons):.2f}]  "
            f"lat=[{min(lats):.2f}, {max(lats):.2f}]")

        # Sozlamalar
        cfg = Settings(
            cloud_cover_max  = CLOUD_MAX,
            roi_cloud_max    = ROI_CLOUD_MAX,
            albedo_method    = ALBEDO_METHOD,
            g_method         = "bastiaanssen",
            zom_method       = "lai",
            lapse_rate       = LAPSE_RATE,
            etrf_bare        = ETRF_BARE,
            timezone_lon_deg = TIMEZONE_LON,
            wrs_path         = path,
            wrs_row          = row,
            etrf_interpolation="article",   # ← 
            rl_down_temp="ts",       # pyMETRIC usuli (default)
            # rl_down_temp="t_cold", # cold piksel Ts
            # rl_down_temp="era5",   # ERA5 Ta (eski usul)
            lai_method="METRIC",    # 11×SAVI³ (METRIC usuli, default)
        )

        # Pipeline
        pipe    = METRICPipeline(settings=cfg)
        results = pipe.run_monthly(geometry, start, end)

        elapsed = (datetime.now() - t0).total_seconds()

        if 'error' in results:
            log(f"\n  ❌ XATO: {results['error']}")
            return None

        monthly   = results['monthly_et']
        daily_res = results.get('daily_results', [])
        etrf_res  = results.get('etrf_results', [])
        et_mean   = float(monthly.get('et_monthly_mean') or 0)

        # ── Kunlik CSV ────────────────────────────────────────────
        daily_csv = os.path.join(out_dir, "daily_et.csv")
        img_dates = {r['date'] for r in etrf_res}
        with open(daily_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(
                f, fieldnames=["date", "etrf", "etr_mm", "crad", "et_mm", "is_image"])
            w.writeheader()
            for d in daily_res:
                w.writerow({
                    "date":     d['date'],
                    "etrf":     f"{d.get('etrf', 0):.4f}",
                    "etr_mm":   f"{d.get('etr',  0):.3f}",
                    "crad":     f"{d.get('crad', 1):.4f}",
                    "et_mm":    f"{d.get('et',   0):.3f}",
                    "is_image": "1" if d['date'] in img_dates else "0",
                })

        # ── Tasvir ETrF CSV ───────────────────────────────────────
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
                    "g_wm2":     f"{s.get('G',  0):.2f}",
                    "h_wm2":     f"{s.get('H',  0):.2f}",
                })

        # ── Summary JSON ──────────────────────────────────────────
        summary = {
            "tag":           tag,
            "path":          path,
            "row":           row,
            "year":          year,
            "month":         month,
            "start":         start,
            "end":           end,
            "run_time_sec":  round(elapsed, 1),
            "n_images":      len(etrf_res),
            "image_dates":   [r['date'] for r in etrf_res],
            "et_monthly_mm": round(et_mean, 2),
            "albedo_method": ALBEDO_METHOD,
            "etrf_images": [
                {
                    "date":      r['date'],
                    "etrf_tile": round(r.get('etrf_mean', 0), 4),
                    "etr_mm_hr": round(r.get('etr_hourly_val', 0), 4),
                    "crad":      round(r.get('crad_mean', 1.0), 4),
                    "le_wm2":    round(r.get('stats', {}).get('LE', 0), 1),
                    "rn_wm2":    round(r.get('stats', {}).get('Rn', 0), 1),
                    "h_wm2":     round(r.get('stats', {}).get('H',  0), 1),
                }
                for r in etrf_res
            ],
        }
        summary_path = os.path.join(out_dir, "summary.json")
        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)

        log(f"\n  ✅  ET={et_mean:.1f} mm  |  {len(etrf_res)} tasvir  "
            f"|  {elapsed/60:.1f} daq")
        log(f"  Lokal: {out_dir}")

        # ── Google Drive Export ───────────────────────────────────
        if EXPORT_TO_DRIVE:
            et_img = monthly.get('et_monthly_image')
            if et_img is not None:
                desc = tag   # METRIC_monthly_ET_2025-04_P040_R030
                try:
                    task = to_drive(
                        image       = et_img.rename('ET_monthly_mm'),
                        description = desc,
                        folder      = EXPORT_FOLDER,
                        geometry    = geometry,
                        scale       = EXPORT_SCALE,
                    )
                    status = task.status()['state']
                    log(f"  📤 Drive: {EXPORT_FOLDER}/{desc}  [{status}]")

                    # Task ID saqlash
                    export_info = {
                        "task_id":     task.id,
                        "description": desc,
                        "folder":      EXPORT_FOLDER,
                        "scale_m":     EXPORT_SCALE,
                        "status":      status,
                    }
                    with open(os.path.join(out_dir, "export_task.json"),
                              "w", encoding="utf-8") as f:
                        json.dump(export_info, f, indent=2)

                    summary['export_task_id'] = task.id
                    summary['export_status']  = status

                except Exception as ex:
                    log(f"  ⚠️  Drive export xato: {ex}")
            else:
                log("  ⚠️  et_monthly_image yo'q — export o'tkazildi.")

        return summary

    except Exception as ex:
        log(f"\n  ❌  {tag}: {ex}")
        import traceback
        log(traceback.format_exc())
        return None

    finally:
        sys.stdout = sys.__stdout__
        tile_log_fh.close()


# ================================================================
# BATCH RUNNER
# ================================================================

def main():
    os.makedirs(OUT_ROOT, exist_ok=True)

    run_ts   = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(OUT_ROOT, f"batch_{run_ts}.log")
    log_fh   = open(log_path, "w", encoding="utf-8")

    def log(msg: str):
        print(msg)
        log_fh.write(msg + "\n")
        log_fh.flush()

    n_total = len(TILES) * len(SEASONAL_MONTHS)
    log("=" * 60)
    log(f"METRIC-GEE Batch  |  {run_ts}")
    log(f"Tillar: {TILES}")
    log(f"Oylar:  {SEASONAL_MONTHS}")
    log(f"Jami:   {n_total} ta kombinatsiya")
    log("=" * 60)

    ee.Initialize()
    log("GEE initialized\n")

    n_done = n_skip = n_err = 0
    all_results = []

    for (path, row) in TILES:
        for (year, month) in SEASONAL_MONTHS:
            tag     = make_tag(path, row, year, month)
            out_dir = os.path.join(OUT_ROOT,
                                   f"{year}-{month:02d}_P{path:03d}_R{row:03d}")

            # Allaqachon bajarilganlarni o'tkazib yuborish
            if SKIP_EXISTING and is_done(out_dir):
                log(f"  ⏩ {tag} — mavjud, o'tkazildi")
                n_skip += 1
                continue

            log(f"\n[{n_done+n_err+1}/{n_total-n_skip}]  {tag}")

            try:
                res = run_tile_month(path, row, year, month, log_fh)
                if res:
                    all_results.append(res)
                    n_done += 1
                else:
                    n_err += 1
            except Exception as ex:
                log(f"  ❌ {tag}: {ex}")
                n_err += 1

    # ── Yakuniy hisobot ───────────────────────────────────────────
    log(f"\n{'='*60}")
    log(f"BATCH YAKUNLANDI  |  {run_ts}")
    log(f"  Jami:             {n_total}")
    log(f"  ✅ Muvaffaqiyatli: {n_done}")
    log(f"  ⏩ O'tkazilgan:    {n_skip}")
    log(f"  ❌ Xatolik:        {n_err}")
    log(f"  Log: {log_path}")
    log("=" * 60)

    if all_results:
        log("\n  Natijalar:")
        for r in all_results:
            log(f"    {r['tag']:50s}  ET={r['et_monthly_mm']:6.1f} mm  "
                f"n={r['n_images']}")

    # Batch summary JSON
    batch_summary = {
        "run_time":  run_ts,
        "n_total":   n_total,
        "n_done":    n_done,
        "n_skipped": n_skip,
        "n_errors":  n_err,
        "results":   all_results,
    }
    batch_json = os.path.join(OUT_ROOT, f"batch_{run_ts}_summary.json")
    with open(batch_json, "w", encoding="utf-8") as f:
        json.dump(batch_summary, f, indent=2, ensure_ascii=False)
    log(f"  Batch summary: {batch_json}")

    log_fh.close()


if __name__ == "__main__":
    main()