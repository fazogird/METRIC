"""
mosaic.py — METRIC GeoTIFF tillarini mean mozaykaga birlashtirish.

Qo'llab-quvvatlanadigan fayl nomlari:
  A) METRIC_monthly_ET_2025-04_P040_R030.tif            ← to'liq fayl
  B) METRIC_monthly_ET_2025-04_P040_R030-0000023296-0000000000.tif  ← GEE chunk
  C) ET_2025_03-0000000000-0000000000.tif               ← sodda GEE chunk

GEE chunklar → bir oyning barcha fayllari (A+B+C) birgalikda
mosaic_mean ga beriladi: overlap → o'rtacha, nodata → boshqa tile qiymati.

Foydalanish:
  python mosaic.py --input D:/METRIC --output D:/Mosaics
  python mosaic.py --input D:/METRIC --output D:/Mosaics --tiles P040_R030 P041_R030
  python mosaic.py --input D:/METRIC --output D:/Mosaics --months 2025-07 2025-08
"""
import os
import re
import sys
import shutil
import argparse
from collections import defaultdict

import numpy as np

try:
    from osgeo import gdal
    gdal.UseExceptions()
except ImportError:
    sys.exit("GDAL o'rnatilmagan:  conda install -c conda-forge gdal")


# ================================================================
# FAYL NOMLARI PATTERN (2 xil format + GEE chunk suffix)
# ================================================================

# Format A/B: METRIC_monthly_ET_2025-04_P040_R030[-\d{10}-\d{10}].tif

# # Format A/B — METRIC to'liq nom (chunk bor yoki yo'q):
# r"METRIC_monthly_ET_(?P<month>\d{4}-\d{2})_(?P<tile>P\d{3}_R\d{3})"
# r"(?:-\d{10}-\d{10})?\.tif$"
# #              └─ bu qism ixtiyoriy (?:...)?

# # Format C — sodda GEE chunk (ET_2025_03-...):
# r"ET_(?P<year>\d{4})_(?P<mm>\d{2})"
# r"(?:-\d{10}-\d{10})?\.tif$"

_METRIC_RE = re.compile(
    r"METRIC_monthly_ET_(?P<month>\d{4}-\d{2})_(?P<tile>P\d{3}_R\d{3})"
    r"(?:-\d{10}-\d{10})?\.tif$",
    re.IGNORECASE,
)

# Format C: ET_2025_03[-\d{10}-\d{10}].tif  (path/row yo'q)
_SIMPLE_RE = re.compile(
    r"ET_(?P<year>\d{4})_(?P<mm>\d{2})"
    r"(?:-\d{10}-\d{10})?\.tif$",
    re.IGNORECASE,
)

_SEBAL_RE = re.compile(
    r"SEBAL_monthly_(?P<product>[a-zA-Z]+)_(?P<month>\d{4}-\d{2})_(?P<tile>P\d{3}_R\d{2,3})"
    r"(?:-\d{10}-\d{10})?\.tif$", 
    re.IGNORECASE
)

def parse_filename(fname: str) -> tuple[str | None, str | None]:
    """
    Fayl nomidan (month, tile) qaytaradi.
    tile = None — agar fayl nomida path/row yo'q bo'lsa (format C).
    """
    m = _METRIC_RE.match(fname)
    if m:
        return m.group("month"), m.group("tile").upper()

    m = _SIMPLE_RE.match(fname)
    if m:
        return f"{m.group('year')}-{m.group('mm')}", None

    m = _SEBAL_RE.match(fname)
    if m:
        return m.group("month"), m.group("tile").upper()

    return None, None


# ================================================================
# MEAN MOZAYKA (GDAL)
# ================================================================

def mosaic_mean(files: list[str], output_path: str, nodata: float = -9999.0):
    """
    Ko'p tiledan mean mozayka.
    GEE chunklar + turli tilalar barchasini birga qabul qiladi.
    Overlap → o'rtacha | NoData → boshqa tiledagi qiymat.
    """
    datasets = [gdal.Open(f) for f in files]

    # ── Umumiy extent ──────────────────────────────────────────────
    min_x, max_y = float('inf'),  float('-inf')
    max_x, min_y = float('-inf'), float('inf')

    for ds in datasets:
        gt = ds.GetGeoTransform()
        w, h = ds.RasterXSize, ds.RasterYSize
        x0, y0 = gt[0], gt[3]
        x1 = x0 + w * gt[1]
        y1 = y0 + h * gt[5]
        min_x = min(min_x, x0, x1)
        max_x = max(max_x, x0, x1)
        min_y = min(min_y, y0, y1)
        max_y = max(max_y, y0, y1)

    gt0 = datasets[0].GetGeoTransform()
    px, py = gt0[1], gt0[5]   # px > 0, py < 0

    cols    = int(np.ceil((max_x - min_x) / px))
    rows    = int(np.ceil((max_y - min_y) / abs(py)))
    n_bands = datasets[0].RasterCount

    print(f"     Extent: {cols}×{rows} px | {n_bands} band | {len(files)} fayl")

    # ── Sum + count ────────────────────────────────────────────────
    sum_arr   = np.zeros((n_bands, rows, cols), dtype=np.float64)
    count_arr = np.zeros((n_bands, rows, cols), dtype=np.int16)

    for i, ds in enumerate(datasets):
        gt      = ds.GetGeoTransform()
        col_off = int(round((gt[0] - min_x) /  px))
        row_off = int(round((gt[3] - max_y) /  py))

        for b in range(n_bands):
            arr   = ds.GetRasterBand(b + 1).ReadAsArray().astype(np.float64)
            valid = (~np.isnan(arr)) & (arr > nodata + 1)

            h_t, w_t = arr.shape
            r0 = max(0, row_off);       c0 = max(0, col_off)
            r1 = min(rows, row_off + h_t)
            c1 = min(cols, col_off + w_t)
            tr0 = r0 - row_off;         tc0 = c0 - col_off
            tr1 = tr0 + (r1 - r0);      tc1 = tc0 + (c1 - c0)

            chunk_data  = arr[tr0:tr1, tc0:tc1]
            chunk_valid = valid[tr0:tr1, tc0:tc1]

            sum_arr  [b, r0:r1, c0:c1] += np.where(chunk_valid, chunk_data, 0)
            count_arr[b, r0:r1, c0:c1] += chunk_valid.astype(np.int16)

        print(f"     [{i+1}/{len(files)}] {os.path.basename(files[i])}")

    # ── O'rtacha ───────────────────────────────────────────────────
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.where(count_arr > 0, sum_arr / count_arr, np.nan)

    # ── GeoTIFF yozish ─────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(
        output_path, cols, rows, n_bands, gdal.GDT_Float32,
        options=['COMPRESS=LZW', 'TILED=YES',
                 'BLOCKXSIZE=512', 'BLOCKYSIZE=512', 'BIGTIFF=YES'],
    )
    out_ds.SetGeoTransform((min_x, px, 0, max_y, 0, py))
    out_ds.SetProjection(datasets[0].GetProjection())

    for b in range(n_bands):
        data = np.where(np.isnan(result[b]), nodata, result[b]).astype(np.float32)
        out_b = out_ds.GetRasterBand(b + 1)
        out_b.WriteArray(data)
        out_b.SetNoDataValue(nodata)

    out_ds.FlushCache()
    out_ds = None
    for ds in datasets:
        ds = None

    valid_vals = result[~np.isnan(result)]
    if valid_vals.size:
        return valid_vals.min(), valid_vals.max(), valid_vals.mean(), cols, rows
    return None, None, None, cols, rows


# ================================================================
# FAYLLARNI GURUHLASH
# ================================================================

def scan_input(input_dir: str,
               tiles:  list[str] | None = None,
               months: list[str] | None = None) -> dict[str, list[str]]:
    """
    Papkadagi barcha tif fayllarni oylar bo'yicha guruhlaydi.
    Bir oy uchun barcha chunklar + barcha tilalar bitta ro'yxatga tushadi.

    Returns:
        {"2025-04": [file1, file2, ...], "2025-05": [...], ...}
    """
    groups: dict[str, list[str]] = defaultdict(list)

    for fname in sorted(os.listdir(input_dir)):
        if not fname.lower().endswith('.tif'):
            continue

        month, tile = parse_filename(fname)
        if month is None:
            continue

        # Filtr: tile (faqat tile nomi bor fayllar uchun)
        if tiles and tile is not None:
            if tile not in [t.upper() for t in tiles]:
                continue

        # Filtr: month
        if months and month not in months:
            continue

        groups[month].append(os.path.join(input_dir, fname))

    return dict(sorted(groups.items()))


# ================================================================
# ASOSIY RUN
# ================================================================

def run(input_dir: str, output_dir: str,
        tiles:  list[str] | None = None,
        months: list[str] | None = None,
        nodata: float = -9999.0):

    print(f"\nInput:  {input_dir}")
    print(f"Output: {output_dir}")
    if tiles:  print(f"Tillar: {tiles}")
    if months: print(f"Oylar:  {months}")

    groups = scan_input(input_dir, tiles=tiles, months=months)

    if not groups:
        sys.exit(
            "Hech qanday .tif fayl topilmadi.\n"
            "Kutilgan fayl formatlari:\n"
            "  METRIC_monthly_ET_YYYY-MM_PXXX_RXXX.tif\n"
            "  METRIC_monthly_ET_YYYY-MM_PXXX_RXXX-⁠\d{10}-\d{10}.tif\n"
            "  ET_YYYY_MM-\d{10}-\d{10}.tif"
        )

    print(f"\n{len(groups)} ta oy topildi:\n")

    ok = err = skip = 0

    for month, files in groups.items():
        out_name = f"METRIC_monthly_ET_{month}_mosaic.tif"
        out_path = os.path.join(output_dir, out_name)

        # Fayl nomlarini chiqarish
        fnames = [os.path.basename(f) for f in files]
        print(f"  {month}  ({len(files)} fayl)")
        for fn in fnames:
            print(f"    ← {fn}")

        if os.path.exists(out_path):
            print(f"    ⏩ mavjud — o'tkazildi  ({out_name})\n")
            skip += 1
            continue

        # 1 ta fayl → oddiy nusxa
        if len(files) == 1:
            os.makedirs(output_dir, exist_ok=True)
            shutil.copy2(files[0], out_path)
            print(f"    ✅ nusxa ko'chirildi → {out_name}\n")
            ok += 1
            continue

        try:
            mn, mx, mean, w, h = mosaic_mean(files, out_path, nodata=nodata)
            if mean is not None:
                print(f"    ✅ {out_name}")
                print(f"       {w}×{h} px  |  ET: {mn:.1f}–{mx:.1f}, o'rta {mean:.1f} mm\n")
            else:
                print(f"    ⚠️  {out_name} — barcha piksellar NoData\n")
            ok += 1
        except Exception as ex:
            print(f"    ❌ XATO: {ex}\n")
            import traceback; traceback.print_exc()
            err += 1

    print("─" * 45)
    print(f"Jami: {len(groups)}  |  ✅ {ok}  |  ⏩ {skip}  |  ❌ {err}")
    print(f"Output: {output_dir}\n")


# ================================================================
# CLI
# ================================================================

def main():
    p = argparse.ArgumentParser(
        description="METRIC GeoTIFF mean mozayka (GEE chunk + multi-tile qo'llab-quvvatlaydi)",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--input",  "-i", required=True,
                   help="GeoTIFF fayllar papkasi")
    p.add_argument("--output", "-o", required=True,
                   help="Mozayka chiqish papkasi")
    p.add_argument("--tiles",  "-t", nargs="+", metavar="PxRy",
                   help="Faqat shu tile(lar): P040_R030 P041_R030")
    p.add_argument("--months", "-m", nargs="+", metavar="YYYY-MM",
                   help="Faqat shu oy(lar): 2025-07 2025-08")
    p.add_argument("--nodata", type=float, default=-9999.0,
                   help="NoData qiymati (default: -9999)")
    args = p.parse_args()

    run(
        input_dir  = args.input,
        output_dir = args.output,
        tiles      = args.tiles,
        months     = args.months,
        nodata     = args.nodata,
    )


if __name__ == "__main__":
    main()
