"""
OpenET → Tile → METRIC: teskari yondashuv.

1. OpenET Idaho uchun mavjud oylarni topadi (2025 yoki oxirgi)
2. Shu oylar uchun eng yaxshi Landsat tile tanlaydi
3. METRIC pipeline ishlatadi
4. OpenET bilan taqqoslaydi
"""
import ee
import calendar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pipeline import METRICPipeline
from metric_gee.config.settings import Settings

ee.Initialize()

# ══════════════════════════════════════════════════════════════════
# 1. OPENET — MAVJUD OYLARNI TOP
# ══════════════════════════════════════════════════════════════════
OPENET_ENSEMBLE = 'projects/openet/assets/ensemble/conus/gridmet/monthly/v2_1'
OPENET_MODELS = {
    'Ensemble': ('projects/openet/assets/ensemble/conus/gridmet/monthly/v2_1', 'et_ensemble_mad'),
    'eeMETRIC': ('projects/openet/assets/eemetric/conus/gridmet/monthly/v2_1', 'et'),
    'disALEXI': ('projects/openet/assets/disalexi/conus/gridmet/monthly/v2_1', 'et'),
    'geeSEBAL': ('projects/openet/assets/geesebal/conus/gridmet/monthly/v2_1', 'et'),
    'PT-JPL':   ('projects/openet/assets/ptjpl/conus/gridmet/monthly/v2_1',    'et'),
    'SIMS':     ('projects/openet/assets/sims/conus/gridmet/monthly/v2_1',     'et'),
    'SSEBop':   ('projects/openet/assets/ssebop/conus/gridmet/monthly/v2_1',   'et'),
}

TARGET_YEAR = 2025   # shu yilni sinab ko'ramiz
MONTHS      = [3, 4, 5, 6, 7, 8, 9, 10]   # mart-oktyabr
CLOUD_MAX    = 50
# reduceRegion uchun katta scale — xotirani tejaydi, aniqlik yetarli
SAMPLE_SCALE = 2000   # 500 → 2000 (16× kam piksel, GEE memory error yo'q)

# Idaho bbox
IDAHO_AOI = ee.Geometry.Rectangle([-113.10, 42.50, -112.80, 42.70])

print('='*60)
print(f'  OpenET {TARGET_YEAR} — Idaho uchun mavjud oylar')
print('='*60)

available_months = []
for mo in MONTHS:
    start = f'{TARGET_YEAR}-{mo:02d}-01'
    end   = f'{TARGET_YEAR}-{mo:02d}-{calendar.monthrange(TARGET_YEAR, mo)[1]:02d}'
    n = (ee.ImageCollection(OPENET_ENSEMBLE)
           .filterDate(start, end)
           .filterBounds(IDAHO_AOI)
           .size().getInfo())
    status = '✅' if n > 0 else '❌'
    print(f'  {TARGET_YEAR}-{mo:02d}: {n} obraz {status}')
    if n > 0:
        available_months.append(f'{TARGET_YEAR}-{mo:02d}')

if not available_months:
    print(f'\n  {TARGET_YEAR} uchun OpenET yo\'q — oxirgi mavjud yil qidirilmoqda...')
    for yr in range(TARGET_YEAR - 1, 1998, -1):
        n = (ee.ImageCollection(OPENET_ENSEMBLE)
               .filterDate(f'{yr}-01-01', f'{yr+1}-01-01')
               .filterBounds(IDAHO_AOI)
               .size().getInfo())
        if n > 0:
            for mo in MONTHS:
                start = f'{yr}-{mo:02d}-01'
                end   = f'{yr}-{mo:02d}-{calendar.monthrange(yr, mo)[1]:02d}'
                n2 = (ee.ImageCollection(OPENET_ENSEMBLE)
                        .filterDate(start, end)
                        .filterBounds(IDAHO_AOI)
                        .size().getInfo())
                if n2 > 0:
                    available_months.append(f'{yr}-{mo:02d}')
            break
    print(f'  Topildi: {available_months[0][:4]} yil, {len(available_months)} oy')

print(f'\n  Taqqoslash oylari: {available_months}')
USE_YEAR = int(available_months[0][:4])

# ══════════════════════════════════════════════════════════════════
# 2. ENG YAXSHI LANDSAT TILE TANLASH
#    Idaho uchun eng ko'p toza tasvirga ega WRS-2 tileni topamiz
# ══════════════════════════════════════════════════════════════════
print('\n' + '='*60)
print(f'  Idaho uchun eng yaxshi Landsat tile qidirilmoqda ({USE_YEAR})')
print('='*60)

# OpenET qamrovini biladigan bitta obraz olish (ixtiyoriy oy)
ref_month = available_months[len(available_months)//2]   # o'rtadagi oy
rm_year, rm_mo = int(ref_month[:4]), int(ref_month[5:])
openet_ref = (ee.ImageCollection(OPENET_ENSEMBLE)
              .filterDate(f'{rm_year}-{rm_mo:02d}-01',
                          f'{rm_year}-{rm_mo:02d}-{calendar.monthrange(rm_year, rm_mo)[1]:02d}')
              .filterBounds(IDAHO_AOI)
              .first())

# Barcha mavjud oylar uchun Landsat L8/L9 tasvirlarini sanash
start_all = f'{USE_YEAR}-{MONTHS[0]:02d}-01'
end_all   = f'{USE_YEAR}-{MONTHS[-1]:02d}-{calendar.monthrange(USE_YEAR, MONTHS[-1])[1]:02d}'

l8 = (ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
        .filterBounds(IDAHO_AOI)
        .filterDate(start_all, end_all)
        .filter(ee.Filter.lt('CLOUD_COVER', CLOUD_MAX)))
l9 = (ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
        .filterBounds(IDAHO_AOI)
        .filterDate(start_all, end_all)
        .filter(ee.Filter.lt('CLOUD_COVER', CLOUD_MAX)))
merged = l8.merge(l9)
count_total = merged.size().getInfo()
print(f'  Jami Landsat tasvirlar (Idaho, {USE_YEAR}): {count_total}')

img_list = merged.toList(count_total)

# PATH/ROW bo'yicha sanash
from collections import Counter
tile_counts = Counter()
for i in range(count_total):
    img = ee.Image(img_list.get(i))
    path = int(ee.Number(img.get('WRS_PATH')).getInfo())
    row  = int(ee.Number(img.get('WRS_ROW')).getInfo())
    tile_counts[(path, row)] += 1

print(f'\n  Top 5 tile (eng ko\'p tasvirli):')
for (path, row), cnt in tile_counts.most_common(5):
    print(f'    P{path:03d}/R{row:03d}: {cnt} tasvir')

best_path, best_row = tile_counts.most_common(1)[0][0]
print(f'\n  Tanlangan tile: P{best_path:03d}/R{best_row:03d}  '
      f'({tile_counts[(best_path, best_row)]} tasvir)')

# ══════════════════════════════════════════════════════════════════
# 3. TILE GEOMETRY — OPENET VA METRIC UCHUN UMUMIY AOI
# ══════════════════════════════════════════════════════════════════
tile_geometry = (
    ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filter(ee.Filter.eq('WRS_PATH', best_path))
    .filter(ee.Filter.eq('WRS_ROW',  best_row))
    .first().geometry()
)

# CDL cropland mask
cdl_year = USE_YEAR
cdl_coll = ee.ImageCollection('USDA/NASS/CDL').filter(
    ee.Filter.calendarRange(cdl_year, cdl_year, 'year'))
if cdl_coll.size().getInfo() > 0:
    crop_mask = cdl_coll.first().select('cropland').lte(60)
    frac = (crop_mask.reduceRegion(
        ee.Reducer.mean(), tile_geometry, 500,
        maxPixels=1e9, bestEffort=True).getInfo().get('cropland', 0))
    print(f'  CDL {cdl_year}: ekin {frac*100:.1f}%')
else:
    crop_mask = None
    print(f'  CDL {cdl_year} topilmadi')

TIMEZONE_LON = -105.0   # Idaho UTC-7

# ══════════════════════════════════════════════════════════════════
# 4. METRIC + OPENET — OY BO'YICHA HISOBLASH
# ══════════════════════════════════════════════════════════════════
def run_metric(month_str):
    year, mo = int(month_str[:4]), int(month_str[5:])
    last = calendar.monthrange(year, mo)[1]
    start, end = f'{year}-{mo:02d}-01', f'{year}-{mo:02d}-{last:02d}'
    cfg = Settings(
        cloud_cover_max=CLOUD_MAX, roi_cloud_max=20,
        albedo_method='olmedo', g_method='bastiaanssen',
        zom_method='lai', lapse_rate=6.5, etrf_bare=0.05,
        timezone_lon_deg=TIMEZONE_LON,
        wrs_path=best_path, wrs_row=best_row,
    )
    try:
        res = METRICPipeline(settings=cfg).run_monthly(tile_geometry, start, end)
    except Exception as e:
        print(f'    Pipeline xato: {e}')
        return None, None
    if 'error' in res:
        return None, None
    monthly = res.get('monthly_et', {})
    et_img  = monthly.get('et_monthly_image')
    et_all  = monthly.get('et_monthly_mean')
    if et_img is None:
        return et_all, None
    img = et_img.updateMask(crop_mask) if crop_mask else et_img
    try:
        val = (img.reduceRegion(
            reducer   = ee.Reducer.mean(),
            geometry  = tile_geometry,
            scale     = SAMPLE_SCALE,
            maxPixels = 1e9,
            tileScale = 8,       # ← muhim: hisoblashni 8 qismga bo'ladi
            bestEffort= True,
        ).getInfo().get('ET_monthly'))
        return (float(val) if val else et_all), et_img
    except Exception as e:
        print(f'    reduceRegion xato ({e.__class__.__name__}), et_all ishlatildi')
        return et_all, et_img


def get_openet(coll_id, band, month_str):
    year, mo = int(month_str[:4]), int(month_str[5:])
    last = calendar.monthrange(year, mo)[1]
    start, end = f'{year}-{mo:02d}-01', f'{year}-{mo:02d}-{last:02d}'
    coll = (ee.ImageCollection(coll_id)
              .filterDate(start, end)
              .filterBounds(tile_geometry))
    if coll.size().getInfo() == 0:
        return None
    img = coll.first().select(band)
    if crop_mask:
        img = img.updateMask(crop_mask)
    val = (img.reduceRegion(
        reducer   = ee.Reducer.mean(),
        geometry  = tile_geometry,
        scale     = SAMPLE_SCALE,
        maxPixels = 1e9,
        tileScale = 4,
        bestEffort= True,
    ).getInfo().get(band))
    return float(val) if val is not None else None


print('\n' + '='*60)
print(f'  Hisoblash: P{best_path:03d}/R{best_row:03d} — {len(available_months)} oy')
print('='*60)

def get_raster_stats(metric_img, coll_id, band, month_str):
    """
    METRIC va OpenET raster o'rtasida piksel bo'yicha taqqoslash.
    Qaytadi: {'r2', 'rmse', 'bias', 'mae'} yoki None.
    """
    year, mo = int(month_str[:4]), int(month_str[5:])
    last = calendar.monthrange(year, mo)[1]
    start, end = f'{year}-{mo:02d}-01', f'{year}-{mo:02d}-{last:02d}'
    coll = (ee.ImageCollection(coll_id)
              .filterDate(start, end)
              .filterBounds(tile_geometry))
    if coll.size().getInfo() == 0:
        return None
    oe_img = coll.first().select(band).rename('oe')
    if crop_mask:
        oe_img = oe_img.updateMask(crop_mask)
    m_img = metric_img.rename('metric')
    if crop_mask:
        m_img = m_img.updateMask(crop_mask)

    # Pearson r (R²)
    try:
        corr = m_img.addBands(oe_img).reduceRegion(
            reducer=ee.Reducer.pearsonsCorrelation(),
            geometry=tile_geometry,
            scale=SAMPLE_SCALE, maxPixels=1e9, tileScale=4, bestEffort=True,
        ).getInfo()
        r_val = corr.get('correlation') if corr else None
        r2 = (r_val ** 2) if r_val is not None else None
    except Exception:
        r2 = None

    # bias, RMSE, MAE
    try:
        diff = m_img.subtract(oe_img)
        stats = (diff.rename('bias')
                 .addBands(diff.pow(2).rename('mse'))
                 .addBands(diff.abs().rename('mae'))
                 .reduceRegion(
                     reducer=ee.Reducer.mean(),
                     geometry=tile_geometry,
                     scale=SAMPLE_SCALE, maxPixels=1e9, tileScale=4, bestEffort=True,
                 ).getInfo()) or {}
        bias = stats.get('bias')
        rmse = (np.sqrt(stats['mse']) if stats.get('mse') is not None else None)
        mae  = stats.get('mae')
    except Exception:
        bias = rmse = mae = None

    return {'r2': r2, 'rmse': rmse, 'bias': bias, 'mae': mae}


records = []
for month_str in available_months:
    print(f'\n  {month_str}')
    metric_val, metric_img = run_metric(month_str)
    print(f'    METRIC: {f"{metric_val:.1f} mm" if metric_val else "xato"}')
    row = {'month': month_str, 'METRIC': metric_val, 'metric_img': metric_img}
    for name, (cid, band) in OPENET_MODELS.items():
        val = get_openet(cid, band, month_str)
        row[name] = val
        print(f'    {name:10s}: {f"{val:.1f} mm" if val else "yoq"}')
    records.append(row)

df_all = pd.DataFrame(records)
df = df_all.drop(columns=['metric_img'], errors='ignore').set_index('month')
print('\n\nNATIJA (skalyar oy o\'rtachalari):')
print(df.round(1).to_string())

# ── Raster taqqoslash statistikasi ───────────────────────────────────
print('\n' + '='*60)
print('  RASTER STATISTIKA  (piksel bo\'yicha, METRIC vs OpenET)')
print('  ' + '-'*56)
print(f'  {"Model":12s}  {"R²":>6}  {"RMSE":>7}  {"Bias":>8}  {"MAE":>7}')
print('  ' + '-'*56)

raster_agg = {}   # {model_name: {'r2':[], 'rmse':[], 'bias':[], 'mae':[]}}
for name in OPENET_MODELS:
    raster_agg[name] = {'r2': [], 'rmse': [], 'bias': [], 'mae': []}

for rec in records:
    m_img = rec.get('metric_img')
    if m_img is None:
        continue
    for name, (cid, band) in OPENET_MODELS.items():
        try:
            st = get_raster_stats(m_img, cid, band, rec['month'])
            if st is None:
                continue
            for k in ('r2', 'rmse', 'bias', 'mae'):
                if st.get(k) is not None:
                    raster_agg[name][k].append(st[k])
        except Exception as e:
            print(f'    ⚠️ {name} {rec["month"]}: {e}')

raster_summary = {}
for name, vals in raster_agg.items():
    r2_m   = float(np.mean(vals['r2']))   if vals['r2']   else np.nan
    rmse_m = float(np.mean(vals['rmse'])) if vals['rmse'] else np.nan
    bias_m = float(np.mean(vals['bias'])) if vals['bias'] else np.nan
    mae_m  = float(np.mean(vals['mae']))  if vals['mae']  else np.nan
    raster_summary[name] = {'r2': r2_m, 'rmse': rmse_m, 'bias': bias_m, 'mae': mae_m}
    r2_s   = f'{r2_m:.3f}'   if not np.isnan(r2_m)   else '—'
    rmse_s = f'{rmse_m:.1f}' if not np.isnan(rmse_m) else '—'
    bias_s = f'{bias_m:+.1f}'if not np.isnan(bias_m) else '—'
    mae_s  = f'{mae_m:.1f}'  if not np.isnan(mae_m)  else '—'
    print(f'  {name:12s}  {r2_s:>6}  {rmse_s:>7}  {bias_s:>8}  {mae_s:>7} mm')
print('  ' + '-'*56)

# ══════════════════════════════════════════════════════════════════
# 5. GRAFIK  (3 qism: bar chart | scatter Ensemble | per-model stats)
# ══════════════════════════════════════════════════════════════════
from scipy import stats as sp_stats

MODEL_COLORS = {
    'Ensemble': '#2c3e50', 'eeMETRIC': '#27ae60', 'disALEXI': '#2980b9',
    'geeSEBAL': '#8e44ad', 'PT-JPL':   '#e67e22', 'SIMS':     '#c0392b',
    'SSEBop':   '#16a085', 'METRIC':   '#f39c12',
}
months_u = df.index.tolist()
cols     = list(OPENET_MODELS.keys()) + ['METRIC']

# ── Panel 1: Bar chart ────────────────────────────────────────────
fig = plt.figure(figsize=(22, 14))
gs  = fig.add_gridspec(3, 4, hspace=0.45, wspace=0.35)

ax_bar = fig.add_subplot(gs[0, :])
n_cols = len(cols)
bar_w  = 0.8 / n_cols
x      = np.arange(len(months_u))

for i, col in enumerate(cols):
    vals  = [df.loc[m, col] if m in df.index else np.nan for m in months_u]
    color = MODEL_COLORS.get(col, 'gray')
    kw    = dict(edgecolor='black', linewidth=1.5, zorder=5) if col == 'METRIC' else {}
    ax_bar.bar(x + i*bar_w, vals, bar_w*0.9, label=col, color=color, alpha=0.75, **kw)

ax_bar.set_xticks(x + bar_w*n_cols/2)
ax_bar.set_xticklabels([m[5:] for m in months_u], rotation=30)
ax_bar.set_ylabel('ET (mm/month)')
ax_bar.set_title(f'P{best_path:03d}/R{best_row:03d} — METRIC vs OpenET ({USE_YEAR})',
                 fontsize=12, fontweight='bold')
ax_bar.legend(ncol=4, fontsize=8)
ax_bar.grid(axis='y', alpha=0.3)

# ── Panel 2–8: Per-model scatter (METRIC vs each OpenET model) ───
model_names = list(OPENET_MODELS.keys())
scatter_axes = [fig.add_subplot(gs[1 + i//4, i % 4]) for i in range(len(model_names))]

for ax_sc, name in zip(scatter_axes, model_names):
    oe_s = df[name].dropna() if name in df.columns else pd.Series(dtype=float)
    me_s = df.loc[oe_s.index, 'METRIC'].dropna()
    common = oe_s.index.intersection(me_s.index)
    ax_sc.set_xlabel(f'OpenET {name} (mm)', fontsize=8)
    ax_sc.set_ylabel('METRIC (mm)', fontsize=8)
    ax_sc.tick_params(labelsize=7)
    ax_sc.set_title(name, fontsize=9)
    ax_sc.grid(alpha=0.3)
    if len(common) < 2:
        ax_sc.text(0.5, 0.5, 'yetarli ma\'lumot yo\'q',
                   ha='center', va='center', transform=ax_sc.transAxes, fontsize=8)
        continue
    x_v = oe_s[common].values.astype(float)
    y_v = me_s[common].values.astype(float)
    ax_sc.scatter(x_v, y_v, s=60, zorder=3,
                  c=MODEL_COLORS.get(name, 'gray'), edgecolors='k', lw=0.8)
    for xi, yi, mo in zip(x_v, y_v, common):
        ax_sc.annotate(mo[5:], (xi, yi), xytext=(3, 2),
                       textcoords='offset points', fontsize=6)
    vmin = min(x_v.min(), y_v.min()) * 0.85
    vmax = max(x_v.max(), y_v.max()) * 1.1
    ax_sc.plot([vmin, vmax], [vmin, vmax], 'k--', lw=1, alpha=0.4)
    sl, ic, r, *_ = sp_stats.linregress(x_v, y_v)
    xr = np.linspace(vmin, vmax, 40)
    ax_sc.plot(xr, sl*xr + ic, 'r-', lw=1.2)
    rmse_v = np.sqrt(np.mean((y_v - x_v)**2))
    bias_v = float(np.mean(y_v - x_v))
    # Raster r² (agar hisoblangan bo'lsa)
    rs = raster_summary.get(name, {})
    rr2 = rs.get('r2', np.nan)
    rr2_s = f'R²_rast={rr2:.2f}' if not np.isnan(rr2) else ''
    ax_sc.text(0.04, 0.97,
               f'R²={r**2:.2f}  RMSE={rmse_v:.1f}\nBias={bias_v:+.1f} mm  {rr2_s}',
               transform=ax_sc.transAxes, fontsize=6.5, va='top',
               bbox=dict(boxstyle='round', fc='white', alpha=0.8))
    ax_sc.set_xlim(vmin, vmax); ax_sc.set_ylim(vmin, vmax)

plt.suptitle(f'Idaho P{best_path:03d}/R{best_row:03d} — {USE_YEAR}  (CDL cropland)',
             fontsize=13, fontweight='bold')

out = f'compare_P{best_path:03d}_R{best_row:03d}_{USE_YEAR}.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f'\nGrafik saqlandi: {out}')
plt.show()

# ── Saqlash ──────────────────────────────────────────────────────
try:
    df.round(2).to_excel(f'compare_P{best_path:03d}_R{best_row:03d}_{USE_YEAR}.xlsx')
    print(f'Excel saqlandi: compare_P{best_path:03d}_R{best_row:03d}_{USE_YEAR}.xlsx')
except ImportError:
    out_csv = f'compare_P{best_path:03d}_R{best_row:03d}_{USE_YEAR}.csv'
    df.round(2).to_csv(out_csv)
    print(f'openpyxl topilmadi — CSV saqlandi: {out_csv}')

# Raster stats CSV
rs_df = pd.DataFrame(raster_summary).T.rename_axis('model')
rs_df.round(3).to_csv(f'raster_stats_P{best_path:03d}_R{best_row:03d}_{USE_YEAR}.csv')
print(f'Raster stats saqlandi: raster_stats_P{best_path:03d}_R{best_row:03d}_{USE_YEAR}.csv')
