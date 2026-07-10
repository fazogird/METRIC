"""
ERA5-Land soatlik meteorologik ma'lumotlar.

METRIC aslida meteostansiya ishlatadi, biz ERA5 bilan almashtiramiz.
Bu global ishlash imkonini beradi — stansiya kerak emas.

ERA5-Land GEE da: ECMWF/ERA5_LAND/HOURLY
Band nomlari:
  temperature_2m        → T2m (K)
  dewpoint_temperature_2m → Td (K)  → ea hisoblash uchun
  u_component_of_wind_10m → u10 (m/s)
  v_component_of_wind_10m → v10 (m/s)
  surface_pressure       → sp (Pa)
  total_precipitation    → tp (m)  → mm ga o'tkazish kerak
  surface_solar_radiation_downwards → ssrd (J/m²)
  surface_net_solar_radiation → ssr (J/m²)
"""
import ee


class ERA5Loader:
    """ERA5-Land soatlik ma'lumotlar yuklagich."""
    
    COLLECTION_ID = "ECMWF/ERA5_LAND/HOURLY"
    
    def __init__(self, settings):
        self.settings = settings
    
    def get_hourly(self, geometry, datetime_str):
        """
        Landsat image vaqtiga eng yaqin ERA5 soatlik ma'lumot.
        
        Args:
            geometry: ee.Geometry
            datetime_str: str — 'YYYY-MM-DDTHH:MM:SS' yoki ee.Date
            
        Returns:
            ee.Image — T2m, Td, u10, v10, sp, tp, ssrd bands
        """
        date = ee.Date(datetime_str)
        
        # ±1 soat oraliqda filtr
        era5 = (ee.ImageCollection(self.COLLECTION_ID)
            .filterBounds(geometry)
            .filterDate(date.advance(-1, 'hour'), date.advance(1, 'hour'))
            .first())
        
        return self._process(era5)
    
    def get_daily(self, geometry, date_str):
        """Kunlik yig'indi/o'rtacha — ET24 va ETr24 uchun."""
        date = ee.Date(date_str)
        
        day_coll = (ee.ImageCollection(self.COLLECTION_ID)
            .filterBounds(geometry)
            .filterDate(date, date.advance(1, 'day')))
        
        # Temperatura va shamol → kunlik o'rtacha
        t2m_mean = day_coll.select('temperature_2m').mean()
        td_mean = day_coll.select('dewpoint_temperature_2m').mean()
        u10_mean = day_coll.select('u_component_of_wind_10m').mean()
        v10_mean = day_coll.select('v_component_of_wind_10m').mean()
        sp_mean = day_coll.select('surface_pressure').mean()
        
        # Yog'in → kunlik summa (m → mm)
        tp_sum = day_coll.select('total_precipitation').sum().multiply(1000)
        
        # Radiatsiya → kunlik summa (J/m² → W/m² o'rtacha uchun /86400)
        ssrd_sum = day_coll.select('surface_solar_radiation_downwards').sum()
        
        return (t2m_mean
            .addBands(td_mean)
            .addBands(u10_mean)
            .addBands(v10_mean)
            .addBands(sp_mean)
            .addBands(tp_sum.rename('tp_daily_mm'))
            .addBands(ssrd_sum.rename('ssrd_daily')))
    
    def get_precip_sum(self, geometry, end_date, days=60):
        """
        Oxirgi N kunlik yog'in summasi — CIMEC Tfac (F.8) uchun.

        Args:
            days: int — default 60 kun

        Returns:
            ee.Image — 'precip_sum_mm'
        """
        end   = ee.Date(end_date)
        start = end.advance(-days, 'day')

        tp = (ee.ImageCollection(self.COLLECTION_ID)
              .filterBounds(geometry)
              .filterDate(start, end)
              .select('total_precipitation')
              .sum()
              .multiply(1000))   # m → mm

        return tp.rename('precip_sum_mm')

    def get_etr_sum(self, geometry, end_date, days=60):
        """
        Oxirgi N kunlik ETr summasi — CIMEC Tfac (F.8) uchun.

        Hargreaves-Samani yaqinlashtirilgan formula:
            ETr_day = 0.0023 × sqrt(Tmax-Tmin) × (Tmean+17.8) × Ra
        bu erda Ra = ssrd / τ_sw  (MJ/m²/day),  τ_sw ≈ 0.75

        Args:
            days: int — default 60 kun

        Returns:
            ee.Image — 'etr_sum_mm'
        """
        end   = ee.Date(end_date)
        start = end.advance(-days, 'day')

        daily_coll = (ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
                      .filterBounds(geometry)
                      .filterDate(start, end))

        def daily_etr(img):
            t_max  = img.select('temperature_2m_max').subtract(273.15)
            t_min  = img.select('temperature_2m_min').subtract(273.15)
            t_mean = t_max.add(t_min).divide(2)
            t_rng  = t_max.subtract(t_min).max(0)

            # Ra: surface ssrd / τ_sw yaqinlashtirilgan (J/m²/day → MJ/m²/day)
            ssrd = img.select('surface_solar_radiation_downwards_sum').divide(1e6)
            Ra   = ssrd.divide(0.75).max(0)

            # Hargreaves-Samani: ETr = 0.0023 × sqrt(ΔT) × (Tmean+17.8) × Ra
            etr = t_rng.sqrt().multiply(t_mean.add(17.8)).multiply(Ra).multiply(0.0023)
            return etr.max(0).rename('ETr_daily')

        return daily_coll.map(daily_etr).sum().rename('etr_sum_mm')
    
    def _process(self, image):
        """ERA5 bandlarni qayta ishlash va nomlash."""
        # Wind speed (m/s): ws = sqrt(u10² + v10²)
        u10 = image.select('u_component_of_wind_10m')
        v10 = image.select('v_component_of_wind_10m')
        ws10 = u10.pow(2).add(v10.pow(2)).sqrt().rename('wind_speed_10m')
        
        # Vapor pressure (kPa): ea = 0.6108 * exp(17.27*Td / (Td+237.3))
        # Td Kelvin da keladi, Celsius ga o'tkazamiz
        td_c = image.select('dewpoint_temperature_2m').subtract(273.15)
        ea = td_c.multiply(17.27).divide(td_c.add(237.3)).exp().multiply(0.6108).rename('ea_kpa')
        
        # Surface pressure Pa → kPa
        sp_kpa = image.select('surface_pressure').divide(1000).rename('sp_kpa')
        
        # Precipitation m → mm
        tp_mm = image.select('total_precipitation').multiply(1000).rename('tp_mm')
        
        return (image
            .addBands(ws10)
            .addBands(ea)
            .addBands(sp_kpa)
            .addBands(tp_mm))


def prepare_era5_for_pm_hourly(era5_img):
    """ERA5 soatlik -> RefETCalculator uchun band tayyorlash."""
    T_air = era5_img.select('temperature_2m').subtract(273.15).rename('T_air')
    P_kPa = era5_img.select('sp_kpa').rename('P_kPa')
    u2 = era5_img.select('wind_speed_10m').multiply(0.748).rename('u2')
    ea = era5_img.select('ea_kpa').rename('ea')
    Rs = era5_img.select('surface_solar_radiation_downwards').divide(1e6).max(0).rename('Rs')
    return (T_air.addBands(P_kPa).addBands(u2).addBands(ea).addBands(Rs)
            .copyProperties(era5_img, era5_img.propertyNames()))


# def prepare_era5_for_pm_daily(geometry, date_str):
#     """ERA5 kunlik -> RefETCalculator daily mode uchun."""
#     date = ee.Date(date_str)
#     coll = (ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")
#             .filterBounds(geometry)
#             .filterDate(date, date.advance(1, 'day')))
#     t2m = coll.select('temperature_2m')
#     T_max = t2m.max().subtract(273.15).rename('T_max')
#     T_min = t2m.min().subtract(273.15).rename('T_min')
#     T_air = t2m.mean().subtract(273.15).rename('T_air')
#     P_kPa = coll.select('surface_pressure').mean().divide(1000).rename('P_kPa')
#     u10 = coll.select('u_component_of_wind_10m').mean()
#     v10 = coll.select('v_component_of_wind_10m').mean()
#     u2 = u10.pow(2).add(v10.pow(2)).sqrt().multiply(0.748).rename('u2')
#     td = coll.select('dewpoint_temperature_2m').mean().subtract(273.15)
#     ea = td.multiply(17.27).divide(td.add(237.3)).exp().multiply(0.6108).rename('ea')
#     # Rs = coll.select('surface_solar_radiation_downwards').max().divide(1e6).max(0).rename('Rs')
    
#     # # ✅ YANGI — ikkala 12-soatlik davr yig'indisi:
#     # ssrd = coll.select('surface_solar_radiation_downwards')

#     # # 1-davr: 00:00-11:59 UTC
#     # period1 = ssrd.filterDate(date, date.advance(12, 'hour')).max()
#     # # 2-davr: 12:00-23:59 UTC  
#     # period2 = ssrd.filterDate(date.advance(12, 'hour'), date.advance(1, 'day')).max()

#     # Rs = period1.add(period2).divide(1e6).max(0).rename('Rs')
    
#     # Eng ishonchli:
#     Rs = coll.select('surface_solar_radiation_downwards').sum().divide(1e6).max(0).rename('Rs')
#     # sum() to'g'ri AGAR GEE dagi ERA5-Land soatlik de-akkumulyativ bo'lsa
#     # GEE docs: "accumulated over the 1 hour ending at validity time" → sum() to'g'ri
    
#     result = (T_air.addBands(T_max).addBands(T_min)
#               .addBands(P_kPa).addBands(u2).addBands(ea).addBands(Rs))
#     return result.set('system:time_start', date.millis())

def prepare_era5_for_pm_daily(geometry, date_str):
    date = ee.Date(date_str)
    
    # Kunlik aggregat — ssrd allaqachon summa
    daily = (ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
             .filterBounds(geometry)
             .filterDate(date, date.advance(1, 'day'))
             .first())
    
    T_max = daily.select('temperature_2m_max').subtract(273.15).rename('T_max')
    T_min = daily.select('temperature_2m_min').subtract(273.15).rename('T_min')
    T_air = T_max.add(T_min).divide(2).rename('T_air')
    
    P_kPa = daily.select('surface_pressure').divide(1000).rename('P_kPa')
    
    u = daily.select('u_component_of_wind_10m')
    v = daily.select('v_component_of_wind_10m')
    u2 = u.pow(2).add(v.pow(2)).sqrt().multiply(0.748).rename('u2')
    
    td = daily.select('dewpoint_temperature_2m').subtract(273.15)
    ea = td.multiply(17.27).divide(td.add(237.3)).exp().multiply(0.6108).rename('ea')
    
    # ✅ Allaqachon kunlik summa (J/m²) → MJ/m²/day
    Rs = daily.select('surface_solar_radiation_downwards_sum').divide(1e6).max(0).rename('Rs')
    
    result = (T_air.addBands(T_max).addBands(T_min)
              .addBands(P_kPa).addBands(u2).addBands(ea).addBands(Rs))
    return result.set('system:time_start', date.millis())
