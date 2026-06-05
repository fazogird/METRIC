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
            ee.Image — precipitation_60 (mm)
        """
        end = ee.Date(end_date)
        start = end.advance(-days, 'day')
        
        tp = (ee.ImageCollection(self.COLLECTION_ID)
            .filterBounds(geometry)
            .filterDate(start, end)
            .select('total_precipitation')
            .sum()
            .multiply(1000))  # m → mm
        
        return tp.rename('precip_sum_mm')
    
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


def prepare_era5_for_pm_daily(geometry, date_str):
    """ERA5 kunlik -> RefETCalculator daily mode uchun."""
    date = ee.Date(date_str)
    coll = (ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")
            .filterBounds(geometry)
            .filterDate(date, date.advance(1, 'day')))
    t2m = coll.select('temperature_2m')
    T_max = t2m.max().subtract(273.15).rename('T_max')
    T_min = t2m.min().subtract(273.15).rename('T_min')
    T_air = t2m.mean().subtract(273.15).rename('T_air')
    P_kPa = coll.select('surface_pressure').mean().divide(1000).rename('P_kPa')
    u10 = coll.select('u_component_of_wind_10m').mean()
    v10 = coll.select('v_component_of_wind_10m').mean()
    u2 = u10.pow(2).add(v10.pow(2)).sqrt().multiply(0.748).rename('u2')
    td = coll.select('dewpoint_temperature_2m').mean().subtract(273.15)
    ea = td.multiply(17.27).divide(td.add(237.3)).exp().multiply(0.6108).rename('ea')
    Rs = coll.select('surface_solar_radiation_downwards').max().divide(1e6).max(0).rename('Rs')
    result = (T_air.addBands(T_max).addBands(T_min)
              .addBands(P_kPa).addBands(u2).addBands(ea).addBands(Rs))
    return result.set('system:time_start', date.millis())
