"""
core/reference_et.py — ASCE Standardized Reference ET (ETr).

To'liq ASCE-EWRI 2005 implementatsiyasi.
Soatlik va kunlik rejimlar. Alfalfa va grass reference.

Kirish bandlari: T_air(°C), P_kPa, u2(m/s), ea(kPa), Rs(MJ/m²/hr yoki /day)
Chiqish: ETr band qo'shiladi.
"""
import ee
import math


class RefETCalculator:
    """ASCE Standardized Penman-Monteith ETr."""

    Gsc = 0.0820       # MJ m⁻² min⁻¹
    albedo = 0.23
    sigma_sb = 4.903e-9  # MJ K⁻⁴ m⁻² day⁻¹

    def __init__(self, timezone_lon_deg=0.0, ref_type='alfalfa'):
        self.Lz = ee.Number(timezone_lon_deg)
        self.ref_type = ref_type
        if ref_type == 'alfalfa':
            self.Cn_hr_day = 66.0
            self.Cn_hr_night = 66.0
            self.Cd_hr_day = 0.25
            self.Cd_hr_night = 1.7
            self.Cn_day = 1600.0
            self.Cd_day = 0.38
        elif ref_type == 'grass':
            self.Cn_hr_day = 37.0
            self.Cn_hr_night = 37.0
            self.Cd_hr_day = 0.24
            self.Cd_hr_night = 0.96
            self.Cn_day = 900.0
            self.Cd_day = 0.34

    # === YORDAMCHI ===

    @staticmethod
    def _esat(T_c):
        return T_c.expression('0.6108 * exp((17.27 * T) / (T + 237.3))', {'T': T_c})

    @staticmethod
    def _delta_svp(T_c, es_kPa):
        return es_kPa.multiply(4098).divide(T_c.add(237.3).pow(2))

    @staticmethod
    def _gamma(P_kPa):
        return P_kPa.multiply(0.000665)

    @staticmethod
    def _doy_hour(img):
        img = ee.Image(img)
        date = ee.Date(img.get('system:time_start'))
        doy = ee.Number(date.getRelative('day', 'year')).add(1)
        hour = ee.Number(date.get('hour'))
        return doy, hour

    @staticmethod
    def _dr_decl(doy):
        J = ee.Number(doy)
        dr = J.multiply(2 * math.pi / 365.0).cos().multiply(0.033).add(1.0)
        dec = J.multiply(2 * math.pi / 365.0).subtract(1.39).sin().multiply(0.409)
        return dr, dec

    @staticmethod
    def _seasonal_correction(doy):
        J = ee.Number(doy)
        b = J.subtract(81).multiply(2 * math.pi / 364.0)
        Sc = (b.multiply(2).sin().multiply(0.1645)
              .subtract(b.cos().multiply(0.1255))
              .subtract(b.sin().multiply(0.025)))
        return Sc

    def _omega_mid_hour(self, lon_deg, doy, hour):
        """
        Soat burchagi hisoblash.

        ASCE-EWRI (2005) formula mahalliy standart vaqtni (LST) talab qiladi.
        GEE image timestampi UTC da keladi → LST ga o'tkazish shart:
          LST = UTC + Lz/15
          Idaho: Lz=-105 → UTC_offset = -7 → LST = UTC - 7

        F.28: t_sol = t_std + 0.06667*(Lm - Lz) + Sc
        omega = (t_sol - 12) * π/12
        """
        Sc = self._seasonal_correction(doy)

        # UTC → Local Standard Time (LST)
        # UTC_offset = Lz / 15  (Idaho: -105/15 = -7 soat)
        t_local_mid = ee.Number(hour).add(self.Lz.divide(15)).add(0.5)

        lon_img = ee.Image(lon_deg)
        # sol_time = t_local + (lon - Lz)*0.06667 + Sc
        sol_time = lon_img.subtract(ee.Image.constant(self.Lz)) \
                          .multiply(0.06667) \
                          .add(t_local_mid) \
                          .add(Sc)
        return sol_time.subtract(12).multiply(math.pi / 12.0)

    # === RADIATSIYA ===

    def Ra_hourly(self, lat_rad, lon_deg, doy, hour, dr, dec):
        omega = self._omega_mid_hour(lon_deg, doy, hour)
        omega1 = omega.subtract(math.pi / 24.0)
        omega2 = omega.add(math.pi / 24.0)
        ws = lat_rad.tan().multiply(-1).multiply(dec.tan()).acos()
        omega1c = omega1.max(ws.multiply(-1)).min(ws)
        omega2c = omega2.max(ws.multiply(-1)).min(ws)
        term1 = omega2c.subtract(omega1c).multiply(lat_rad.sin()).multiply(dec.sin())
        term2 = lat_rad.cos().multiply(dec.cos()).multiply(omega2c.sin().subtract(omega1c.sin()))
        Ra = term1.add(term2).multiply(dr).multiply(12 * 60 / math.pi * self.Gsc)
        return Ra.max(0)

    def Ra_daily(self, lat_rad, doy, dr, dec):
        ws = lat_rad.tan().multiply(-1).multiply(dec.tan()).acos()
        term1 = ws.multiply(lat_rad.sin()).multiply(dec.sin())
        term2 = lat_rad.cos().multiply(dec.cos()).multiply(ws.sin())
        Ra = term1.add(term2).multiply(dr).multiply(24 * 60 / math.pi * self.Gsc)
        return Ra

    @staticmethod
    def Rso(Ra, z_m):
        return z_m.multiply(2e-5).add(0.75).multiply(Ra)

    @staticmethod
    def Rnl_hourly(T_k, ea_kPa, Rs, Rso):
        Rso_safe = Rso.max(1e-6)
        rs_rso = Rs.divide(Rso_safe).clamp(0.0, 1.0)
        fcd_day = rs_rso.multiply(1.35).subtract(0.35)
        fcd_night = ee.Image.constant(1.0)
        is_day = Rs.gt(0.01)
        fcd = fcd_day.where(is_day.Not(), fcd_night)
        emiss = ea_kPa.sqrt().multiply(-0.14).add(0.34)
        sb = T_k.pow(4).multiply(2.043e-10)
        return sb.multiply(emiss).multiply(fcd)

    @staticmethod
    def Rnl_daily(Tmax_k, Tmin_k, ea_kPa, Rs, Rso):
        rs_rso = Rs.divide(Rso.max(1e-6)).clamp(0.0, 1.0)
        fcd = rs_rso.multiply(1.35).subtract(0.35)
        emiss = ea_kPa.sqrt().multiply(-0.14).add(0.34)
        tmean4 = Tmax_k.pow(4).add(Tmin_k.pow(4)).divide(2.0)
        sb = ee.Image.constant(4.903e-9)
        return sb.multiply(tmean4).multiply(emiss).multiply(fcd)

    def soil_heat_flux(self, Rn, mode):
        if mode == 'daily':
            return ee.Image.constant(0)
        if self.ref_type == 'grass':
            g_day = Rn.multiply(0.1)
            g_night = Rn.multiply(0.5)
        else:
            g_day = Rn.multiply(0.04)
            g_night = Rn.multiply(0.20)
        return g_day.where(Rn.lt(0), g_night)

    # === ASOSIY HISOBLASH ===

    def calculate(self, img, dem, mode='hourly'):
        """
        ASCE PM ETr hisoblash.
        Kirish: img (T_air, P_kPa, u2, ea, Rs bandlari), dem
        Chiqish: img + ETr band
        """
        img = ee.Image(img)
        dem = ee.Image(dem)

        ll = ee.Image.pixelLonLat()
        lat_rad = ll.select('latitude').multiply(math.pi / 180.0)
        lon_deg = ll.select('longitude')

        doy, hour = self._doy_hour(img)
        dr, dec = self._dr_decl(doy)

        if mode == 'hourly':
            Ra = self.Ra_hourly(lat_rad, lon_deg, doy, hour, dr, dec)
        else:
            Ra = self.Ra_daily(lat_rad, doy, dr, dec)

        Rso_img = self.Rso(Ra, dem)

        T = img.select('T_air')
        u2 = img.select('u2')
        P = img.select('P_kPa')
        ea = img.select('ea')
        Rs = img.select('Rs')

        if mode == 'daily':
            has_tmax = img.bandNames().contains('T_max')
            has_tmin = img.bandNames().contains('T_min')
            T_max_c = ee.Image(ee.Algorithms.If(has_tmax, img.select('T_max'), T))
            T_min_c = ee.Image(ee.Algorithms.If(has_tmin, img.select('T_min'), T))
            es = self._esat(T_max_c).add(self._esat(T_min_c)).divide(2).rename('es')
        else:
            es = self._esat(T).rename('es')

        Delta = self._delta_svp(T, es).rename('Delta')
        Gamma = self._gamma(P).rename('Gamma')

        Rns = Rs.multiply(1.0 - self.albedo)

        if mode == 'hourly':
            T_k = T.add(273.15)
            Rnl = self.Rnl_hourly(T_k, ea, Rs, Rso_img)
        else:
            Tmax_k = ee.Image(ee.Algorithms.If(
                img.bandNames().contains('T_max'),
                img.select('T_max').add(273.15), T.add(273.15)))
            Tmin_k = ee.Image(ee.Algorithms.If(
                img.bandNames().contains('T_min'),
                img.select('T_min').add(273.15), T.add(273.15)))
            Rnl = self.Rnl_daily(Tmax_k, Tmin_k, ea, Rs, Rso_img)

        Rn = Rns.subtract(Rnl).rename('Rn_ref')
        G = self.soil_heat_flux(Rn, mode).rename('G_ref')

        if mode == 'hourly':
            num1_day = Delta.multiply(Rn.subtract(G)).multiply(0.408)
            num2_day = (Gamma.multiply(u2).multiply(es.subtract(ea))
                        .multiply(ee.Image.constant(self.Cn_hr_day).divide(T.add(273.15))))
            den_day = Delta.add(Gamma.multiply(u2.multiply(self.Cd_hr_day).add(1.0)))
            ETr_day = num1_day.add(num2_day).divide(den_day)

            num1_night = Delta.multiply(Rn.subtract(G)).multiply(0.408)
            num2_night = (Gamma.multiply(u2).multiply(es.subtract(ea))
                          .multiply(ee.Image.constant(self.Cn_hr_night).divide(T.add(273.15))))
            den_night = Delta.add(Gamma.multiply(u2.multiply(self.Cd_hr_night).add(1.0)))
            ETr_night = num1_night.add(num2_night).divide(den_night)

            ETr = ETr_day.where(Rn.lt(0), ETr_night).max(0).rename('ETr')
        else:
            if mode == 'daily':
                T_max_c2 = ee.Image(ee.Algorithms.If(
                    img.bandNames().contains('T_max'), img.select('T_max'), T))
                T_min_c2 = ee.Image(ee.Algorithms.If(
                    img.bandNames().contains('T_min'), img.select('T_min'), T))
                es = self._esat(T_max_c2).add(self._esat(T_min_c2)).divide(2).rename('es')
                Delta = self._delta_svp(T, es).rename('Delta')

            num1 = Delta.multiply(Rn.subtract(G)).multiply(0.408)
            num2 = (Gamma.multiply(u2).multiply(es.subtract(ea))
                    .multiply(ee.Image.constant(self.Cn_day).divide(T.add(273.15))))
            den = Delta.add(Gamma.multiply(u2.multiply(self.Cd_day).add(1.0)))
            ETr = num1.add(num2).divide(den).max(0).rename('ETr')

        return img.addBands([es, Delta, Gamma, Ra.rename('Ra'),
                             Rso_img.rename('Rso'), Rnl.rename('Rnl'),
                             Rn, G, ETr])
