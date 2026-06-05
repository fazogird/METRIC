"""
Quyosh geometriyasi — F.7, F.8, F.9.
DOY, d², delta, omega, cos(theta_hor), cos(theta_rel).
"""
import ee
import math


class SolarGeometry:
    """Quyosh pozitsiyasi hisoblash."""

    def calc_solar_params(self, image_date, geometry, sun_elevation=None):
        """
        Barcha quyosh parametrlari.
        
        Returns:
            dict: d2, cos_theta_hor, cos_theta_rel, sun_elevation
        """
        date = ee.Date(image_date)
        doy = date.getRelative('day', 'year').add(1)

        # F.9: d² = 1 / (1 + 0.033·cos(DOY·2π/365))
        d2 = (doy.multiply(2 * math.pi / 365).cos()
              .multiply(0.033).add(1).pow(-1))

        # Declination (rad) — Cooper equation
        delta = (doy.subtract(1).multiply(2 * math.pi / 365)
                 .subtract(ee.Number(math.pi * 80 / 365))
                 .sin().multiply(0.4093))

        # Hour angle — Landsat overpass ~10:30 local
        omega = ee.Number(-0.3)

        # F.8: cos(theta_hor) — tekis yer
        lat_deg = ee.Number(geometry.centroid().coordinates().get(1))
        lat_rad = lat_deg.multiply(math.pi / 180)

        cos_theta = (delta.sin().multiply(lat_rad.sin())
                     .add(delta.cos().multiply(lat_rad.cos())
                          .multiply(omega.cos())))

        if sun_elevation is not None:
            sun_elev = ee.Number(sun_elevation)
        else:
            sun_elev = cos_theta.acos().multiply(-1).add(math.pi / 2)
            sun_elev = sun_elev.multiply(180 / math.pi)

        return {
            'doy': doy,
            'd2': d2,
            'delta': delta,
            'omega': omega,
            'cos_theta_hor': cos_theta,
            'cos_theta_rel': cos_theta,
            'sun_elevation': sun_elev,
        }
