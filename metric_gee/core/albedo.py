"""
Module 2: Yer yuzasi albedosi — F.15.

L2 C2 Surface Reflectance allaqachon atmosfera tuzatilgan →
F.10-14 (TOA → surface) SKIP.

Mavjud usullar:
  olmedo  — Olmedo et al. (2012), L8/9 OLI SR uchun
  liang   — Liang (2001), umumiy Landsat, barcha vaznlar musbat
  ke      — Ke et al. (2016), L8 uchun, B_UB kanal qo'shilgan
  tasumi  — Tasumi et al. (2008), METRIC Idaho, F.15 ruhida
  avg3    — olmedo + liang + ke o'rtachasi (barqarorroq)

run.py da foydalanuvchi qo'lda tanlaydi:
  cfg.albedo_method = 'liang'   # yoki 'olmedo' / 'ke' / 'tasumi' / 'avg3'

Albedo diapazoni: [0, 1]
  Qor: ~0.80  |  Suv: ~0.06  |  O'simlik: 0.10-0.25  |  Quruq tuproq: 0.25-0.45
"""
import ee


class AlbedoCalculator:
    """Surface albedo hisoblash — L8/9 C2 L2 SR uchun."""

    def __init__(self, settings):
        self.settings = settings

    def calc_albedo(self, image, method=None):
        """
        F.15: αs = Σ(ρs,b × wb) — surface albedo.

        Args:
            image  : ee.Image — scaled SR bands (B_UB, B_BLUE, B_GREEN, ...)
            method : str — 'olmedo' | 'liang' | 'ke' | 'tasumi' | 'avg3'
                     None → settings.albedo_method

        Returns:
            ee.Image — 'albedo' band, clamp [0, 1]
        """
        method = method or self.settings.albedo_method

        # Umumiy bandlar — barcha usullarda ishlatiladi
        blue  = image.select('B_BLUE')
        green = image.select('B_GREEN')
        red   = image.select('B_RED')
        nir   = image.select('B_NIR')
        swir1 = image.select('B_SWIR1')
        swir2 = image.select('B_SWIR2')

        if method == 'olmedo':
            # Olmedo et al. (2012) — L8/9 OLI SR
            # Manba: Remote Sensing of Environment, 118, 320-331
            albedo = (blue.multiply(0.4739)
                     .add(green.multiply(-0.4372))
                     .add(red.multiply(0.1652))
                     .add(nir.multiply(0.2831))
                     .add(swir1.multiply(0.1072))
                     .add(swir2.multiply(0.1029))
                     .add(0.0366))

        elif method == 'liang':
            # Liang (2001) — Landsat uchun universal, barcha vaznlar musbat
            # Manba: Remote Sensing of Environment, 76, 213-238
            albedo = (blue.multiply(0.356)
                     .add(red.multiply(0.130))
                     .add(nir.multiply(0.373))
                     .add(swir1.multiply(0.085))
                     .add(swir2.multiply(0.072))
                     .subtract(0.0018))

        elif method == 'ke':
            # Ke et al. (2016) — L8 OLI, B_UB (Band 1) bilan
            # Manba: Int. J. Remote Sensing, 37, 4490-4507
            ub = image.select('B_UB')
            albedo = (ub.multiply(0.130)
                     .add(blue.multiply(0.115))
                     .add(green.multiply(0.143))
                     .add(red.multiply(0.180))
                     .add(nir.multiply(0.281))
                     .add(swir1.multiply(0.108))
                     .add(swir2.multiply(0.042)))

        elif method == 'tasumi':
            # Tasumi et al. (2008) — METRIC Idaho, F.15 ruhida
            # Allen 2007 maqolasidagi vaznlarga yaqin (L5 uchun asl)
            # Manba: ASCE Journal of Irrigation and Drainage Engineering
            albedo = (blue.multiply(0.300)
                     .add(green.multiply(0.277))
                     .add(red.multiply(0.233))
                     .add(nir.multiply(0.143))
                     .add(swir1.multiply(0.036))
                     .add(swir2.multiply(0.012)))

        elif method == 'avg3':
            # Olmedo + Liang + Ke o'rtachasi — uchta usul bilan barqarorroq
            ub = image.select('B_UB')

            a_olmedo = (blue.multiply(0.4739)
                       .add(green.multiply(-0.4372))
                       .add(red.multiply(0.1652))
                       .add(nir.multiply(0.2831))
                       .add(swir1.multiply(0.1072))
                       .add(swir2.multiply(0.1029))
                       .add(0.0366))

            a_liang = (blue.multiply(0.356)
                      .add(red.multiply(0.130))
                      .add(nir.multiply(0.373))
                      .add(swir1.multiply(0.085))
                      .add(swir2.multiply(0.072))
                      .subtract(0.0018))

            a_ke = (ub.multiply(0.130)
                   .add(blue.multiply(0.115))
                   .add(green.multiply(0.143))
                   .add(red.multiply(0.180))
                   .add(nir.multiply(0.281))
                   .add(swir1.multiply(0.108))
                   .add(swir2.multiply(0.042)))

            albedo = a_olmedo.add(a_liang).add(a_ke).divide(3)

        else:
            raise ValueError(
                f"Noma'lum albedo usuli: '{method}'. "
                f"Mavjudlar: olmedo | liang | ke | tasumi | avg3"
            )

        return albedo.clamp(0, 1).rename('albedo')
