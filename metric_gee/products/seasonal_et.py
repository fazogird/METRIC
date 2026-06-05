"""
Module 12: Mavsumiy ET.
F.58: ETperiod = sum(ETrFi * ETr24i)
F.59: ETrF_period = sum(ETrFi * ETr24i) / sum(ETr24i)
"""
import ee


class SeasonalET:
    """Mavsumiy/oylik kumulyativ ET."""

    def calc_et_period(self, etrf_images, etr24_images):
        """F.58: Davr bo'yicha kumulyativ ET (mm)."""
        total = ee.Image.constant(0)
        for etrf, etr24 in zip(etrf_images, etr24_images):
            total = total.add(etrf.multiply(etr24))
        return total.rename('ET_period')

    def calc_etrf_period(self, etrf_images, etr24_images):
        """F.59: ETr vazn berilgan o'rtacha ETrF."""
        numer = ee.Image.constant(0)
        denom = ee.Image.constant(0)
        for etrf, etr24 in zip(etrf_images, etr24_images):
            numer = numer.add(etrf.multiply(etr24))
            denom = denom.add(etr24)
        return numer.divide(denom.max(0.01)).rename('ETrF_period')
