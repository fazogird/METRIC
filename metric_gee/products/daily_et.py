"""
Module 11: Kunlik ET.
F.54: ETrF = ETinst / ETr
F.55: ET24 = Crad * ETrF * ETr24
F.56: Crad — tekis yer: 1.0, tog'li: slope correction
"""
import ee


class DailyET:
    """Kunlik bug'lanish hisoblash."""
    

    def calc_etrf(self, et_inst, etr):
        """F.54: ETrF = ETinst / ETr — crop coefficient analog."""
        
        etr_matched = (
            etr
            .resample('bilinear')
            .reproject(crs=et_inst.projection())
            .max(0.01)      # division by zero himoya
        )

        return (et_inst
                .divide(etr_matched.max(0.01))
                .clamp(0, 1.05)
                .rename('ETrF'))

    def calc_et24(self, etrf, etr_24, crad=None):
        """F.55: ET24 = Crad * ETrF * ETr24 (mm/day)."""
        if crad is None:
            crad = ee.Image.constant(1.0)
        return (crad
                .multiply(etrf)
                .multiply(etr_24)
                .rename('ET24'))
