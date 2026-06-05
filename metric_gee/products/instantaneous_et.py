"""
Module 10: Lahzali ET.
F.1:  LE = Rn - G - H
F.52: ETinst (mm/hr) = 3600 * LE / lambda
F.53: lambda = (2.501 - 0.00236*(Ts-273.15)) * 1e6

BIRLIK TAHLILI:
  LE       [W/m²] = [J / (m² · s)]
  lambda   [J/kg]
  rho_w    [kg/m³] = 1000

  ETinst [m/s]   = LE / (lambda × rho_w)
  ETinst [m/hr]  = LE × 3600 / (lambda × rho_w)
  ETinst [mm/hr] = LE × 3600 / (lambda × rho_w) × 1000

  rho_w = 1000 va ×1000 mm/m → bir-birini bekor qiladi:
  ETinst [mm/hr] = LE × 3600 / lambda   ← soddalashtirilgan formula
"""
import ee


class InstantaneousET:
    """Lahzali bug'lanish hisoblash."""

    def calc(self, rn, g, h, ts):
        """
        LE → ETinst (mm/hr).

        Returns:
            tuple: (le, et_inst) — ee.Image
        """
        # F.1: LE = Rn - G - H
        le = rn.subtract(g).subtract(h).rename('LE')

        # F.53: lambda (J/kg) = (2.501 - 0.00236*(Ts-273.15)) * 1e6
        lam = (ts.subtract(273.15)
               .multiply(-0.00236)
               .add(2.501)
               .multiply(1e6))

        # F.52: ETinst [mm/hr] = LE × 3600 / lambda
        # (rho_w=1000 va mm/m=1000 bir-birini bekor qiladi)
        et_inst = (le
                   .multiply(3600)
                   .divide(lam)
                   .rename('ETinst'))

        return le, et_inst
