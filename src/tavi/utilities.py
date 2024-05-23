import numpy as np

# --------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------

# import scipy.constants as co
# ksq2E = (co.Planck / co.elementary_charge / 2.0 / np.pi) ** 2.0 * co.elementary_charge / 2.0 / co.neutron_mass * 1e23
ksq2E = 2.072124855  # calculated with scipy.constants using the formula above

sig2fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0))
cm2A = 1e8
min2rad = 1.0 / 60.0 / 180.0 * np.pi
rad2deg = 180.0 / np.pi
# --------------------------------------------------------------------------
# d_spacing table from Shirane Appendix 3, in units of Angstrom, from
# --------------------------------------------------------------------------
mono_ana_xtal = {
    "PG002": 3.35416,
    "Pg002": 3.35416,
    "PG004": 1.67708,
    "Cu111": 2.08717,
    "Cu220": 1.27813,
    "Ge111": 3.26627,
    "Ge220": 2.00018,
    "Ge311": 1.70576,
    "Ge331": 1.29789,
    "Be002": 1.79160,
    "Be110": 1.14280,
    "Heusler": 3.435,  # Cu2MnAl(111)
}
