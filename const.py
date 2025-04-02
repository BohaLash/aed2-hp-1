import math

R = 8.314462  # J^1 mol^-1 K^-1
avogadro = 6.02214076e23

molar_mass_air = 0.02896  # kg^1 mol^-1
R_air = R / molar_mass_air  # J^1 kg^-1 K^-1

Cp_air = 1003.5  # J^1 kg^-1 K^-1
Cv_air = 716.5  # J^1 kg^-1 K^-1

heat_capacity_ratio_air = Cp_air / Cv_air


def temperature_air(altitude: float):
    if altitude <= 11000:
        return 288.19 - 0.00649 * altitude
    elif altitude <= 25000:
        return 216.69
    else:
        return 141.94 + 0.00299 * altitude


def pressure_air(altitude: float):
    if altitude <= 11000:
        atm = 101290
        return atm * ((temperature_air(altitude) / 288.08) ** 5.256)
    elif altitude <= 25000:
        atm = 22650
        return atm * (math.e ** (1.73 - 0.000157 * altitude))
    else:
        atm = 2488
        return atm * ((temperature_air(altitude) / 216.6) ** -11.388)


def density_air(altitude: float):
    return pressure_air(altitude) / (286.9 * temperature_air(altitude))
