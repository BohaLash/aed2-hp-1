import math
import numpy as np
import matplotlib.pyplot as plt

from const import *


def speed_of_sound(
    temperature: float,
    heat_capacity_ratio: float,
    gas_constant: float,
):
    return math.sqrt(heat_capacity_ratio * gas_constant * temperature)


def speed_of_sound_air(temperature: float):
    return speed_of_sound(temperature, heat_capacity_ratio_air, R_air)


def mach_number(velocity: float, altitude: float):
    return velocity / speed_of_sound_air(temperature_air(altitude))


def total_temperature_ratio(
    heat_capacity_ratio: float,
    mach_number: float,
):
    return 1 + 0.5 * (heat_capacity_ratio - 1) * mach_number * mach_number


def total_pressure_ratio(
    heat_capacity_ratio: float,
    mach_number: float,
):
    exponent = heat_capacity_ratio / (heat_capacity_ratio - 1)
    return 1 + (0.5 * (heat_capacity_ratio - 1) * mach_number**2) ** exponent


atmosphere = np.genfromtxt("atmosphere.csv", delimiter=",", names=True)

plt.plot(
    atmosphere["altitude"],
    atmosphere["density"],
    label="csv",
    color="blue",
    linestyle="--",
)

altitude = np.linspace(0, 20000, 200)
density = np.array([density_air(h) for h in altitude])

plt.plot(
    altitude,
    density,
    label="func",
    color="red",
    linestyle=":",
)

plt.title("Density [kg/m3]")
plt.xlabel("altitude")
plt.ylabel("density")
plt.legend()
plt.grid(True)
plt.show()
