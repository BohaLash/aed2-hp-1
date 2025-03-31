import math
import numpy as np  # Import numpy for ln

# -----------------------------------------------------------------------------
# Task A: Data Collection
# -----------------------------------------------------------------------------

print("--- Task A: Data Collection ---")

# A.1: Physical Constants and Standard Atmosphere (ISA) Data
# Source: Standard definitions and Google Search results
g = 9.80665  # Gravity (m/s^2)
R_air = 287.053  # Specific gas constant for dry air (J/kg*K)
gamma_air = 1.4  # Ratio of specific heats for air

# ISA Sea Level Conditions
T0_isa = 288.15  # Temperature (K)
P0_isa = 101325  # Pressure (Pa)
rho0_isa = 1.225  # Density (kg/m^3)

# ISA Lapse Rate (Troposphere)
L_isa = 0.0065  # Temperature lapse rate (K/m)


# Function to calculate ISA properties in the Troposphere (h < 11km)
def calculate_isa_troposphere(h_m):
    """Calculates ISA T, P, rho, a for a given altitude h_m (meters) in the troposphere."""
    T = T0_isa - L_isa * h_m
    P = P0_isa * (T / T0_isa) ** (g / (L_isa * R_air))
    rho = P / (R_air * T)
    a = math.sqrt(gamma_air * R_air * T)
    return T, P, rho, a


# Calculate/Define conditions at relevant altitudes
# Takeoff Altitude: 1.6 km = 1600 m
h_takeoff_m = 1600
T_isa_takeoff, P_isa_takeoff, rho_isa_takeoff, a_isa_takeoff = (
    calculate_isa_troposphere(h_takeoff_m)
)

# Single-Engine Cruise Altitude: 5 km = 5000 m
h_sec_m = 5000
T_isa_sec, P_isa_sec, rho_isa_sec, a_isa_sec = calculate_isa_troposphere(h_sec_m)

# Initial Cruise Altitude: 11 km = 11000 m (Tropopause)
# Source: Google Search result (Wikipedia) & calculated 'a'
h_cruise_m = 11000
T_isa_cruise = 216.65  # Temperature (K)
P_isa_cruise = 22632  # Pressure (Pa)
rho_isa_cruise = 0.3639  # Density (kg/m^3)
a_isa_cruise = math.sqrt(gamma_air * R_air * T_isa_cruise)  # Speed of sound (m/s)

print(f"ISA Conditions at Takeoff Altitude ({h_takeoff_m} m):")
print(f"  T = {T_isa_takeoff:.2f} K")
print(f"  P = {P_isa_takeoff:.2f} Pa")
print(f"  rho = {rho_isa_takeoff:.4f} kg/m^3")
print(f"  a = {a_isa_takeoff:.2f} m/s")
print(f"ISA Conditions at Single-Engine Altitude ({h_sec_m} m):")
print(f"  T = {T_isa_sec:.2f} K")
print(f"  P = {P_isa_sec:.2f} Pa")
print(f"  rho = {rho_isa_sec:.4f} kg/m^3")
print(f"  a = {a_isa_sec:.2f} m/s")
print(f"ISA Conditions at Cruise Altitude ({h_cruise_m} m):")
print(f"  T = {T_isa_cruise:.2f} K")
print(f"  P = {P_isa_cruise:.2f} Pa")
print(f"  rho = {rho_isa_cruise:.4f} kg/m^3")
print(f"  a = {a_isa_cruise:.2f} m/s\n")

# A.2: HP-1 Aircraft and Mission Parameters
# Source: manual-md.txt (Task A description)
M0_cruise = 0.83  # Cruise Mach number
h_takeoff_km = 1.6  # Takeoff altitude (km)
T_day_takeoff_C = 38  # Takeoff temperature on hot day (deg C)
T_day_takeoff_K = T_day_takeoff_C + 273.15  # Takeoff temperature on hot day (K)
runway_m = 3650  # Runway length (m)
climb_gradient_oei = 0.024  # Single-engine climb gradient (%)
payload_kg = 22770  # Payload (253 passengers * 90 kg)
Range_km = 11120  # Design range (km)
Range_m = Range_km * 1000  # Design range (m)
loiter_min = 30  # Reserve fuel loiter time (min)
h_cruise_km = 11  # Initial cruise altitude (km)
Ps_cruise = 1.5  # Rate of climb at initial cruise (m/s)
h_sec_km = 5  # Single-engine cruise altitude (km)
M0_sec = 0.45  # Single-engine cruise Mach number
Ps_sec = 1.5  # Rate of climb at single-engine cruise (m/s)
N_eng = 2  # Number of engines

print("HP-1 Aircraft and Mission Parameters:")
print(f"  Cruise Mach: {M0_cruise}")
print(f"  Cruise Altitude: {h_cruise_km} km")
print(f"  Payload: {payload_kg} kg")
print(f"  Design Range: {Range_km} km")
print(f"  Number of Engines: {N_eng}\n")

# A.3: Data from Example 1.2 (Mattingly Book) - **ASSUMED VALUES**
# Source: manual-md.txt references Example 1.2. Actual values needed from book file.
# !! These are HYPOTHETICAL values for demonstration !!
W_TO_kg = 150000.0  # Assumed Max Takeoff Weight (kg)
W_empty_kg = 70000.0  # Assumed Empty Weight (kg)
LD_cruise = 18.0  # Assumed Lift-to-Drag ratio at cruise
LD_loiter = 15.0  # Assumed Lift-to-Drag ratio at loiter

print("Data assumed from Example 1.2 (Mattingly Book):")
print(f"  Assumed W_TO: {W_TO_kg} kg")
print(f"  Assumed W_empty: {W_empty_kg} kg")
print(f"  Assumed L/D (Cruise): {LD_cruise}")
print(f"  Assumed L/D (Loiter): {LD_loiter}\n")

# A.4: Mission Fuel Fractions from Table 1 - **ASSUMED VALUES**
# Source: manual-md.txt mentions Table 1. Actual values needed from manual file.
# !! These are HYPOTHETICAL values for demonstration !!
# W_i / W_{i-1} fuel fractions for mission segments
fuel_frac = {
    "start_warmup": 0.995,
    "taxi": 0.997,
    "takeoff": 0.998,
    "climb": 0.980,
    "cruise": 0.820,  # This fraction significantly impacts range/TSFC calculations
    "descent": 0.990,
    "loiter": 0.985,  # Represents reserve fuel
    "landing_taxi": 0.997,
}

print("Assumed Fuel Fractions from Table 1 (Manual):")
for segment, fraction in fuel_frac.items():
    print(f"  {segment}: {fraction}")
print("\n")


# -----------------------------------------------------------------------------
# Task B: Thrust and Fuel Consumption Requirements
# -----------------------------------------------------------------------------

print("--- Task B: Requirement Calculations ---")

# B.1: Calculate Overall Fuel Fraction and Fuel Weight
# Product of all segment fractions gives W_final / W_initial (W_landing / W_TO)
W_ratio_overall = 1.0
for segment in fuel_frac:
    W_ratio_overall *= fuel_frac[segment]

# Product of fractions for mission *excluding* loiter (reserve)
W_ratio_mission = (
    fuel_frac["start_warmup"]
    * fuel_frac["taxi"]
    * fuel_frac["takeoff"]
    * fuel_frac["climb"]
    * fuel_frac["cruise"]
    * fuel_frac["descent"]
    * fuel_frac["landing_taxi"]
)

W_landing_kg = W_TO_kg * W_ratio_overall
W_fuel_total_kg = W_TO_kg * (1.0 - W_ratio_overall)
W_fuel_mission_kg = W_TO_kg * (
    1.0 - W_ratio_mission / fuel_frac["loiter"]
)  # Approx mission fuel
W_fuel_reserve_kg = (
    W_TO_kg * (1.0 - fuel_frac["loiter"]) * W_ratio_mission
)  # Approx reserve fuel

print("Fuel Consumption Requirements (based on assumed fractions):")
print(f"  Overall Weight Ratio (W_landing / W_TO): {W_ratio_overall:.4f}")
print(f"  Total Fuel Weight: {W_fuel_total_kg:.2f} kg")
print(
    f"  Mission Fuel Fraction (W_fuel_mission / W_TO): {(W_fuel_mission_kg / W_TO_kg):.4f}"
)
print(
    f"  Reserve Fuel Fraction (W_fuel_reserve / W_TO): {(W_fuel_reserve_kg / W_TO_kg):.4f}\n"
)

# B.2: Calculate Weight at Start of Cruise
W_start_climb = (
    W_TO_kg * fuel_frac["start_warmup"] * fuel_frac["taxi"] * fuel_frac["takeoff"]
)
W_start_cruise_kg = W_start_climb * fuel_frac["climb"]

print(f"Weight at Start of Cruise: {W_start_cruise_kg:.2f} kg\n")

# B.3: Calculate Thrust Requirement at Start of Cruise
# Assuming Level Flight (Thrust = Drag) at the beginning of the cruise segment
# Drag = Weight / (L/D)
Thrust_total_cruise_N = (W_start_cruise_kg * g) / LD_cruise
Thrust_per_engine_cruise_N = Thrust_total_cruise_N / N_eng

print("Thrust Requirements (Start of Cruise, Level Flight):")
print(f"  Total Thrust Required: {Thrust_total_cruise_N:.2f} N")
print(f"  Thrust per Engine: {Thrust_per_engine_cruise_N:.2f} N\n")

# B.4: Implied Cruise TSFC (Informational)
# Calculate TSFC (S) implied by the assumed cruise fuel fraction and L/D
# Breguet Range: R = (V / (g*S)) * (L/D) * ln(W_start / W_end)
# S = (V / (g*R)) * (L/D) * ln(W_start / W_end)
V_cruise_mps = M0_cruise * a_isa_cruise
W_end_cruise_kg = W_start_cruise_kg * fuel_frac["cruise"]

# Check for valid weights and cruise fraction before calculation
if (
    W_start_cruise_kg > W_end_cruise_kg
    and W_end_cruise_kg > 0
    and fuel_frac["cruise"] < 1.0
):
    # Use the design range (Range_m) and cruise segment weight ratio
    S_implied_cruise_sec_kg_N = (
        (V_cruise_mps / (g * Range_m))
        * LD_cruise
        * np.log(W_start_cruise_kg / W_end_cruise_kg)
    )
    S_implied_cruise_hr_kg_N = (
        S_implied_cruise_sec_kg_N * 3600
    )  # Convert from per sec to per hour
    # Convert to lbm/(hr*lbf) = (kg/(hr*N)) * (1/0.45359 kg/lbm) / (1/4.4482 N/lbf)
    S_implied_cruise_hr_lbm_lbf = (
        S_implied_cruise_hr_kg_N * (1 / 0.45359) / (1 / 4.4482)
    )

    print("Implied Cruise TSFC (based on assumed W_TO, L/D, and fuel fractions):")
    print(f"  V_cruise: {V_cruise_mps:.2f} m/s")
    print(f"  Cruise Fuel Ratio (W_end / W_start): {fuel_frac['cruise']:.4f}")
    print(f"  Implied S: {S_implied_cruise_sec_kg_N:.4e} kg/(N*s)")
    print(f"  Implied S: {S_implied_cruise_hr_kg_N:.4f} kg/(N*hr)")
    print(f"  Implied S: {S_implied_cruise_hr_lbm_lbf:.4f} lbm/(lbf*hr)")
    print(
        "  Note: This TSFC is highly sensitive to the assumed cruise fuel fraction and L/D."
    )
else:
    print(
        "Could not calculate implied cruise TSFC due to invalid weight or fuel fraction values."
    )

print("\n--- End of Script ---")
