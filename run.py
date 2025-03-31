import math

# --------------------------------------------------------------------------
# Constants and Atmospheric Data at 11 km (ISA Model)
# --------------------------------------------------------------------------
GAMMA = 1.4  # Ratio of specific heats for air (Replace if different)
R = 287.0  # Specific gas constant for air (J/kg·K) (Replace if different)
ALTITUDE = 11000  # meters
MACH_CRUISE = 0.83
PS = 1.5  # Rate of Climb (m/s)
PHI_LOSS = 0.02  # Inlet + Nozzle loss factor

# ISA Conditions at 11km
T_11km = 216.65  # Temperature (K)
P_11km = 22632  # Pressure (Pa)
rho_11km = 0.3639  # Density (kg/m^3)

# Calculate Speed of Sound and Cruise Velocity
a_11km = math.sqrt(GAMMA * R * T_11km)
V_cruise = MACH_CRUISE * a_11km

print(f"--- Conditions at {ALTITUDE} m ---")
print(f"Temperature (T): {T_11km:.2f} K")
print(f"Pressure (P): {P_11km:.2f} Pa")
print(f"Density (rho): {rho_11km:.4f} kg/m^3")
print(f"Speed of Sound (a): {a_11km:.2f} m/s")
print(f"Cruise Velocity (V0): {V_cruise:.2f} m/s")
print("-" * 30)

# --------------------------------------------------------------------------
# User Input Required Below
# --------------------------------------------------------------------------

# 1. Aircraft Weight and Drag/LD data
W_cruise = None  # <<<--- ENTER Aircraft Weight at cruise (Newtons)
LD_cruise = None  # <<<--- ENTER Lift-to-Drag Ratio at cruise (dimensionless)
# OR
# D_cruise = None # <<<--- ENTER Drag at cruise (Newtons) if L/D is not known

if W_cruise is None or (LD_cruise is None and D_cruise is None):
    print("ERROR: Please enter Aircraft Weight and L/D ratio (or Drag).")
    exit()

# Calculate Drag if L/D is given
# Assuming Lift approx equals Weight in cruise (or small climb angle)
if D_cruise is None:
    # Note: For climbing flight L = W * cos(climb_angle), but angle is small for Ps=1.5m/s
    # Theta = arcsin(Ps/V_cruise) approx Ps/V_cruise
    # For Ps=1.5, V0 approx 245 m/s, theta is very small, cos(theta) approx 1
    # If a more precise calculation is needed, adjust lift: L = W * cos(math.asin(PS / V_cruise))
    Lift_approx = W_cruise
    D_cruise = Lift_approx / LD_cruise
    print(f"Calculated Drag (D): {D_cruise:.2f} N (assuming L=W)")


# 2. Equation (1.28) Implementation
def calculate_installed_thrust(Drag, Weight, Velocity, Ps_climb):
    """
    Calculates required installed thrust based on Eq. 1.28.
    *** YOU MUST REPLACE THE FORMULA BELOW WITH THE ACTUAL EQ 1.28 ***
    """
    print(
        "\nWARNING: Using placeholder formula for Eq 1.28. Replace with actual equation."
    )
    # Placeholder formula: Thrust = Drag + Weight * (Ps / Velocity)
    # This is the basic thrust required for steady level flight + climb component
    # *** REPLACE THIS LINE with the calculation from Eq 1.28 from Mattingly ***
    installed_thrust_eq128 = Drag + Weight * (Ps_climb / Velocity)
    # ***********************************************************************

    if installed_thrust_eq128 is None:
        print(
            "ERROR: Implement Equation 1.28 in the 'calculate_installed_thrust' function."
        )
        exit()
    return installed_thrust_eq128


# 3. Maximum Inlet Mass Flow Equation Implementation
def calculate_max_mass_flow(
    diameter, rho_air, velocity_air, mach_air, temp_air, press_air
):
    """
    Calculates maximum inlet mass flow based on equation from Mattingly Ch 1.
    *** YOU MUST REPLACE THE FORMULA BELOW WITH THE ACTUAL EQUATION ***
    """
    print(
        "WARNING: Using placeholder formula for Max Mass Flow. Replace with actual equation."
    )
    inlet_area = math.pi * (diameter / 2) ** 2
    # Placeholder formula: Basic ram air capture: mdot = rho * Area * Velocity
    # *** REPLACE THIS LINE with the calculation from Mattingly Ch 1 ***
    # Example might involve stagnation properties or corrected flow parameters.
    # E.g., mdot = A * P0 * sqrt(gamma / (R * T0)) * M * (1 + (gamma-1)/2 * M^2)^(-(gamma+1)/(2*(gamma-1)))
    # where P0, T0 are stagnation properties, but use the specific one from the book.
    mdot_max_placeholder = rho_air * inlet_area * velocity_air
    # ***********************************************************************

    if mdot_max_placeholder is None:
        print(
            "ERROR: Implement Mass Flow equation in 'calculate_max_mass_flow' function."
        )
        exit()
    return mdot_max_placeholder


# --------------------------------------------------------------------------
# Calculations for Section C
# --------------------------------------------------------------------------
print("\n--- Section C Calculations ---")

inlet_diameters = [2.2, 2.5, 2.75, 3.0, 3.25, 3.5]  # meters

results = {}

# Step B1: Calculate Required Thrust (Done once as it doesn't depend on inlet diameter)
try:
    T_installed = calculate_installed_thrust(D_cruise, W_cruise, V_cruise, PS)
    T_uninstalled = T_installed / (1.0 - PHI_LOSS)
    print(f"\nStep B1 Results:")
    print(f"Required Installed Thrust (T_inst): {T_installed:.2f} N")
    print(f"Required Uninstalled Thrust (T_uninst): {T_uninstalled:.2f} N")
except Exception as e:
    print(f"\nError in Step B1 calculation: {e}")
    print("Cannot proceed without valid Step B1 results.")
    exit()

print("-" * 30)

# Steps B2 & B3: Iterate through diameters
for D_inlet in inlet_diameters:
    print(f"\nCalculating for Inlet Diameter: {D_inlet} m")
    inlet_area = math.pi * (D_inlet / 2) ** 2
    print(f"Inlet Area (A0): {inlet_area:.4f} m^2")

    try:
        # Step B2: Calculate Max Mass Flow
        mdot_max = calculate_max_mass_flow(
            D_inlet, rho_11km, V_cruise, MACH_CRUISE, T_11km, P_11km
        )
        print(f"Step B2: Max Inlet Mass Flow (mdot_max): {mdot_max:.3f} kg/s")

        # Step B3: Calculate Minimum Uninstalled Specific Thrust
        F_sp_uninst = T_uninstalled / mdot_max  # N / (kg/s) or m/s
        print(
            f"Step B3: Min Uninstalled Specific Thrust (Fsp_uninst): {F_sp_uninst:.2f} N/(kg/s)"
        )

        results[D_inlet] = F_sp_uninst

    except Exception as e:
        print(f"Error calculating for diameter {D_inlet} m: {e}")
        results[D_inlet] = None

print("\n--- Summary: Minimum Uninstalled Specific Thrust vs Diameter ---")
for diameter, f_sp in results.items():
    if f_sp is not None:
        print(f"Diameter: {diameter:.2f} m  ->  Min Fsp_uninst: {f_sp:.2f} N/(kg/s)")
    else:
        print(f"Diameter: {diameter:.2f} m  ->  Calculation Error")
print("-" * 60)
