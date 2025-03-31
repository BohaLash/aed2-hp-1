import numpy as np

# --- Constants ---
g0 = 9.80665  # Standard gravity (m/s^2)
R_air = 287.05  # Specific gas constant for air (J/kg-K)
gamma_air = 1.4 # Ratio of specific heats for air

# --- Standard Atmosphere Functions (ISA Model) ---
def get_isa_temp(alt_m):
    """Calculates ISA temperature (Kelvin) at a given altitude (meters)."""
    if alt_m < 11000:
        # Troposphere
        return 288.15 - 0.0065 * alt_m
    elif alt_m <= 20000:
        # Lower Stratosphere
        return 216.65
    else:
        # Simplified: Assume constant temp above 20km for this scope
        # For more accuracy, extend the ISA model
        # print("Warning: Altitude > 20km, ISA model simplified.")
        return 216.65 # Placeholder for higher altitudes if needed

def get_isa_pressure(alt_m):
    """Calculates ISA pressure (Pascals) at a given altitude (meters). Corrected for recursion."""
    T0 = 288.15 # Sea level temperature (K)
    P0 = 101325 # Sea level pressure (Pa)
    lapse_rate = 0.0065 # K/m
    
    # Calculate conditions at the tropopause (11km) once
    T11k = 288.15 - lapse_rate * 11000 # Temp at 11km
    P11k = P0 * (T11k / T0)**(g0 / (lapse_rate * R_air)) # Pressure at 11km using troposphere formula
    
    if alt_m < 11000:
        # Troposphere
        T = get_isa_temp(alt_m) # Or use T = T0 - lapse_rate * alt_m
        # Ensure T is not zero or negative which might happen for very high incorrect altitudes
        if T <= 0: 
             raise ValueError("Temperature calculated to be zero or negative in Troposphere model.")
        return P0 * (T / T0)**(g0 / (lapse_rate * R_air))
    elif alt_m <= 20000:
        # Lower Stratosphere (Isothermal layer)
        # Use the pre-calculated P11k and T11k (which is constant 216.65 K)
        T_strat = 216.65 # Constant temperature in lower stratosphere
        # Correct calculation uses P(h) = P_base * exp(-g0*(h-h_base)/(R*T_base))
        return P11k * np.exp(-g0 * (alt_m - 11000) / (R_air * T_strat))
    else:
        # Simplified: Use 20km value above
         P20k = get_isa_pressure(20000) # Calculate pressure at 20km using the stratosphere formula above
         T20k = 216.65 # Temp is constant up to 20km
         # Placeholder for altitudes above 20km. A more complete ISA model is needed.
         # For this problem focusing on 11km cruise, this simplification is likely acceptable.
         print("Warning: Altitude > 20km, ISA model simplified. Returning 20km pressure.")
         # Example: Extend using the next layer's lapse rate if needed
         # lapse_rate_upper_strat = 0.001 # K/m (Example)
         # T = T20k + lapse_rate_upper_strat * (alt_m - 20000)
         # P = P20k * (T / T20k)**(-g0 / (lapse_rate_upper_strat * R_air)) # Check formula for positive lapse rate
         return P20k # Return 20km pressure as a fallback

def get_isa_density(alt_m):
    """Calculates air density (kg/m^3) at a given altitude (meters)."""
    P = get_isa_pressure(alt_m)
    T = get_isa_temp(alt_m)
    if T == 0:
        raise ValueError("ISA Temperature is zero, cannot calculate density.")
    return P / (R_air * T)

def get_speed_of_sound(alt_m):
    """Calculates speed of sound (m/s) at a given altitude (meters)."""
    T = get_isa_temp(alt_m)
    if T < 0:
        raise ValueError("ISA Temperature is negative, cannot calculate speed of sound.")
    return np.sqrt(gamma_air * R_air * T)

# --- HP-1 Aircraft and Mission Parameters ---
# These values should be cross-referenced with manual.md.txt and book.md.txt (Example 1.2)
cruise_mach = 0.83
cruise_alt_m = 11000  # meters
range_km = 11120
range_m = range_km * 1000

# !!! Placeholder: Aircraft weight at cruise start/average (kg).
# This needs to be determined from mission analysis (payload, empty weight, fuel burn).
# Refer to Table 1 fuel fractions and Example 1.2 data.
W_cruise_kg = 150000 # Example Value - **MUST BE UPDATED**

# !!! Placeholder: Lift-to-Drag ratio at cruise.
# Check Example 1.2 in book.md.txt. Typical values might be 15-20 for transports.
L_D_cruise = 17.0 # Example Value - **MUST BE UPDATED**

# !!! Placeholder: Mission fuel fractions (Wi+1 / Wi).
# These should come from Table 1 in manual.md.txt or detailed mission analysis.
# The keys represent mission segments, values are weight fractions.
mission_fuel_fractions = {
    "warmup_takeoff": 0.990, # Example Value - **MUST BE UPDATED**
    "climb": 0.980,          # Example Value - **MUST BE UPDATED**
    # "cruise": This is what we often solve for or relate to range/TSFC
    "descent": 0.990,        # Example Value - **MUST BE UPDATED**
    "loiter": 0.985          # Example Value - **MUST BE UPDATED** (Based on 30 min loiter requirement)
}

# Calculate W_cruise_frac = W_cruise_end / W_cruise_start
# This calculation depends heavily on how Table 1 data is presented.
# Assuming Table 1 gives fuel consumed as fraction of W_TO (Wf/W0),
# Let W_payload = 22770 kg. Need W_empty from Ex 1.2. W0 = W_empty + W_payload + W_fuel_total
# Let's assume Table 1 gives W_{i+1}/W_i for non-cruise segments and Wf_cruise/W0.
# Then W_cruise_start = W0 * frac_warmup * frac_climb
# W_cruise_end = W_cruise_start - Wf_cruise = W_cruise_start - (Wf_cruise/W0)*W0
# W_cruise_frac = W_cruise_end / W_cruise_start = 1 - (Wf_cruise/W0) / (frac_warmup * frac_climb)
# !!! Placeholder - This calculation MUST be adapted based on actual Table 1 structure in manual.md.txt
# Example value, assuming Wf_cruise / W0 = 0.35 (hypothetical)
Wf_cruise_over_W0 = 0.35 # Example Value - **MUST BE UPDATED from Table 1**
if mission_fuel_fractions['warmup_takeoff'] * mission_fuel_fractions['climb'] == 0:
     raise ValueError("Warmup or Climb fuel fraction is zero, cannot calculate cruise fraction.")
W_cruise_frac = 1 - Wf_cruise_over_W0 / (mission_fuel_fractions['warmup_takeoff'] * mission_fuel_fractions['climb']) # Example calc

# Installation losses (Task C)
install_loss_factor = 0.02 # Phi_inlet + Phi_noz

# --- Task B: Max Allowable Uninstalled TSFC (S) ---
def calculate_max_tsfc(range_m, cruise_mach, cruise_alt_m, L_D, W_cruise_frac):
    """
    Calculates the maximum allowable TSFC (S) in mg/(N*s) using the Breguet Range Equation.
    W_cruise_frac = W_end_cruise / W_start_cruise
    """
    print("\n--- Task B: Maximum Allowable TSFC ---")
    if W_cruise_frac <= 0 or W_cruise_frac >= 1:
        print(f"Error: Invalid W_cruise_frac ({W_cruise_frac:.4f}) calculated. Check mission fuel fraction data and calculation.")
        return None

    a = get_speed_of_sound(cruise_alt_m)
    V = cruise_mach * a # Cruise speed in m/s

    # Breguet Range Eq: R = (V / (g0 * S)) * (L/D) * ln(W_start / W_end)
    # Where S is TSFC in kg/N-s (mass-based specific fuel consumption)
    # Or R = (V / c_t) * (L/D) * ln(W_start / W_end), where c_t = S*g0 is in 1/s (weight-based)
    # Let's calculate S (TSFC in kg/N-s) directly:
    # S = (V / (g0 * R)) * (L/D) * ln(W_start / W_end)
    # S = (V / (g0 * R)) * (L/D) * ln(1 / W_cruise_frac)
    
    tsfc_kg_N_s = (V / (g0 * range_m)) * L_D * np.log(1.0 / W_cruise_frac) 
    
    # Convert to milligrams per Newton-second (mg/Ns)
    tsfc_mg_N_s = tsfc_kg_N_s * 1000 * 1000 # kg to mg

    print(f"Cruise Speed (V): {V:.2f} m/s")
    print(f"Cruise Fuel Weight Fraction (W_end / W_start): {W_cruise_frac:.4f}")
    print(f"Maximum Allowable Uninstalled TSFC (S): {tsfc_kg_N_s:.4e} kg/(N*s)")
    print(f"Maximum Allowable Uninstalled TSFC (S): {tsfc_mg_N_s:.2f} mg/(N*s)")
    
    return tsfc_mg_N_s

# --- Task C: Minimum Uninstalled Specific Thrust ---
def calculate_min_specific_thrust(cruise_alt_m, cruise_mach, W_cruise_kg, L_D, loss_factor, inlet_dia_m):
    """
    Calculates the minimum required uninstalled specific thrust (F/m_dot) in (N*s)/kg.
    Uses steady level cruise assumption (Thrust = Drag).
    """
    print(f"\n--- Task C: Minimum Specific Thrust (Inlet Diameter: {inlet_dia_m} m) ---")
    
    # Atmospheric conditions at cruise
    rho = get_isa_density(cruise_alt_m)
    a = get_speed_of_sound(cruise_alt_m)
    T = get_isa_temp(cruise_alt_m)
    P = get_isa_pressure(cruise_alt_m)
    print(f"Cruise Alt Conditions: rho={rho:.4f} kg/m^3, a={a:.2f} m/s, T={T:.2f} K, P={P:.2f} Pa")

    # Cruise speed
    V = cruise_mach * a
    print(f"Cruise Speed (V): {V:.2f} m/s")
    
    # Aircraft weight in Newtons
    W_cruise_N = W_cruise_kg * g0
    print(f"Assumed Cruise Weight: {W_cruise_kg:.0f} kg ({W_cruise_N:.0f} N)")

    # Required INSTALLED thrust for steady level cruise (Thrust = Drag)
    # Drag = Weight / (L/D)
    # Check manual.md.txt Task C.1: "Determine the required installed thrust to *attain* the cruise condition, using Eq. (1.28)."
    # Check manual.md.txt Task A.3: "It attains an initial altitude of 11 km at beginning of cruise (Ps=1.5 m/s)."
    # If Eq (1.28) incorporates climb rate (Ps), then Thrust = Drag + W*Ps/V
    Ps_initial_cruise = 1.5 # m/s (from Task A.3)
    # Thrust_installed_N = (W_cruise_N / L_D) # Level cruise T=D
    Thrust_installed_N = (W_cruise_N / L_D) + (W_cruise_N * Ps_initial_cruise / V) # Thrust for initial climb T=D+W*Ps/V
    print(f"Required Installed Thrust (T=D + W*Ps/V for Ps={Ps_initial_cruise} m/s): {Thrust_installed_N:.0f} N")
    
    # Required UNINSTALLED thrust (accounting for losses)
    # T_installed = T_uninstalled * (1 - loss_factor) => T_uninstalled = T_installed / (1 - loss_factor)
    thrust_uninstalled_N = Thrust_installed_N / (1.0 - loss_factor)
    print(f"Required Uninstalled Thrust (Total for Aircraft): {thrust_uninstalled_N:.0f} N")
    
    # Calculate max inlet mass flow rate (m_dot_max = rho * V * A_inlet)
    # This is total flow for the specified inlet diameter. If diameter is per engine, adjust A_inlet.
    # Assuming inlet_dia_m refers to a single engine inlet (consistent with F/mdot per engine)
    A_inlet_per_engine = np.pi * (inlet_dia_m / 2.0)**2
    m_dot_max_per_engine_kgs = rho * V * A_inlet_per_engine
    print(f"Inlet Area per Engine (A_inlet): {A_inlet_per_engine:.3f} m^2")
    print(f"Maximum Inlet Mass Flow per Engine (m_dot_air): {m_dot_max_per_engine_kgs:.2f} kg/s")
    
    # Calculate minimum uninstalled specific thrust (F / m_dot) per engine
    num_engines = 2
    thrust_uninstalled_per_engine_N = thrust_uninstalled_N / num_engines
    
    if m_dot_max_per_engine_kgs > 0:
        specific_thrust_N_s_kg = thrust_uninstalled_per_engine_N / m_dot_max_per_engine_kgs
        print(f"Required Uninstalled Thrust per Engine: {thrust_uninstalled_per_engine_N:.0f} N")
        # print(f"Maximum Inlet Mass Flow per Engine: {m_dot_max_per_engine_kgs:.2f} kg/s") # Already printed above
        print(f"Minimum Uninstalled Specific Thrust (F_uninstalled / m_dot_air): {specific_thrust_N_s_kg:.2f} N*s/kg")
        return specific_thrust_N_s_kg
    else:
        print("Error: Mass flow calculation resulted in zero or negative value.")
        return None

# --- Main Execution ---
if __name__ == "__main__":
    print("HP-1 Engine Preliminary Analysis")
    print("="*30)
    print(f"Inputs: Cruise Mach={cruise_mach}, Alt={cruise_alt_m}m, Range={range_km}km")
    print(f"Assumed Cruise L/D = {L_D_cruise} (Update Needed!)")
    print(f"Assumed Cruise Weight = {W_cruise_kg} kg (Update Needed!)")
    print(f"Assumed Cruise Fuel Fraction (W_end/W_start) = {W_cruise_frac:.4f} (Update Needed based on Table 1!)")
    print(f"Installation Loss Factor = {install_loss_factor}")
    
    try:
        # Execute Task B
        max_tsfc = calculate_max_tsfc(range_m, cruise_mach, cruise_alt_m, L_D_cruise, W_cruise_frac)
        
        # Execute Task C for baseline and other diameters
        base_inlet_diameter_m = 2.2
        min_spec_thrust_base = calculate_min_specific_thrust(cruise_alt_m, cruise_mach, W_cruise_kg, L_D_cruise, install_loss_factor, base_inlet_diameter_m)
        
        other_diameters_m = [2.5, 2.75, 3.0, 3.25, 3.5]
        for dia in other_diameters_m:
            calculate_min_specific_thrust(cruise_alt_m, cruise_mach, W_cruise_kg, L_D_cruise, install_loss_factor, dia)

    except ValueError as e:
        print(f"\n!!! Calculation Error: {e}")
    except Exception as e:
        print(f"\n!!! An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()


    print("\nScript Finished.")
    print("="*30)
    print("REMINDER: Ensure placeholder values (W_cruise_kg, L_D_cruise, Wf_cruise_over_W0/mission_fuel_fractions)")
    print("are updated based on data from manual.md.txt (Table 1) and book.md.txt (Example 1.2).")
    print("The calculation for W_cruise_frac needs careful implementation based on Table 1 format.")
    print("Installed Thrust calculation now includes initial cruise climb rate (Ps=1.5 m/s) as per Task A.3 & C.1 suggestion.")