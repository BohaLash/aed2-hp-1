import numpy as np
from scipy.optimize import newton

# Constants
g = 9.81  # m/s^2
R = 287  # J/(kg·K) for air
gamma = 1.4  # ratio of specific heats
cp = 1.004  # kJ/(kg·K)
hPR = 42800  # kJ/kg, fuel heating value

# Aircraft parameters from the problem
W_start_cruise = 1577940  # N, starting cruise weight
M_cruise = 0.83
cruise_alt_start = 11000  # m, 11 km
cruise_distance = 10650  # km
payload = 22770  # kg
fuel_used_cruise = 50240  # kg
fuel_used_total = 61200  # kg
S = 427.8  # m^2, wing area (from Example 1.2)
engine_thrust_takeoff = 214000  # N per engine (from problem)


# Atmosphere model functions
def atm_properties(h):
    """Calculate atmospheric properties at altitude h (m)"""
    if h <= 11000:
        # Troposphere
        T0 = 288.16  # K
        P0 = 101325  # Pa
        a = -0.0065  # K/m
        T = T0 + a * h
        P = P0 * (T / T0) ** (-g / (a * R))
    else:
        # Lower stratosphere
        T = 216.66  # K
        P = 22632 * np.exp(-g * (h - 11000) / (R * T))

    rho = P / (R * T)
    a_sound = np.sqrt(gamma * R * T)
    return T, P, rho, a_sound


# Task B1: Cruise TSFC calculations
def cruise_tsfc_calculations():
    """Calculate allowable TSFC for cruise conditions"""
    # Part (a): Cruise climb (constant CD/CL)
    # Using Breguet range equation: R = (V/SFC)*(L/D)*ln(Wi/Wf)
    V = M_cruise * atm_properties(cruise_alt_start)[3]  # cruise speed m/s
    W_end_cruise = W_start_cruise - fuel_used_cruise * g  # N

    # Convert cruise distance to meters
    R_meters = cruise_distance * 1000

    # Calculate average L/D (assuming constant CD/CL means constant L/D)
    L_D_avg = (R_meters * (cp * V)) / (V * np.log(W_start_cruise / W_end_cruise))

    # Allowable TSFC for cruise climb
    TSFC_cruise_climb = (V / (L_D_avg * R_meters)) * np.log(
        W_start_cruise / W_end_cruise
    )

    # Calculate final altitude for cruise climb (constant CL)
    # At constant CL, W/rho = constant
    rho_start = atm_properties(cruise_alt_start)[2]
    rho_end = rho_start * (W_end_cruise / W_start_cruise)

    # Solve for altitude with this density
    def alt_eqn(h):
        return atm_properties(h)[2] - rho_end

    h_end = newton(alt_eqn, cruise_alt_start)

    # Part (b): Constant altitude cruise
    CL_start = W_start_cruise / (0.5 * rho_start * V**2 * S)
    CL_end = W_end_cruise / (0.5 * rho_start * V**2 * S)

    # Assuming parabolic drag polar: CD = CD0 + K*CL^2
    # From Example 1.2: CD0 = 0.015, K = 0.04
    CD0 = 0.015
    K = 0.04

    CD_start = CD0 + K * CL_start**2
    CD_end = CD0 + K * CL_end**2

    L_D_start = CL_start / CD_start
    L_D_end = CL_end / CD_end
    L_D_avg_const_alt = (L_D_start + L_D_end) / 2

    TSFC_const_alt = (V / (L_D_avg_const_alt * R_meters)) * np.log(
        W_start_cruise / W_end_cruise
    )

    return {
        "cruise_climb": {
            "TSFC": TSFC_cruise_climb * 1e5,  # convert to kg/(N·s) from kg/(N·m)
            "end_altitude": h_end,
            "L_D": L_D_avg,
        },
        "constant_altitude": {
            "TSFC": TSFC_const_alt * 1e5,
            "L_D_start": L_D_start,
            "L_D_end": L_D_end,
            "L_D_avg": L_D_avg_const_alt,
        },
    }


# Task B2: Loiter Mach numbers
def loiter_mach_numbers():
    """Calculate loiter Mach numbers at different altitudes"""
    altitudes = [10000, 9000, 8000, 7000, 6000]  # m
    W_loiter = 0.64 * (W_start_cruise + fuel_used_cruise * g)  # 64% of WTO

    mach_numbers = []
    for h in altitudes:
        T, P, rho, a = atm_properties(h)

        # For loiter, we want maximum endurance (max L/D)
        # Maximum L/D occurs when CD0 = K*CL^2
        CL_opt = np.sqrt(CD0 / K)
        CD_opt = 2 * CD0
        L_D_max = CL_opt / CD_opt

        # Required velocity for this CL
        V = np.sqrt(2 * W_loiter / (rho * S * CL_opt))
        M = V / a
        mach_numbers.append(M)

    return dict(zip(altitudes, mach_numbers))


# Task B3: Aircraft drag calculations
def calculate_drag_points():
    """Calculate drag at various flight points"""
    drag_points = {}

    # (a) Takeoff, M=0.23, sea level
    M = 0.23
    T, P, rho, a = atm_properties(0)
    V = M * a
    W_TO = W_start_cruise + fuel_used_total * g  # approx takeoff weight
    CL = W_TO / (0.5 * rho * V**2 * S)
    CD = CD0 + K * CL**2
    D = 0.5 * rho * V**2 * S * CD
    drag_points["takeoff"] = D

    # (b) Start of cruise, M=0.83, 11 km
    M = 0.83
    T, P, rho, a = atm_properties(11000)
    V = M * a
    CL = W_start_cruise / (0.5 * rho * V**2 * S)
    CD = CD0 + K * CL**2
    D = 0.5 * rho * V**2 * S * CD
    drag_points["start_cruise"] = D

    # (c) End of cruise climb - use results from B1a
    h_end = cruise_tsfc_calculations()["cruise_climb"]["end_altitude"]
    T, P, rho, a = atm_properties(h_end)
    V = M * a
    W_end = W_start_cruise - fuel_used_cruise * g
    CL = W_end / (0.5 * rho * V**2 * S)
    CD = CD0 + K * CL**2
    D = 0.5 * rho * V**2 * S * CD
    drag_points["end_cruise_climb"] = D

    # (d) End of 11-km cruise
    # Same as start but with lower weight
    W_end_constant = W_start_cruise - fuel_used_cruise * g
    CL = W_end_constant / (0.5 * rho * V**2 * S)
    CD = CD0 + K * CL**2
    D = 0.5 * rho * V**2 * S * CD
    drag_points["end_constant_cruise"] = D

    # (e) Engine out (88% of WTO), M=0.45, 5 km
    M = 0.45
    T, P, rho, a = atm_properties(5000)
    V = M * a
    W_engine_out = 0.88 * (W_start_cruise + fuel_used_cruise * g)  # 88% of WTO
    CL = W_engine_out / (0.5 * rho * V**2 * S)
    CD = CD0 + K * CL**2
    D = 0.5 * rho * V**2 * S * CD
    drag_points["engine_out"] = D

    return drag_points


if __name__ == "__main__":
    print("HP-1 Aircraft Engine Design Calculations\n")

    # Task B1 calculations
    print("Task B1: Cruise TSFC Calculations")
    tsfc_results = cruise_tsfc_calculations()
    print(f"(a) Cruise climb:")
    print(f"  Allowable TSFC: {tsfc_results['cruise_climb']['TSFC']:.2e} kg/(N·s)")
    print(f"  End altitude: {tsfc_results['cruise_climb']['end_altitude']/1000:.1f} km")
    print(f"  Average L/D: {tsfc_results['cruise_climb']['L_D']:.1f}")

    print(f"\n(b) Constant altitude cruise:")
    print(f"  Allowable TSFC: {tsfc_results['constant_altitude']['TSFC']:.2e} kg/(N·s)")
    print(f"  L/D start: {tsfc_results['constant_altitude']['L_D_start']:.1f}")
    print(f"  L/D end: {tsfc_results['constant_altitude']['L_D_end']:.1f}")
    print(f"  Average L/D: {tsfc_results['constant_altitude']['L_D_avg']:.1f}")

    # Task B2 calculations
    print("\nTask B2: Loiter Mach Numbers at 64% WTO")
    loiter_mach = loiter_mach_numbers()
    for alt, M in loiter_mach.items():
        print(f"  {alt/1000:.0f} km: M = {M:.2f}")

    # Task B3 calculations
    print("\nTask B3: Aircraft Drag at Key Points")
    drag_points = calculate_drag_points()
    for point, drag in drag_points.items():
        print(f"  {point.replace('_', ' ')}: {drag/1000:.1f} kN")
