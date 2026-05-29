# materials.py
# Thermophysical properties for aerospace and engineering alloys.
# Sources: ASM Aerospace Specification Metals, Carpenter Technology,
#          MatWeb, MMPDS (Metallic Materials Properties Development and
#          Standardisation Handbook)

MATERIALS = {
    "Inconel 718": {
        "k":    11.4,
        "rho":  8190,
        "cp":   435,
        "E":    200e9,      # Young's modulus [Pa]
        "alpha_t": 13.0e-6, # thermal expansion coefficient [1/K]
        "T_max": 1200,
        "T_melt": 1336,
        "category": "Nickel Superalloy",
        "use_case": "HPT blades, combustor liners (GE9X, LEAP)",
        "color": "#4A90D9"
    },
    "Ti-6Al-4V": {
        "k":    7.2,
        "rho":  4430,
        "cp":   560,
        "E":    114e9,
        "alpha_t": 8.6e-6,
        "T_max": 600,
        "T_melt": 1660,
        "category": "Titanium Alloy",
        "use_case": "Fan blades, compressor stages (CFM56, GEnx)",
        "color": "#E8A838"
    },
    "Aluminium 6061-T6": {
        "k":    167.0,
        "rho":  2700,
        "cp":   896,
        "E":    68.9e9,
        "alpha_t": 23.6e-6,
        "T_max": 300,
        "T_melt": 652,
        "category": "Aluminium Alloy",
        "use_case": "Nacelle structures, non-hot-section components",
        "color": "#A8D8A8"
    },
    "Stainless Steel 316L": {
        "k":    16.3,
        "rho":  8000,
        "cp":   500,
        "E":    193e9,
        "alpha_t": 16.0e-6,
        "T_max": 870,
        "T_melt": 1400,
        "category": "Stainless Steel",
        "use_case": "Exhaust components, low-pressure turbine casings",
        "color": "#888888"
    },
    "Inconel 625": {
        "k":    9.8,
        "rho":  8440,
        "cp":   410,
        "E":    207e9,
        "alpha_t": 12.8e-6,
        "T_max": 1050,
        "T_melt": 1290,
        "category": "Nickel Superalloy",
        "use_case": "Combustor transition ducts, exhaust ducting",
        "color": "#7B68EE"
    },
    "Hastelloy X": {
        "k":    9.1,
        "rho":  8220,
        "cp":   461,
        "E":    197e9,
        "alpha_t": 13.9e-6,
        "T_max": 1175,
        "T_melt": 1315,
        "category": "Nickel Superalloy",
        "use_case": "Combustion chambers, afterburner components",
        "color": "#FF6B6B"
    },
    "Titanium CP Grade 2": {
        "k":    16.4,
        "rho":  4510,
        "cp":   520,
        "E":    105e9,
        "alpha_t": 8.9e-6,
        "T_max": 400,
        "T_melt": 1668,
        "category": "Titanium Alloy",
        "use_case": "Airframe brackets, low-temperature ducting",
        "color": "#FFA07A"
    },
}

def get_material(name: str) -> dict:
    if name not in MATERIALS:
        raise ValueError(f"Unknown material '{name}'. "
                         f"Available: {list(MATERIALS.keys())}")
    return MATERIALS[name]

def thermal_diffusivity(mat: dict) -> float:
    """α = k / (ρ · Cp)  [m²/s] — governs transient heat spread rate"""
    return mat["k"] / (mat["rho"] * mat["cp"])

def thermal_stress_index(mat: dict, delta_T: float) -> float:
    """
    Approximate thermal stress σ = E · α · ΔT  [Pa]
    Higher values indicate greater risk of thermal fatigue cracking.
    """
    return mat["E"] * mat["alpha_t"] * delta_T

def biot_number(mat: dict, h: float = 1000.0, L: float = 0.05) -> float:
    """
    Bi = hL/k — ratio of convective to conductive resistance.
    h: convective heat transfer coefficient [W/m²K], default: film cooling
    L: characteristic length [m], default: 5cm chord
    """
    return (h * L) / mat["k"]