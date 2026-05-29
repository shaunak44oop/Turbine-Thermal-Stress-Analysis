# app.py
# Streamlit web application — Turbine Blade Thermal & Stress Simulator
# Run: streamlit run "Turbine Heat Simulation/app.py"

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import streamlit as st
from solver import solve_temperature, compute_heat_flux
from materials import (get_material, MATERIALS,
                       thermal_diffusivity, thermal_stress_index, biot_number)

THERMAL_CMAP = LinearSegmentedColormap.from_list(
    "thermal",
    ["#0d0221", "#1a1a8c", "#c0392b", "#e67e22", "#f1c40f", "#ffffff"],
    N=256
)

st.set_page_config(
    page_title="Turbine Blade Thermal & Stress Simulator",
    page_icon="𖣐", layout="wide"
)

st.title("𖣐 Turbine Blade Thermal & Stress Simulator")
st.markdown("""
**2D steady-state heat conduction and thermal stress analysis across a NACA aerofoil.**  
Modelled after real HPT operating conditions for GE Aerospace / Safran Pratt & Whitney components.  
`Solver: Finite Difference Method (Gauss-Seidel)`
`Geometry: NACA 4-digit aerofoil`
`Physics: Fourier heat conduction + thermal stress indexing`
""")

# ── Sidebar ──────────────────────────────────────────────────────────────────
st.sidebar.header("⚙️ Simulation Parameters")

material_name = st.sidebar.selectbox("Blade material", list(MATERIALS.keys()))
mat = get_material(material_name)

st.sidebar.markdown(f"> **{mat['category']}** — {mat['use_case']}")

T_hot = st.sidebar.slider("Hot gas temp °C — leading edge",
                           800, 1600, 1200, 50)
T_cool = st.sidebar.slider("Cooling air temp °C — blade surfaces",
                            100, 600, 400, 50)
h_conv = st.sidebar.slider("Film cooling coefficient h (W/m²K)",
                            500, 5000, 1000, 100)
naca = st.sidebar.selectbox("Aerofoil profile (NACA)", ["0012","2412","4412"])

st.sidebar.markdown("---")
show_flux    = st.sidebar.checkbox("Show heat flux vectors")
show_compare = st.sidebar.checkbox("Compare all materials")
show_stress  = st.sidebar.checkbox("Show thermal stress analysis")

# Warn if material near its service limit
if T_hot > mat["T_max"] * 0.9:
    st.sidebar.warning(
        f"⚠️ T_hot approaches {material_name}'s service limit "
        f"({mat['T_max']}°C). Consider a higher-temp alloy."
    )

# ── Solve ────────────────────────────────────────────────────────────────────
with st.spinner(f"Solving heat equation for {material_name}..."):


   X, Y, T, mask, iters = solve_temperature(
       mat, T_hot, T_cool, nx=80, ny=40, naca=naca, h_conv=h_conv
   )

T_max_val  = T[mask].max()
T_min_val  = T[mask].min()
delta_T    = T_max_val - T_min_val
alpha      = thermal_diffusivity(mat)
sigma_idx  = thermal_stress_index(mat, delta_T)
bi         = biot_number(mat, h=h_conv)
T_plot     = np.where(mask, T, np.nan)

# ── Metrics row ──────────────────────────────────────────────────────────────
c1,c2,c3,c4,c5,c6 = st.columns(6)
c1.metric("T_max",            f"{T_max_val:.0f}°C")
c2.metric("T_min",            f"{T_min_val:.0f}°C")
c3.metric("ΔT across blade",  f"{delta_T:.0f}°C")
c4.metric("Thermal diffusivity α", f"{alpha*1e6:.2f} ×10⁻⁶ m²/s")
c5.metric("Thermal stress index",  f"{sigma_idx/1e6:.0f} MPa")
c6.metric("Biot number Bi",   f"{bi:.3f}")

st.markdown("---")

# ── Main heat map ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 4.5))
cf = ax.contourf(X, Y, T_plot, levels=80, cmap=THERMAL_CMAP,
                 vmin=T_cool, vmax=T_hot)
ax.contour(X, Y, mask.astype(float), levels=[0.5],
           colors='white', linewidths=1.0, alpha=0.5)

if show_flux:
    qx, qy = compute_heat_flux(T, mat["k"])
    step = 4
    Xq,Yq = X[::step,::step], Y[::step,::step]
    Uq = np.where(mask[::step,::step], qx[::step,::step], np.nan)
    Vq = np.where(mask[::step,::step], qy[::step,::step], np.nan)
    ax.quiver(Xq, Yq, Uq, Vq, color='white', alpha=0.45,
              scale=5e6, width=0.002)

cb = plt.colorbar(cf, ax=ax, label='Temperature (°C)', pad=0.01)
cb.ax.yaxis.label.set_color('white')
cb.ax.tick_params(colors='white')
ax.set_aspect('equal')
ax.set_title(f"{material_name} ({mat['category']}) — NACA {naca} | "
             f"T_hot={T_hot}°C  T_cool={T_cool}°C  Bi={bi:.3f}",
             color='white')
ax.set_xlabel("Chord (normalised)", color='white')
ax.set_ylabel("Span", color='white')
ax.tick_params(colors='white')
ax.set_facecolor('#0d0221')
fig.patch.set_facecolor('#0d0221')
st.pyplot(fig)
plt.close()

# ── Thermal stress analysis ───────────────────────────────────────────────────
if show_stress:
    st.markdown("---")
    st.subheader("🔩 Thermal Stress Analysis")
    st.markdown("""
    Thermal stress arises when differential heating causes constrained expansion.
    Estimated using: **σ = E · α · ΔT** (linear elastic approximation)
    """)

    names  = list(MATERIALS.keys())
    sigmas = []
    alphas = []
    for n in names:
        m = get_material(n)

        _, _, T2, mask2, _ = solve_temperature(m, T_hot, T_cool, nx=80, ny=40, h_conv=h_conv)
        dT = T2[mask2].max() - T2[mask2].min()
        sigmas.append(thermal_stress_index(m, dT) / 1e6)
        alphas.append(thermal_diffusivity(m) * 1e6)

    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))
    fig2.patch.set_facecolor('#0d0221')

    colors = [MATERIALS[n]["color"] for n in names]

    # Thermal stress bar chart
    bars = ax1.barh(names, sigmas, color=colors, edgecolor='white', linewidth=0.5)
    ax1.set_xlabel("Thermal Stress Index (MPa)", color='white')
    ax1.set_title("Thermal Stress Index — σ = E·α·ΔT", color='white')
    ax1.set_facecolor('#0d0221')
    ax1.tick_params(colors='white')
    ax1.spines['bottom'].set_color('white')
    ax1.spines['left'].set_color('white')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    # Highlight current material
    idx = names.index(material_name)
    bars[idx].set_edgecolor('yellow')
    bars[idx].set_linewidth(2.5)

    # Thermal diffusivity bar chart
    bars2 = ax2.barh(names, alphas, color=colors, edgecolor='white', linewidth=0.5)
    ax2.set_xlabel("Thermal Diffusivity α (×10⁻⁶ m²/s)", color='white')
    ax2.set_title("Thermal Diffusivity — α = k/ρCp", color='white')
    ax2.set_facecolor('#0d0221')
    ax2.tick_params(colors='white')
    ax2.spines['bottom'].set_color('white')
    ax2.spines['left'].set_color('white')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    bars2[idx].set_edgecolor('yellow')
    bars2[idx].set_linewidth(2.5)

    plt.tight_layout()
    st.pyplot(fig2)
    plt.close()

    # Material suitability table
    st.markdown("#### Material Suitability Summary")
    rows = []
    for n in names:
        m = get_material(n)
        suitable = "✅" if T_hot <= m["T_max"] else "❌ Exceeds limit"
        bi_val = biot_number(m, h=h_conv)
        rows.append({
            "Material": n,
            "Category": m["category"],
            "Max Service Temp": f"{m['T_max']}°C",
            "Suitable at T_hot": suitable,
            "Biot Number": f"{bi_val:.3f}",
            "Use Case": m["use_case"]
        })
    st.table(rows)

# ── Material comparison heat maps ─────────────────────────────────────────────
if show_compare:
    st.markdown("---")
    st.subheader("📊 Full Material Comparison — Thermal Distribution")
    names = list(MATERIALS.keys())
    fig3, axes = plt.subplots(len(names), 1, figsize=(13, 3.5*len(names)))
    fig3.patch.set_facecolor('#0d0221')

    for ax3, name in zip(axes, names):
        m = get_material(name)

        _, _, T3, mask3, iters3 = solve_temperature(
            m, T_hot, T_cool, nx=80, ny=40, naca=naca, h_conv=h_conv)
        T3_plot = np.where(mask3, T3, np.nan)
        cf3 = ax3.contourf(X, Y, T3_plot, levels=60, cmap=THERMAL_CMAP,
                           vmin=T_cool, vmax=T_hot)
        ax3.contour(X, Y, mask3.astype(float), levels=[0.5],
                    colors='white', linewidths=0.8, alpha=0.4)
        plt.colorbar(cf3, ax=ax3, label='°C').ax.tick_params(colors='white')
        ax3.set_aspect('equal')
        ax3.set_facecolor('#0d0221')
        dT3 = T3[mask3].max() - T3[mask3].min()
        si  = thermal_stress_index(m, dT3) / 1e6
        ax3.set_title(
            f"{name} | k={m['k']} W/m·K | "
            f"T_max={T3[mask3].max():.0f}°C | "
            f"σ_idx={si:.0f} MPa | iters={iters3}",
            color='white', fontsize=9)
        ax3.tick_params(colors='white')

    plt.tight_layout()
    st.pyplot(fig3)
    plt.close()