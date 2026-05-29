# visualiser.py
# Publication-quality visualisation for turbine blade thermal analysis.
# Produces: heat maps, material comparisons, and heat flux vector overlays.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from solver import solve_temperature, compute_heat_flux
from materials import get_material, MATERIALS

# Custom colourmap mimicking real thermal imaging cameras
THERMAL_CMAP = LinearSegmentedColormap.from_list(
    "thermal",
    ["#0d0221", "#1a1a8c", "#c0392b", "#e67e22", "#f1c40f", "#ffffff"],
    N=256
)

def plot_single(material_name: str, T_hot: float = 1200, T_cool: float = 400,
                ax=None, show: bool = True):
    """Single blade heat map for a given material."""
    mat = get_material(material_name)
    X, Y, T, mask, iters = solve_temperature(mat, T_hot, T_cool)
    T_plot = np.where(mask, T, np.nan)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.get_figure()

    cf = ax.contourf(X, Y, T_plot, levels=60, cmap=THERMAL_CMAP,
                     vmin=T_cool, vmax=T_hot)
    ax.contour(X, Y, mask.astype(float), levels=[0.5],
               colors='white', linewidths=1.0, alpha=0.6)
    plt.colorbar(cf, ax=ax, label='Temperature (°C)', shrink=0.8)
    ax.set_aspect('equal')
    ax.set_title(f"{material_name}  |  T_hot={T_hot}°C  T_cool={T_cool}°C  "
                 f"|  converged: {iters} iterations", fontsize=10)
    ax.set_xlabel("Chord (normalised)")
    ax.set_ylabel("Span")
    ax.set_facecolor('#000000')

    if show:
        plt.tight_layout()
        plt.show(block=True)

    return T, mask, iters


def plot_comparison(T_hot: float = 1200, T_cool: float = 400):
    """
    Side-by-side comparison of all materials in the library.
    The headline output for the Streamlit app and GitHub README.
    """
    names = list(MATERIALS.keys())
    fig, axes = plt.subplots(len(names), 1,
                             figsize=(12, 4 * len(names)))
    fig.suptitle(
        f"Turbine Blade Thermal Analysis — Material Comparison\n"
        f"T_hot = {T_hot}°C  |  T_cool = {T_cool}°C  |  NACA 0012",
        fontsize=13, fontweight='bold', y=1.01
    )

    results = {}
    for ax, name in zip(axes, names):
        print(f"  Solving for {name}...")
        T, mask, iters = plot_single(name, T_hot, T_cool, ax=ax, show=False)
        results[name] = {"T_max": T[mask].max(), "T_min": T[mask].min(),
                         "iters": iters}

    plt.tight_layout()
    plt.savefig("material_comparison.png", dpi=150, bbox_inches='tight')
    print("Saved: material_comparison.png")
    plt.show(block=True)
    return results


def plot_heat_flux(material_name: str, T_hot: float = 1200, T_cool: float = 400):
    """
    Heat map with heat flux vector arrows overlaid.
    Shows direction and magnitude of heat flow through the blade.
    """
    mat = get_material(material_name)
    X, Y, T, mask, iters = solve_temperature(mat, T_hot, T_cool)
    qx, qy = compute_heat_flux(T, mat["k"])

    T_plot = np.where(mask, T, np.nan)

    # Subsample arrows for readability
    step = 4
    Xq = X[::step, ::step]
    Yq = Y[::step, ::step]
    Uq = np.where(mask[::step, ::step], qx[::step, ::step], np.nan)
    Vq = np.where(mask[::step, ::step], qy[::step, ::step], np.nan)

    fig, ax = plt.subplots(figsize=(12, 5))
    cf = ax.contourf(X, Y, T_plot, levels=60, cmap=THERMAL_CMAP,
                     vmin=T_cool, vmax=T_hot)
    ax.quiver(Xq, Yq, Uq, Vq, color='white', alpha=0.5,
              scale=5e6, width=0.002)
    ax.contour(X, Y, mask.astype(float), levels=[0.5],
               colors='cyan', linewidths=1.0, alpha=0.5)
    plt.colorbar(cf, ax=ax, label='Temperature (°C)')
    ax.set_aspect('equal')
    ax.set_title(f"Heat flux vectors — {material_name}  (arrows show heat flow direction)",
                 fontsize=11)
    ax.set_xlabel("Chord (normalised)")
    ax.set_facecolor('#0d0221')
    plt.tight_layout()
    plt.show(block=True)


if __name__ == "__main__":
    print("=== Turbine Blade Thermal Visualiser ===\n")

    print("1. Single blade — Inconel 718")
    plot_single("Inconel 718", T_hot=1200, T_cool=400)

    print("\n2. Material comparison (all alloys)...")
    results = plot_comparison(T_hot=1200, T_cool=400)
    print("\nSummary:")
    for name, r in results.items():
        print(f"  {name:25s}  T_max={r['T_max']:.1f}°C  iters={r['iters']}")

    print("\n3. Heat flux overlay — Inconel 718")
    plot_heat_flux("Inconel 718")