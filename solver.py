# solver.py
# 2D steady-state heat conduction solver using the finite difference method.
# Solves the Laplace equation: d²T/dx² + d²T/dy² = 0
# Boundary conditions mimic aeroengine operating environment:
#   - Leading edge (hot gas impingement): T_hot
#   - Trailing edge and surfaces: T_cool (cooling air)
#   - Interior: solved iteratively via Gauss-Seidel relaxation

import numpy as np
from geometry import build_mesh, blade_mask

def solve_temperature(
        material: dict,
        T_hot: float = 1200.0,  # °C — hot gas temperature at leading edge
        T_cool: float = 400.0,  # °C — cooling air temperature
        nx: int = 80,
        ny: int = 40,
        naca: str = "0012",
        max_iter: int = 5000,
        tol: float = 1e-4,
        h_conv: float = 1000.0  # Added to pass the UI's cooling coefficient
) -> tuple:
    """
    Solves 2D steady-state heat conduction across a turbine blade cross-section
    using material-dependent convective boundary conditions.
    """
    X, Y = build_mesh(nx, ny)
    mask = blade_mask(X, Y, naca)

    k = material["k"]
    # Calculate real spatial steps based on mesh boundaries (x: [0,1], y: [-0.15, 0.15])
    dx = 1.0 / (nx - 1)
    dy = 0.30 / (ny - 1)

    # Initialise temperature field
    T = np.full((ny, nx), T_cool)

    # Fixed boundary conditions for baseline constraints
    T[:, 0] = T_hot  # Leading edge impingement
    T[:, -1] = T_cool  # Trailing edge drop

    # Gauss-Seidel iterative solver with Convective (Robin) Boundary Conditions
    iterations = 0
    for iteration in range(max_iter):
        T_old = T.copy()

        for j in range(1, nx - 1):
            col = mask[:, j]
            if not col.any():
                continue

            # Find the true top and bottom boundary cells of the blade profile
            idx_bottom = np.argmax(col)
            idx_top = ny - 1 - np.argmax(col[::-1])

            # 1. Update strictly interior cells using the 5-point stencil
            for i in range(idx_bottom + 1, idx_top):
                T[i, j] = 0.25 * (
                        T[i + 1, j] + T[i - 1, j] +
                        T[i, j + 1] + T[i, j - 1]
                )

            # 2. Dynamically calculate surface temperatures based on material conductivity (k)
            # Top surface node balanced against convective cooling fluid
            T[idx_top, j] = (k * T[idx_top - 1, j] + h_conv * dy * T_cool) / (k + h_conv * dy)

            # Bottom surface node balanced against convective cooling fluid
            T[idx_bottom, j] = (k * T[idx_bottom + 1, j] + h_conv * dy * T_cool) / (k + h_conv * dy)

        iterations = iteration + 1
        delta = np.max(np.abs(T - T_old))
        if delta < tol:
            break

    return X, Y, T, mask, iterations


def compute_heat_flux(T: np.ndarray, k: float, dx: float = 1/80, dy: float = 1/40) -> tuple:
    """
    Computes heat flux vectors from the temperature gradient.
    q = -k * ∇T

    Returns qx, qy — heat flux components [W/m²]
    """
    dTdy, dTdx = np.gradient(T, dy, dx)
    qx = -k * dTdx
    qy = -k * dTdy
    return qx, qy

if __name__ == "__main__":
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from materials import get_material

    mat = get_material("Inconel 718")

    print("Running solver — this may take 10–20 seconds...")
    X, Y, T, mask, iters = solve_temperature(
        material=mat,
        T_hot=1200,
        T_cool=400,
        nx=80,
        ny=40
    )
    print(f"Converged in {iters} iterations")
    print(f"Max temperature: {T[mask].max():.1f}°C")
    print(f"Min temperature: {T[mask].min():.1f}°C")

    T_plot = np.where(mask, T, np.nan)

    fig, ax = plt.subplots(figsize=(10, 5))
    cf = ax.contourf(X, Y, T_plot, levels=50, cmap='hot')
    ax.contour(X, Y, mask.astype(float), levels=[0.5], colors='white', linewidths=1.5)
    plt.colorbar(cf, ax=ax, label='Temperature (°C)')
    ax.set_aspect('equal')
    ax.set_title("Turbine blade — 2D steady-state thermal distribution (Inconel 718)")
    ax.set_xlabel("Chord (normalised)")
    ax.set_ylabel("Span (normalised)")
    plt.tight_layout()
    plt.show(block=True)

