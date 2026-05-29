# geometry.py
# Generates a NACA 4-digit aerofoil profile and a 2D Cartesian mesh.

import numpy as np

def naca4_profile(number: str = "0012", n_points: int = 100) -> tuple:
    m = int(number[0]) / 100
    p = int(number[1]) / 10
    t = int(number[2:]) / 100

    x = np.linspace(0, 1, n_points)

    yt = 5 * t * (0.2969*np.sqrt(x)
                - 0.1260*x
                - 0.3516*x**2
                + 0.2843*x**3
                - 0.1015*x**4)

    if m == 0 or p == 0:
        yc = np.zeros_like(x)
        dyc_dx = np.zeros_like(x)
    else:
        yc = np.where(
            x < p,
            (m / p**2) * (2*p*x - x**2),
            (m / (1-p)**2) * ((1-2*p) + 2*p*x - x**2)
        )
        dyc_dx = np.where(
            x < p,
            (2*m / p**2) * (p - x),
            (2*m / (1-p)**2) * (p - x)
        )

    theta = np.arctan(dyc_dx)

    xu = x  - yt * np.sin(theta)
    xl = x  + yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    yl = yc - yt * np.cos(theta)

    return xu, xl, yu, yl, yc


def build_mesh(nx: int = 80, ny: int = 40) -> tuple:
    """
    Returns (X, Y) meshgrid for a unit chord blade cross-section.
    The mesh spans x=[0,1], y=[-0.15, 0.15] — covers typical NACA profiles.
    """
    x = np.linspace(0, 1, nx)
    y = np.linspace(-0.15, 0.15, ny)
    return np.meshgrid(x, y)


def blade_mask(X: np.ndarray, Y: np.ndarray, naca: str = "0012") -> np.ndarray:
    """
    Returns a boolean mask — True inside the blade, False outside.
    Used by the solver to apply boundary conditions correctly.
    """
    _, _, yu, yl, _ = naca4_profile(naca, n_points=X.shape[1])
    x_vals = np.linspace(0, 1, X.shape[1])

    yu_interp = np.interp(X[0], x_vals, yu)
    yl_interp = np.interp(X[0], x_vals, yl)

    mask = np.zeros_like(X, dtype=bool)
    for j in range(X.shape[1]):
        mask[:, j] = (Y[:, j] >= yl_interp[j]) & (Y[:, j] <= yu_interp[j])

    return mask


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from geometry import naca4_profile, build_mesh, blade_mask

    # Plot the aerofoil profile
    xu, xl, yu, yl, yc = naca4_profile("0012")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Left: aerofoil outline
    axes[0].plot(xu, yu, 'b-', linewidth=2, label='Upper surface')
    axes[0].plot(xl, yl, 'r-', linewidth=2, label='Lower surface')
    axes[0].plot(xu, yc, 'g--', linewidth=1, label='Camber line')
    axes[0].set_aspect('equal')
    axes[0].legend()
    axes[0].set_title("NACA 0012 — aerofoil profile")
    axes[0].set_xlabel("Chord (normalised)")
    axes[0].grid(True, alpha=0.3)

    # Right: mesh mask (shows which cells are inside the blade)
    X, Y = build_mesh(80, 40)
    mask = blade_mask(X, Y, "0012")
    axes[1].contourf(X, Y, mask.astype(float), cmap='Blues')
    axes[1].set_aspect('equal')
    axes[1].set_title("2D mesh — blade interior (blue)")
    axes[1].set_xlabel("Chord (normalised)")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()