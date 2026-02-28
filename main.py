import numpy as np

from src.scaling import log_linear_regression, plot_scaling
from src.thermodynamics import plot_C_V, plot_prob


def main() -> None:
    """Run scaling and thermodynamics analyses and save plots."""
    # ---------- Scaling (SAW) ----------
    N = np.array([20, 40, 80, 160], dtype=float)
    r_ee = np.sqrt(np.array([66.7, 193.0, 549.0, 1555.0], dtype=float))  # r_ee/a

    lin_fit = log_linear_regression(N, r_ee)
    plot_scaling(*lin_fit)

    # ---------- Thermodynamics (density of states) ----------
    k = 1.0
    eps = 1.0
    E_min = -13.0 * eps

    # Energies in units of eps (E/eps)
    E_eps = np.array([0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13], dtype=float)
    g_E = np.array(
        [18671059783.5, 15687265041, 5351538782, 1222946058, 234326487, 40339545,
         5824861, 710407, 77535, 9046, 645, 86, 0, 1],
        dtype=float
    )

    T = np.linspace(0.1, 1.0, 100, dtype=float)

    plot_C_V(E_eps, g_E, T, k, eps, E_min)
    plot_prob(E_eps, g_E, T, k, eps, E_min)


if __name__ == "__main__":
    main()