import numpy as np
import matplotlib.pyplot as plt
from src.scaling import log_linear_regression, plot_scaling
from src.thermodynamics import plot_prob, exp_E, exp_E_squared, C_v_calc, plot_C_V, prob_calc

# Simulation data for the end-effector radius (r_ee) as a function of the number of particles (N)
N = np.array([20, 40, 80, 160])
r_ee = np.sqrt(np.array([66.7, 193, 549, 1555]))



def main():
#initializing data for the scaling analysis
    N = np.array([20, 40, 80, 160])
    r_ee = np.sqrt(np.array([66.7, 193, 549, 1555]))

#initializing data for the thermodynamics analysis
    k = eps = 1
    E_min = -13*eps
    E_eps = np.array([0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13])
    g_E = np.array([18671059783.5, 15687265041,5351538782, 1222946058, 234326487, 40339545, 5824861, 710407, 77535, 9046, 645, 86, 0, 1])
    T = np.linspace(0.1, 1, 100
                    )
#performing the scaling analysis
    lin_fit = log_linear_regression(N, r_ee)
    plot_scaling(*lin_fit)

#performing the thermodynamics analysis
    plot_C_V(E_eps, g_E, T, k, eps, E_min)
    plot_prob(E_eps, g_E, T, k, eps, E_min)


if __name__ == "__main__":
    main()