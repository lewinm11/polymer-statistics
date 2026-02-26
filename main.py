import numpy as np
import matplotlib.pyplot as plt
from src.scaling import log_linear_regression, plot_scaling

# Simulation data for the end-effector radius (r_ee) as a function of the number of particles (N)
N = np.array([20, 40, 80, 160])
r_ee = np.sqrt(np.array([66.7, 193, 549, 1555]))



def main():
    N = np.array([20, 40, 80, 160])
    r_ee = np.sqrt(np.array([66.7, 193, 549, 1555]))
    lin_fit = log_linear_regression(N, r_ee)
    plot_scaling(*lin_fit)


if __name__ == "__main__":
    main()