import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def log_linear_regression(N, r_ee):
    """Return slope, intercept, r_value, std_err and log-transformed data."""
    log_N = np.log(N)
    log_r_ee = np.log(r_ee)

    slope, intercept, r_value, _, std_err = stats.linregress(log_N, log_r_ee)
    return slope, intercept, r_value, std_err, log_N, log_r_ee


def plot_scaling(slope, intercept, r_value, std_err, log_N, log_r_ee):
    """Plot log-log scaling fit and save figure."""
    plt.figure()
    plt.scatter(log_N, log_r_ee)
    plt.plot(log_N, slope * log_N + intercept, color="red",
             label=f"Fit: slope = {slope:.3f}")
    plt.xlabel("log(N)")
    plt.ylabel("log(r_ee/a)")
    plt.title("Scaling Exponent Estimation")
    plt.legend()
    plt.tight_layout()
    plt.savefig("scaling.png")
    plt.close()