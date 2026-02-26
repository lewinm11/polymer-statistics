import numpy as np
import matplotlib.pyplot as plt
import typing
from scipy import stats

def log_linear_regression(N, r_ee):
    log_N = np.log(N)
    log_r_ee = np.log(r_ee)
    # Perform linear regression on the log-log data
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_N, log_r_ee)
    return slope, intercept, r_value, std_err, log_N, log_r_ee

def plot_scaling(slope, intercept, r_value, std_err, log_N, log_r_ee):
    # Plot the data and the fitted line
    plt.scatter(log_N, log_r_ee)
    plt.plot(log_N, slope * log_N + intercept, color='red', label=f'Fit: slope={slope:.2f}')
    plt.xlabel('log(N)')
    plt.ylabel('log(r_ee/a)')
    plt.title('Estimation of Scaling Exponent')
    plt.legend()
    plt.savefig('scaling.png')
    plt.show()
    print(slope, intercept, r_value, std_err)