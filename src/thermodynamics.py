import numpy as np
import matplotlib.pyplot as plt
import typing
from scipy import stats

def exp_E (E_eps, g_E, T, k, eps, E_min):
    return (np.sum(g_E * E_eps * np.exp(-(E_eps - E_min) / (k * T)))) / np.sum(g_E * np.exp(-(E_eps-E_min) / (k * T)))

def exp_E_squared (E_eps, g_E, T, k, eps, E_min):
    return (np.sum(g_E * E_eps**2 * np.exp(-(E_eps - E_min) / (k * T)))) / np.sum(g_E * np.exp(-(E_eps-E_min) / (k * T)))

def C_v_calc (E_eps, g_E, T, k, eps, E_min):
    return 1/(k * T**2) * (exp_E_squared(E_eps, g_E, T, k, eps, E_min) - exp_E(E_eps, g_E, T, k, eps, E_min)**2)

def prob_calc(E_eps, g_E, T, k, eps, E_min):
    return 1/np.sum(g_E * np.exp(-(E_eps-E_min)/ (k * T)))

def plot_C_V(E_eps, g_E, T, k, eps, E_min):
# Calculate C_v for each temperature and find T_max
    C_v_values = [C_v_calc(E_eps, g_E, T_i, k, eps, E_min) for T_i in T]
    T_max = T[np.argmax(C_v_values)]
    print(f"Temperature T_max = {T_max}")
# Plot C_v vs T and mark T_max
    plt.plot(T, C_v_values)
    plt.vlines(T_max, min(C_v_values), max(C_v_values), colors='red', linestyles='dashed', label=f'T_max = {T_max:.2f}')
    plt.legend()
    plt.xlabel('Temperature (T)')
    plt.ylabel('Heat Capacity (C_v)')
    plt.title('Heat Capacity vs Temperature')
    plt.savefig('thermodynamics.png')
    plt.legend()
    plt.show()

def plot_prob(E_eps, g_E, T, k, eps, E_min):
    prob_values = [prob_calc(E_eps, g_E, T_i, k, eps, E_min) for T_i in T]
    plt.plot(T, prob_values)
    plt.xlabel('Temperature (T)')
    plt.ylabel('Probability of Ground State')
    plt.title('Ground State Probability vs Temperature')
    plt.savefig('probability.png')
    plt.show()