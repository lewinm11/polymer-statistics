import numpy as np
import matplotlib.pyplot as plt


def exp_E(E_eps, g_E, T, k, eps, E_min):
    """Return ⟨E⟩ at temperature T using shifted Boltzmann weights."""
    weights = g_E * np.exp(-(E_eps - E_min) / (k * T))
    return np.sum(E_eps * weights) / np.sum(weights)


def exp_E_squared(E_eps, g_E, T, k, eps, E_min):
    """Return ⟨E²⟩ at temperature T."""
    weights = g_E * np.exp(-(E_eps - E_min) / (k * T))
    return np.sum(E_eps**2 * weights) / np.sum(weights)


def C_v_calc(E_eps, g_E, T, k, eps, E_min):
    """Return heat capacity C_V at temperature T."""
    mean_E = exp_E(E_eps, g_E, T, k, eps, E_min)
    mean_E2 = exp_E_squared(E_eps, g_E, T, k, eps, E_min)
    return (mean_E2 - mean_E**2) / (k * T**2)


def prob_calc(E_eps, g_E, T, k, eps, E_min):
    """Return ground-state probability at temperature T."""
    return 1.0 / np.sum(g_E * np.exp(-(E_eps - E_min) / (k * T)))


def plot_C_V(E_eps, g_E, T, k, eps, E_min):
    """Plot heat capacity as a function of temperature."""
    C_v_values = [C_v_calc(E_eps, g_E, T_i, k, eps, E_min) for T_i in T]
    T_max = T[np.argmax(C_v_values)]

    plt.figure()
    plt.plot(T, C_v_values)
    plt.axvline(T_max, linestyle="--", label=f"T_max = {T_max:.2f}")
    plt.xlabel("Temperature (T)")
    plt.ylabel("Heat Capacity (C_v)")
    plt.title("Heat Capacity vs Temperature")
    plt.legend()
    plt.tight_layout()
    plt.savefig("thermodynamics.png")
    plt.close()


def plot_prob(E_eps, g_E, T, k, eps, E_min):
    """Plot ground-state probability as a function of temperature."""
    prob_values = [prob_calc(E_eps, g_E, T_i, k, eps, E_min) for T_i in T]

    plt.figure()
    plt.plot(T, prob_values)
    plt.xlabel("Temperature (T)")
    plt.ylabel("Probability of Ground State")
    plt.title("Ground State Probability vs Temperature")
    plt.tight_layout()
    plt.savefig("probability.png")
    plt.close()