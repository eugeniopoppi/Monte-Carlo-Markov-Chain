#!/usr/bin/env python3
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Modello originale: tau = a * exp(z*N)
def exp_model(N, a, z):
    return a * np.exp(z * N)

# Percorso del file (non directory!)
path = "/home/eugenio/Documenti/c/qm_topo/corr_local/fit_results.txt"

# Carico i dati: N, A, dA, tau, dtau, chi2
N, A, dA, tau, dtau, chi2 = np.loadtxt(path, unpack=True)

# ------------------------------------------------------------
# FIT LOGARITMICO (stabilizzato numericamente)
# log(tau) = log(a) + z*N
# ------------------------------------------------------------
def lin_model(N, loga, z):
    return loga + z * N

# stime iniziali robuste
loga0 = np.log(tau[0])
z0 = (np.log(tau[-1]) - np.log(tau[0])) / (N[-1] - N[0])

popt_lin, pcov_lin = curve_fit(
    lin_model,
    N,
    np.log(tau),
    sigma=dtau / tau,
    p0=[loga0, z0],
    absolute_sigma=False
)

loga, z = popt_lin
dloga, dz = np.sqrt(np.diag(pcov_lin))

# Ricostruisco i parametri del modello originale
a = np.exp(loga)
da = a * dloga

print(f"a = {a}  {da}")
print(f"z = {z}  {dz}")

residuals = tau - exp_model(N, a, z)
chi2_val = np.sum((residuals / dtau)**2)
ndof = len(tau) - len([a, z])
chi2_reduced = chi2_val / ndof
print(f"chi2_red = {chi2_reduced}")

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
N_fit = np.linspace(min(N), max(N), 300)
tau_fit = exp_model(N_fit, a, z)

plt.errorbar(N, tau, yerr=dtau, fmt='o')
plt.plot(N_fit, tau_fit, '-', label=rf'$a \approx {a:.3g}, \,  z \approx {z:.3g}$')

plt.xlabel(r"$N$", fontsize = 14)
plt.ylabel(r"$\tau_{exp}(N)$", fontsize = 14)
plt.yscale("log")
plt.legend(loc='upper left', fontsize = 14)
# plt.grid(True, which='both')
plt.title(r"Exponential Fit $\tau_{exp}(N) = a · exp(z·N)$", fontsize = 14)

plt.show()
