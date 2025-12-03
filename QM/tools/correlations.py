#!/usr/bin/env python3
import numpy as np
import sys, os, re
from corr import correlation
from scipy.optimize import curve_fit

# Modello esponenziale decrescente
def exp_decay(n, A, tau):
    return A * np.exp(-n / tau)

# Fit
def fit_exponential(corr_array):
    n = np.arange(len(corr_array))
    y = np.array(corr_array)

    sigma = np.ones_like(y)

    popt, pcov = curve_fit(exp_decay, n, y, p0=[1.0, 10.0], sigma=sigma, absolute_sigma=False)

    A, tau = popt
    dA, dtau = np.sqrt(np.diag(pcov))

    residuals = y - exp_decay(n, *popt)
    chi2 = np.sum((residuals / sigma) ** 2)
    ndof = len(y) - len(popt)
    chi2_reduced = chi2 / ndof if ndof > 0 else np.nan

    return A, dA, tau, dtau, chi2_reduced


# Funzione per estrarre il numero N dal nome del file
def extract_N(filename):
    numbers = re.findall(r'\d+', filename)
    if len(numbers) == 0:
        return None
    return int(numbers[0])   # Usa il primo intero trovato


# Read directory
dir_path = input("Write the file directory (e.g. /home/eugenio/.../file):\n").strip()

if not os.path.isdir(dir_path):
    print(f"Directory '{dir_path}' not found.")
    exit(1)

# Lista file .txt
file_paths = [
    os.path.join(dir_path, file)
    for file in os.listdir(dir_path)
    if file.endswith(".txt")
]

if len(file_paths) == 0:
    print("No .txt files found in the specified directory.")
    exit(1)

entries = []

# Loop su ogni file
for path in file_paths:
    filename = os.path.basename(path)
    N_value = extract_N(filename)

    if N_value is None:
        print(f"Warning: no integer found in filename '{filename}', skipping.")
        continue

    data = np.loadtxt(path)
    data = data*data
    corr_values = correlation(data, 5000, 50000)

    try:
        A, dA, tau, dtau, chi2 = fit_exponential(corr_values)
    except Exception as e:
        print(f"Fit failed for {path}: {e}")
        continue

    entries.append((N_value, filename, corr_values, A, dA, tau, dtau, chi2))

# Ordino per N
entries.sort(key=lambda x: x[0])


#faccio i grafici
import matplotlib.pyplot as plt

plt.figure(figsize=(8,6))

colors = plt.cm.viridis(np.linspace(0, 1, len(entries)))

for idx, entry in enumerate(entries):
    N, filename, corr_values, A, dA, tau, dtau, chi2 = entry

    n = np.arange(len(corr_values))
    yerr = np.zeros_like(corr_values)

    plt.errorbar(
        n,
        corr_values,
        yerr=yerr,
        fmt="o",
        markersize=3,
        color=colors[idx],
        label=f"N = {N}"
    )


plt.xlabel("$n$", fontsize = 14)
plt.ylabel(r"$C_{Q^2}(n)$", fontsize = 14)
plt.title("Correlazioni per diversi N", fontsize = 14)
plt.legend(loc='upper right',fontsize = 14)
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()

# Salvo risultati
output_file = os.path.join(dir_path, "fit_results.txt")

with open(output_file, "w") as f:
    f.write("# N  A   dA   tau   dtau   chi2_reduced\n")
    for entry in entries:
        N, filename, corr_values, A, dA, tau, dtau, chi2 = entry
        f.write(f"{N} {A:.6g} {dA:.6g} {tau:.6g} {dtau:.6g} {chi2:.6g}\n")


print(f"Fit results saved to: {output_file}")
