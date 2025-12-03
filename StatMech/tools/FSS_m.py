#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math
import scipy
from scipy.optimize import curve_fit


print(f"-----FIT OF THE FINITE SIZE SCALING OF THE MAGNETIZATION-----")
print("-----------------------<|m|> vs L------------------------------------")
file_path = '/home/eugenio/Documenti/c/O_N/tools/m_vs_L.txt'
L,m, dm =   np.loadtxt(file_path, unpack=True)


#define the fit function
def line(x,a,beta_nu):
    return a*x**(-beta_nu) #beta_nu = beta/nu

#fit
popt,pcov = curve_fit(line,L,m,sigma=dm,absolute_sigma = True)
a,beta_nu = popt
# check pcov validità
if np.any(~np.isfinite(pcov)):
    print("Attenzione: la matrice di covarianza contiene NaN o inf. Il fit potrebbe essere instabile.")
# estrai varianze e deviazioni
var = np.diag(pcov)
da,dbeta_nu = np.sqrt(var)

# dof e chi^2
dof = L.size - popt.size
chi = (m - line(L, *popt)) / dm
chisq = np.sum(chi**2)

print(f"a = {a:.6g} ± {da:.2g}")
print(f"beta_nu= {beta_nu:.6g} ± {dbeta_nu:.2g}")
beta = beta_nu * 0.70  # nu = 0.71 +- 0.25 da fit precedente
dbeta = np.sqrt((dbeta_nu*0.80)**2+(0.31*beta_nu)**2)
print(f"beta= {beta:.6g} ± {dbeta:.2g}")
print(pcov)
print(f"Chi2 = {chisq:.3g} con dof = {dof}\n")


plt.errorbar(L,m,dm,fmt=".")
x = np.linspace(L.min(),L.max(),1000)
plt.plot(
    x,
    line(x, *popt),
    #label = rf"$\nu = {nu:.6g} \pm {dnu:.6g}$" "\n" rf"$\beta_c = {b:.6g} \pm {db:.6g}$",
    color = "green"
)
plt.ylabel(r"$\langle |\mathbf{m}| \rangle (L)$", fontsize=12)
plt.xlabel(r"$L$", fontsize=12)

#plt.legend(loc='upper left', fontsize=12)
plt.show()
