#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math
import scipy
from scipy.optimize import curve_fit


print(f"-----FIT OF THE FINITE SIZE SCALING OF THE IMPROVED SUSCEPTIBILITY-----")
print("-----------------------beta_pc vs L------------------------------------")
file_path = '/home/eugenio/Documenti/c/O_N/tools/chi1_done.txt'
L,beta_pc, dbeta_pc, Xmax, dXmax =   np.loadtxt(file_path, unpack=True)



#define the fit function
def line(x,a,nu,b):
    return a*x**(-1/nu)+b

#fit
popt,pcov = curve_fit(line,L,beta_pc,sigma=dbeta_pc,absolute_sigma = True)
a, nu,b = popt
# check pcov validità
if np.any(~np.isfinite(pcov)):
    print("Attenzione: la matrice di covarianza contiene NaN o inf. Il fit potrebbe essere instabile.")
# estrai varianze e deviazioni
var = np.diag(pcov)
da,dnu,db = np.sqrt(var)

# dof e chi^2
dof = L.size - popt.size
chi = (beta_pc - line(L, *popt)) / dbeta_pc
chisq = np.sum(chi**2)
reduced_chi2 = chisq / dof if dof>0 else np.nan

print(f"a = {a:.6g} ± {da:.2g}")
print(f"beta_c= {b:.6g} ± {db:.2g}")
print(f"nu = {nu:.6g} ± {dnu:.2g}")
print(pcov)
print(f"Chi2 = {chisq:.3g} con dof = {dof}\n") #, chi2/dof = {reduced_chi2:.3g}


plt.errorbar(L,beta_pc,dbeta_pc,fmt=".")
x = np.linspace(L.min(),L.max(),1000)
plt.plot(
    x,
    line(x, *popt),
    #label = rf"$\nu = {nu:.6g} \pm {dnu:.6g}$" "\n" rf"$\beta_c = {b:.6g} \pm {db:.6g}$",
    color = "green"
)

plt.ylabel(r"$\beta_{p.c.}(L)$", fontsize=12)
plt.xlabel(r"$L$", fontsize=12)
#plt.legend(loc='upper left', fontsize=12)
plt.show()
