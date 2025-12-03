
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy
from scipy.optimize import curve_fit


print(f"-----FIT OF THE FINITE SIZE SCALING OF THE IMPROVED SUSCEPTIBILITY-----")
print("-----------------chi'_{max} vs L-------------------")


file_path = '/home/eugenio/Documenti/py/Chi1_analyzed.txt'
L,beta_pc, dbeta_pc, Xmax, dXmax =   np.loadtxt(file_path, unpack=True)
dXmax = dXmax #*2
#define the fit function
def curve(x,a,g_n): #g_n sarebbe gamma/nu
    return a*x**g_n

#fit
popt,pcov = curve_fit(curve,L,Xmax,sigma=dXmax,absolute_sigma = True)
a,g_n= popt
# check pcov validità
if np.any(~np.isfinite(pcov)):
    print("Attenzione: la matrice di covarianza contiene NaN o inf. Fit instabile.")
# estrai varianze e deviazioni
var = np.diag(pcov)
da,dg_n= np.sqrt(var)

# dof e chi^2
dof = L.size - popt.size
chi = (Xmax - curve(L, *popt)) / dXmax
chisq = np.sum(chi**2)
reduced_chi2 = chisq / dof if dof>0 else np.nan

print(f"a = {a:.6g} ± {da:.2g}")
print(f"gamma/nu = {g_n:.6g} ± {dg_n:.2g}")
print(pcov)
print(f"Chi2 = {chisq:.3g} con dof = {dof}\n") #chi2/dof = {reduced_chi2:.3g}


plt.errorbar(L,Xmax,dXmax,fmt=".")
x = np.linspace(L.min()-1,L.max()+1,1000)
gamma = g_n * 0.71  # nu = 0.71 +- 0.25 da fit precedente
dgamma = np.sqrt((dg_n*0.71)**2+(0.25*g_n)**2)


plt.plot(x,curve(x,*popt),label = f"gamma = {gamma:.6g} ± {dgamma:.6g}",color="green")
#label = f"gamma = {gamma:.6g} ± {dgamma:.6g}",

plt.ylabel(r"$\chi'_{max}(L)$", fontsize = 12)
plt.xlabel('L', fontsize = 12)
# Posiziona la legenda in alto a destra
plt.legend(loc='upper left', fontsize = 12)

plt.show()
