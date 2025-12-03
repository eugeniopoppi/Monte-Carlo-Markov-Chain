import matplotlib.pyplot as plt
import numpy as np
import math
import scipy
from scipy.optimize import curve_fit

# L =  [8,12,16,20,24,28,32]
# m =  [4.45,9.3,15.4,24.0,37.5,40.0,51.0]
# dm = [0.08,0.2,0.2,0.3,0.5,1,1]
#
# L =  [8,12,16,20,24,28,32,36]
# m =  [5.6,8.9,15.5,22.6,37.4,37,51.0,58]
# dm = [0.3,0.1,0.3,0.3,0.6,1,1,1]

L =  [8,12,16,20,24,28,32,36]
# m =  [-0.6245326694482924,-1.045830901077977,-1.4855321614385264,-2.710386229089189,-2.766241236036998,-5.076629939809792,
#      -5.949876510608418,-8.121754409413285]
m =  [0.6245326694482924,1.045830901077977,1.4855321614385264,2.710386229089189,2.766241236036998,5.076629939809792,
     5.949876510608418,8.121754409413285]

dm = [0.08561266486403339,0.06944786777792852,0.09740277485339834,0.15488372549394042, 0.10985003850926374,0.48699137942916493,
     0.408119309582087,0.4362281673600251]
dm = np.array(dm)
dm *= 2.8
def line(x,a,nu):
    return a*x**(1/nu)

#fit
popt,pcov = curve_fit(line,L,m,sigma=dm,absolute_sigma = True)
a,nu= popt
da,dnu = np.sqrt(pcov.diagonal())
# dof e chi^2
dof = len(L)- popt.size
chi = (m - line(L, *popt)) / dm
chisq = np.sum(chi**2)
print(f"a = {a}+-{da},\n nu = {nu}+-{dnu},\n")
print(pcov)
print(f"Chi2 = {chisq:.3g} con dof = {dof}\n")


plt.errorbar(L,m,dm,fmt=".")
x = np.linspace(7,37,100)
plt.plot(x,line(x,*popt),color="green")
#plt.plot(x,line(x,a,0.6717,),color = "blue")
plt.ylabel('Binder Slope (L)', fontsize = 14)
plt.xlabel('L', fontsize = 14)
# Posiziona la legenda in alto a destra
#plt.legend(loc='upper left', fontsize = 14)

plt.show()
