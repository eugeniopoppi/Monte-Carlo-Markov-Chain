#!/usr/bin/env python3
# Computation of observables (defined below) with jackknife and blocking

import matplotlib.pyplot as plt
import numpy as np
import math
import os

print("This code computes the following observables (and relative errors) through jackknife:\n")
print("- Binder Cumulant U")
print("- Susceptibility Chi")
print("- A different definition of Susceptibility Chi'")
print("- Specific heat C")

# Read the directory from the prompt line
dir_path = input("Write the file directory (e.g. /home/eugenio/.../file):\n").strip()

# Check if directory exists
if not os.path.isdir(dir_path):
    print(f"Directory '{dir_path}' not found.")
    exit(1)

# Create a list of file paths
file_paths = []
for file in os.listdir(dir_path):
    if file.endswith(".txt"):
        full_path = os.path.join(dir_path, file)
        file_paths.append(full_path)

# Check if any .txt files were found
if len(file_paths) == 0:
    print("No .txt files found in the specified directory.")
    exit(1)

# Read the output directory from the prompt line
output_file = input("The file to save the data on (e.g. /home/eugenio/.../output.txt):\n").strip()
with open(f"{output_file}","w") as f:
    f.write("#beta U dU Chi dChi Chi1 dChi1 C dC \n")

#remove the termaliazion phase
therm = int(input("Insert the thermalization to be excluded: \n"))
block = int(input('Insert block length: \n'))

#now the computations
for idx, path in enumerate(file_paths):
    #unpack the file
    print(f"open file_path_{idx} = {path}")
    mx,my,e = np.loadtxt(path, unpack=True)
    mx = mx[therm:]
    my = my[therm:]
    e = e[therm:]

    m = np.sqrt(my*my + mx*mx)
    #observables to compute
    # U = []       #<m^4>/<m^2>^2
    # CHI = []     #<m^2>
    # CHI1 = []     #<m^2>-<|m|>^2
    # C = []       #<e^2>-<e>^2
    k = block
    N = len(m)
    blocknumb = N // k
    m = m[:blocknumb * k]
    e = e[:blocknumb * k]

    M = []
    M2 = []
    M4 = []
    E = []
    E2 = []
    for j in range(blocknumb):
        Mj = m[(j*k):((j+1)*k)]
        Ej = e[(j*k):((j+1)*k)]
        Mj = np.mean(Mj)
        Ej = np.mean(Ej)
        M.append(Mj)
        M2.append(Mj*Mj)
        M4.append(Mj*Mj*Mj*Mj)
        E.append(Ej)
        E2.append(Ej*Ej)
    M = np.array(M)
    M2 = np.array(M2)
    M4 = np.array(M4)
    E = np.array(E)
    E2 = np.array(E2)

    #jackknife
    S = sum(M)
    S2 = sum(M2)
    S4 = sum(M4)
    SE = sum(E)
    SE2 = sum(E2)

    U= []
    CHI= []
    CHI1 = []
    C = []

    for i in range(blocknumb):
        S_i  = (S-M[i])/(blocknumb-1)
        S2_i = (S2-M2[i])/(blocknumb-1)
        S4_i = (S4-M4[i])/(blocknumb-1)
        SE_i = (SE-E[i])/(blocknumb-1)
        SE2_i = (SE2-E2[i])/(blocknumb-1)
        U_i = S4_i/(S2_i)**2
        CHI_i = S2_i
        CHI1_i = S2_i-S_i**2
        C_i = SE2_i-SE_i**2
        U.append(U_i)
        CHI.append(CHI_i)
        CHI1.append(CHI1_i)
        C.append(C_i)

    #error bars
    dU = 0
    dCHI = 0
    dCHI1= 0
    dC = 0

    U_mean = np.mean(U)
    CHI_mean = np.mean(CHI)
    CHI1_mean = np.mean(CHI1)
    C_mean = np.mean(C)
    for k in range(len(U)):
        dU += (U[k] - U_mean)*(U[k] - U_mean)
        dCHI += (CHI[k]-CHI_mean)*(CHI[k]-CHI_mean)
        dCHI1 += (CHI1[k]-CHI1_mean)*(CHI1[k]-CHI1_mean)
        dC += (C[k]-C_mean)*(C[k]-C_mean)
    dU=  np.sqrt(((blocknumb - 1) / blocknumb) *dU)
    dCHI =  np.sqrt(((blocknumb - 1) / blocknumb) *dCHI)
    dCHI1 =  np.sqrt(((blocknumb - 1) / blocknumb) *dCHI1)
    dC = np.sqrt(((blocknumb - 1) / blocknumb) *dC)


    #normalize chi chi1 and c
    #read beta and L from the file path ......./xyL40B0.4547.txt
    beta = float(path[-10:-4])
    L = int(path[-13:-11])
    #L = int(path[-12])
    CHI_mean = L*L*L*CHI_mean
    dCHI = L*L*L*dCHI
    CHI1_mean = L*L*L*CHI1_mean
    dCHI1 = L*L*L*dCHI1
    C_mean = L*L*L*C_mean
    dC = L*L*L*dC

    #for each iteration (different beta) append a row
    with open(f"{output_file}","a") as f:
        f.write(f"{beta} {U_mean} {dU} {CHI_mean} {dCHI} {CHI1_mean} {dCHI1} {C_mean} {dC} \n")


print("Processing complete.")
