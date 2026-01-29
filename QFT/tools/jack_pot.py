#!/usr/bin/env python3
import numpy as np
import sys, os

def jackknife_for_primary(data,k):
    N = len(data)
    blocknumb = N // k
    x = data[:blocknumb * k]

    #blocking
    X = []
    for j in range(blocknumb):
        X.append(np.mean(x[j*k:(j+1)*k]))
    X = np.array(X)

    #jackknife
    S = sum(X)
    S_iter = []
    for i in range(blocknumb):
        S_i  = (S-X[i])/(blocknumb-1)
        S_iter.append(S_i)

    #error bars
    dS_iter = 0
    S_mean = np.mean(S_iter)
    for l in range(len(S_iter)):
        dS_iter += (S_iter[l]-S_mean)*(S_iter[l]-S_mean)
    dS_iter =  np.sqrt(((blocknumb - 1) / blocknumb) *(dS_iter))

    return S_mean ,dS_iter


def jackknife_for_potential(data,k,wt): #-log(data)
    N = len(data)
    blocknumb = N // k
    x = data[:blocknumb * k]

    #blocking
    X = []
    for j in range(blocknumb):
        X.append(np.mean(x[j*k:(j+1)*k]))
    X = np.array(X)

    #jackknife
    S = sum(X)
    pot_iter = []
    for i in range(blocknumb):
        S_i  = (S-X[i])/(blocknumb-1)
        pot_i = -np.log(S_i)/wt
        pot_iter.append(pot_i)

    #error bars
    dpot_iter = 0
    pot_mean = np.mean(pot_iter)
    for l in range(len(pot_iter)):
        dpot_iter += (pot_iter[l]-pot_mean)*(pot_iter[l]-pot_mean)
    dpot_iter =  np.sqrt(((blocknumb - 1) / blocknumb) *(dpot_iter))

    return pot_mean,dpot_iter


# main
if __name__=="__main__":

  try:
    infile =sys.argv[1]
  except:
    print("USE: %s file_name block_size wt\n" % sys.argv[0])
    sys.exit(1)

  try:
    blocksize =int(sys.argv[2])
  except:
    print("USE: %s file_name block_size wt\n" % sys.argv[0])
    sys.exit(1)
  try:
    wt =int(sys.argv[3])
  except:
    print("USE: %s file_name block_size wt\n" % sys.argv[0])
    sys.exit(1)
  if not os.path.isfile(infile):
    print("ERROR: file %s does not exists\n" % infile)
    sys.exit(1)
  # data acquisition
  indata=np.loadtxt(infile, skiprows=1, dtype=float) #np.float
  #jackknife for the potential
  ris, err = jackknife_for_primary(indata,blocksize)
  print(ris, err, end=' ')
  ris, err = jackknife_for_potential(indata,blocksize,wt)
  print(ris, err, end='\n')
