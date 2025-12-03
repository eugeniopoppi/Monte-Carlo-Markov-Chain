#!/usr/bin/env python3
import numpy as np
import sys, os

def jackknife_for_x(data,k):
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


def jackknife_for_x2(data,k):
    N = len(data)
    blocknumb = N // k
    x = data[:blocknumb * k]

    #blocking on x^2
    X2 = []
    for j in range(blocknumb):
        X2.append(np.mean(x[j*k:(j+1)*k]**2))
    X2 = np.array(X2)

    #jackknife
    S2 = np.mean(X2)
    S2_iter = []
    for i in range(blocknumb):
        S2_i = (np.sum(X2) - X2[i]) / (blocknumb - 1)
        S2_iter.append(S2_i)
    S2_iter = np.array(S2_iter)

    #error bars
    dS2 = np.sqrt((blocknumb - 1) / blocknumb * np.sum((S2_iter - np.mean(S2_iter))**2))

    return S2, dS2

# main
if __name__=="__main__":

  try:
    infile =sys.argv[1]
  except:
    print("USE: %s file_name block_size\n" % sys.argv[0])
    sys.exit(1)

  try:
    blocksize =int(sys.argv[2])
  except:
    print("USE: %s file_name block_size\n" % sys.argv[0])
    sys.exit(1)

  if not os.path.isfile(infile):
    print("ERROR: file %s does not exists\n" % infile)
    sys.exit(1)
  # data acquisition
  indata=np.loadtxt(infile, skiprows=1, dtype=float) #np.float
  #jackknife for the potential
  ris, err = jackknife_for_x(indata,blocksize)
  print(ris, err, end=' ')
  ris, err = jackknife_for_x2(indata,blocksize)
  print(ris, err, end='\n')
