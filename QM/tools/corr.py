#!/usr/bin/env python3
import numpy as np
import sys, os


#funzione che calcola le correlazioni prendendo in pasto una colonna  e massima dist a cui osservare le corr
def correlation(x,n_max,therm):

    x = x - np.mean(x)

    if therm != 0 :
        x = x[therm:]

    N = len(x)
    if N == 0:
        raise ValueError("Thermalization is too large, no data left!")

    n_max = min(n_max, N - 1)  # Avoid out-of-bounds

    x2 = np.mean(x**2)
    if x2 == 0:
        raise ValueError("Variance is zero, autocorrelation not defined!")

    # Vectorized autocorrelation
    Cx = [np.mean(x[:N-n] * x[n:N]) / x2 for n in range(n_max)]

    return Cx

# main
if __name__=="__main__":

  try:
    infile =sys.argv[1]
  except:
    print("USE: %s input_file n_max therm output_file \n" % sys.argv[0])
    sys.exit(1)

  try:
    n_max =int(sys.argv[2])
  except:
    print("USE: %s input_file n_max therm output_file \n" % sys.argv[0])
    sys.exit(1)

  try:
    therm = int(sys.argv[3])
  except:
    print("USE: %s input_file n_max therm output_file \n" % sys.argv[0])
    sys.exit(1)

  try:
    outfile = sys.argv[4]
  except:
    print("USE: %s input_file n_max therm output_file \n" % sys.argv[0])
    sys.exit(1)

  if not os.path.isfile(infile):
    print("ERROR: file %s does not exists\n" % infile)
    sys.exit(1)


  # data acquisition
  indata=np.loadtxt(infile, skiprows=0, dtype=float)
  indata = indata*indata
  #jackknife for the potential
  Cx = correlation(indata,n_max,therm)

  # Save results
  with open(outfile, "w") as f:
      f.write("#n C_n \n")
      for n, Cn in enumerate(Cx):
          f.write(f"{n} {Cn}\n")

  print(f"Autocorrelation saved to {outfile}")
