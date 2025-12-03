#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<complex.h>
#include"/home/eugenio/Documenti/c/header/random.h"
#define PI 3.14159265358979323846
#define numb_of_upd 1000000

//distance on a circle
double d(double tmp1, double tmp2){
  double dx = tmp1 - tmp2;
  if (fabs(dx) <= 0.5) return dx;
  else if (dx > 0.5)   return dx - 1.0;
  else                 return dx + 1.0;
}

//winding number
double top_charge(double *lattice,long int *nnp,long int N){
  double Q = 0.0;
  long int j,tmp;

  for (j=0;j<N;j++){
    tmp = nnp[j];
    Q += d(lattice[tmp],lattice[j]);
  }
  //Q = llround(Q);
  return Q;
}

//metropolis
int metropolis(int i,double delta,double eta,
               double *lattice,long int *nnp,long int *nnn){
  double dE,x,xt,xn,xp;
  long int idx;
  int acc=0;

  x = lattice[i];      //the site to be updated
  xt = lattice[i] + delta*(1.0-2.0*myrand()); //trial state, keep it in (0,1)
  if (xt < 0) xt += 1.0;
  if (xt >= 1.0) xt -= 1.0;
  idx = nnp[i];
  xn = lattice[idx];  //next site
  idx = nnn[i];
  xp = lattice[idx];  //previus site

  //energy variation
  dE = -d(xn,x)*d(xn,x)-d(x,xp)*d(x,xp)+d(xn,xt)*d(xn,xt)+d(xt,xp)*d(xt,xp);
  dE = dE/(2.0*eta);

  //accept-reject metropolis
  if(dE<0){
    lattice[i] = xt;
    acc = 1;
  }
  else{
    if (myrand()<exp(-dE)){
      lattice[i]=  xt;
      acc = 1;
    }
  }
  return acc;
}

//main
int main(int argc, char **argv){
  double *lattice,etaN,eta,delta,Q;
  int i,k,n,meas_every;
  long int *nnp,*nnn,N,acc = 0;
  const unsigned long int seed1=(unsigned long int) time(NULL);
  const unsigned long int seed2=seed1+127;
  char datafile[100];
  FILE *fp;
  myrand_init(seed1, seed2);
  time_t start_time = time(NULL);

  if(argc != 5)
    {
    fprintf(stdout, "How to use this program:\n");
    fprintf(stdout, "  %s etaN N measure_every out_file \n", argv[0]);
    // fprintf(stdout, "  eta*N = double\n");
    // fprintf(stdout, "  N = long int\n");
    // fprintf(stdout, "  delta = double, metro step\n");
    // fprintf(stdout, "  measure_every = int\n");
    // fprintf(stdout, "  datafile = name of the file on which to write the data\n\n");
    fprintf(stdout, "Output: Q measured every -measure_every- complete update of the lattice\n");

    return EXIT_SUCCESS;
    }
  else
    {
    // read input values
    etaN=atof(argv[1]);
    N=atol(argv[2]);
    eta = etaN/N;
    delta= 2.5*sqrt(eta); //atof(argv[3]);
    meas_every=atoi(argv[3]);
    strcpy(datafile, argv[4]);
    }

  //allocate the lattice, and the next neighbor in positive and negative directions
  lattice=(double*)malloc((unsigned long int)(N)*sizeof(double));
  if(lattice == NULL)
    {
    fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }
  nnp=(long int *)malloc((unsigned long int)(N)*sizeof(long int));
  if(nnp == NULL)
    {
    fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }
  nnn=(long int *)malloc((unsigned long int)(N)*sizeof(long int));
  if(nnn == NULL)
    {
    fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }
  // open data file
  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    return EXIT_FAILURE;
    }
  //initialize the vectors
  for(i=0; i<N; i++){
    lattice[i] = 0.0;  //ordered start, = double*rand() also valid
    nnp[i] = (i+1)%N;  //next neighbor in the positive direction
    nnn[i] = (i-1+N)%N; //next neighbor in the negative direction
  }

  //update
  for(i=0; i<numb_of_upd; i++){
    for(k=0; k<N; k++){
      n = rand()%N;
      acc+=metropolis(n,delta,eta,lattice,nnp,nnn);
    }
    if(i%meas_every==0){
      Q = top_charge(lattice,nnp,N);
      fprintf(fp, "%lf \n",Q);
    }
  }

  // free the memory
  free(lattice);
  free(nnp);
  free(nnn);

  //close the file
  fclose(fp);
  //running time
  time_t end_time = time(NULL);
  double time_spent = difftime(end_time,start_time);
  printf("Execution time = %lf \n",time_spent);
  printf("acceptance rate (with delta=%lf) = %lf\n",delta,(double)acc / (N * numb_of_upd));
  return EXIT_SUCCESS;
  }
