#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<complex.h>
#include"/home/eugenio/Documenti/c/header/random.h"
#define PI 3.14159265358979323846
#define numb_of_upd 1000000


//structure of each lattice site,
//u0 is the temporal direction and u1 the spatial direction
typedef struct site {
  double complex u0;
  double complex u1;
} site;


//staple in the TEMPORAL DIRECTOIN
double complex staple0(site **lattice,int i, int j,int Nt, int Ns){
  double complex S,Stmp;
  long int tp,tn,xp,xn;

  //next neighbors
  tp = (i+1)%Nt;
  tn = (i-1+Nt)%Nt;
  xp = (j+1)%Ns;
  xn = (j-1+Ns)%Ns;

  //staple in positive direction
  S = lattice[tp][j].u1;
  S *= conj(lattice[i][xp].u0);
  S *= conj(lattice[i][j].u0);

  //staple in negative direction
  Stmp = conj(lattice[tp][xn].u1);
  Stmp *= lattice[i][xn].u0;
  Stmp *= conj(lattice[i][xn].u1);
  S+= Stmp;

  return S;
}

//staple in the SPATIAL DIRECTION
double complex staple1(site **lattice,int i, int j,int Nt, int Ns){
  double complex S,Stmp;
  long int tp,tn,xp,xn;

  //next neighbors
  tp = (i+1)%Nt;
  tn = (i-1+Nt)%Nt;
  xp = (j+1)%Ns;
  xn = (j-1+Ns)%Ns;

  //staple in positive direction
  S = conj(lattice[i][j].u0);
  S *= conj(lattice[tp][j].u1);
  S *= lattice[i][j].u0;

  //staple in negative direction
  Stmp = lattice[tn][j].u0;
  Stmp *= conj(lattice[tn][j].u1);
  Stmp *= conj(lattice[tn][xp].u0);
  S += Stmp;

  return S;
}

//metropolis
int metropolis(site **lattice, int i, int j,double eps, double beta, int Nt, int Ns){
  double complex Ut,Eold,Enew,Stmp,dE;
  double theta,eold,enew;
  int acc=0;

  //lets update temporal direction
  theta = eps*(2*myrand()-1);
  //Ut = (1+I*theta)/sqrt(1+theta*theta);
  Ut = cexp(-I*theta);
  Stmp = staple0(lattice,i,j,Nt,Ns);

  //old energy
  Eold = lattice[i][j].u0*Stmp;
  eold = -beta*creal(Eold);

  //new energy
  Enew = Ut*Stmp;
  enew = -beta*creal(Enew);

  //energy variation
  dE = enew-eold;

  if((double)(dE)<0){
    lattice[i][j].u0 =  Ut;
    acc+=1;
  }
  //accept-reject with exp{-beta dE}
  else{
    if (myrand()<exp(-dE)){
      lattice[i][j].u0 =  Ut;
      acc+=1;
    }
  }

  //lets update spatial direction
  theta = eps*(2*myrand()-1);
  Ut = (1+I*theta)/sqrt(1+theta*theta);

  Stmp = staple1(lattice,i,j,Nt,Ns);

  //old energy
  Eold = lattice[i][j].u1*Stmp;
  eold = -beta*creal(Eold);
  //new energy
  Enew = Ut*Stmp;
  enew = -beta*creal(Enew);
  //energy variation
  dE = enew-eold;


  if((double)(dE)<0){
    lattice[i][j].u1 =  Ut;
    acc+=1;
  }
  //accept-reject with exp{-beta dE}
  else{
    if (myrand()<exp(-dE)){
      lattice[i][j].u1 =  Ut;
      acc+=1;
    }
  }

  return acc;  //0 se non ha accettato ne lo step temporale che quello spaziale
               //1 se accetta o uno o l'altro
               //2 se li accetta entrambi
}

//unitariziation
void project(site  **lattice,int Nt,int Ns){
  double complex tmp,U0,U1;
  long int i,j;
  for (i=0; i<Nt; i++){
    for (j=0;j<Ns;j++){
      U0= lattice[i][j].u0;
      tmp = U0/cabs(U0);
      lattice[i][j].u0 = tmp;
      U1= lattice[i][j].u1;
      tmp = U1/cabs(U1);
      lattice[i][j].u1 = tmp;
    }
  }
}

//initiazation of the lattice
void init_lattice(site **lattice,int Nt,int Ns){
  int k,l;
  for (k=0;k<Nt;k++){
    for(l=0;l<Ns;l++){
      lattice[k][l].u0 = 1+I*0;
      lattice[k][l].u1 = 1+I*0;
    }
  }
}
//main
int main(int argc, char **argv){
  site **lattice;
  double beta,eps;
  int i,j,k,l,Nt,Ns;
  long int acc = 0,iter;
  const unsigned long int seed1=(unsigned long int) time(NULL);
  const unsigned long int seed2=seed1+127;
  char datafile[100];
  FILE *fp;
  myrand_init(seed1, seed2);
  time_t start_time = time(NULL);

  if(argc != 6)
    {
    fprintf(stdout, "How to use this program:\n");
    fprintf(stdout, "  %s Nt Ns beta epsilon out_file \n", argv[0]);
    fprintf(stdout, "Output: \n");

    return EXIT_SUCCESS;
    }
  else
    {
    // read input values
    Nt=atoi(argv[1]);
    Ns=atoi(argv[2]);
    beta=atof(argv[3]);
    eps=atof(argv[4]);
    strcpy(datafile, argv[5]);
    }

    // open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

      lattice = malloc(Nt * sizeof(site*));
      if (lattice == NULL) {
          fprintf(stderr, "allocation problem at (%s,%d)\n", __FILE__, __LINE__);
          return EXIT_FAILURE;
      }

      for (i = 0; i < Nt; i++) {
          lattice[i] = malloc(Ns * sizeof(site));
          if (lattice[i] == NULL) {
              fprintf(stderr, "allocation problem at row %d\n", i);
              return EXIT_FAILURE;
          }
      }

    //initialize the lattice
    init_lattice(lattice,Nt,Ns);

    //update
    for(iter=0;iter<numb_of_upd;iter++){
      for(i=0; i<Nt; i++){
        for(j=0; j<Ns; j++){
          acc+=metropolis(lattice,i,j,eps,beta,Nt,Ns);
        }
      }
    }

    //free the lattice
    for (i = 0; i < Nt; i++)
        free(lattice[i]);

    free(lattice);


    //close the file
    fclose(fp);
    //running time
    time_t end_time = time(NULL);
    double time_spent = difftime(end_time,start_time);
    printf("Execution time = %lf \n",time_spent);
    printf("acceptance rate (with eps=%lf) = %lf\n",eps,(double)acc / (2*Nt*Ns*numb_of_upd));
    return EXIT_SUCCESS;
    }
