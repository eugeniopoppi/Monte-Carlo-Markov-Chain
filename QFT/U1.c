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
  double complex S,Sn,Sp;
  int tp,xp,xn;

  //next neighbors
  tp = (i+1)%Nt;
  xp = (j+1)%Ns;
  xn = (j-1+Ns)%Ns;

  //staple in positive direction
  Sp = lattice[tp][j].u1
     * conj(lattice[i][xp].u0)
     * conj(lattice[i][j].u1);

  //staple in negative direction
  Sn = conj(lattice[tp][xn].u1)
     * conj(lattice[i][xn].u0)
     * lattice[i][xn].u1;
  S = Sn+Sp;

  return S;
}

//staple in the SPATIAL DIRECTION
double complex staple1(site **lattice,int i, int j,int Nt, int Ns){
  double complex S,Sp,Sn;
  int tp,tn,xp;

  //next neighbors
  tp = (i+1)%Nt;
  tn = (i-1+Nt)%Nt;
  xp = (j+1)%Ns;

  //staple in positive direction
  Sp = conj(lattice[tn][xp].u0)
     * conj(lattice[tn][j].u1)
     *  lattice[tn][j].u0;

  //staple in negative direction
  Sn = lattice[i][xp].u0
     * conj(lattice[tp][j].u1)
     * conj(lattice[i][j].u0);
  S = Sn + Sp;

  return S;
}

//metropolis
int metropolis(site **lattice, int i, int j,double eps, double beta, int Nt, int Ns){
  double complex Ut,Eold,Enew,Stmp;
  double theta,eold,enew,dE;
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

  if(dE<0){
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
  Ut = cexp(-I*theta);
  //Ut = (1+I*theta)/sqrt(1+theta*theta);

  Stmp = staple1(lattice,i,j,Nt,Ns);

  //old energy
  Eold = lattice[i][j].u1*Stmp;
  eold = -beta*creal(Eold);
  //new energy
  Enew = Ut*Stmp;
  enew = -beta*creal(Enew);
  //energy variation
  dE = enew-eold;


  if(dE<0){
    lattice[i][j].u1 =  Ut;
    acc+=1;
  }
  //accept-reject with exp{-beta dS}
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


//microcanonical update
void micro(site **lattice, int i, int j, int Nt, int Ns){
  double complex Uold,Unew,S;
  double absS;


  //update the temporal direction
  Uold = lattice[i][j].u0;
  S = staple0(lattice,i,j,Nt,Ns);
  absS = cabs(S);
  if(absS>0.000000000001){
    //Unew = conj(Uold*S*S);
    //Unew/=(absS*absS);
    Unew = conj(Uold)*(conj(S)/absS)*(conj(S)/absS);
    lattice[i][j].u0 = Unew;
  }

  //update the spacial direction
  Uold = lattice[i][j].u1;
  S = staple1(lattice,i,j,Nt,Ns);
  absS = cabs(S);
  if(absS> 1e-16){
    Unew = conj(Uold)*(conj(S)/absS)*(conj(S)/absS);
    lattice[i][j].u1 = Unew;
  }
}

//wilson loop
double wilson(site **lattice, int wt, int ws,int Nt,int Ns){
  double complex w,W;
  int t0,x0,t,x;

  W = 0+0*I;
  for(t0=0;t0<Nt;t0++){
    for(x0=0;x0<Ns;x0++){
        w = lattice[t0][x0].u0;
        for(t=1;t<wt;t++){
          w*= lattice[(t0+t)%Nt][x0].u0;
        }
        for(x=0;x<ws;x++){
          w*= lattice[(t0+wt)%Nt][(x0+x)%Ns].u1;
        }
        for(t=0;t<wt;t++){
          w*= conj(lattice[(t0+wt-1-t+Nt)%Nt][(x0+ws)%Ns].u0);
        }
        for(x=0;x<ws;x++){
          w*= conj(lattice[t0][(x0+ws-1-x+Ns)%Ns].u1);
        }
        W+=w;
      }
    }
    W/=(Nt*Ns);
  return creal(W);
}


//unitariziation
void project(site  **lattice,int Nt,int Ns){
  double complex U0,U1;
  double r;
  int i,j;
  for (i=0; i<Nt; i++){
    for (j=0;j<Ns;j++){
      U0= lattice[i][j].u0;
      r = cabs(U0);
      if(r > 1e-16) lattice[i][j].u0 = U0 / r;
      else lattice[i][j].u0 = 1.0 + I*0; // o altra scelta di fallback

      U1= lattice[i][j].u1;
      r = cabs(U1);
      if(r > 1e-16) lattice[i][j].u1 = U1 / r;
      else lattice[i][j].u1 = 1.0 + I*0; // o altra scelta di fallback
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
  double beta,eps,p,W,perc_of_metro=0.2;
  //double complex W;
  int i,j,Nt,Ns,numbofmet=0,wt,ws;
  long int acc = 0,iter;
  const unsigned long int seed1=(unsigned long int) time(NULL);
  const unsigned long int seed2=seed1+127;
  char datafile[100];
  FILE *fp;
  myrand_init(seed1, seed2);
  time_t start_time = time(NULL);

  if(argc != 8)
    {
    fprintf(stdout, "How to use this program:\n");
    fprintf(stdout, "  %s Nt Ns beta epsilon wt ws out_file \n", argv[0]);
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
    wt =atoi(argv[5]);
    ws = atoi(argv[6]);
    strcpy(datafile, argv[7]);
    }

    // open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    fprintf(fp, "#real(wilson) beta=%lf Nt/Ns=%d/%d wt/ws=%d/%d \n",beta,Nt,Ns,wt,ws);

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
      if(iter%1000==0){
        project(lattice,Nt,Ns);
      }
      p=myrand();
      if (p<perc_of_metro){
        for(i=0; i<Nt; i++){
          for(j=0; j<Ns; j++){
            acc+=metropolis(lattice,i,j,eps,beta,Nt,Ns);
            numbofmet+=2;    //+2 since metro updates two link variables
          }
        }
      }
      else{
          for(i=0; i<Nt; i++){
            for(j=0; j<Ns; j++){
              micro(lattice,i,j,Nt,Ns);
            }
           }

      }
      if(iter%50==0&&iter>10000){            //10000 of thermalization
        W=wilson(lattice,wt,ws,Nt,Ns)        //every 50 to take uncorrelated value
        fprintf(fp, "%lf \n", W);
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
    printf("acceptance rate (with eps=%lf) = %lf\n",eps,(double)acc / (numbofmet));
    return EXIT_SUCCESS;
    }
