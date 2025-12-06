#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<complex.h>
#include"/home/eugenio/Documenti/c/header/random.h"
#define PI 3.14159265358979323846
#define MULTIHIT
#define THERM 10000
#define MEAS_EVERY 50
#define PROJECT_EVERY 1000
#define METRO_PERC 0.2

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

  //staple in positive spatial direction
  Sp = lattice[tp][j].u1
   * conj(lattice[i][xp].u0)
   * conj(lattice[i][j].u1);


  //staple in negative spatial direction
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

  //staple in negative temporal direction
  Sn = conj(lattice[tn][xp].u0)
     * conj(lattice[tn][j].u1)
     *  lattice[tn][j].u0;

  //staple in positive temporal direction
  Sp = lattice[i][xp].u0
     * conj(lattice[tp][j].u1)
     * conj(lattice[i][j].u0);

  S = Sn + Sp;

  return S;
}

//metropolis
int metropolis(site **lattice, int i, int j,double eps, double beta, int Nt, int Ns){
  double complex Ut,Eold,Enew,S;
  double theta,dE;
  int acc=0;

  //lets update TEMPORAL DIRECTION
  theta = eps*(2*myrand()-1);
  //Ut = lattice[i][j].u0*(1+I*theta)/sqrt(1+theta*theta);
  Ut = lattice[i][j].u0*cexp(-I*theta);
  S = staple0(lattice,i,j,Nt,Ns);

  //energy variation
  Eold = lattice[i][j].u0*S;
  Enew = Ut*S;
  dE = -beta*creal(Enew-Eold);

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

  //lets update SPATIAL DIRECTION
  theta = eps*(2*myrand()-1);
  Ut = lattice[i][j].u1*cexp(-I*theta);
  //Ut = lattice[i][j].u1*(1+I*theta)/sqrt(1+theta*theta);
  S = staple1(lattice,i,j,Nt,Ns);

  //energy variation
  Eold = lattice[i][j].u1*S;
  Enew = Ut*S;
  dE = -beta*creal(Enew-Eold) ;
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

  return acc;  //0 no step accepted
               //1 only one of the two spetes accepted
               //2 sboth accepted
}


//microcanonical update
void micro(site **lattice, int i, int j, int Nt, int Ns){
  double complex Uold,Unew,S;
  double absS;

  //update the temporal direction
  Uold = lattice[i][j].u0;
  S = staple0(lattice,i,j,Nt,Ns);
  absS = cabs(S);
  if(absS> 1e-16){
    Unew = conj(Uold)*(conj(S)/absS)*(conj(S)/absS);
    lattice[i][j].u0 = Unew;
  }

  //update the spatial direction
  Uold = lattice[i][j].u1;
  S = staple1(lattice,i,j,Nt,Ns);
  absS = cabs(S);
  if(absS> 1e-16){
    Unew = conj(Uold)*(conj(S)/absS)*(conj(S)/absS);
    lattice[i][j].u1 = Unew;
  }
}


//multihit
double complex multihit(site **lattice,int i,int j,int dir,double eps,double beta,int mh_upd,int Nt,int Ns){
  double complex U=0+0*I ,Uaux = 0+0*I,S,Utmp,Uprop,Eold,Enew,var;
  double theta,dE,absS;
  int hit=0,idx;

  //updates for the temporal link
  if(dir==0){
    Utmp = lattice[i][j].u0;
    S = staple0(lattice,i,j,Nt,Ns);
    absS  = cabs(S);
    for(idx=0;idx<mh_upd;idx++){
      //metropolis
      theta = eps*(2*myrand()-1);
      //variation to be proposed in the Utmp
      var = cexp(-I*theta); //* (1 + I*theta)/sqrt(1 + theta*theta);
      //enery variation
      Eold = Utmp * S;
      Enew = Utmp*var* S;
      dE = -beta * creal(Enew - Eold);

      if (dE < 0 || myrand() < exp(-dE)) {
          Utmp *= var ;
      }
      Uaux += Utmp;
      hit+=1;

     //microcanonical
     if(absS> 1e-10){
       Uprop = conj(Utmp)*(conj(S)/absS)*(conj(S)/absS);
       Uprop/=cabs(Uprop);
       Uaux+=Uprop;
       hit+=1;
     }
     else{
       Uprop = Utmp;
       hit+=1;
     }
     Utmp = Uprop;

   }
   U+= Uaux;
  }

  //updates for the spatial link
  if(dir==1){
    Utmp = lattice[i][j].u1;
    S = staple1(lattice,i,j,Nt,Ns);
    absS = cabs(S);
    for(idx=0;idx<mh_upd;idx++){
      //metropolis
      theta = eps*(2*myrand()-1);
      //variation proposed
      var = cexp(-I*theta);// * (1 + I*theta)/sqrt(1 + theta*theta);
      //enery veriation
      Eold = Utmp * S;
      Enew = Utmp*var* S;
      dE = -beta * creal(Enew - Eold);

      if (dE < 0 || myrand() < exp(-dE)) {
          Utmp*= var;
      }
      Uaux += Utmp;
      hit++;
     //microcanonical
     if(absS> 1e-10){
       Uprop = conj(Utmp)*(conj(S)/absS)*(conj(S)/absS);
       Uprop/=cabs(Uprop);
       Uaux+=Uprop;
       hit+=1;
     }
     else{
       Uprop = Utmp;
       hit+=1;
     }
     Utmp = Uprop;
   }
   U+= Uaux;
  }

  U/=hit;
  return U;
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


double wilson_multihit(site **lattice, int wt, int ws,double eps,double beta,int mh_upd,int Nt,int Ns){
  double complex w,W,U;
  int t0,x0,t,x,i,j;
  W = 0+0*I;
  if(wt==1 || ws==1){
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
  }
  else{
    for(t0=0;t0<Nt;t0++){
      for(x0=0;x0<Ns;x0++){
          w = 1+0*I;
          for(t=0;t<wt;t++){
            i = (t0+t)%Nt;
            j = x0;
            U = multihit(lattice,i,j,0,eps,beta,mh_upd,Nt,Ns);
            w*= U;
          }
          for(x=0;x<ws;x++){
            if(x==0 || x==ws-1){
                w*= lattice[(t0+wt)%Nt][(x0+x)%Ns].u1;
            }
            else{
              i = (t0+wt)%Nt;
              j = (x0+x)%Ns;
              U=multihit(lattice,i,j,1,eps,beta,mh_upd,Nt,Ns);
              w*=U;
            }
          }
          for(t=0;t<wt;t++){
            i = (t0+wt-1-t+Nt)%Nt;
            j = (x0+ws)%Ns;
            U = multihit(lattice,i,j,0,eps,beta,mh_upd,Nt,Ns);
            w*= conj(U);
          }
          for(x=0;x<ws;x++){
            if(x==0 || x==ws-1){
              w*= conj(lattice[t0][(x0+ws-1-x+Ns)%Ns].u1);
            }
            else{
              i=t0;
              j=(x0+ws-1-x+Ns)%Ns;
              U = multihit(lattice,i,j,1,eps,beta,mh_upd,Nt,Ns);
              w*= conj(U);
            }
          }
          W+=w;
        }
      }
      W/=(Nt*Ns);

  }
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
      if(r > 1e-13) lattice[i][j].u0 = U0 / r;
      else lattice[i][j].u0 = 1.0 + I*0;

      U1= lattice[i][j].u1;
      r = cabs(U1);
      if(r > 1e-13) lattice[i][j].u1 = U1 / r;
      else lattice[i][j].u1 = 1.0 + I*0;
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
  double beta,eps,p,W;
  //double complex W;
  int i,j,Nt,Ns,numbofmet=0,wt,ws,mh_upd;
  long int acc = 0,iter,numb_of_upd;
  const unsigned long int seed1=(unsigned long int) time(NULL);
  const unsigned long int seed2=seed1+127;
  char datafile[100];
  FILE *fp;
  myrand_init(seed1, seed2);
  time_t start_time = time(NULL);

  if(argc != 10)
    {
    fprintf(stdout, "How to use this program:\n");
    fprintf(stdout, "  %s Nt Ns beta epsilon wt ws numb_of_upd multihit_upd out_file \n", argv[0]);
    fprintf(stdout, " Nt Ns  the temporal and spatial extension of the lattice\n");
    fprintf(stdout, " wt ws the temporal and spatial extension of the silson loop\n");
    fprintf(stdout, "Output: Re(wilson(wt,ws) measured every %d\n",MEAS_EVERY);

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
    numb_of_upd = atoi(argv[7]);
    mh_upd = atoi(argv[8]);
    strcpy(datafile, argv[9]);
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
      if(iter%PROJECT_EVERY==0){
        project(lattice,Nt,Ns);
      }
      p=myrand();
      if (p<METRO_PERC){
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
      if(iter%MEAS_EVERY==0&&iter>THERM){
        #ifdef MULTIHIT
        W=wilson_multihit(lattice,wt,ws,eps,beta,mh_upd,Nt,Ns);
        #else
        W=wilson(lattice,wt,ws,Nt,Ns);
        #endif
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
 
