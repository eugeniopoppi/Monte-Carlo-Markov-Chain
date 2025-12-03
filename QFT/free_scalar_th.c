#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include"/home/eugenio/Documenti/c/header/random.h"
#define STRING_LENGTH 50
#define Nx 60
#define Nt 15


//computation of the observables epsilon/T^2 = O1+O2-O3
double obs1(double lattice[Nx][Nt],double hatm)
            {
            int x,t,n_new;
            double ris= 0.0;
            for(x=0; x<Nx; x++)
                {
                for(t=0; t<Nt; t++)
                  {
                  ris+=(lattice[x][t]*lattice[x][t]);
                  }
               }
             ris /= (double) Nt*Nx;
             ris *= (double) hatm*hatm;
             return ris;
            }
double obs2(double lattice[Nx][Nt])
            {
            int x,t,n_new;
            double ris= 0.0;
            double phi_old,phi_new;
            for(x=0; x<Nx; x++)
               {
                for(t=0; t<Nt; t++)
                  {
                  phi_old = lattice[x][t];
                  n_new =(x+1)%Nx;
                  phi_new=lattice[n_new][t];
                  ris += (phi_new - phi_old)*(phi_new - phi_old);
                  }
            }
            ris /= (double) Nt*Nx;
            return ris;
        }
double obs3(double lattice[Nx][Nt])
            {
            int x,t,n_new;
            double ris = 0.0;
            double phi_old,phi_new;
            for(x=0; x<Nx; x++)
               {
                for(t=0; t<Nt; t++)
                  {
                  phi_old = lattice[x][t];
                  n_new =(t+1)%Nt;
                  phi_new=lattice[x][n_new];
                  ris += (phi_new - phi_old)*(phi_new - phi_old);
                  }
            }
            ris /= (double) Nt*Nx;
            return ris;
        }

//define microcanonic
//S_n = sum_directions hat_phi_(n+direction)
void microcanonic(double lattice[Nx][Nt],
               double hatm,
               int x,
               int t)
    {int n_new;
     double phi_test,Sn;

     //phi_n^NEW = 2(S_n/(hatm^2+2D)) - phi_n
     //compute of Sn
     Sn=0.0;
     n_new =(x+1)%Nx;
     Sn+=lattice[n_new][t];
     n_new = (x - 1 + Nx) % Nx;
     Sn+=lattice[n_new][t];
     n_new=(t+1)%Nt;
     Sn+=lattice[x][n_new];
     n_new = (t - 1 + Nt) % Nt;
     Sn+=lattice[x][n_new];
     //now modify the lattice site
     phi_test = (2*Sn)/((hatm*hatm)+4) - lattice[x][t];
     lattice[x][t] = phi_test;
  }


//define metropolis
//S_L = 1/2 sum_n (hatm^2+ 2D) hat_phi^2 - 2 hat_phi S_n
int metropolis(double lattice[Nx][Nt],
               double hatm,
               int x,
               int t,
               double Delta)
  {int acc=0,n_new;
   double phi_test,Sn,dS,shift;
   shift = Delta*(myrand()*2-1);
   //define the trial configuration
   phi_test = lattice[x][t]+shift;
   //compute of Sn
   Sn=0.0;
   n_new =(x+1)%Nx;
   Sn+=lattice[n_new][t];
   n_new = (x - 1 + Nx) % Nx;
   Sn+=lattice[n_new][t];
   n_new=(t+1)%Nt;
   Sn+=lattice[x][n_new];
   n_new = (t - 1 + Nt) % Nt;
   Sn+=lattice[x][n_new];
   //now define dS
   dS = 0.5*(((hatm*hatm)+4)*((phi_test*phi_test)-(lattice[x][t]*lattice[x][t]))-2*Sn*shift);
   //accept-reject with exp{-DS}
   if ((double)(dS)<0){
     lattice[x][t] =  phi_test;
     acc = 1;
   }
   else{
     if (myrand()<exp(-dS)){
        lattice[x][t] =  phi_test;
        acc = 1;
        }
    else{
      lattice[x][t] = lattice[x][t];
    }

   }
  return acc;
  }


int main()
    {
    int numbofmet = 0, acce = 0,iter,x,t,N=1000000;
    double acc_perc ,hatm,shift,epst2;
    double Delta = 1,r,rho=0.2,lattice[Nx][Nt];
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;
    //char free_scalar;
    FILE *fp1,*fp2;
    myrand_init(seed1, seed2);
    hatm = 1/Nt;
    // initialize lattice to ordered start
    printf("Osserviamo la denistà di energia su T^2 per un campo scalare libero \n");
    printf("simulazione con Nx = %d e Nt = %d \n ",Nx,Nt);
    for(x=0; x<Nx; x++)
       {
       for(t=0; t<Nt; t++)
          {
          lattice[x][t]=0;
          }
       }

    // open data file for writing
    fp1=fopen("free_scalar.txt", "w");
    if(fp1==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n","free_scalar.txt" , __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    fp2=fopen("free_scalar_th_finiteT.txt", "w");
    if(fp2==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n","free_scalar_th.txt" , __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

      for (iter=0;iter<N;iter++){
        r = myrand();
        for (x=0;x<Nx;x++){
          for(t=0;t<Nt;t++){
              if(r<rho){
                numbofmet += 1;
                acce += metropolis(lattice,hatm,x,t,Delta);
              }
              else{
                microcanonic(lattice,hatm,x,t);
              }
          }
        }
        epst2 = obs1(lattice,hatm)+obs2(lattice)-obs3(lattice);
        fprintf(fp2, "%lf \n" ,epst2);
    }
    for(t=0;t<Nt;t++){
      for(x=0;x<Nx;x++){
     fprintf(fp1, "%lf " ,lattice[x][t]);
      }
    fprintf(fp1,"\n");
    }
    printf("shift = %lf \n",shift);
    printf("acc = %d \n",acce);
    printf("numbofmet = %d \n",numbofmet);
    acc_perc = acce/(double)numbofmet;
    printf("acc perc = %lf , w/ delta = %lf \n",acc_perc,Delta);
    // close the file
    fclose(fp1);
    fclose(fp2);

    return EXIT_SUCCESS;
    }
