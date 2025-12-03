#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<complex.h>
#include"/home/eugenio/Documenti/c/header/random.h"
#define L 5
#define PI 3.14159265358979323846
#define N 1000000
//#define B 0.3

// magnetization per site
double complex magn(double complex lattice[L][L][L])
  {double complex sum;
   int rx,ry,rz;

  sum=0+0*I;
  for(rx=0; rx<L;rx++){
    for(ry=0; ry<L;ry++){
      for(rz=0;rz<L;rz++){
        sum+=lattice[rx][ry][rz];
      }
    }
  }

 return  sum / (L*L*L);
 }
 // // energy per site
double energy(double complex lattice[L][L][L]){
 long int rx, ry,rz, tmp, r_new;
 double sum;
 sum=0;
 for(rx=0;rx<L;rx++){
   for(ry=0;ry<L;ry++){
     for(rz=0;rz<L;rz++){
       r_new=(rx+1)%L;
       sum+=creal(lattice[rx][ry][rz])*creal(lattice[r_new][ry][rz]);
       sum+=cimag(lattice[rx][ry][rz])*cimag(lattice[r_new][ry][rz]);

       r_new=(ry+1)%L;
       sum+=creal(lattice[rx][ry][rz])*creal(lattice[rx][r_new][rz]);
       sum+=cimag(lattice[rx][ry][rz])*cimag(lattice[rx][r_new][rz]);

       r_new=(rz+1)%L;
       sum+=creal(lattice[rx][ry][rz])*creal(lattice[rx][ry][r_new]);
       sum+=cimag(lattice[rx][ry][rz])*cimag(lattice[rx][ry][r_new]);
      }
    }
  }
 return  -sum / (double) (L*L*L);
 }


 //microcanonic:has to flip wrt Sr
 void microcanonic(double complex lattice[L][L][L],
                int rx,
                int ry,
                int rz)
   {double complex sr,Sr,St;
    double modSr;
    int r_new;
    //compute of first closests' energy sum
    sr = lattice[rx][ry][rz];
    Sr=0;
    r_new=(rx+1)%L;
    Sr+=lattice[r_new][ry][rz];
    r_new = (rx - 1 + L) % L;
    Sr+=lattice[r_new][ry][rz];
    r_new=(ry+1)%L;
    Sr+=lattice[rx][r_new][rz];
    r_new = (ry - 1 + L) % L;
    Sr+=lattice[rx][r_new][rz];
    r_new=(rz+1)%L;
    Sr+=lattice[rx][ry][r_new];
    r_new = (rz - 1 + L) % L;
    Sr+=lattice[rx][ry][r_new];
    modSr = cabs(Sr);
    if(modSr<0.00000000001){
      St = sr;
    }
    else{
      St = (2*(creal(sr)*creal(Sr)+cimag(sr)*cimag(Sr))*Sr)/(modSr*modSr)-sr;
      lattice[rx][ry][rz] = St;
    }
  }

// metropolis update at site r
// return 1 if accepted, else 0
int metropolis(double complex lattice[L][L][L],
               int rx,
               int ry,
               int rz,
               double theta,
               double B)
  {int acc=0;
   double complex Sr,St;
   double complex ds;
   int r_new;

   //define the trial configuration
   St = lattice[rx][ry][rz]*cexp(-I*theta);
   ds = lattice[rx][ry][rz]*(cexp(-I*theta)-1);
   //compute of first closests' energy sum
   Sr=0;
   r_new=(rx+1)%L;
   Sr+=lattice[r_new][ry][rz];
   r_new = (rx - 1 + L) % L;
   Sr+=lattice[r_new][ry][rz];
   r_new=(ry+1)%L;
   Sr+=lattice[rx][r_new][rz];
   r_new = (ry - 1 + L) % L;
   Sr+=lattice[rx][r_new][rz];
   r_new=(rz+1)%L;
   Sr+=lattice[rx][ry][r_new];
   r_new = (rz - 1 + L) % L;
   Sr+=lattice[rx][ry][r_new];
   Sr = (creal(ds)*creal(Sr)+cimag(ds)*cimag(Sr));
   Sr = -Sr;
   if((double)(Sr)<0){
     lattice[rx][ry][rz] =  St;
     acc = 1;
   }
   //accept-reject with exp{-beta DE}
   else{
     if (myrand()<exp(-B*Sr)){
       lattice[rx][ry][rz] =  St;
       acc = 1;
     }
   }

  return acc;
  }


int main(){
  double complex lattice[L][L][L];
  double complex M;
  double r,rho = 0.2, eps = PI,theta = 0,B,E;
  int i,j,k,acc,numbofmet;
  int rx,ry,rz,iter;
  const unsigned long int seed1=(unsigned long int) time(NULL);
  const unsigned long int seed2=seed1+127;
  char datafile[50];
  FILE *fp;
  long long V;
  V = L*L*L* (long long) N;
  time_t start_time = time(NULL);

  //initialize random number generator
  myrand_init(seed1, seed2);
  //initialize the  complex lattice
  for (i=0;i<L;i++){
    for(j=0;j<L;j++){
      for(k=0;k<L;k++){
        lattice[i][j][k] = 1 ; //
      }
    }
  }
for (B=0.400;B<0.422;B += 0.002){
  //open the file to save the result on
  sprintf(datafile, "xyL%dB%.4f.txt",L,B);
  fp=fopen(datafile,"w");
  if (fp == NULL)
     {
      perror("Errore di apertura del file");
      exit(1);
     }
  fprintf(fp, "#real(magnetization) imag(magnetization) energy \n");
  acc = 0;
  numbofmet = 0;
  for (iter=0;iter<N;iter++){
    r = myrand();
    for (i=0;i<L;i++){
      for(j=0;j<L;j++){
        for(k=0;k<L;k++){
          rx=i;
          ry=j;
          rz=k;
          if(r<rho){
            theta = (myrand()*2*eps)-eps;
            acc += metropolis(lattice,rx,ry,rz,theta,B);
            numbofmet += 1;
          }
          else{
            microcanonic(lattice,rx,ry,rz);
          }
        }
      }
    }
    M = magn(lattice);
    E = energy(lattice);
    fprintf(fp, "%lf %lf %lf \n", creal(M),cimag(M),E);
  }
fprintf(fp,"#microcanonic-metropolis percentage = %lf-%lf \n",((V-numbofmet)/(double)V),(numbofmet/(double)V));
fprintf(fp,"#acc = %lf\n",acc/(double)numbofmet); //(double)
fclose(fp);
}
  printf("last microcanonic-metropolis percentage = %lf-%lf \n",((V-numbofmet)/(double)V),(numbofmet/(double)V));
  printf("last acc = %lf",acc/(double)numbofmet); //(double)
  time_t end_time = time(NULL);
  double time_spent = difftime(end_time,start_time);
  printf("Execution time = %lf \n",time_spent);
  return EXIT_SUCCESS;
}
