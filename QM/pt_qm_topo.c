#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include"/home/eugenio/Documenti/c/header/random.h"
#define PI 3.14159265358979323846
#define numb_of_upd 1000000

// stucture for the configuration
typedef struct conf {
  double *lattice;
  double  tau;
  long int Nt;
  double eta;

} conf;


//distance on a circle
double d(double tmp1, double tmp2){
  double dx = tmp1 - tmp2;
  if (fabs(dx) <= 0.5) return dx;
  else if (dx > 0.5)   return dx - 1.0;
  else                 return dx + 1.0;
}

//winding number
double top_charge(double *lattice,long int *nnp,long int Nt){
  double Q = 0.0;
  long int j,tmp;

  for (j=0;j<Nt;j++){
    tmp = nnp[j];
    Q += d(lattice[tmp],lattice[j]);
  }

  return round(Q);
}

//metropolis
long int metropolis(conf *configuration,long int *nnp,long int *nnn){
  double dE,x,xt,xn,xp,delta;
  long int s,idx,acc=0;

  delta = 2.5*sqrt(configuration->eta);
  for(s=0;s<configuration->Nt;s++){
    x = configuration->lattice[s];      //the site to be updated
    xt = configuration->lattice[s] + delta*(1.0-2.0*myrand()); //trial state, keep it in (0,1)
    if (xt < 0) xt += 1.0;
    if (xt >= 1.0) xt -= 1.0;
    idx = nnp[s];
    xn = configuration->lattice[idx];  //next site
    idx = nnn[s];
    xp = configuration->lattice[idx];  //previus site
    //energy variation
    dE = -d(xn,x)*d(xn,x)-d(x,xp)*d(x,xp)+d(xn,xt)*d(xn,xt)+d(xt,xp)*d(xt,xp);
    dE = dE/(2.0*(configuration->eta));
    //accept-reject metropolis
    if(dE<0){
      configuration->lattice[s] = xt;
      acc++;
    }
    else{
      if (myrand()<exp(-dE)){
        configuration->lattice[s]=  xt;
        acc++;
      }
    }
  }
  return acc;
}


// energy of the configuration
double calc_energy(conf *configuration,
                   long int *nnp)
  {
  long int r;
  double aux, ris;

  ris=0.0;
  for(r=0; r<configuration->Nt; r++)
     {
     aux=d(configuration->lattice[nnp[r]], configuration->lattice[r]);
     ris+=aux*aux;
     }
  ris/=(2.0*configuration->eta);

  return ris;
  }


// propose the swap of two configurations and return 1 if the swap is accpted
int swap_configuration(conf *config1, conf *config2, long int *nnp)
  {
  double energy1, energy2, energy1new, energy2new;
  double prob, *tmp;
  int itmp;

  energy1=calc_energy(config1, nnp);
  energy2=calc_energy(config2, nnp);

  energy2new=energy1*(config1->eta)/(config2->eta);
  energy1new=energy2*(config2->eta)/(config1->eta);

  prob=exp(-(energy1new-energy1)-(energy2new-energy2));
  if(myrand() < prob)
    {
    tmp=config2->lattice;
    config2->lattice=config1->lattice;
    config1->lattice=tmp;
    return 1;
    }

  return 0;
  }


//initialize the configurations
void init_configuration(conf *configuration,double tau,long int Nt){
  long int r;

  configuration->lattice=(double*)malloc((unsigned long int)(Nt)*sizeof(double));
  if(configuration->lattice == NULL)
    {
    fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<Nt; r++)
     {
     configuration->lattice[r]=0.5;
     }
  configuration->tau=tau;
  configuration->Nt=Nt;
  configuration->eta=tau/(double) Nt;
}


//main
int main(int argc, char **argv){
  conf *configuration;
  double tau,eta,delta,Q,tmptau,max_tau,dtau;  //tau=eta*N
  int i,j,l,k,n,rep,meas_every,numb_of_conf,swap_every=5;
  long int *nnp,*nnn,Nt,acc_local=0,acc_swap=0;
  const unsigned long int seed1=(unsigned long int) time(NULL);
  const unsigned long int seed2=seed1+127;
  myrand_init(seed1, seed2);
  char datafile[100];
  FILE *fp;
  time_t start_time = time(NULL);

  if(argc != 7)
    {
    fprintf(stdout, "How to use this program:\n");
    fprintf(stdout, "  %s tau N numb_of_conf max_tau measure_every out_file \n", argv[0]);
    fprintf(stdout, "Output: Q measured every -measure_every- complete update of the lattice\n");

    return EXIT_SUCCESS;
    }
  else
    {
    // read input values
    tau=atof(argv[1]);
    Nt=atol(argv[2]);
    eta = tau/Nt;
    numb_of_conf=atoi(argv[3]);
    max_tau = atof(argv[4]);
    meas_every=atoi(argv[5]);
    strcpy(datafile, argv[6]);
    }

  //allocate the array of configuations
  configuration=(conf*)malloc((unsigned long int)(numb_of_conf)*sizeof(conf));
  if(configuration == NULL)
    {
    fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }
  //allocate the next neighbor in positive and negative directions
  nnp=(long int *)malloc((unsigned long int)(Nt)*sizeof(long int));
  if(nnp == NULL)
    {
    fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }
  nnn=(long int *)malloc((unsigned long int)(Nt)*sizeof(long int));
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

  //initialize the next neighbors
  for(i=0; i<Nt; i++){
    nnp[i] = (i+1)%Nt;  //next neighbor in the positive direction
    nnn[i] = (i-1+Nt)%Nt; //next neighbor in the negative direction
  }

  //initialize the congigurations
  for(l=0;l<numb_of_conf;l++){
    //dtau  = 0.005 *pow(2, l);
    tmptau = tau*pow(max_tau/tau, l/((double)numb_of_conf-1.0));
    // tmptau = tau*pow(max_tau/((double) tau*(numb_of_conf-1.0)), l);
    init_configuration(&(configuration[l]),tmptau,Nt);
  }

  //update
  for(i=0; i<numb_of_upd; i++){
    for(n=0;n<numb_of_conf;n++){
      if(n==0){
        acc_local+=metropolis(&configuration[n],nnp,nnn); //metropolis(conf configuration,long int *nnp,long int *nnn)
      }
      else{
        (void) metropolis(&configuration[n],nnp,nnn);
      }
    }
    if(i%meas_every==0){
      Q = top_charge(configuration[0].lattice,nnp,Nt);
      fprintf(fp, "%lf \n",Q);
    }
    if(i%swap_every==0){
      if(myrand()>0.5)
        {
        for(rep=0; rep<numb_of_conf-1; rep++)
           {
           acc_swap+=swap_configuration(&(configuration[rep]), &(configuration[rep+1]), nnp);
           }
        }
      else
        {
        for(rep=numb_of_conf-1; rep>0; rep--)
           {
           acc_swap+=swap_configuration(&(configuration[rep]), &(configuration[rep-1]), nnp);
           }
        }
    }
  }

  // free the memory
  for(j=0; j<numb_of_conf; j++)
     {
     free(configuration[j].lattice);
     }
  free(nnp);
  free(nnn);

  //close the file
  fclose(fp);
  //running time
  time_t end_time = time(NULL);
  double time_spent = difftime(end_time,start_time);
  printf("Execution time = %lf \n",time_spent);
  printf("local acceptance rate = %lf\n",(double)acc_local / (Nt*numb_of_upd));
  // printf("" acceptance rate = %lf\n",(double)acc_swap / (swap_every*numb_of_upd*numb_of_conf));
  printf("swap acceptance rate = %lf\n",(double)acc_swap / ((numb_of_upd / swap_every) * (numb_of_conf - 1)));
  return EXIT_SUCCESS;
  }
