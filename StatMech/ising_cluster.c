#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "/home/eugenio/Documenti/c/header/random.h"
#define L 4
#define STRING_LENGTH 50

// Funzione per resettare il reticolo ausiliario a zero
void reset_lattice(int auxlatt[L][L]) {
    memset(auxlatt, 0, sizeof(auxlatt[0][0]) * L * L);
}

// Funzione per la crescita del cluster
int cluster(int latt[L][L], int auxlatt[L][L], int x, int y, double p) {
    int i, xt, yt;
    int dir_x[4] = {0, 1, 0, -1}; // Direzioni in x
    int dir_y[4] = {1, 0, -1, 0}; // Direzioni in y
    int clustdim = 1;
    auxlatt[x][y] = 1;
    for (i = 0; i < 4; i++) {
        xt = (x + dir_x[i] + L) % L; // Bordo periodico in x
        yt = (y + dir_y[i] + L) % L; // Bordo periodico in y
        if (latt[x][y] == latt[xt][yt] && auxlatt[xt][yt] == 0) {
            if (myrand() <= p) {
                auxlatt[xt][yt] = 1;
                clustdim += cluster(latt, auxlatt, xt, yt, p); // Chiamata ricorsiva
            }
        }
    }
return clustdim;
}

// Funzione per calcolare la magnetizzazione
double magn(int lattice[L][L]) {
    int i, j, sum = 0;

    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            sum += lattice[i][j];
        }
    }
    return (double) sum / (double) (L * L);
}

int main(int argc, char **argv) {
    int i, j,k, rx, ry, clustdim = 0, latt[L][L], auxlatt[L][L];
    long int N, iter;
    double B, M, seed1, seed2, p;
    char datafile[STRING_LENGTH];
    FILE *fp;

    if (argc != 4) {
        fprintf(stdout, "How to use this program:\n");
        fprintf(stdout, "  %s beta iterations filename \n\n", argv[0]);
        fprintf(stdout, "  beta = double : 1/(K_b*T) \n");
        fprintf(stdout, "  iterations = long int : the numbers of draws to be generated\n");
        fprintf(stdout, "  filename = name of the file on which the result is written\n");
        return EXIT_SUCCESS;
    } else {
        // Leggi i valori di input
        B = atof(argv[1]);
        N = atol(argv[2]);

        // Controlla la lunghezza del nome del file
        if (strlen(argv[3]) >= STRING_LENGTH) {
            fprintf(stderr, "File name too long. Increase STRING_LENGTH or shorten the name (%s, %d)\n", __FILE__, __LINE__);
            return EXIT_FAILURE;
        } else {
            strcpy(datafile, argv[3]);
        }
    }

    // Apri il file per scrivere i risultati
    fp = fopen(datafile, "w");
    if (fp == NULL) {
        perror("Errore di apertura del file");
        exit(1);
    }
    //per calcolare quanto tempo impiega a girare il codice
    time_t start_time = time(NULL);
    // Inizializza il seme per myrand
    seed1 = time(NULL);
    seed2 = seed1 + 127;
    myrand_init(seed1, seed2);

    // Definisci la probabilità di accettazione
    p = 1 - exp(-2 * B);

    // Inizializza il reticolo
    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            latt[i][j] = 1;
        }
    }

    // Inizializza il reticolo ausiliario
    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            auxlatt[i][j] = 0;
        }
    }

    // Ciclo principale delle iterazioni
    for (iter = 0; iter < N; iter++) {
          rx=(int)((double)L * myrand());
          ry=(int)((double)L * myrand());
          clustdim = cluster(latt, auxlatt, rx, ry, p);

          for (i = 0; i < L; i++) {
              for (j = 0; j < L; j++) {
                if (1 == auxlatt[i][j]){
                  latt[i][j]  = -latt[i][j];
                }    
              }
          }
          // Resetta il reticolo ausiliario
          //reset_lattice(auxlatt);
          for (i = 0; i < L; i++) {
              for (j = 0; j < L; j++) {
                  auxlatt[i][j] = 0;
              }
           }

        M = magn(latt);
        fprintf(fp, "%lf \n", M);

    }

    // Chiudi il file dei dati
    fclose(fp);

    // Mostra il tempo di esecuzione
    time_t end_time = time(NULL);
    double time_spent = difftime(end_time, start_time);
    printf("last cluster dimension = %d \n", clustdim);
    printf("Execution time = %lf \n", time_spent);

    return EXIT_SUCCESS;
}
