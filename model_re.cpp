#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string>
#include<fstream>

#define DT 1.0e-2
#define OUTDT 1.0e+1
#define TMAX 1.0e+2

#define NUM_PARTICLE  10
#define L0 1.0
#define P1 1.0e-3
#define rt 1.0e-4
#define r0rt 10e+4
#define DT 1.0e-2
#define run_max 100
#define ITMAX ((int)((TMAX)/(DT)))

double Uniform(){
  return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double randGauss(double value1, double value2, double sigma){
  return sigma*sqrt(-2.0*log(value1))*sin(2.0*M_PI*value2);
}

typedef struct{
  double x;
}S_DIS;

double Force_Calc(S_DIS par[NUM_PARTICLE], int i_number){
  double ll = par[i_number + 0].x - par[i_number - 1].x;
  double lr = par[i_number + 1].x - par[i_number + 0].x;
  return (L0 - ll) + (lr - L0);
}


int main(void){

  srand((unsigned)time(NULL));

  S_DIS par[NUM_PARTICLE];
  double F[NUM_PARTICLE];
  double V[NUM_PARTICLE];
  double sigma[ITMAX][run_max];
  double ave[ITMAX][run_max];


  FILE* fp0;
  char filename[256];
  sprintf(filename, "n%d_r0rt%.0f.dat", NUM_PARTICLE,r0rt);
  if((fp0 = fopen(filename, "w")) == NULL){printf("FAILED TO OPEN FILE.\n"); exit(1);};

  for(int run=0; run<run_max; run++){

    for(int ini=0; ini<NUM_PARTICLE; ini++){
      
      par[ini].x = (double)ini*L0;
      //printf("%f\t",par[ini].x);

      if(ini != 0 and ini != NUM_PARTICLE - 1){
        par[ini].x += r0rt*randGauss(Uniform(), Uniform(), rt*L0);
      }
      //printf("%f\n",par[ini].x);
    }


    for(int iter=0; iter<ITMAX; iter++){ 
     //printf("%d",iter); 
        for(int ini=1; ini<NUM_PARTICLE-1; ini++){
        F[ini] = Force_Calc(par, ini);
        V[ini] = F[ini];
        }

        for(int num=1; num<NUM_PARTICLE - 1; num++){
          par[num].x += V[num]*DT + randGauss(Uniform(), Uniform(), rt*L0)*sqrt(DT);
        }

        double sum_l = 0.0;
        double sum_l2 = 0.0;
        for(int num=0; num<NUM_PARTICLE - 1; num++){
          double dl = par[num + 1].x - par[num + 0].x - L0;
          sum_l += dl;
          sum_l2 += dl*dl;
          //pintf("%d\n",dl);
        }
        ave[iter][run] = sum_l/(double)(NUM_PARTICLE - 1);
        sigma[iter][run] = sqrt(sum_l2/(double)(NUM_PARTICLE - 1) - ave[iter][run]*ave[iter][run]);
       /* 
        for(int par_num=0; par_num<NUM_PARTICLE; par_num++){
          printf("n=%d\t%f\n",par_num,par[par_num]);
        }
    */
    }  
  }
  
  for(int t=0; t<ITMAX; t++){
    double SD_ave = 0.0;
    double SD_ave2 = 0.0; 
    double SD_SD;
    for(int r=0; r<run_max; r++){
      SD_ave += sigma[t][r];
      SD_ave2 += sigma[t][r]*sigma[t][r];
    }
    double SD = SD_ave/(double)run_max;
    SD_SD = sqrt(SD_ave2/(double)run_max - SD*SD);

    double t_real = (double)t*DT;
    //printf("t=%f\tave=%f\tSD=%15f\n",t_real,SD,SD_SD);
    fprintf(fp0, "%f\t%f\t%f\n",t_real,SD,SD_SD);
    
  }
  fclose(fp0);

  return 0;
}
  
