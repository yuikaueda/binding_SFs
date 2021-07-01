#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define DT              1.0e-2
#define OUTDT           1.0e-1
#define TMAX            1.0e+1
#define MAXITER         ((int)((TMAX)/(DT)))
#define OUTGAP          ((int)((OUTDT)/(DT)))

#define NUM_PARTICLE    10
#define L0              1.0
//#define ETA             1.0
//#define K               1.0e-1
#define NUM_HISTOGRAM   25
#define P1		1.0e-3
//#define P2		0.001

double randRange(double min, double max){return (rand()/(RAND_MAX/(max - min)) + min);}
double randBasic(){return (rand()/(double)(RAND_MAX));}
double randGauss(double value1, double value2, double sigma){return sigma*sqrt(-2.0*log(value1))*sin(2.0*M_PI*value2);}

class Fiber
{
  public:
    Fiber();
    ~Fiber();

  public:
    void    CalcForce();      // calculate force
    void    CalcVelocity();   // calculate velocity
    void    Update();         // update position
    void    Output(int iter); // output
    double  mu_, sigma_, *x_; 
    int     numParticle_;
    double  l0_;  // unit length
  private:
    double   *f_, *v_;//, ;*x_;
    double  k_;   // spring constant
    double  lt_;  // total length
    double  eta_; // viscosity
//    double  l0_;  // unit length

};

Fiber::Fiber()
{
  numParticle_ = NUM_PARTICLE; 

  l0_  = L0;
  lt_  = l0_*(numParticle_ - 1);
//  eta_ = ETA;
//  k_   = K;

  x_ = new double[numParticle_];
  f_ = new double[numParticle_];
  v_ = new double[numParticle_];

  srand((unsigned) time(NULL));
  for(int idParticle = 0; idParticle < numParticle_; idParticle++){
//    x_[idParticle] = idParticle*l0_; //Nosyokinoise
    
//    if(idParticle != 0 and idParticle != numParticle_ - 1) // if not wall
//     x_[idParticle] += randGauss(randBasic(), randBasic(), 0.1*l0_);

    if(idParticle != 0 and idParticle != numParticle_ - 1) // if not wall
     x_[idParticle] += 1.0e+2*P1*randGauss(randBasic(), randBasic(), 1)*sqrt(DT);
  }
}

Fiber::~Fiber()
{
  delete [] x_;
  delete [] f_;
}

void Fiber::CalcForce()
{
  for(int idParticle = 1; idParticle < numParticle_ - 1; idParticle++){ // loop excluding wall
    double  ll = x_[idParticle + 0] - x_[idParticle - 1];  
    double  lr = x_[idParticle + 1] - x_[idParticle + 0];  

//     f_[idParticle] = k_*((l0_ - ll) + (lr - l0_)); 
      f_[idParticle] = ((l0_ - ll) + (lr - l0_));
  }
}

void Fiber::CalcVelocity()
{
  for(int idParticle = 1; idParticle < numParticle_ - 1; idParticle++){ // loop excluding wall
   // v_[idParticle] = f_[idParticle]/eta_;
    v_[idParticle] = f_[idParticle];
  }
}

void Fiber::Update()
{
  for(int idParticle = 1; idParticle < numParticle_ - 1; idParticle++){ // loop excluding wall
    //x_[idParticle] += v_[idParticle]*DT + randGauss(randBasic(), randBasic(), 0.00*l0_);
    x_[idParticle] += v_[idParticle]*DT + P1*randGauss(randBasic(), randBasic(), 1.0)*sqrt(DT);
  }
}

void Fiber::Output(int iter)
{
  if(iter%OUTGAP != 0) return;

  FILE *fp;
  char filename[256];

  sprintf(filename, "result/result%05d.dat", iter/OUTGAP);
  if((fp = fopen(filename, "w")) == NULL){printf("FAILED TO OPEN FILE.\n"); exit(1);};
  for(int idParticle = 0; idParticle < numParticle_; idParticle++)
    fprintf(fp, "%15e %15e\n", x_[idParticle], x_[idParticle]/lt_);
  fclose(fp);

  double lmax = l0_*1.05;
  double lmin = l0_*0.95;
  double dl   = (lmax - lmin)/NUM_HISTOGRAM;
  int    histogram[NUM_HISTOGRAM] = {0};

  for(int idParticle = 0; idParticle < numParticle_ - 1; idParticle++){ // loop for spring
    double  length = x_[idParticle + 1] - x_[idParticle + 0];

    for(int idHistogram = 0; idHistogram < NUM_HISTOGRAM; idHistogram++){
      double t_lmin = dl*(idHistogram + 0) + lmin; 
      double t_lmax = dl*(idHistogram + 1) + lmin; 

      if(t_lmin < length and length < t_lmax) histogram[idHistogram]++;
    }
  }






  sprintf(filename, "result/count%05d.dat", iter/OUTGAP);
  if((fp = fopen(filename, "w")) == NULL){printf("FAILED TO OPEN FILE.\n"); exit(1);};

  for(int idHistogram = 0; idHistogram < NUM_HISTOGRAM; idHistogram++){
    double t_lmean = dl*((double)idHistogram + 0.5) + lmin; 

    fprintf(fp, "%15e %15e\n", t_lmean/l0_, histogram[idHistogram]/(double)(numParticle_ - 1));
  }

  fclose(fp);

}

int main()
{
  Fiber fiber;
  FILE *fp;
  char filename[256];
  sprintf(filename, "2sigma-10^2.%d.dat", NUM_PARTICLE);
  if((fp = fopen(filename, "w")) == NULL){printf("FAILED TO OPEN FILE.\n"); exit(1);};


  for(int iter = 0; iter <= MAXITER; iter++){
    fiber.CalcForce();
    fiber.CalcVelocity();

    fiber.Output(iter);

    fiber.Update();

    if(iter%OUTGAP == 0){
      double sum_l  = 0.0;
      double sum_l2 = 0.0;
      for(int idParticle = 0; idParticle < fiber.numParticle_ - 1; idParticle++){
        double dl = fiber.x_[idParticle + 1] - fiber.x_[idParticle + 0] - fiber.l0_;
        sum_l += dl;
        sum_l2 += dl*dl;
        }

      fiber.mu_    = sum_l/ (double) (fiber.numParticle_ - 1)  ;
      fiber.sigma_ = sqrt(sum_l2/(double) (fiber.numParticle_ - 1) - fiber.mu_*fiber.mu_) ;

      double t = (double)iter *DT;
      fprintf(fp, "%15e %15e\n", t, fiber.sigma_);
      }  
     }
  fclose(fp);


  return 0;
}
