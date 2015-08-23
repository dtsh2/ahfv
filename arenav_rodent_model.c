/* pomp model file: seir arena */

  #include <R.h>
  #include <Rmath.h>
  #include <R_ext/Rdynload.h>
  #include <C:/Users/David Hayman/Documents/R/win-library/3.0/pomp/include/pomp.h>

// define parameters
// TO DO CHANGE NAMES

  #define KAPPA	    (p[parindex[0]]) 			// birth rate (peak size)
  #define S	    	(p[parindex[1]]) 			// synchrony
  #define OMEGA	    (p[parindex[2]]) 			// pulse per yr
  #define PHI	    (p[parindex[3]]) 			// timing during year
  #define BETAAJ    (p[parindex[4]]) 			// transmission rate ad-juv
  #define BETAJJ    (p[parindex[5]]) 			// transmission rate juv-juv
  #define BETAC     (p[parindex[6]]) 			// transmission rate persistently infected
  #define BETAAA    (p[parindex[7]]) 			// transmission rate ad-ad
  #define BETAJA    (p[parindex[7]]) 			// transmission rate juv-ad
  #define GAMMA	    (p[parindex[8]]) 			// aging from juv to adult
  #define DELTA     (p[parindex[9]]) 			// juvenile mortality rate
  #define MU        (p[parindex[10]]) 			// adult death rate
  #define TAU	    (p[parindex[11]]) 			// seroconversion rate

// define states

  #define SUSJ       (x[stateindex[0]]) // number of susceptible juveniles
  #define CARJ	     (x[stateindex[1]]) // number of exposed to infectious juveniles
  #define INFJ       (x[stateindex[2]]) // number of infected juveniles
  #define RECJ       (x[stateindex[3]]) // number of recovered juveniles
  #define CARA	     (x[stateindex[4]]) // number of exposed to infectious adults
  #define SUSA       (x[stateindex[5]]) // number of susceptibles adults
  #define INFA       (x[stateindex[6]]) // number of infected adults
  #define RECA       (x[stateindex[7]]) // number of recovered adults

// the process model:
// an SIR model with Euler-multinomial step,

  void sir_euler_simulator (double *x, const double *p, 
			  const int *stateindex, const int *parindex, const int *covindex,
			  int covdim, const double *covar, 
			  double t, double dt)
{

  int nrate = 16; 					// number of rates
  double rate[nrate];					// transition rates
  double trans[nrate];					// transition numbers
  double NCa = x[1]+x[4];			// Carrier population size
  double NIa = x[4]+x[6];			// Infected adult population size
  double NAd = x[5]+x[7];			// Breeding adult population size (non-infected)
  void (*reulmult)(int,double,double*,double,double*);

// to evaluate the basis functions and compute the transmission rate, use some of 
// pomp's C-level eulermultinomial simulator

  reulmult = (void (*)(int,double,double*,double,double*)) R_GetCCallable("pomp","reulermultinom");

// define pi

  const double pi = 3.14159265359;

// in C --- pow(a,b) to do a^b 

// compute the transition rates

  rate[0] = (KAPPA*(1/sqrt((1/S)*pi)*exp(-pow((cos(pi*OMEGA*t-PHI)),2)/(1/S))))*(NAd);	// approx delta function birth into susceptible class
  rate[1] = GAMMA;			// aging from juv to ad
  rate[2] = BETAJA*INFA+BETAJJ*INFJ+BETAC*NCa;		// sus to infectious
  rate[3] = DELTA;			// mortality - Juvenile   
  rate[4] = (KAPPA*(1/sqrt((1/S)*pi)*exp(-pow((cos(pi*OMEGA*t-PHI)),2)/(1/S))))*(NIa);	// approx delta function birth into carrier class
  rate[5] = GAMMA;			// aging from juv to ad
  rate[6] = DELTA;			// mortality - Juvenile 
  rate[7] = TAU; 			// seroconversion rate 
  rate[8] = DELTA;			// mortality - Juvenile
  rate[9] = GAMMA;			// aging from juv to ad
  rate[10] = DELTA;			// mortality - Juvenile 
  rate[11] = MU;			// mortality - Adult
  rate[12] = BETAAJ*INFJ+BETAAA*INFA+BETAC*NCa;		// sus to exposed- infectious route
  rate[13] = MU;			// mortality - Adult
  rate[14] = TAU; 			// seroconversion rate 
  rate[15] = MU;			// mortality - adult
  rate[16] = MU;			// mortality - adult
    
// compute the transition numbers  // in reulmult, first # is transitions, state name, rate # trans start, trans name at start

  trans[0] = rpois(rate[0]*dt);	               // births are Poisson
  (*reulmult)(3,SUSJ,&rate[1],dt,&trans[1]);   // euler-multinomial exits from SUSJ class
  trans[4] = rpois(rate[4]*dt);	               // births are Poisson
  (*reulmult)(2,CARJ,&rate[5],dt,&trans[5]);   // euler-multinomial exits from CARJ class
  (*reulmult)(2,INFJ,&rate[7],dt,&trans[7]);   // euler-multinomial exits from INFJ class
  (*reulmult)(2,RECJ,&rate[9],dt,&trans[9]);   // euler-multinomial exits from RECJ class
  (*reulmult)(1,CARA,&rate[11],dt,&trans[11]); // euler-multinomial exits from CARA class
  (*reulmult)(2,SUSA,&rate[12],dt,&trans[12]); // euler-multinomial exits from SUSA class
  (*reulmult)(2,INFA,&rate[14],dt,&trans[14]); // euler-multinomial exits from INFA class
  (*reulmult)(1,RECA,&rate[16],dt,&trans[16]); // euler-multinomial exits from RECA class

// balance the equations

  SUSJ += trans[0]-trans[1]-trans[2]-trans[3];  	// IN births; OUT juv mort, aging, inf
  CARJ += trans[4]-trans[5]-trans[6];  	// IN births; OUT aging, juv mortality
  INFJ += trans[2]-trans[7]-trans[8];  	// IN inf j / OUT seroconversion, juv mortality
  RECJ += trans[7]-trans[9]-trans[10]; 			// IN sero j /OUT  aging, juv mortality
  CARA += trans[5]-trans[11]-trans[12]; 	// IN aging /OUT ad mortality
  SUSA += trans[1]-trans[12]-trans[13]; 			// IN aging / OUT inf, ad mortality
  INFA += trans[14]-trans[14]-trans[15]; 	// IN inf / OUT aging, juv mortality
  RECA += trans[14]+trans[9]-trans[16]; 			// IN from serconverstion & aging; OUT from death

}
