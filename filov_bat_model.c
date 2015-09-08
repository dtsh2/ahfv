/* pomp model file: seir filo */

  #include <R.h>
  #include <Rmath.h>
  #include <R_ext/Rdynload.h>
  #include <C:/Users/David Hayman/Documents/R/win-library/3.0/pomp/include/pomp.h>

// define parameters

  #define BETA      (p[parindex[0]]) 			// transmission rate
  #define MU        (p[parindex[1]]) 			// adult death rate
  #define DELTA     (p[parindex[2]]) 			// juvenile mortality rate
  #define SIGMA	    (p[parindex[3]]) 			// incubation period
  #define K	    	(p[parindex[4]]) 			// carrying capacity 
  #define EPSILON   (p[parindex[5]]) 			// rate of juvenile aging
  #define TAU	    (p[parindex[6]]) 			// seroconversion
  #define KAPPA	    (p[parindex[7]]) 			// birth rate (peak size)
  #define S	    	(p[parindex[8]]) 			// synchrony
  #define OMEGA	    (p[parindex[9]]) 			// pulse per yr
  #define PHI	    (p[parindex[10]]) 			// timing during year

// define states

  #define SUSJ       (x[stateindex[0]]) // number of susceptible juveniles
  #define EXPJ	     (x[stateindex[1]]) // number of exposed to infectious juveniles
  #define INFJ       (x[stateindex[2]]) // number of infected juveniles
  #define RECJ       (x[stateindex[3]]) // number of recovered juveniles
  #define SUSA       (x[stateindex[4]]) // number of susceptibles adults
  #define EXPA	     (x[stateindex[5]]) // number of exposed to infectious adults
  #define INFA       (x[stateindex[6]]) // number of infected adults
  #define RECA       (x[stateindex[7]]) // number of recovered adults

// the process model:
// an SIR model with Euler-multinomial step,

  void sir_euler_simulator (double *x, const double *p, 
			  const int *stateindex, const int *parindex, const int *covindex,
			  int covdim, const double *covar, 
			  double t, double dt)
{

  int nrate = 19; 										// number of rates
  double rate[nrate];									// transition rates
  double trans[nrate];									// transition numbers
  double N = x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7];	// population size
  double NAd = x[4]+x[5]+x[6]+x[7];						// Adult population size
  double NJu = x[0]+x[1]+x[2]+x[3];						// Juvenile population size
  void (*reulmult)(int,double,double*,double,double*);

// to evaluate the basis functions and compute the transmission rate, use some of 
// pomp's C-level eulermultinomial simulator

  reulmult = (void (*)(int,double,double*,double,double*)) R_GetCCallable("pomp","reulermultinom");

// define pi

  const double pi = 3.14159265359;

// in C --- pow(a,b) to do a^b 

// compute the transition rates

  rate[0] = (KAPPA*(1/sqrt((1/S)*pi)*exp(-pow((cos(pi*OMEGA*t-PHI)),2)/(1/S))))*(NAd);	// approx delta function birth into susceptible class
  rate[1] = EPSILON;				// aging from juv to ad
  rate[2] = BETA*(INFJ+INFA);		// sus to exposed 
  rate[3] = DELTA*(N/K);			// density dept mortality - Juvenile   
  rate[4] = SIGMA;					// incubation rate
  rate[5] = DELTA*(N/K);			// density dept mortality - Juvenile 
  rate[6] = EPSILON;				// aging from juv to ad
  rate[7] = TAU; 					// seroconversion rate 
  rate[8] = EPSILON;				// aging from juv to ad
  rate[9] = DELTA*(N/K);			// density dept mortality - Juvenile
  rate[10] = EPSILON;				// aging from juv to ad
  rate[11] = MU;					// mortality - Adult
  rate[12] = BETA*(INFJ+INFA);		// sus to exposed
  rate[13] = MU;					// mortality - Adult
  rate[14] = SIGMA;					// incubation rate
  rate[15] = MU;					// mortality - adult
  rate[16] = TAU; 					// seroconversion rate 
  rate[17] = MU;					// mortality - adult
  rate[18] = MU;					// mortality - adult
    
// compute the transition numbers  // in reulmult, first # is transitions, state name, rate # trans start, trans name at start

  trans[0] = rpois(rate[0]*dt);	               // births are Poisson
  (*reulmult)(3,SUSJ,&rate[1],dt,&trans[1]);   // euler-multinomial exits from SUSJ class
  (*reulmult)(3,EXPJ,&rate[4],dt,&trans[4]);   // euler-multinomial exits from EXPJ class
  (*reulmult)(3,INFJ,&rate[7],dt,&trans[7]);   // euler-multinomial exits from INFJ class
  (*reulmult)(2,RECJ,&rate[10],dt,&trans[10]); // euler-multinomial exits from RECJ class
  (*reulmult)(2,SUSA,&rate[12],dt,&trans[12]); // euler-multinomial exits from SUSA class
  (*reulmult)(2,EXPA,&rate[14],dt,&trans[14]); // euler-multinomial exits from EXPA class
  (*reulmult)(2,INFA,&rate[16],dt,&trans[16]); // euler-multinomial exits from INFA class
  (*reulmult)(1,RECA,&rate[18],dt,&trans[18]); // euler-multinomial exits from RECA class

// balance the equations

  SUSJ += trans[0]-trans[1]-trans[2]-trans[3];  	// IN births / OUT juv mortality, aging, inf exp
  EXPJ += trans[2]-trans[4]-trans[5]-trans[6];  	// IN inf exp / OUT incubation, juv mortality, aging
  INFJ += trans[4]-trans[7]-trans[8]-trans[9];  	// IN inc / OUT juv mortality, aging, seroconversion
  RECJ += trans[7]-trans[10]-trans[11]; 			// IN serocon /OUT juv mortality, aging
  SUSA += trans[1]-trans[12]-trans[13]; 			// IN aging / OUT mortality, inf exp
  EXPA += trans[12]+trans[6]-trans[14]-trans[15]; 	// IN inf exp, aging / OUT incubation, mortality
  INFA += trans[14]+trans[8]-trans[16]-trans[17]; 	// IN inc, aging / OUT seroconv, mortality
  RECA += trans[16]+trans[10]-trans[18]; 			// IN from serconverstion, aging / OUT from death

}
