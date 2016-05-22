# arenavirus
pomp(
  data = data.frame(
    time=seq(from=0,to=365*20,by=1),  # time for simulations to run
    X = NA # dummy variables
  ),
  times='time',
  t0=0,
  rprocess=euler.sim(step.fun="sir_euler_simulator",delta.t=1),
  rmeasure="binomial_rmeasure",
  dmeasure="binomial_dmeasure",
  statenames=c(
    "SUSJ","CARJ","INFJ", "RECJ", "CARA", "SUSA","INFA","RECA"),
  states=c(
    "SUSJ","CARJ","INFJ", "RECJ", "CARA", "SUSA","INFA","RECA"),
  ## the order of the parameters assumed in the native routines:
  paramnames=c(
    'KAPPA',
    'S',
    "OMEGA",
    "PHI",
    'BETAAJ',
    'BETAJJ',
    'BETAC',
    'BETAAA',
    'BETAJA',
    'GAMMA',
    "DELTA",
    "MU",
    "TAU",
    "SUSJ.0","CARJ.0","INFJ.0", "RECJ.0", "CARA.0", "SUSA.0","INFA.0","RECA.0"),
  initializer=function(params,t0,states,...){
    x0<-params[paste(states,".0",sep="")]
    names(x0)<-states
    return(x0)
  }
) -> seir

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


params <- c(
  KAPPA=15/365,
  S=100,
  OMEGA=1/365,
  PHI=0.5,
  BETAAJ=0.001,
  BETAJJ=0.001,
  BETAC=0.001,
  BETAAA=0.001,
  BETAJA=0.001,
  GAMMA=0.1,
  DELTA=0.002,
  MU=0.05,
  TAU=0.001,
  SUSJ.0=2000,CARJ.0=1000,INFJ.0=1000,RECJ.0=1000,
  CARA.0=2000, SUSA.0=1000,INFA.0=1000, RECA.0=1000) # this adds to the initial conditions given the state variables

sim <- simulate(seir,params=c(params),seed=3493895L,
                nsim=10,states=T,obs=F,as.data.frame=T) # 
class(seir) # pomp object
class(sim) # data frame - even if I remove "as.data.frame" in the above code (sim)
#sim <- simulate(sir,params=c(params),nsim=1,states=T,obs=F)#,as.data.frame=T) # saves as an array
matplot(sim[1:7301,1:8],type='l')

