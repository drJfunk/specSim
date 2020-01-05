#include "/Users/jburgess/Research/specSim/physPulse.h"

double physPulse(double t, double ene, double tf, double u0, double etaT, double etaR, double etaW)
{

  double eMin = 8.; // keV
  double eMax = 40000.; // keV


  double c = 2.99E10;

  double Gamma = 300.;
  double z=1.;
  //double tf = 1.;


  double beta = sqrt(1.-1./(Gamma*Gamma));

  double muJ = cos(4./Gamma);

  double deltaTprime = 2.*Gamma*etaT*tf/(1.+z);
  double deltaRprime = 2.*Gamma*etaW*c*tf/(1.+z);
  double deltaR = deltaRprime/Gamma;

  double dL = 2.02E28;



  double r0 = 2.* etaR * Gamma * Gamma * c * tf/(1.+z);


  double norm = c*u0/(6.*dL*dL*deltaRprime);
  

  double ruThis, rlThis;
  double term;
  double factor1 = 0.0;
  double factor2 = 0.0;
  // Integration params
  double nj=100000.;
  double mu ,dmu=(1.-muJ)/nj;


  




  //  for (ene=eMin; ene<=eMax; ene+=dEne )
  //  {

  
      for(mu=muJ; mu<=1.; mu+=dmu)
	{
	  
	  ruThis = ru(mu,t, beta, Gamma, r0, deltaTprime, deltaR, z);
	  rlThis = rl(mu,t, r0, beta, z);
	  
	  if(ruThis <= rlThis )
	    {
	      term = 0.0;
	    }
	  else
	    {
	      term = pow(Doppler(mu, Gamma, beta),3.)*(pow(ruThis,3.) - pow(rlThis,3.) )*spectrum(ene,mu,z,Gamma);
	      
	    }
	  
	  factor1 += dmu*term;
	  
	}

      factor1*=norm;

      //  factor2 += factor1*dEne;

      
      // }


      return factor1/(ene*1.62e-9);


}


double H(double x)
{

  if(x<0.)
    {
      return 0.;
    }
  else
    {
      return 1.;
    }

}


double Doppler(double mu, double Gamma, double beta)
{
  double val;

  val = 1./(Gamma*(1-beta*mu));

  return val;

}

double spectrum(double ene, double  mu, double z, double Gamma)
{
  double a=4./3.;
  double b=-1./2.;

  double beta = sqrt(1.-1./(Gamma*Gamma));

  double ep0 = 200.;
  double epPrime = (1.+z)*ep0/(2.*Gamma);
  double x = (ene/epPrime) * (1.+z)/Doppler(mu,Gamma,beta);


  double val;

  val = pow(x,a)*H(1.-x) + pow(x,b)*H(x-1.);

  return val;

}



double rl(double mu, double t, double r0, double beta, double z)
{
  double tz = t/(1.+z);
  double a,b;
  double c = 2.99E10;

  a = beta*c*tz/(1.-beta*mu);
  b = (r0/beta - c*tz)/mu;

  if (a>b)
    {
      return a;
    }
  else
    {
      return b;
    }


}

double ru(double mu, double t, double beta, double Gamma, double r0, double deltaTprime, double deltaR, double z)
{

  double tz = t/(1.+z);
  double a,b;
  double c = 2.99E10;

  a=((beta*c*tz) + deltaR)/(1. - beta*mu);
  b=(r0/beta - c*tz + c*Gamma*deltaTprime)/mu;

  if (a<b)
    {
      return a;
    }
  else
    {
      return b;
    }


}

