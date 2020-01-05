from numpy import *




def KRL(t,tmax,c,r,d,fmax):

        f = (fmax*(power((((t+c)/(tmax+c))),r)/power(((d+(r*power((((t+c)/(tmax+c))),(1+r))))/(d+r)),((d+r)/(1+r)))))
        return f



def Norris(t,ts,A,tr,td):

       f = A*exp(2*(tr/ td)**(1/2) ) * exp( -tr / (t - ts) - (t - ts) / td )
       return f
