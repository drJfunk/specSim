from scipy.integrate import quad, quadrature
from numpy import dot, array, zeros

from numba import autojit

def matConv(vec,mat):

    outVec= zeros(128)
    
    
    for i in range(128):
        for j in range(140):
            outVec[i]+= vec[j]*mat[j,i]


    return outVec 

autoDot = autojit(matConv)


class model2counts(object):



    def __init__(self,drm,model,intType='quad'):

        self.model = model
        if intType =='quad':
            self.q=quad
        else:
            self.q=quadrature
        
        self.drm = drm

        self._CreatePhotonVector()
        self._CreateCounts()


    def _CreatePhotonVector(self):

        
        phtVec = map(lambda x: self.q(self.model, x[0],x[1])[0], self.drm.phtBins )
        self.phtVec = array(phtVec)

        


    def _CreateCounts(self):

        self.counts = array(dot(self.phtVec,self.drm.drm))[0]#/(self.drm.binWidths)
        #counts = autoDot(self.phtVec, self.drm.drm)
        #self.counts = counts
        
        

            
