from scipy.integrate import quad, quadrature
from numpy import dot, array

class model2counts(object):



    def __init__(self,drm,model):

        self.model = model

        self.drm = drm

        self._CreatePhotonVector()
        self._CreateCounts()


    def _CreatePhotonVector(self):

        
        phtVec = map(lambda x: quad(self.model, x[0],x[1])[0], self.drm.phtBins )
        self.phtVec = array(phtVec)
        


    def _CreateCounts(self):

        self.counts = array(dot(self.phtVec,self.drm.drm))[0]
        
        

            
