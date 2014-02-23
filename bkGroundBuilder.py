from photonBuilder import photonBuilder


def PowerLaw(x, A, index):

    return A*x**(index)

class bkGroundBuilder(photonBuilder):




    def SetSpectrum(self,A,index):
        
        
        self.paramEvo = [[A,index]]*len(self.timeRange)
        self.SetModel(PowerLaw)
        
        


    


        
