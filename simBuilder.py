from rsp import rsp


class simBuilder(object):


    def __inti__(self):


        pass


    def _ReadRSP(self):
        '''
        This method imports the RSPs needed for the simulation
        All fourteen detectors need to be present so that the simulation
        can read calculate the effective area properly


        The effective area will be calculated and then the fractional 
        area will be used to augment the rate for each spectrum


        '''
        pass


    def _GenerateBackgrounds(self):
        '''
        For each detector a background will be generated from the master
        background specified




        '''

        pass



    def _GenerateSignals(self):
        '''
        For each detector generate the pulse from the master parameters
        '''


        pass


    def _CreateTTE(self):

        pass
