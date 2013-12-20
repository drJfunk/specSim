import astropy.io.fits as pf
from numpy import load as npload



class tteBuilder(object):

    def __init__(self,fileName):
        
        self.fileName =fileName
        self.bgo = False
        self.nai=True


    def _CreatePrimaryHeader(self):


        primaryHeader = pf.Header()

        primaryHeader["CREATOR"] = "specSim"
        primaryHeader["FILETYPE"] = 'GBM PHOTON LIST'
        primaryHeader["FILE-VER"] = '1.0.0'
        primaryHeader["TELESCOP"] = 'GLAST   '
        primaryHeader["INSTRUME"] = 'GBM     '
        primaryHeader["DETNAM"] = 'NAI_03  '
        primaryHeader["OBSERVER"] = 'Meegan  '
        primaryHeader["ORIGIN"] = 'GIOC    '
        primaryHeader["DATE"] = '2009-05-19T18:49:32'
        primaryHeader["DATE-OBS"] = '2008-09-16T00:12:20'
        primaryHeader["DATE-END"] = '2008-09-16T00:17:47'
        primaryHeader["TIMESYS"] = 'TT      '
        primaryHeader["TIMEUNIT"] = 's       '
        primaryHeader["MJDREFI"] = 51910
        primaryHeader["MJDREFF"] = 7.428703703703703E-4
        primaryHeader["TSTART"] = 243216740.670344
        primaryHeader["TSTOP"] = 243217067.272056
        primaryHeader["FILENAME"] = self.fileName
        primaryHeader["DATATYPE"] = 'TTE     '
        primaryHeader["TRIGTIME"] = 243216766.613542
        primaryHeader["OBJECT"] = 'GRB000000000'
        primaryHeader["RADECSYS"] = 'FK5     '
        primaryHeader["EQUINOX"] = 2000.0
        primaryHeader["RA_OBJ"] = 30.0000
        primaryHeader["DEC_OBJ"] = -15.000
        primaryHeader["ERR_RAD"] = 3.000
        primaryHeader["INFILE01"] = 'glg_lutcs_nai_080915695_v00.fit'
        primaryHeader["INFILE02"] = 'GLAST_2008260_074300_VC09_GBTTE.0.00'
        primaryHeader["CHECKSUM"] = 'PZaEQWTCPWZCPWZC'
        primaryHeader["DATASUM"] = '         0'

        self.prihdu = pf.PrimaryHDU(header=primaryHeader)


    def _MakeEboundsEXT(self):

        if self.bgo:
            ebounds = npload('ebounds_bgo.npy')
        if self.nai:
            ebounds = npload('ebounds_nai.npy')


        channel = pf.Column(name='CHANNEL ', format='1I      ', unit='none', array=ebounds[:,0])
        emin =  pf.Column(name='E_MIN   ', format='1E      ', unit='keV     ', array=ebounds[:,1])
        emax =  pf.Column(name='E_MAX   ', format='1E      ', unit='keV     ', array=ebounds[:,2])


        cols = pf.ColDefs([channel,emin,emax])

        tbhdu=pf.new_table(cols,tbtype='BinTableHDU')
        
        header=tbhdu.header

        header["TELESCOP"] = 'GLAST   '
        header["INSTRUME"] = 'GBM     '
        header["DETNAM"] = 'NAI_03  '
        header["OBSERVER"] = 'Meegan  '
        header["ORIGIN"] = 'GIOC    '
        header["DATE"] = '2009-05-19T18:49:32'
        header["DATE-OBS"] = '2008-09-16T00:12:20'
        header["DATE-END"] = '2008-09-16T00:17:47'
        header["TIMESYS"] = 'TT      '
        header["TIMEUNIT"] = 's       '
        header["MJDREFI"] = 51910
        header["MJDREFF"] = 7.428703703703703E-4
        header["TSTART"] = 243216740.670344
        header["TSTOP"] = 243217067.272056
        header["FILENAME"] = self.fileName
        header["DATATYPE"] = 'TTE     '
        header["TRIGTIME"] = 243216766.613542
        header["OBJECT"] = 'GRB000000000'
        header["RADECSYS"] = 'FK5     '
        header["EQUINOX"] = 2000.0
        header["RA_OBJ"] = 30.0000
        header["DEC_OBJ"] = -15.000
        header["ERR_RAD"] = 3.000
        header["INFILE01"] = 'glg_lutcs_nai_080915695_v00.fit'
        header["INFILE02"] = 'GLAST_2008260_074300_VC09_GBTTE.0.00'
        header["CHECKSUM"] = 'PZaEQWTCPWZCPWZC'
        header["DATASUM"] = '         0'




        self.eboundshdu = tbhdu


    def _MakeEventsEXT(self,events,channels):

        if len(events) != len(channels):
            print "channels and events are not the same length!!!"
            return

       


        pha  = pf.Column(name= 'PHA     ', format='1I      ', unit='none', array=channels)
        time = pf.Column(name= 'TIME    ', format='1D      ', unit='s       ', array=events)
      


        cols = pf.ColDefs([time,pha])

        tbhdu=pf.new_table(cols,tbtype='BinTableHDU')
        
        header=tbhdu.header
        header["TZERO1"]=243216766.613542
        header["EXTNAME"] = 'EVENTS  '
        header["TELESCOP"] = 'GLAST   '
        header["INSTRUME"] = 'GBM     '
        header["DETNAM"] = 'NAI_03  '
        header["OBSERVER"] = 'Meegan  '
        header["ORIGIN"] = 'GIOC    '
        header["DATE"] = '2009-05-19T18:49:32'
        header["DATE-OBS"] = '2008-09-16T00:12:20'
        header["DATE-END"] = '2008-09-16T00:17:47'
        header["TIMESYS"] = 'TT      '
        header["TIMEUNIT"] = 's       '
        header["MJDREFI"] = 51910
        header["MJDREFF"] = 7.428703703703703E-4
        header["TSTART"] = 243216740.670344
        header["TSTOP"] = 243217067.272056
        header["FILENAME"] = self.fileName
        header["DATATYPE"] = 'TTE     '
        header["TRIGTIME"] = 243216766.613542
        header["OBJECT"] = 'GRB000000000'
        header["RADECSYS"] = 'FK5     '
        header["EQUINOX"] = 2000.0
        header["RA_OBJ"] = 30.0000
        header["DEC_OBJ"] = -15.000
        header["ERR_RAD"] = 3.000
        header["RESPFILE"]  = 'none    '
        header["EVT_DEAD"]  = 2.6000E-06
        header["EVTDEDHI"]  = 1.04170E-05
        header["DETCHANS"]  = 128
        header["HDUCLASS"]  = 'OGIP    '
        header["HDUCLAS1"] = 'EVENTS  '
        header["CHECKSUM"] = 'PZaEQWTCPWZCPWZC'
        header["DATASUM"] = '         0'




        self.eventshdu = tbhdu




    def _MakeHDUList(self):

        hduList = pf.HDUList([self.prihdu, self.eboundshdu,self.eventshdu])
        hduList.writeto(self.fileName)
