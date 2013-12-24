import astropy.io.fits as pf
from numpy import load as npload



eboundDict = {'NAI_00': 'n0' ,'NAI_01': 'n1'  , 'NAI_02': 'n2'  , 'NAI_03': 'n3' , 'NAI_04': 'n4' , 'NAI_05': 'n5' , 'NAI_06': 'n6' , 'NAI_07': 'n7' , 'NAI_08': 'n8' , 'NAI_09': 'n9' , 'NAI_10': 'na'  , 'NAI_11': 'nb' , 'BGO_00': 'b0'  , 'BGO_01': 'b1'}



class tteBuilder(object):

    def __init__(self,fileName,evts,chans,det,simInfo):
        
        self.fileName =fileName
        self.det = det
        self.evts = evts
        self.chans = chans
        self.tz, self.tStart, self.tStop = simInfo
        self.trigT = self.tz
        self.ebounds = npload(eboundDict[det]+'_ebounds.npy')


        self._CreatePrimaryHeader()
        self._MakeEboundsEXT()
        self._MakeEventsEXT()
        self._MakeGTI()
        self._MakeHDUList()
        


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
        primaryHeader["TSTART"] = self.tStart
        primaryHeader["TSTOP"] = self.tStop
        primaryHeader["FILENAME"] = self.fileName
        primaryHeader["DATATYPE"] = 'TTE     '
        primaryHeader["TRIGTIME"] = self.trigT
        primaryHeader["OBJECT"] = 'GRB000000000'
        primaryHeader["RADECSYS"] = 'FK5     '
        primaryHeader["EQUINOX"] = 2000.0
        primaryHeader["RA_OBJ"] = 30.0000
        primaryHeader["DEC_OBJ"] = -15.000
        primaryHeader["ERR_RAD"] = 3.000
        primaryHeader["INFILE01"] = 'glg_lutcs_nai_080915695_v00.fit'
        primaryHeader["INFILE02"] = 'GLAST_2008260_074300_VC09_GBTTE.0.00'
      #  primaryHeader["CHECKSUM"] = 'PZaEQWTCPWZCPWZC'
      #  primaryHeader["DATASUM"] = '         0'

        self.prihdu = pf.PrimaryHDU(header=primaryHeader)


    def _MakeEboundsEXT(self):

        

        channel = pf.Column(name='CHANNEL ', format='1I', unit='none', array=self.ebounds[:,0])
        emin =  pf.Column(name='E_MIN   ', format='1E', unit='keV     ', array=self.ebounds[:,1])
        emax =  pf.Column(name='E_MAX   ', format='1E', unit='keV     ', array=self.ebounds[:,2])


        cols = pf.ColDefs([channel,emin,emax])

        tbhdu=pf.new_table(cols,tbtype='BinTableHDU')
        
        header=tbhdu.header
        header["EXTNAME"] = 'EBOUNDS '
        header["TELESCOP"] = 'GLAST   '
        header["INSTRUME"] = 'GBM     '
        header["DETNAM"] = self.det
        header["OBSERVER"] = 'Meegan  '
        header["ORIGIN"] = 'GIOC    '
        header["DATE"] = '2009-05-19T18:49:32'
        header["DATE-OBS"] = '2008-09-16T00:12:20'
        header["DATE-END"] = '2008-09-16T00:17:47'
        header["TIMESYS"] = 'TT      '
        header["TIMEUNIT"] = 's       '
        header["MJDREFI"] = 51910
        header["MJDREFF"] = 7.428703703703703E-4
        header["TSTART"] = self.tStart
        header["TSTOP"] = self.tStop
       
        header["TRIGTIME"] = self.trigT
        header["OBJECT"] = 'GRB000000000'
        header["RADECSYS"] = 'FK5     '
        header["EQUINOX"] = 2000.0
        header["RA_OBJ"] = 30.0000
        header["DEC_OBJ"] = -15.000
        header["ERR_RAD"] = 3.000
        header["CHANTYPE"] = 'PHA     '
        header["DETCHANS"] = 128
        header["EXTVER"] = 1
        header["CH2E_VER"] ='SPLINE 2.0'
        header["GAIN_COR"] = 1.0
      #  header["CHECKSUM"] = 'PZaEQWTCPWZCPWZC'
      #  header["DATASUM"] = '         0'




        self.eboundshdu = tbhdu


    def _MakeEventsEXT(self):

        if len(self.evts) != len(self.chans):
            print "channels and events are not the same length!!!"
            return

       
        

        pha  = pf.Column(name= 'PHA     ', format='1I', unit='none', array=self.chans)
        time = pf.Column(name= 'TIME    ', format='1D', unit='s       ', array=self.evts-self.tz)
      
        

        cols = pf.ColDefs([time,pha])

        tbhdu=pf.new_table(cols,tbtype='BinTableHDU')
        
        header=tbhdu.header
        header["TZERO1"]=self.tz
        header["EXTNAME"] = 'EVENTS  '
        header["TELESCOP"] = 'GLAST   '
        header["INSTRUME"] = 'GBM     '
        header["DETNAM"] = self.det
        header["OBSERVER"] = 'Meegan  '
        header["ORIGIN"] = 'GIOC    '
        header["DATE"] = '2009-05-19T18:49:32'
        header["DATE-OBS"] = '2008-09-16T00:12:20'
        header["DATE-END"] = '2008-09-16T00:17:47'
        header["TIMESYS"] = 'TT      '
        header["TIMEUNIT"] = 's       '
        header["MJDREFI"] = 51910
        header["MJDREFF"] = 7.428703703703703E-4
        
        header["TSTART"] = self.tStart
        header["TSTOP"] = self.tStop
        header["DATATYPE"] = 'TTE     '
        header["TRIGTIME"] = self.trigT
        
        
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
        header["EXTVER"] = 1
        header["FIFO_END"] = self.tStart
        header["FIFO_END"] = self.tStart
        header["CHECKSUM"] = 'PZaEQWTCPWZCPWZC'
        header["DATASUM"] = '         0'




        self.eventshdu = tbhdu



    def _MakeGTI(self):
        

        

        tstart = pf.Column(name='START   ', format='1D', unit='s', array=[self.tStart])
        tstop =  pf.Column(name='STOP    ', format='1D', unit='s     ', array=[self.tStop])
       


        cols = pf.ColDefs([tstart,tstop])

        tbhdu=pf.new_table(cols,tbtype='BinTableHDU')
        
        header=tbhdu.header
        header["EXTNAME"] = 'GTI     ' 
        header["TELESCOP"] = 'GLAST   '
        header["INSTRUME"] = 'GBM     '
        header["DETNAM"] = self.det
        header["OBSERVER"] = 'Meegan  '
        header["ORIGIN"] = 'GIOC    '
        header["DATE"] = '2009-05-19T18:49:32'
        header["DATE-OBS"] = '2008-09-16T00:12:20'
        header["DATE-END"] = '2008-09-16T00:17:47'
        header["TIMESYS"] = 'TT      '
        header["TIMEUNIT"] = 's       '
        header["MJDREFI"] = 51910
        header["MJDREFF"] = 7.428703703703703E-4
        header["TSTART"] = self.tStart
        header["TSTOP"] = self.tStop
        
        header["TRIGTIME"] = self.trigT
       
        header["TZERO1"]=self.tz
        header["TZERO2"]=self.tz
       
        header["OBJECT"] = 'GRB000000000'
        header["RADECSYS"] = 'FK5     '
        header["EQUINOX"] = 2000.0
        header["RA_OBJ"] = 30.0000
        header["DEC_OBJ"] = -15.000
        header["ERR_RAD"] = 3.000
        header["CHANTYPE"] = 'PHA     '
        header["DETCHANS"] = 128
        header["EXTVER"] = 1
        header["HDUCLASS"]  = 'OGIP    '
        header["HDUCLAS1"] = 'GTI     '
        header["HDUVERS"] = '1.2.0'
        header["EXTVER"] = 1
       
       # header["CHECKSUM"] = 'PZaEQWTCPWZCPWZC'
       # header["DATASUM"] = '         0'




        self.gtihdu = tbhdu
        


    def _MakeHDUList(self):

        self.hduList = pf.HDUList([self.prihdu, self.eboundshdu,self.eventshdu, self.gtihdu])
        self.hduList.writeto(self.fileName,clobber=True)
        
