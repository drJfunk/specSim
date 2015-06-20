import astropy.io.fits as fits
import os, sys
import shutil
import numpy as np



def MakeDiagRSP(rspFile):
    "Create a perfect RSP matrix"

    # First we duplicate the file
    origFile = rspFile
    idx = origFile.rfind(".rsp")

    perfectFile = origFile[:idx]+"_PERFECT"+origFile[idx:]

    origRSP = fits.open(origFile)
    


    primHeader = origRSP["PRIMARY"]
    ebounds    = origRSP["EBOUNDS"]
    mat        = origRSP["SPECRESP MATRIX"]










    eLO = origRSP["EBOUNDS"].data["E_MIN"]
    eHI = origRSP["EBOUNDS"].data["E_MAX"]

    nChannels = len(eLO)

    # Ideally, one would create the columns this way,
    # but pyfits is fucking retarted at the moment and
    # you have to do this in a special way because of the matrix.
    # One day of my life gone away

    
    # loSide = fits.Column(name='ENERG_LO', format='1E', unit='keV', array=eLO)
    # hiSide = fits.Column(name='ENERG_HI', format='1E', unit='keV', array=eHI)
    # ngrp = fits.Column(name='N_GRP', format='1I',  array=np.ones(nChannels))
    # fchan = fits.Column(name='F_CHAN', format='PI', array=np.zeros(nChannels))
    # nchan = fits.Column(name='N_CHAN', format='PI',  array=np.zeros(nChannels)+nChannels)
    # drm = fits.Column(name='MATRIX', format='E', unit='cm**2', array=np.identity(nChannels))
   

    # cols = fits.ColDefs([drm,loSide,hiSide,ngrp,fchan,nchan])

    # First construct a record array to hold all the info
    tmp = []
    for a,b,c,d,e,f in zip(eLO,eHI,np.ones(nChannels),np.zeros(nChannels),np.zeros(nChannels)+nChannels,np.identity(nChannels)):
        tmp.append([a,b,c,d,e,f.tolist()])

    data=np.rec.array(tmp,dtype=[('ENERG_LO', '>f4'), ('ENERG_HI', '>f4'), ('N_GRP', '>f4'), ('F_CHAN', '>f4'),('N_CHAN', '>f4'), ('MATRIX', '>f4', (nChannels, nChannels))])

    tbhdu=fits.BinTableHDU.from_columns(data)    
    header=tbhdu.header
    #header["TZERO1"]=self.tz
    header["EXTNAME"] = 'SPECRESP MATRIX'
    header["TELESCOP"] = 'GLAST   '
    header["INSTRUME"] = 'GBM     '
    header["DETNAM"] = mat.header["DETNAM"]
    header["OBSERVER"] = 'Meegan  '
    header["ORIGIN"] = 'GIOC    '
    header["DATE"] = mat.header["DATE"]
    header["DATE-OBS"] = mat.header["DATE-OBS"]
    header["DATE-END"] =mat.header["DATE-END"]# '2008-09-16T00:17:47'
    header["TIMESYS"] = 'TT      '
    header["TIMEUNIT"] = 's       '
    header["MJDREFI"] = 51910
    header["MJDREFF"] = 7.428703703703703E-4

    header["TSTART"] = mat.header["TSTART"]#self.tStart
    header["TSTOP"] = mat.header["TSTOP"]#self.tStop
    #header["DATATYPE"] = 'TTE     '
    header["TRIGTIME"] = mat.header["TRIGTIME"]#self.trigT
    header["DETCHANS"]  = nChannels
    header["NUMEBINS"]  = nChannels
    header["MAT_TYPE"] = mat.header["MAT_TYPE"]
    header["RSP_NUM"] = mat.header["RSP_NUM"]
    header["OBJECT"] = mat.header["OBJECT"]#'GRB000000000'
    header["RADECSYS"] = 'FK5     '
    header["EQUINOX"] = mat.header["EQUINOX"]
    header["RA_OBJ"] = mat.header["RA_OBJ"]
    header["DEC_OBJ"] = mat.header["DEC_OBJ"]
    header["HDUCLASS"]  = 'OGIP    '
    header["HDUCLAS1"] = 'RESPONSE'
    header["HDUCLAS2"] = 'RSP_MATRIX'
    header["EXTVER"] = 1
        

    hduList = fits.HDUList([primHeader, ebounds, tbhdu])
    hduList.writeto(perfectFile,clobber=True)


    origRSP.close()    


    
