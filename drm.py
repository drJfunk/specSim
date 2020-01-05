#!/usr/bin/env python

import math
import random

import numpy as np
import astropy.io.fits as pf

class RSP:
	def __init__(self, drm, photonE, chanE, detector, angle = False):
		self.drm = drm
		self.photonE = photonE
		self.chanE = chanE
		self.detector = detector
		self.beta = angle
		if self.beta:
			self.calculateProb()
		else:
			self.prob = False		
		# get centre of bins for plotting purposes
		self.photonCent = (photonE[1]-photonE[0])/2. + photonE[0]
		self.chanCent = (chanE[1]-chanE[0])/2. + chanE[0]
	def calculateProb(self):
		'''
		Convert the drm into a probablility. To do so, we need the angle 
		relative to the normal of the detector. For the NaI the 
		angle should be relative to the face of the cylinder. For the BGO
		it is given relative to the side of the detector.
		'''

		area = lambda r, h, b: math.fabs(math.pi *(r**2) * math.cos(math.radians(b))) + math.fabs(2 * r * h * math.sin(math.radians(b)))

		if self.detector[0].lower() == "n":
			# NaI's are cylinders of diameter 12.7 cm and thickness 1.27 cm
			geoArea = area(12.7/2., 1.27, self.beta)
		else:
			# BGO's are cylinders of diameter 12.7 cm and thickness 12.7 cm.
			# For the BGO, we have to add an offset of 90 deg to the angle. This
			# is because for the NaI the angles are calculate relative to the face
			# of the cylinder whereas for the BGO they are calculated relative to the 
			# cylindrical face.
			geoArea = area(12.7/2., 12.7, self.beta + 90.)
		self.geoArea = geoArea
		self.prob = self.drm/geoArea
		# print "Det: %s, Beta: %.1f, Area: %.1f cm^2" %(self.detector, self.geoArea, self.beta)	
	def findNearestPhotonBin(self, e):
		''' 
		Take an energy (e) and find return the 
		corresponding column from the DRM which 
		has been normalised by the geometric area
		'''
		idx=(np.abs(self.photonE[0]-e)).argmin()
		return self.prob[idx,:]

def decompressDRM(fo):
	'''
	Decompress a GBM DRM. Takes in a python pyfits object as an argument,
	e.g. fo = pyfits.open("someResponseFile.rsp")
	'''

	nPhotonBins = fo[2].header["NAXIS2"]

	nChan = fo[2].data["N_CHAN"]
	fChan = fo[2].data["F_CHAN"]

	nGrp = fo[2].data["N_GRP"]
	matrix = fo[2].data["MATRIX"]

	eLow = fo[2].data["ENERG_LO"]
	eHi  = fo[2].data["ENERG_HI"]

	eChanMin = fo[1].data["E_MIN"]
	eChanMax = fo[1].data["E_MAX"]
	nChanBins = eChanMax.size

	drm = np.zeros((nPhotonBins, nChanBins ))
	for i, prob in enumerate(matrix):
		# print nGrp[i], matrix[i].size
		if nChan[i] != 128:
			row = np.zeros(fChan[i]-1)
			row = np.concatenate((row, matrix[i]))
		else:
			row = matrix[i]
		drm[i,:] = row

	return drm

def readRSP(fname, angle):
	data = pf.open(fname)
	drm = decompressDRM(data)
	det = data[0].header["DETNAM"]
	photonE = (data[2].data["ENERG_LO"], data[2].data["ENERG_HI"])
	channelE = (data[1].data["E_MIN"], data[1].data["E_MAX"])		
	rsp = RSP(drm, photonE, channelE, det, angle = angle)
	return rsp

def multiplyResponse(t,e, rsp):
	'''
	Loop over input energies and multiply photons
	by the response and obtain a probablility
	distribution in observed energy for that photon.
	'''
	chans = []
	ts = []
	c = 0
	
	tCount, eCount = [], []
	
	for i,j in zip(t,e):
		prob = rsp.findNearestPhotonBin(j)
		pTot = prob.sum()
		
		if pTot == 0:
			continue

		normProb = prob/pTot	

		while pTot > 0:
			r = random.random()
			if r < pTot:
				pCum = np.cumsum(normProb)
				idx=(np.abs(pCum-r)).argmin() 
				tCount.append(i)	
				eCount.append(idx)
			pTot -= 1
		
	return tCount, eCount


def getDetAngles(met, ra, dec, posfile):
	'''
	Take in a MET, source RA & Dec, a relevant
	posfile and obtain the source-detector angles.
	Returns a dictionary of angles indexed by 
	detector name
	'''
	fo = pf.open(posfile)

	dtorad = 180./math.acos(-1.)

	sc_time = fo[1].data.SCLK_UTC
	ind = [(sc_time >= met) & (sc_time <= met +1.)]
	sc_time = fo[1].data.SCLK_UTC[ind]
	
	sc_quat   = np.zeros((1,4), float)
	sc_pos    = np.zeros((1,3), float)
			
	sc_quat[:,0] = fo[1].data.QSJ_1[ind]
	sc_quat[:,1] = fo[1].data.QSJ_2[ind]
	sc_quat[:,2] = fo[1].data.QSJ_3[ind]
	sc_quat[:,3] = fo[1].data.QSJ_4[ind]
	sc_pos[:,0] =  fo[1].data.POS_X[ind]
	sc_pos[:,1] =  fo[1].data.POS_Y[ind]
	sc_pos[:,2] =  fo[1].data.POS_Z[ind]
	
	dZ, dGeo, dDet = __calcDetAngles(sc_time, sc_pos, sc_quat, ra, dec)

	dets = ["n0","n1","n2","n3","n4","n5", 
			"n6","n7","n8","n9","na","nb",
			"b0","b1"]

	angles = {}
	for i, det in zip(dDet[0], dets):
		angles.update({det:i})
	return angles

def parDtFilter(t, pha = np.empty(0)):
	'''
	Read in an array of times and apply GBM deadtime, includes paralyzable & 
	non-paralyzable effects. Can optionally read in a array of channel values 
	(pha). For each count, if the channel corresponds to the overflow (127), 
	then the deadtime is taken as 10.6 us. 
	Returns the filtered time array and the filtered pha array.
	
	The PHA accepts a count if a local max is found and the readings from the 
	next 4 samples (sample ~ 0.1 us) are lower than the max. It then applies a
	deadtime of 2.6 us from the peak. If a photon arrives during the 4 sample 
	window then the first count will be lost and the deadtime will be applied 
	from the arrival time of the second count.
	
	For a given count, there are three possibilities:
	1) The separation time of the count relative to the previous count is
		greater than 2.6 us. In this case the count is kept.
	2) The separation time of the count relative to the previous count is 
		less than 2.6 us but greater than 0.5 us. In this case the count is
		discarded.
	3) The separation time of the count relative to the previous count is
		less than 0.5 us. In this case the previous count is removed and this
		count is kept.
	
	'''


	if not pha.size:
		pha = np.ones(t.size)
	deadTime = 2.6e-6
	ofDeadTime = 10.6e-6    
	parDeadTime = 0.5e-6
	
	filteredPha = [pha[0]]
	filteredT = [t[0]]  
	  
	if pha[0] == 127:
		tEnd = t[0] + ofDeadTime
	else:
		tEnd = t[0] + deadTime
	tEndPar = t[0] + parDeadTime
		
	for i, chan in zip(t[1:], pha[1:]):
		# Case (1): no deadtime
		#print "%f %f %f" %(i, tEnd, tEndPar)
		if i > tEnd:
			# Add current event to filtered data
			filteredT.extend([i])
			filteredPha.extend([chan])
			# Get new end time
			if chan == 127:
				tEnd = i + ofDeadTime
			else:
				tEnd = i + deadTime
			tEndPar = i + parDeadTime                
		# Case (3): par deadtime
		elif i < tEndPar:
			#print "paralyzable Deadtime"
			# Remove most recent entry from filtered data
			filteredT.pop()
			filteredPha.pop()
			# Add current event to filtered data
			filteredT.extend([i])
			filteredPha.extend([chan])
			# Get new end time
			if chan == 127:
				tEnd = i + ofDeadTime
			else:
				tEnd = i + deadTime
			tEndPar = i + parDeadTime
	
	return np.array(filteredT), np.array(pha)

def main():
	angle = 67.
	rsp = readRSP("glg_cspec_b1_bn081001392_v200.rsp", angle)
	
	t = np.arange(10)
	e = np.random.randint(100,2000,size=10)

	tCount, eChannel = multiplyResponse(t,e,rsp)

	# deadtime filter
	tCountFilt, eChannelFilt = parDtFilter(tCount, np.array(eChannel))


if __name__ == "__main__":
	 main()

