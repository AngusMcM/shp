import numpy as np
import math
from astroquery.jplhorizons import Horizons

# function reads catalogue and returns a list where each element is a plate
# replaces spaces with zeros to avoid errors
def getData():
	cat = open('cat.txt', 'r')
	cat = cat.readlines()
	for i in range(0,len(cat)):
		cat[i] = cat[i].replace(" ","0")
	return cat

# converts the catalogue from string list to data array
def mkArray(cat, lng):
	data = np.zeros((len(cat), 7))
	for i in range(0, len(cat)):
		data[i][0] = cat[i][2:7]
		data[i][1] = getJD(cat[i][30:40],lng)
		data[i][2] , data[i][3] = RaDecRad(cat[i][20:30])
		data[i][4] = cat[i][52:56]
	return data

# LEGACY gets the epoch of each plate number
def getEpochs(cat, lng):
	epochs = np.zeros(len(cat), 2)
	for i in range(0,len(cat)):
		JD = getJD(cat[i][30:40],lng)
		epochs[i][0] = cat[i][2:7]
		epochs[i][1] = JD
	return epochs

# determines which plates are to be discarded
def discards():
	return

# calculates JDate (currently only sensible years)
def getJDate(y,m,d):
	if m <=2:
		y_ =y-1
		m_ = m+12
	else:
		y_ =y
		m_ =m
	A = math.trunc(y_/100.)
	B = 2 - A + math.trunc(A/4.)
	C = math.trunc(365.25*y_)
	D = math.trunc(30.6001*(m_+1))
	JDate = B+C+D+d+1720994.5
	return JDate

#converts JD at 0h and GST to UT
def GSTUT(JD, GST):
	S =JD-2451545.0
	T = S/36525.0
	T0 = 6.697374558+(2400.051336*T)+(0.000025862*T**2)
	if T0<0:
		while T0<0: T0=T0+24
	elif T0>=24:
		while T0>24: T0=T0-24
	else: T0=T0
	A = GST-T0
	if A<0:
		while A<0: A=A+24
	elif A>=24:
		while A>=24: A=A-24
	else: A=A
	UT = A*0.9972675663
	return UT

# takes date and time string and longitude in degrees and gives JD
def getJD(Raw, lng):
	if int(Raw[0:2])>=50: year = 1900+float(Raw[0:2])
	else: year = 2000+float(Raw[0:2])
	month = float(Raw[2:4])
	day = float(Raw[4:6])
	hour = Raw[6:8]
	hour = float(hour.replace(" ","0"))
	minute = float(Raw[8:10])
	# converting LST to fractional hours
	LST = hour+(minute/60.)
	lng = lng/15.
	GST = LST-lng
	#GST in hours
	if GST > 24: GST=GST-24
	elif GST <= 0: GST=GST+24
	else: GST=GST
	#finding the Julian Date
	JDate = getJDate(year,month,day)
	if LST<lng: JDate=JDate-1
	else: JDate=JDate
	UT = GSTUT(JDate, GST)
	JD = JDate + UT/24.
	return JD

# takes radec as it appears in the cat and returns them seperately and in radians
def RaDecRad(radec):
	dmin = int(radec[8:10])/60.
	dec = int(radec[6:8])+dmin
	dec = math.pi*dec/180.
	if radec[5]=='+':
		Dec=dec
	else: Dec = -dec
	ramin = int(radec[2:5])/600.
	ra = int(radec[0:2])+ ramin
	Ra = ra*15*math.pi/180.
	return Ra, Dec

# converts from spherical to tangent plane coordinates wrt a reference point z
# returns a value for xi eta in radians
"""
reference provided code
"""
def ds2tp(raz,decz,ra,dec):
	radif = ra-raz
	denom = math.sin(dec)*math.sin(decz)+math.cos(dec)*math.cos(decz)*math.cos(radif)
	xi = math.cos(dec)*math.sin(radif)/denom
	eta = (math.sin(dec)*math.cos(decz)-math.cos(dec)*math.sin(decz)*math.cos(radif)/denom)
	return xi , eta

# checks if an object will appear on a plate
def match(xi,eta,thresh):
	if abs(xi)<=thresh and abs(eta)<=thresh:
		return 1
	else:
		return 0

# returns x-y position of object on plate in actual size
def onPlate(xi,eta):
	x = (xiz-xi+xithresh)*xsize/(xithresh*2)
	y = (etaz-eta+etathresh)*ysize/(etathresh*2)
	return x , y

# changes the JD to be half way through the exposure time
def midExposure(cat):
	for i in range(0, len(cat)):
		extime = (cat[i][4]/10.)*60
		extimedays = extime/86400.
		print extimedays
		cat[i][1] += extimedays/2.
	return cat

# submits horizons queries and returns radec of object for every epoch
# places the values inside the data array
def horizons(location, target, cat):
	Rarray = np.empty(1)
	Decarray = np.empty(1)
	for i in range(0,len(cat),350):
		print 'loop ' + str(i)
		epoch = cat[i:i+350]
		print len(epoch)
		epochs = []
		for j in range(0,len(epoch)):
			epochs.append(epoch[j][1])
		print epochs
		epochs = np.array(epochs)
		obj = Horizons(id=target, location=location, id_type='majorbody', epochs=epochs)
		eph = obj.ephemerides()
		newRA = np.array(eph['RA'])
		print eph['targetname']
		newDec = np.array(eph['DEC'])
		Rarray = np.concatenate((Rarray, newRA))
		Decarray = np.concatenate((Decarray, newDec))
	for i in range(0, len(cat)):
		cat[i][5] = Rarray[i+1]
		cat[i][6] = Decarray[i+1]
	return cat

# gives a list of plate numbers with potential hits on them
def hitList(cat, thresh):
	hitList = []
	hitposx = []
	hitposy = []
	for i in range(0, len(cat)):
		raz = cat[i][5]
		decz = cat[i][6]
		ra = cat[i][2]
		ra = ra*math.pi/180.
		dec = cat[i][3]
		dec = dec*math.pi/180.
		xi , eta = ds2tp(raz, decz, ra, dec)
		hit = match(xi, eta, thresh)
		if hit == 1:
			hitList.append(cat[i][0])
			print cat[i][0]
			hitposx.append(xi)
			hitposy.append(eta)
	hitpos = np.column_stack((hitList,hitposx,hitposy))
	return hitpos

def main():
	lng = 149.0661
	thresh = 6*math.pi/180.
	data = getData()
	cat = mkArray(data,lng)
	print 'len(cat) = ' + str(len(cat))
	print 'catalogue'
	print cat[1:5]
#	cat = midExposure(cat)
	cat = horizons(260,499,cat)
	np.savetxt('check.txt',cat)
	hitpos = hitList(cat, thresh)
	print hitpos
	"""
	#cat = discards(cat)
	raz = 0
	decz = 0
	ra = 0
	dec = 0.15
	xi, eta = ds2tp(raz,decz,ra,dec)
	print 'xi = ' +str(xi)
	print 'eta = ' +str(eta)
	thresh = 6*math.pi/180.
	mat = match(xi,eta,thresh)
	print mat
	"""

main()