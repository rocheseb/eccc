#!/usr/bin/env python2.7
 # -*- coding: ascii -*-

####################
# Code Descirption #
####################

'''
Reads TCCON and GEM-MACH-GHG data.
Homogenize TCCON a priori and GEM-MACH-GHG profiles.
Computes TCCON a priori columns.
Computes GEM-MACH-GHG columns.
Use TCCON averaging kernel to smooth the data.
Save the data to be read by EC_SMOOTH_PLOTS.py for statistics and plots.

This version of the code uses the GGG 2014 data output/format for the TCCON measurements.
'''

#############
# Libraries #
#############

# manipulate paths
import os.path

# some math functions
from math import log,exp

# plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# interaction with netcdf files
import netCDF4

# time handling
import calendar
import time
#import solar

# special lists with special functions
import numpy as np

# easily read/write python objects in a file
import pickle

# commandline interactions
import sys

#############################################################################################################
#############################################################################################################
####################################### MODIFY HERE #########################################################
#############################################################################################################
#############################################################################################################
TCCON_path ='/data/02/sroche/EC-CAS & TCCON data comparison/main/Data/TCCON7'                        ######## path to the TCCON netcdf files
model_path ='/data/02/sroche/EC-CAS & TCCON data comparison/main/Data/model/mn117b_at_TCCON.nc'      ######## path to the model file
save_path='/data/02/sroche/EC-CAS & TCCON data comparison/main/SMOOTH100/mn117b'                     ######## path to save the output files
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################

#########
# SETUP #
#########

#To keep track of time while the code runs
codetsart=time.time()

#Used to divide lists in even chunks when possible, if not the first/last chunk will have 1-2 elements more or less than the rest.
#The profile have 71 LEVEL values, those are interpolated to 351 LEVEL values.
#chunkIt will divide the profile in chunks with 4-5 level values each.
#Those level values within each layer are averaged to obtain the layer mean profile, 
def chunkIt(l,n):
 	n = max(1, n)
 	return [l[i:i + n] for i in range(0, len(l), n)]

#a fancy loadbar to be displayed in the prompt while executing a time consuming loop.
def progress(i,tot,bar_length=20):
	percent=float(i)/tot
	hashes='#' * int(round(percent*bar_length))
	spaces=' ' * (bar_length - len(hashes))
	sys.stdout.write("\rPercent:[{0}] {1}%".format(hashes + spaces, int(round(percent * 100)))+"    "+str(i)+"/"+str(tot))
	sys.stdout.flush()

## Some constants
Av=6.0221415E+23
R=8.314
#model gravity constant m/s
g=9.80616
#standard mean ocean water
SMOW=3.10693E-4
mass_H2O=18.01534E-03
mass_dry_air=28.9644E-03

files=os.listdir(TCCON_path)

##############
# READ FILES #
##############

#Read the Lamont kernels, used for sites without averaging kernels
f=netCDF4.Dataset('/data/02/sroche/EC-CAS & TCCON data comparison/main/Data/Avk/Lam.public.nc','r')
lam_ak_co2=f.variables['ak_co2'][:]
lam_ak_ch4=f.variables['ak_ch4'][:]
lam_akP=f.variables['ak_P_hPa'][:]
lam_akZ=f.variables['ak_zenith'][:]
f.close

## Arrays the program will output with apriori/EC/EC-smooth/TCCON xCO2 column amounts for each site
save_data_co2_path=os.path.join(save_path,'DMF_co2.txt')
save_data_co2=[[] for i in range(len(files)*4)]

save_data_ch4_path=os.path.join(save_path,'DMF_ch4.txt')
save_data_ch4=[[] for i in range(len(files)*4)]

save_errors_path=os.path.join(save_path,'errors.txt')
save_errors=[[] for i in range(len(files)*4)]

#time associated with each column
save_time_path=os.path.join(save_path,'time.txt')
save_time=[[] for i in range(len(files)*2)]

##### READ THE MODEL ###########
f=netCDF4.Dataset(model_path,'r')

t=calendar.timegm((1980,1,1,0,0,0))
Tepoch1=[i+t for i in f.variables['time'][:]]
print 'Model file:',model_path.split('/')[-1]
print '\nModel run start:', time.strftime('%m %d %Y %H:%M:%S',time.gmtime(Tepoch1[0]))
print 'Model run end:', time.strftime('%m %d %Y %H:%M:%S',time.gmtime(Tepoch1[-1])),'\n'

#CO2 mixing ratio
mr_co2=f.variables['CO2'][:]
#CH4 mixing ratio
mr_ch4=f.variables['CH4'][:]
#pressure levels
pr=f.variables['air_pressure'][:]
#convert specific humidity from EC-CAS to H2O mixing ratio
#wr=[(1.0/((1.0/i)-1))*28.9644/18.01534 for i in f.variables['specific_humidity'][:]]
#specific humidity
wr=[i for i in f.variables['specific_humidity'][:]]

#names of sites in the EC simulation
names=[]
for station in f.variables['station_name'][:]:
	name=''
	for letter in station:
		if letter==',':
			break
		if str(letter)=='--':
			break
		name+=letter
	names.append(name)
f.close()

#names of TCCON sites
Tnames=[]
for station in files:
	Tname=''
	for letter in station:
		try:
			if float(letter)/float(letter)==True:
				break
		except ValueError:
			pass
		Tname+=letter
	Tnames.append(Tname)

print 'TCCON sites:'
for i in Tnames:
	print '\t',i

##MODIFY HERE
#t0 is the January 1st 2009 at 00:00:00
#tf is the January 1st 2011 at 00:00:00
t0=calendar.timegm((2009,1,1,0,0,0,0))
tf=calendar.timegm((2011,1,1,0,0,0,0))
################################

#Loop over TCCON sites
for nameoffile in files:

	#nameoffile=files[-1]

	completeName=os.path.join(TCCON_path,nameoffile)

	#READ TCCON netCDF file.
	f=netCDF4.Dataset(completeName,'r')

	t=calendar.timegm((1970,1,1,0,0,0)) #epoch, this should be useless (t=0) on most computer, but it's there anyway just in case
	TCCONepoch=[i*24*3600+t for i in f.variables['time'][:]] #time since epoch in seconds for each measurement
	TCCONepoch=np.array(TCCONepoch)

	## Select the data that overlaps in time between EC-CAS and TCCON, if there is no overlap it will do nothing and the 
	## TCCON site loop will move to the next site
	overlap=True
	if TCCONepoch[0]>Tepoch1[-1] or TCCONepoch[-1]<Tepoch1[0]:
		overlap=False

	if ((Tnames[files.index(nameoffile)] in names)==True or Tnames[files.index(nameoffile)]=='LLauder') and overlap==True:
		Tepoch=[i for i in Tepoch1]
		nicetime=[time.strftime('%b %d %Y %H:%M:%S',time.gmtime(i)) for i in Tepoch]

		#nameoffile=files[3]

		print '\n', nameoffile

		#to keep track of time
		timestart=time.time()

		#select the model data of the corresponding station
		for name in names:
			if (name==Tnames[files.index(nameoffile)]) or (name=='Caltech' and Tnames[files.index(nameoffile)]=='Pasadena') or (name=='Lauder' and Tnames[files.index(nameoffile)]=='LLauder'):
				model_co2=[i[names.index(name)] for i in mr_co2]
				model_ch4=[i[names.index(name)] for i in mr_ch4]
				model_pressure=[i[names.index(name)] for i in pr]
				model_water=[i[names.index(name)] for i in wr]

		#######################################################

		#TCCON dry-air mole fraction column, associated error and solar zenith angle
		x_co2=f.variables['xco2_ppm'][:]
		xco2_error=f.variables['xco2_ppm_error'][:]

		x_ch4=[i*1000.0 for i in f.variables['xch4_ppm'][:]] #get xch4 in ppb instead of ppm
		xch4_error=[i*1000.0 for i in f.variables['xch4_ppm_error'][:]] # same for errors
		
		TCCON_sza=f.variables['asza_deg'][:]

		# Lauder switched instrument to 125HR in 2010, append the 125 HR 2010 data to the 2009 data (they are in different netCDF files)
		# the 120 and 125 measurements overlap but we do not mix the data, there is only 120HR data until the 125HR set starts, then there is only 125HR data
		# https://tccon-wiki.caltech.edu/@api/deki/files/1579/=ll20100202_20121231.txt
		# taken out on 17/02/2016 because the TCCON time series would be too low starting from july 2010 for an unknown reason
		# I thus separate the treament of Lauder 120HR and 125HR data
		'''		
		if Tnames[files.index(nameoffile)]=='Lauder':
			g=netCDF4.Dataset('/data/02/sroche/EC-CAS & TCCON data comparison/main/Data/125HR7/Lauder12520100202_20150101.public.nc','r')
			g_x_co2=g.variables['xco2_ppm'][:]
			g_xco2_error=g.variables['xco2_ppm_error'][:]
			g_TCCON_sza=g.variables['asza_deg'][:]
			g_TCCONepoch=[i*24*3600+t for i in g.variables['time'][:]]
			g.close()

			TCCONepoch=TCCONepoch[TCCONepoch<g_TCCONepoch[0]]
			
			x_co2=np.array(list(x_co2)+list(g_x_co2))
			xco2_error=np.array(list(xco2_error)+list(g_xco2_error))
			TCCON_sza=np.array(list(TCCON_sza)+list(g_TCCON_sza))
			TCCONepoch=np.array(list(TCCONepoch)+list(g_TCCONepoch))
			
			print "start of 125HR is", time.strftime('%b %d %Y %H:%M:%S',time.gmtime(g_TCCONepoch[0]))
			print "Index:",list(TCCONepoch).index(g_TCCONepoch[0])
		
		'''
		#to get a nice time output ! (month day year hr:min:sec)
		TCCONnicetime=[time.strftime('%b %d %Y %H:%M:%S',time.gmtime(i)) for i in TCCONepoch]

		'''
		TCCON_nicetime=[]
		TCCON_epoch=[]
		TCCON_SZA=[]
		xco2=[]
		xch4=[]
		#only select TCCON data from 2009 and 2010
		
		for a in TCCONnicetime:
			if list(a)[10]=='9' or (list(a)[9]+list(a)[10])=='10':
				TCCON_nicetime.append(a)
				TCCON_epoch.append(TCCONepoch[TCCONnicetime.index(a)])
				TCCON_SZA.append(TCCON_sza[TCCONnicetime.index(a)])
				xco2.append(x_co2[TCCONnicetime.index(a)])
		'''
		#this is MUCH faster !!! (in general list comprehension is faster than calls to dot functions like list.append())
		TCCON_epoch=TCCONepoch[(TCCONepoch>t0) & (TCCONepoch<tf)]
		T_ID=np.nonzero(np.in1d(TCCONepoch,TCCON_epoch))[0]
		xco2=[x_co2[i] for i in T_ID]
		xch4=[x_ch4[i] for i in T_ID]
		TCCON_SZA=[TCCON_sza[i] for i in T_ID]
		TCCON_nicetime=[TCCONnicetime[i] for i in T_ID]

		#TCCON a priori profiles / time / pressure / density
		prio_time=[i*24*3600+t for i in f.variables['prior_date'][:]]

		prionicetime=[time.strftime('%m %d %Y %H:%M:%S',time.gmtime(i)) for i in prio_time]

		##### MODIFY HERE
		#here i select the a priori time to have all apriori in 2009-2010
		prio_nicetime=[]
		prioID=[]
		for i in range(len(prionicetime)):
			if list(prionicetime[i])[9]=='9' or (list(prionicetime[i])[8]+list(prionicetime[i])[9])=='10': # very lousy selection of the years !
				prio_nicetime.append(prionicetime[i])
				prioID.append(prionicetime.index(prionicetime[i]))

		prioco2=f.variables['prior_co2'][:]
		prioch4=f.variables['prior_ch4'][:]
		prioh2o=f.variables['prior_h2o'][:]
		priohdo=f.variables['prior_hdo'][:]
		priodensity=f.variables['prior_Density'][:]
		priopressure=f.variables['prior_Pressure'][:]
		priogravity=f.variables['prior_gravity'][:]
		#altitude of the station, needed for the column integration to start from the station altitude instead of 0
		zmin=f.variables['zobs_km'][:]

		#Averaging kernels
		try:
			ak_co2=f.variables['ak_co2'][:]
			ak_ch4=f.variables['ak_ch4'][:]
			ak_P=f.variables['ak_P_hPa'][:]
			ak_zenith=f.variables['ak_zenith'][:]
		except KeyError:
			print "Using Lamont Averaging Kernels"
			ak_co2=lam_ak_co2
			ak_ch4=lam_ak_ch4
			ak_P=lam_akP
			ak_zenith=lam_akZ

		f.close()

		############################################
		########   TCCON A PRIORI COLUMNS   ########
		############################################ 

		print 'readnetCDF DONE',time.time()-timestart
		
		#remove the model data earlier than 1hr before the first TCCON time
		while TCCON_epoch[0]-Tepoch[0]>3600:
			if TCCON_epoch[0]-Tepoch[0]>3600:
				del Tepoch[0], nicetime[0], model_co2[0], model_ch4[0], model_pressure[0], model_water[0]

		#remove the model data later than 1hr after the last TCCON time
		while TCCON_epoch[-1]-Tepoch[-1]<(-3600):
			if TCCON_epoch[-1]-Tepoch[-1]<(-3600):
				del Tepoch[-1], nicetime[-1], model_co2[-1], model_ch4[-1], model_pressure[-1], model_water[-1]
		

		temp=[[] for i in range(t0,tf,3600)]
		temp1=[[] for i in range(t0,tf,3600)]
		#I divide the time between t0 and tf into intervals of 1 hour
		Tepoch=np.array(Tepoch) #make sure that Tepoch is a numpy.array() object
		print TCCON_nicetime[0],' to', TCCON_nicetime[-1]
		a=0
		for i in range(t0,tf,3600):
			progress(1+(i-t0)/3600,1+(tf-t0)/3600)
			'''
			# I leave this there so I won't forget how bad that loop is
			# this loop is SUPER TIME CONSUMING; it checks how many measurements and simulation there are in each hour
			for j in Tepoch:
				if int(j)>=int(i) and int(j)<int(i+3600):
					temp[a].append(Tepoch.index(j))
			for k in TCCON_epoch:
				if int(k)>=int(i) and int(k)<int(i+3600):
					temp1[a].append(TCCON_epoch.index(k))
			'''
			#EDIT 16/02/2016, the new code is several orders of magnitudes faster, this loops takes now a few seconds when it could take more than an hour for some sites before....
			times=TCCON_epoch[(TCCON_epoch>=i) & (TCCON_epoch<(i+3600))]
			times_ID=np.nonzero(np.in1d(TCCON_epoch,times))[0]
			temp1[a]=[j for j in times_ID]
			if len(temp1[a])>0:
				times=Tepoch[(Tepoch>=i) & (Tepoch<(i+3600))]
				times_ID=np.nonzero(np.in1d(Tepoch,times))[0]
				temp[a]=[j for j in times_ID]
			a=a+1

		join = lambda it: (y for x in it for y in x) #I end up not using it

		##Get the indices of hours where measurements and simulations coincide
		#indices of the model profiles used
		tp=[i for i in temp if i!=[] and temp1[temp.index(i)]!=[]]
		#tp=[i for i in temp if i!=[]]
		
		#indices of the TCCON columns used
		tp1=[i for i in temp1 if i!=[] and temp[temp1.index(i)]!=[]]
		#tp1=[i for i in temp1 if i!=[]]

		if len(tp)!=len(tp1):
			print 'Length of matches is different: tp1=',len(tp1),'; tp=',len(tp)

		#average model profiles for slices of 1 hour
		resAVG_co2=[[sum([model_co2[j][k] for j in i])/len(i) for k in range(len(model_co2[0]))] for i in tp]
		resAVG_ch4=[[sum([model_ch4[j][k] for j in i])/len(i) for k in range(len(model_ch4[0]))] for i in tp]
		pressureAVG=[[sum([model_pressure[j][k] for j in i])/len(i) for k in range(len(model_co2[0]))] for i in tp]
		waterAVG=[[sum([model_water[j][k] for j in i])/len(i) for k in range(len(model_co2[0]))] for i in tp]
		TepochAVG=[sum([Tepoch[j] for j in i])/len(i) for i in tp]

		#averaged TCCON columns for slices of 1 hour
		xco2AVG=[sum([xco2[j] for j in i])/len(i) for i in tp1]
		xch4AVG=[sum([xch4[j] for j in i])/len(i) for i in tp1]
		TCCON_szaAVG=[sum([TCCON_SZA[j] for j in i])/len(i) for i in tp1]
		xtimeAVG=[sum([TCCON_epoch[j] for j in i])/len(i) for i in tp1]

		
		#select an a priori for each slice. There is at best 1 a priori a day, so no need to be concerned with hourly average here.
		prio_date=[time.strftime('%m %d %y',time.gmtime(prio_time[i])) for i in prioID]

		model_date=[time.strftime('%m %d %y',time.gmtime(Tepoch[i[0]])) for i in tp]

		#good_profiles=[prio_date.index(i) for i in model_date if i==model_date[prio_date.index(i)]]
		good_profiles=[]
		gp=set(prio_date)
		for i in model_date:
			if i in gp:
				good_profiles.append(prioID[prio_date.index(i)])
			else:
				good_profiles.append('0')

		save='0'
		inte=0
		while save=='0':
			save=good_profiles[inte]
			inte=inte+1

		counter=0
		for v in good_profiles:
			if type(v)==int:
				save=v
			if type(v)==str:
				good_profiles[counter]=save
			counter=counter+1

		#when there is no a priori profile, use the last one available
		prio_co2=[prioco2[i] for i in good_profiles]
		prio_ch4=[prioch4[i] for i in good_profiles]
		prio_h2o=[prioh2o[i] for i in good_profiles]
		prio_hdo=[priohdo[i] for i in good_profiles]
		prio_density=[priodensity[i] for i in good_profiles]
		prio_pressure=[priopressure[i] for i in good_profiles]
		prio_gravity=[priogravity[i] for i in good_profiles]
		z_min=[round(zmin[i],2) for i in good_profiles]

		#declare the pressure weighting functions
		model_VC_co2=[0 for i in tp]
		priori_VC_co2=[0 for i in tp]
		model_VC_ch4=[0 for i in tp]
		priori_VC_ch4=[0 for i in tp]
		dryair_VC=[0 for i in tp]

		print '\nAverage per hour and a priori selection DONE', time.time()-timestart

		print 'N= ', len(tp)
		print 'COMPUTATIONS start',time.time()-timestart
		########################################################################################################################################################
		###########################################################     COMPUTATION     ########################################################################
		########################################################################################################################################################
		# good_profiles = list of indices of the a priori profiles to use from the full set of a priori profiles
		# prio_pressure = pressure levels of the TCCON a priori, two indices, prio_pressure[good_profiles #][level #]
		# len(tp) = the number of coincidences for the comparison

		##### Get a new pressure grid to interpolate the profiles on, it will then be averaged to layers.

		########################
		######   TCCON   #######
		######################## 

		height=np.arange(0,71,1) #original 71 levels with 1 km spacing
		new_height=np.arange(0,70.2,0.2) # intermediate 351 levels with 0.2 km spacing

		new_height=[j for j in new_height] #go from numpy array to simple list
		cutID=[[new_height.index(j) for j in new_height if j<=i] for i in z_min] # ID of levels below the station altitude

		#cut the levels below the station altitude
		new_height=[[round(i,1) for i in new_height] for j in good_profiles]
		for i in range(len(good_profiles)):
			for j in cutID[i]:
				del  (new_height[i][0])
			new_height[i].insert(0,z_min[i])

		#linear interpolation of log(pressure) levels
		log_pressure=[[log(i) for i in a] for a in prio_pressure] # log(prio_pressure)
		new_log_pressure=[np.interp(new_height[a],height,log_pressure[a]) for a in range(len(log_pressure))] # log(prio_pressure) on the new 351 level grid
		newp=[[round(exp(i),4) for i in a] for a in new_log_pressure] # prio_pressure on the new 351 levels grid

		# this loop determines the apropriate length of chunks to define layers.
		# i have to do this because the chunkIt will divide the profile in even chunks but the rest will either be in
		# the first (surface) chunk or the last (top) chunk. 
		chunks=5
		chunked=[[[0] for j in range(chunks)] for i in newp]
		while len(chunked[0][-1])<2:
			chunked=[chunkIt(i,chunks) for i in newp]
			if len(chunked[0][-1])==1:
				chunks=chunks-1
			print 'chunks=',chunks,len(chunked[0][-1])

		#### This part is time consuming

		# the pressure thickness of each chunk [first....(j=0)....last] dp(0)=first-last; [...(j-1)..last]GAP[...(j)..last], dp(j)=last(j-1)-last(j) to not loose the GAP contribution
		newdp=[[0 for j in i] for i in chunked]
		for i in range(len(chunked)):
			newdp[i][0]=abs(chunked[i][0][0]-chunked[i][0][-1])
			for j in range(1,len(chunked[i])):
				newdp[i][j]=abs(chunked[i][j-1][-1]-chunked[i][j][-1])

		###########################
		########   MODEL   ########
		########################### 

		log_pmodel=[[log(a) for a in i[::-1]] for i in pressureAVG] # log(pressure) for the model profile levels
		
		#Linearly interpolate height with log(pressure) to get the height levels of the model. Then do a linear regression to extrapolate model pressure down to the surface.
		height_model=[np.interp(log_pmodel[i][::-1],log_pressure[i][::-1],height[::-1])[::-1] for i in range(len(tp))]
		linfit=[np.polyfit(height_model[i],log_pmodel[i],1) for i in range(len(tp))]
		
		#extrapolation on both sides to go above 70km and below 0 km, this is necessary to be able to interpolate between 0 and 70 km.
		for a in range(len(tp)):
			#extrapolate above 70 km
			nextheight=height_model[a][-1]
			step=(70-nextheight)/5.0
			while nextheight<=71:
				nextheight=nextheight+step
				if nextheight<=71:
					log_pmodel[a]=np.append(log_pmodel[a],linfit[a][0]*nextheight+linfit[a][1])
					height_model[a]=np.append(height_model[a],nextheight)

			#extrapolate below 0 km	
			prevheight=height_model[a][0]
			if prevheight!=0.0:
				step=prevheight/5.0
			else:
				step=(-100.0)

			while prevheight>(-step):
				prevheight=prevheight-step
				if prevheight>(-step):
					log_pmodel[a]=np.insert(log_pmodel[a],0,linfit[a][0]*prevheight+linfit[a][1])
					height_model[a]=np.insert(height_model[a],0,prevheight)

		new_model_pressure=[[exp(i) for i in a][::-1] for a in log_pmodel] # model pressure on interpolatable levels 

		###########################################
		########   AVERAGING KERNELS    ###########
		###########################################

		# interpolate averaging kernel pressure onto newp levels (351 minus the levels below zmin)
		akP_1d_co2=[[] for i in range(len(tp))]
		akP_1d_ch4=[[] for i in range(len(tp))]
		for i in range(len(tp)):
			for j in range(len(ak_co2.T)):
				akP_1d_co2[i].append(np.interp(newp[i][::-1],ak_P[::-1],ak_co2.T[j][::-1])[::-1])
				akP_1d_ch4[i].append(np.interp(newp[i][::-1],ak_P[::-1],ak_ch4.T[j][::-1])[::-1])

		#averaging to get the averaging kernels on chunks LAYERS

		#co2
		chunked=[[chunkIt(j,chunks) for j in i] for i in akP_1d_co2]
		akPavg_co2=[[[0 for k in j] for j in i] for i in chunked]
		for i in range(len(chunked)):
			for j in range(len(chunked[i])):
				#special case of the first and last layers
				akPavg_co2[i][j][0]=(sum(chunked[i][j][0])-chunked[i][j][0][-1]/2.0)/(len(chunked[i][j][0])-1/2.0)
				akPavg_co2[i][j][-1]=(sum(chunked[i][j][-1])+chunked[i][j][-2][-1]/2.0)/(len(chunked[i][j][-1])+1/2.0)
				#general case for the rest of the layers
				for k in range(1,len(chunked[i][j])-1):
					akPavg_co2[i][j][k]=(sum(chunked[i][j][k])+(chunked[i][j][k-1][-1]-chunked[i][j][k][-1])/2.0)/(len(chunked[i][j][k]))

		#ch4
		chunked=[[chunkIt(j,chunks) for j in i] for i in akP_1d_ch4]
		akPavg_ch4=[[[0 for k in j] for j in i] for i in chunked]
		for i in range(len(chunked)):
			for j in range(len(chunked[i])):
				#special case of the first and last layers
				akPavg_ch4[i][j][0]=(sum(chunked[i][j][0])-chunked[i][j][0][-1]/2.0)/(len(chunked[i][j][0])-1/2.0)
				akPavg_ch4[i][j][-1]=(sum(chunked[i][j][-1])+chunked[i][j][-2][-1]/2.0)/(len(chunked[i][j][-1])+1/2.0)
				#general case for the rest of the layers
				for k in range(1,len(chunked[i][j])-1):
					akPavg_ch4[i][j][k]=(sum(chunked[i][j][k])+(chunked[i][j][k-1][-1]-chunked[i][j][k][-1])/2.0)/(len(chunked[i][j][k]))

		# interpolate averaging kernels onto measurements SZA
		#co2
		final_ak_co2=[[] for i in range(len(tp))]
		for i in range(len(akPavg_co2)):
			akP_T=np.array(akPavg_co2[i]).T
			for j in range(len(akP_T)):
				final_ak_co2[i].append(np.interp(sorted(TCCON_szaAVG),ak_zenith,akP_T[j])) #to interpolate, the range of TCCON SZA must be sorted with increasing values.

		#ch4
		final_ak_ch4=[[] for i in range(len(tp))]
		for i in range(len(akPavg_ch4)):
			akP_T=np.array(akPavg_ch4[i]).T
			for j in range(len(akP_T)):
				final_ak_ch4[i].append(np.interp(sorted(TCCON_szaAVG),ak_zenith,akP_T[j])) #to interpolate, the range of TCCON SZA must be sorted with increasing values.


		print 'Averaging kernels, DONE',time.time()-timestart

		##############################
		########   A PRIORI   ########
		############################## 

		##convert mixing ratios to mole fractions:
		apriori_co2=[]
		apriori_ch4=[]
		apriori_water=[]
		for i in range(len(tp)):
				apriori_co2.append(1/(1+1/(prio_co2[i]*10**(-6)))) # a priori CO2 (parts)
				apriori_ch4.append(1/(1+1/(prio_ch4[i]*10**(-9)))) # a priori ch4 (parts)
				apriori_water.append((1/(1+1/(prio_h2o[i])))+(1/(1+1/(SMOW*prio_hdo[i])))) #a priori H2O + HDO (parts)

		##linear interpolation to get mole fractions on newp levels (351 minus the levels below zmin)
		fit_co2=[]
		fit_ch4=[]
		fit_water=[]
		#also use gravity
		fit_g=[]
		for i in range(len(tp)):
			fit_co2.append(np.interp(newp[i][::-1],prio_pressure[i][::-1],apriori_co2[i][::-1])[::-1])
			fit_ch4.append(np.interp(newp[i][::-1],prio_pressure[i][::-1],apriori_ch4[i][::-1])[::-1])
			fit_water.append(np.interp(newp[i][::-1],prio_pressure[i][::-1],apriori_water[i][::-1])[::-1])
			fit_g.append(np.interp(newp[i][::-1],prio_pressure[i][::-1],prio_gravity[i][::-1])[::-1])

		##averaging
		#co2 profiles
		chunked=[chunkIt(i,chunks) for i in fit_co2]

		fit_avg_co2=[[0 for j in i] for i in chunked]
		for i in range(len(chunked)):
			fit_avg_co2[i][0]=(sum(chunked[i][0])-chunked[i][0][-1]/2.0)/(len(chunked[i][0])-1/2.0)
			fit_avg_co2[i][-1]=(sum(chunked[i][-1])+chunked[i][-2][-1]/2.0)/(len(chunked[i][-1])+1/2.0)
			for j in range(1,len(chunked[i])-1):
				fit_avg_co2[i][j]=(sum(chunked[i][j])+(chunked[i][j-1][-1]-chunked[i][j][-1])/2.0)/(len(chunked[i][j]))

		#ch4 profiles
		chunked=[chunkIt(i,chunks) for i in fit_ch4]

		fit_avg_ch4=[[0 for j in i] for i in chunked]
		for i in range(len(chunked)):
			fit_avg_ch4[i][0]=(sum(chunked[i][0])-chunked[i][0][-1]/2.0)/(len(chunked[i][0])-1/2.0)
			fit_avg_ch4[i][-1]=(sum(chunked[i][-1])+chunked[i][-2][-1]/2.0)/(len(chunked[i][-1])+1/2.0)
			for j in range(1,len(chunked[i])-1):
				fit_avg_ch4[i][j]=(sum(chunked[i][j])+(chunked[i][j-1][-1]-chunked[i][j][-1])/2.0)/(len(chunked[i][j]))

		#water profiles
		chunked=[chunkIt(i,chunks) for i in fit_water]

		fit_avg_water=[[0 for j in i] for i in chunked]
		for i in range(len(chunked)):
			fit_avg_water[i][0]=(sum(chunked[i][0])-chunked[i][0][-1]/2.0)/(len(chunked[i][0])-1/2.0)
			fit_avg_water[i][-1]=(sum(chunked[i][-1])+chunked[i][-2][-1]/2.0)/(len(chunked[i][-1])+1/2.0)
			for j in range(1,len(chunked[i])-1):
				fit_avg_water[i][j]=(sum(chunked[i][j])+(chunked[i][j-1][-1]-chunked[i][j][-1])/2.0)/(len(chunked[i][j]))

		#wgravity profiles
		chunked=[chunkIt(i,chunks) for i in fit_g]

		fit_avg_g=[[0 for j in i] for i in chunked]
		for i in range(len(chunked)):
			fit_avg_g[i][0]=(sum(chunked[i][0])-chunked[i][0][-1]/2.0)/(len(chunked[i][0])-1/2.0)
			fit_avg_g[i][-1]=(sum(chunked[i][-1])+chunked[i][-2][-1]/2.0)/(len(chunked[i][-1])+1/2.0)
			for j in range(1,len(chunked[i])-1):
				fit_avg_g[i][j]=(sum(chunked[i][j])+(chunked[i][j-1][-1]-chunked[i][j][-1])/2.0)/(len(chunked[i][j]))
		
		####################################################
		########   COMPUTE TCCON A PRIORI COLUMNS   ########
		#################################################### 
		TCCON_column_co2=[0 for i in tp]
		TCCON_column_ch4=[0 for i in tp]
		TCCON_dryair=[0 for i in tp]
		for i in range(len(tp)):
			for j in range(len(newdp[0])):
				TCCON_column_co2[i]=TCCON_column_co2[i]+(fit_avg_co2[i][j]*newdp[i][j]*100)/(fit_avg_g[i][j]*(mass_H2O*fit_avg_water[i][j]+mass_dry_air*(1-fit_avg_water[i][j])))
				TCCON_column_ch4[i]=TCCON_column_ch4[i]+(fit_avg_ch4[i][j]*newdp[i][j]*100)/(fit_avg_g[i][j]*(mass_H2O*fit_avg_water[i][j]+mass_dry_air*(1-fit_avg_water[i][j])))
				TCCON_dryair[i]=TCCON_dryair[i]+((1-fit_avg_water[i][j])*newdp[i][j]*100)/(fit_avg_g[i][j]*(mass_H2O*fit_avg_water[i][j]+mass_dry_air*(1-fit_avg_water[i][j])))

		TCCON_col_co2=[]
		TCCON_col_ch4=[]
		for i in range(len(tp)):
			TCCON_col_co2.append((10**6)*TCCON_column_co2[i]/TCCON_dryair[i])
			TCCON_col_ch4.append((10**9)*TCCON_column_ch4[i]/TCCON_dryair[i])

		print 'A priori, DONE',time.time()-timestart

		###########################
		########   MODEL   ########
		########################### 

		model_co2=[]
		model_ch4=[]
		model_water=[]
		for i in range(len(tp)):
			model_co2.append(np.array(resAVG_co2[i])*10**(-6)) #EC-CAS CO2 (parts)
			model_ch4.append(np.array(resAVG_ch4[i])*10**(-9)) #EC-CAS ch4 (parts)
			model_water.append(np.array(waterAVG[i])) #EC-CAS specific humidity (kgH2O/kgAIR)

		### interpolations

		#interpolation to get the new_model_co2
		new_model_co2=[np.interp(new_model_pressure[i],pressureAVG[i],model_co2[i]) for i in range(len(tp))] # EC-CAS CO2 on the 71 TCCON levels
		new_model_ch4=[np.interp(new_model_pressure[i],pressureAVG[i],model_ch4[i]) for i in range(len(tp))] # EC-CAS ch4 on the 71 TCCON levels
		new_model_water=[np.interp(new_model_pressure[i],pressureAVG[i],model_water[i]) for i in range(len(tp))] # EC-CAS specific humidity on the 71 TCCON levels

		newpres=[[i/100.0 for i in j] for j in new_model_pressure] # 71 pressure levels in hPa

		# EC-CAS profiles interpolation to newp:

		model_1d_co2=[]
		model_1d_ch4=[]
		model_1d_water=[]
		# the plots are to check the quality of the interpolations steps for the profiles, they should all be overlapped
		for i in range(len(tp)):
			model_1d_co2.append(np.interp(newp[i][::-1],newpres[i],new_model_co2[i])[::-1]) # EC-CAS CO2 on the 351 levels
			model_1d_ch4.append(np.interp(newp[i][::-1],newpres[i],new_model_ch4[i])[::-1]) # EC-CAS ch4 on the 351 levels
			model_1d_water.append(np.interp(newp[i][::-1],newpres[i],new_model_water[i])[::-1]) # EC-CAS specific humidity on the 351 levels

		testp=[i/100.0 for i in model_pressure[tp[0][0]]]
		test_co2=[i*10**(-6) for i in model_co2[tp[0][0]]]
		test_ch4=[i*10**(-6) for i in model_ch4[tp[0][0]]]
		test_water=[i for i in model_water[tp[0][0]]]
		testAVG_co2=[i*10**(-6) for i in resAVG_co2[0]]
		testAVG_ch4=[i*10**(-6) for i in resAVG_ch4[0]]
		testwaterAVG=[i for i in waterAVG[0]]
		testpAVG=[i/100.0 for i in pressureAVG[0]]
		fig=plt.Figure()
		fig.set_canvas(plt.gcf().canvas)
		plt.plot(model_1d_co2[0],newp[0],c='red',label='351 lvls') # interpolation on new 351 levels
		plt.plot(new_model_co2[0],newpres[0],c='blue',label='71 lvls') # interpolation on TCCON 71 levels
		plt.plot(testAVG_co2,testpAVG,c='purple',label='original') # original averaged
		plt.plot(test_co2,testp,c='green') #original
		plt.xlabel('$CO_2$ (parts)')
		plt.ylabel('Pressure (hPa)')
		plt.title(Tnames[files.index(nameoffile)])
		plt.xlim([0.00036,0.000415])
		plt.legend()
		plt.gca().invert_yaxis()
		fig.savefig(os.path.join(save_path,''.join([list(nameoffile)[0],list(nameoffile)[1],list(nameoffile)[2]])+'CO2-interp.png'))
		plt.clf()
		plt.plot(model_1d_water[0],newp[0],c='red',label='351 lvls')  # interpolation on new 351 levels
		plt.plot(new_model_water[0],newpres[0],c='blue',label='71 lvls') # interpolation on TCCON 71 levels
		plt.plot(testwaterAVG,testpAVG,c='purple',label='original') # original averaged
		plt.plot(test_water,testp,c='green')  #original
		plt.xlabel('$H_{2}O$ (parts)')
		plt.ylabel('Pressure (hPa)')
		plt.title(Tnames[files.index(nameoffile)])
		plt.legend()
		plt.gca().invert_yaxis()
		fig.savefig(os.path.join(save_path,''.join([list(nameoffile)[0],list(nameoffile)[1],list(nameoffile)[2]])+'H2O-interp.png'))
		plt.clf()
		plt.plot(model_1d_ch4[0],newp[0],c='red',label='351 lvls') # interpolation on new 351 levels
		plt.plot(new_model_ch4[0],newpres[0],c='blue',label='71 lvls') # interpolation on TCCON 71 levels
		plt.plot(testAVG_ch4,testpAVG,c='purple',label='original') # original averaged
		plt.plot(test_ch4,testp,c='green') #original
		plt.xlabel('$CH_4$ (parts)')
		plt.ylabel('Pressure (hPa)')
		plt.title(Tnames[files.index(nameoffile)])
		plt.legend()
		plt.gca().invert_yaxis()
		fig.savefig(os.path.join(save_path,''.join([list(nameoffile)[0],list(nameoffile)[1],list(nameoffile)[2]])+'CH4-interp.png'))
		plt.clf()
		#break

		#averaging to get mixing ratios on chunks LAYERS
		#co2 profile
		chunked=[chunkIt(i,chunks) for i in model_1d_co2]
		fit_avg_co2_m=[[0 for j in i] for i in chunked]
		for i in range(len(chunked)):
			fit_avg_co2_m[i][0]=(sum(chunked[i][0])-chunked[i][0][-1]/2.0)/(len(chunked[i][0])-1/2.0)
			fit_avg_co2_m[i][-1]=(sum(chunked[i][-1])+chunked[i][-2][-1]/2.0)/(len(chunked[i][-1])+1/2.0)
			for j in range(1,len(chunked[i])-1):
				fit_avg_co2_m[i][j]=(sum(chunked[i][j])+(chunked[i][j-1][-1]-chunked[i][j][-1])/2.0)/(len(chunked[i][j]))

		#ch4 profile
		chunked=[chunkIt(i,chunks) for i in model_1d_ch4]
		fit_avg_ch4_m=[[0 for j in i] for i in chunked]
		for i in range(len(chunked)):
			fit_avg_ch4_m[i][0]=(sum(chunked[i][0])-chunked[i][0][-1]/2.0)/(len(chunked[i][0])-1/2.0)
			fit_avg_ch4_m[i][-1]=(sum(chunked[i][-1])+chunked[i][-2][-1]/2.0)/(len(chunked[i][-1])+1/2.0)
			for j in range(1,len(chunked[i])-1):
				fit_avg_ch4_m[i][j]=(sum(chunked[i][j])+(chunked[i][j-1][-1]-chunked[i][j][-1])/2.0)/(len(chunked[i][j]))

		#water profile
		chunked=[chunkIt(i,chunks) for i in model_1d_water]
		fit_avg_water_m=[[0 for j in i] for i in chunked]
		for i in range(len(chunked)):
			fit_avg_water_m[i][0]=(sum(chunked[i][0])-chunked[i][0][-1]/2.0)/(len(chunked[i][0])-1/2.0)
			fit_avg_water_m[i][-1]=(sum(chunked[i][-1])+chunked[i][-2][-1]/2.0)/(len(chunked[i][-1])+1/2.0)
			for j in range(1,len(chunked[i])-1):
				fit_avg_water_m[i][j]=(sum(chunked[i][j])+(chunked[i][j-1][-1]-chunked[i][j][-1])/2.0)/(len(chunked[i][j]))

		###########################################
		########   COMPUTE MODEL COLUMNS   ########
		########################################### 
		model_column_co2=[0 for i in tp]
		model_column_ch4=[0 for i in tp]
		model_dryair=[0 for i in tp]
		for i in range(len(tp)):
			for j in range(len(newdp[0])):
				model_column_co2[i]=model_column_co2[i]+((1-fit_avg_water_m[i][j])*fit_avg_co2_m[i][j]*newdp[i][j]*100)/(fit_avg_g[i][j]*mass_dry_air)
				model_column_ch4[i]=model_column_ch4[i]+((1-fit_avg_water_m[i][j])*fit_avg_ch4_m[i][j]*newdp[i][j]*100)/(fit_avg_g[i][j]*mass_dry_air)
				model_dryair[i]=model_dryair[i]+((1-fit_avg_water_m[i][j])*newdp[i][j]*100)/(fit_avg_g[i][j]*mass_dry_air)

		model_DMF_co2=[]
		model_DMF_ch4=[]
		for i in range(len(model_column_co2)):
			model_DMF_co2.append((10**6)*model_column_co2[i]/model_dryair[i])
			model_DMF_ch4.append((10**9)*model_column_ch4[i]/model_dryair[i])

		print 'EC-CAS, DONE',time.time()-timestart

		sortSZA=sorted(TCCON_szaAVG)
		## final_ak_co2[i][j][sortSZA.index(TCCON_szaAVG[i])] 
		# i loops over the # of coincidences
		# j loops over the # of pressure layers
		# recall the final column averaging kernels final_ak_co2 were sorted with respect to SZA
		# the ith TCCON measurement has a SZA of TCCON_szaAVG[i]
		# to get the averaging kernel in final_ak_co2 at the jth level associated with the ith TCCON measurement, you must write final_ak_co2[i][j][sortSZA.index(TCCON_szaAVG[i])]

		###############################
		########   SMOOTHING   ########
		###############################
		SMOOTH_co2=[]
		SMOOTH_co2g=[]
		Dif_co2=[]
		Dif_co2g=[]
		nDif_co2=[]

		SMOOTH_ch4=[]
		SMOOTH_ch4g=[]
		Dif_ch4=[]
		Dif_ch4g=[]
		nDif_ch4=[]

		for i in range(len(tp)):
			progress(i+1,len(tp))
			for j in range(len(newdp[0])):
				model_VC_co2[i]=model_VC_co2[i]+(final_ak_co2[i][j][sortSZA.index(TCCON_szaAVG[i])]*(1-fit_avg_water_m[i][j])*fit_avg_co2_m[i][j]*newdp[i][j]*100)/(fit_avg_g[i][j]*mass_dry_air)
				priori_VC_co2[i]=priori_VC_co2[i]+(final_ak_co2[i][j][sortSZA.index(TCCON_szaAVG[i])]*fit_avg_co2[i][j]*newdp[i][j]*100)/(fit_avg_g[i][j]*(mass_H2O*fit_avg_water[i][j]+mass_dry_air*(1-fit_avg_water[i][j])))

				model_VC_ch4[i]=model_VC_ch4[i]+(final_ak_ch4[i][j][sortSZA.index(TCCON_szaAVG[i])]*(1-fit_avg_water_m[i][j])*fit_avg_ch4_m[i][j]*newdp[i][j]*100)/(fit_avg_g[i][j]*mass_dry_air)
				priori_VC_ch4[i]=priori_VC_ch4[i]+(final_ak_ch4[i][j][sortSZA.index(TCCON_szaAVG[i])]*fit_avg_ch4[i][j]*newdp[i][j]*100)/(fit_avg_g[i][j]*(mass_H2O*fit_avg_water[i][j]+mass_dry_air*(1-fit_avg_water[i][j])))


			SMOOTH_co2.append(TCCON_col_co2[i]+(10**6)*(model_VC_co2[i]-priori_VC_co2[i])/TCCON_dryair[i])
			SMOOTH_co2g.append(0.989*TCCON_col_co2[i]+(10**6)*(model_VC_co2[i]-0.989*priori_VC_co2[i])/TCCON_dryair[i])

			SMOOTH_ch4.append(TCCON_col_ch4[i]+(10**9)*(model_VC_ch4[i]-priori_VC_ch4[i])/TCCON_dryair[i])
			SMOOTH_ch4g.append(0.978*TCCON_col_ch4[i]+(10**9)*(model_VC_ch4[i]-0.978*priori_VC_ch4[i])/TCCON_dryair[i])

			Dif_co2.append(SMOOTH_co2[i]-xco2AVG[i])
			Dif_co2g.append(SMOOTH_co2g[i]-xco2AVG[i])
			nDif_co2.append(model_DMF_co2[i]-xco2AVG[i])

			Dif_ch4.append(SMOOTH_ch4[i]-xch4AVG[i])
			Dif_ch4g.append(SMOOTH_ch4g[i]-xch4AVG[i])
			nDif_ch4.append(model_DMF_ch4[i]-xch4AVG[i])

		print '\nSmoothing, DONE',time.time()-timestart

		print 'CO2:'
		print '\tBias is:', sum(Dif_co2)/len(Dif_co2)
		print '\tBias is:', sum(Dif_co2g)/len(Dif_co2g), 'with scale factor'
		print '\tBias is:', sum(nDif_co2)/len(nDif_co2), 'without smoothing\n'

		print 'CH4:'
		print '\tBias is:', sum(Dif_ch4)/len(Dif_ch4)
		print '\tBias is:', sum(Dif_ch4g)/len(Dif_ch4g), 'with scale factor'
		print '\tBias is:', sum(nDif_ch4)/len(nDif_ch4), 'without smoothing\n'

		save_data_co2[files.index(nameoffile)]=TCCON_col_co2
		save_data_co2[files.index(nameoffile)+len(files)]=xco2AVG
		save_data_co2[files.index(nameoffile)+2*len(files)]=model_DMF_co2
		save_data_co2[files.index(nameoffile)+3*len(files)]=SMOOTH_co2

		save_data_ch4[files.index(nameoffile)]=TCCON_col_ch4
		save_data_ch4[files.index(nameoffile)+len(files)]=xch4AVG
		save_data_ch4[files.index(nameoffile)+2*len(files)]=model_DMF_ch4
		save_data_ch4[files.index(nameoffile)+3*len(files)]=SMOOTH_ch4

		save_time[files.index(nameoffile)]=xtimeAVG
		save_time[files.index(nameoffile)+len(files)]=TepochAVG
		print 'DONE',time.time()-timestart

		#break


########################################
########   WRITE OUTPUT FILES   ########
######################################## 

outfile=open(save_time_path,'w')
pickle.dump(save_time,outfile)
outfile.close()

outfile=open(save_data_co2_path,'w')
pickle.dump(save_data_co2,outfile)
outfile.close()

outfile=open(save_data_ch4_path,'w')
pickle.dump(save_data_ch4,outfile)
outfile.close()

print 'DONE',time.time()-codetsart
