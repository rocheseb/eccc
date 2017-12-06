#!/usr/bin/env python2.7
 # -*- coding: ascii -*-

####################
# Code Descirption #
####################

'''
Reads the output of EC_SMOOTH100.py and creates plots of the data
'''

#############
# Libraries #
#############

import os.path
import pickle
from math import sqrt

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import gridspec

from scipy.stats import pearsonr

import numpy as np

import time
import calendar
from datetime import datetime

# interactive plots
from bokeh.plotting import figure, output_file
from bokeh.models import Legend, Panel, Tabs, CustomJS, ColumnDataSource, CheckboxGroup, RadioGroup, Button, VBox, Range1d
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.layouts import gridplot, widgetbox
from bokeh.resources import CDN
from bokeh.embed import file_html

def duplicates(lst,item):
	return [i for i, x in enumerate(lst) if x == item]

savefix=0

################################################
#####    Hex size for correlation plots    #####
################################################
hexsize=20                                 #####
################################################

#######################
#####    PATHS    #####
#######################

sites_path ='/data/02/sroche/EC-CAS & TCCON data comparison/main/Data/TCCON7'
file_path_co2='/data/02/sroche/EC-CAS & TCCON data comparison/main/SMOOTH100/mn117b/DMF_co2.txt'
file_path_ch4='/data/02/sroche/EC-CAS & TCCON data comparison/main/SMOOTH100/mn117b/DMF_ch4.txt'
file_path2='/data/02/sroche/EC-CAS & TCCON data comparison/main/SMOOTH100/mn117b/time.txt'
save_path='/data/02/sroche/EC-CAS & TCCON data comparison/main/SMOOTH100/mn117b'

#######################
#####    SETUP    #####
#######################

######## /!\ important time handling to make sure times don't get shifted from UTC due to computer environment variables when using datetime objects ###########
os.environ['TZ'] = 'UTC'
time.tzset()
################################################################################################################################################################

kelly_colors = {	
					'vivid_yellow':(255, 179, 0),
					'strong_purple':(128, 62, 117),
					'vivid_orange':(255, 104, 0),
					'very_light_blue':(166, 189, 215),
					'vivid_red':(193, 0, 32),
					#'grayish_yellow':(206, 162, 98),
					#'medium_gray':(129, 112, 102),

					# these aren't good for people with defective color vision:
					'vivid_green':(0, 125, 52),
					'strong_purplish_pink':(246, 118, 142),
					'strong_blue':(0, 83, 138),
					'strong_yellowish_pink':(255, 122, 92),
					'strong_violet':(83, 55, 122),
					'vivid_orange_yellow':(255, 142, 0),
					'strong_purplish_red':(179, 40, 81),
					'vivid_greenish_yellow':(244, 200, 0),
					'strong_reddish_brown':(127, 24, 13),
					'vivid_yellowish_green':(147, 170, 0),
					'deep_yellowish_brown':(89, 51, 21),
					'vivid_reddish_orange':(241, 58, 19),
					#'dark_olive_green':(35, 44, 22),

					#filling in for the last sites
					'gold':'gold',
					'chartreuse':'chartreuse',
					'cyan':'cyan',
					'firebrick':'firebrick',
					'lightsalmon':'lightsalmon',
					'peru':'peru',
					'goldenrod':'goldenrod',
					'navy':'navy',
					'green':'green'
				}

tabs = []	

months=['January','February','March','April','May','June','July','August','September','October','November','December']

DATA = {}
DATA['co2']={}
DATA['ch4']={}

### pick the co2 columns data
infile=open(file_path_co2,'r')
DMF_list=pickle.load(infile)
infile.close()

DATA['co2']['priori']=[[] for i in range(len(DMF_list)/4)]
DATA['co2']['TCCON']=[[] for i in range(len(DMF_list)/4)]
DATA['co2']['CAS']=[[] for i in range(len(DMF_list)/4)]
DATA['co2']['smooth']=[[] for i in range(len(DMF_list)/4)]

#priori=[[] for i in range(len(DMF_list)/4)]
#TCCON=[[] for i in range(len(DMF_list)/4)]
#CAS=[[] for i in range(len(DMF_list)/4)]
#smooth=[[] for i in range(len(DMF_list)/4)]

for i in range(0,len(DMF_list)/4):
	DATA['co2']['priori'][i]=[item for item in DMF_list[i]]
	DATA['co2']['TCCON'][i]=[item for item in DMF_list[i+(len(DMF_list)/4)]]
	DATA['co2']['CAS'][i]=[item for item in DMF_list[i+2*(len(DMF_list)/4)]]
	DATA['co2']['smooth'][i]=[item for item in DMF_list[i+3*(len(DMF_list)/4)]]
	#priori[i]=[item for item in DMF_list[i]]
	#TCCON[i]=[item for item in DMF_list[i+(len(DMF_list)/4)]]
	#CAS[i]=[item for item in DMF_list[i+2*(len(DMF_list)/4)]]
	#smooth[i]=[item for item in DMF_list[i+3*(len(DMF_list)/4)]]

### pick the ch4 columns data
infile=open(file_path_ch4,'r')
DMF_list=pickle.load(infile)
infile.close()

DATA['ch4']['priori']=[[] for i in range(len(DMF_list)/4)]
DATA['ch4']['TCCON']=[[] for i in range(len(DMF_list)/4)]
DATA['ch4']['CAS']=[[] for i in range(len(DMF_list)/4)]
DATA['ch4']['smooth']=[[] for i in range(len(DMF_list)/4)]

for i in range(0,len(DMF_list)/4):
	DATA['ch4']['priori'][i]=[item for item in DMF_list[i]]
	DATA['ch4']['TCCON'][i]=[item for item in DMF_list[i+(len(DMF_list)/4)]]
	DATA['ch4']['CAS'][i]=[item for item in DMF_list[i+2*(len(DMF_list)/4)]]
	DATA['ch4']['smooth'][i]=[item for item in DMF_list[i+3*(len(DMF_list)/4)]]

### pick the time data
infile=open(file_path2,'r')
time_list=pickle.load(infile)
infile.close()

TCCON_time=[[] for i in range(len(time_list)/2)]
DMF_time=[[] for i in range(len(time_list)/2)]

for i in range(0,len(time_list)/2):
	TCCON_time[i]=[item for item in time_list[i]]
	DMF_time[i]=[item for item in time_list[i+(len(time_list)/2)]]

#Station labels with 3 first letters. I then append a blank at the beginning and end so that there will
#be some space on both ends of the plot before the first labels.
datafiles=os.listdir(sites_path)
Label_station=[''.join([list(nameoffile)[0],list(nameoffile)[1],list(nameoffile)[2]]) for nameoffile in datafiles]
Label_station.insert(0,'')
Label_station.append('')


### MODIFY when the set of sites changes
Tnames=[]
for station in datafiles:
	Tname=''
	for letter in station:
		try:
			if float(letter)/float(letter)==True:
				break
		except ValueError:
			pass
		Tname+=letter
	if Tname=='ParkFalls':
		Tname='Park Falls'
	if Tname=='ReunionIsland':
		Tname='Reunion Island'
	Tnames.append(Tname)

###########################
#####    USER INPUT   #####
###########################
B='0'
while B!='y' and B!='n':
	B=raw_input('Use statistics with smoothed columns? (\"y\" / \"n\") \n')
	if B!='y' and B!='n':
		print 'Wrong entry: you have to type either "y" or "n" (no quotes)'

for k in range(1,len(Label_station)-1):
	print k-1,': ',Label_station[k]
try:
	custom_station=int(raw_input('Select site (enter the number) or type a letter/word to do all sites: '))
	all_check=0
except ValueError:
	all_check='all'
	custom_station=0

print 'Select starting day:'
custom_year=int(raw_input('Year: '))

if all_check==0:
	m_check=False
	while m_check==False:
		custom_month=int(raw_input('Month (1 for January, 12 for December): '))

		#Get the indices of pairs separated for each month.
		month_tags=[int(time.strftime('%m',time.gmtime(i))) for i in TCCON_time[custom_station]]
		m_idx=[]
		for i in range(1,13):
			m_idx.append(duplicates(month_tags,i))
			if len(m_idx)==0:
				print '  No measurements in ', months[i-1]
			else:
				m_check=True

	custom_check=False
	while custom_check==False:
		custom_day=int(raw_input('Day: '))
		custom_select=(custom_year,custom_month,custom_day,0,0,0)
		span=int(raw_input('Span (in days): '))
		custom_next=(custom_year,custom_month,custom_day+span,0,0,0)

		#Get the indices of pairs separated for each month.
		month_tags=[int(time.strftime('%m',time.gmtime(i))) for i in TCCON_time[custom_station]]
		m_idx=[]
		for i in range(1,13):
			m_idx.append(duplicates(month_tags,i))
			if len(m_idx)==0:
				print '	No measurements in ', months[i-1]

		custom_tags=[]
		for i in TCCON_time[custom_station]:
			if i>=calendar.timegm(custom_select) and i<calendar.timegm(custom_next):
				custom_tags.append(TCCON_time[custom_station].index(i))
		if len(custom_tags)!=0:
			custom_check=True
		else:
			print '  No measurements for: ',custom_year,'-',custom_month,'-',custom_day,'\n'
else:
	custom_month=int(raw_input('Month (1 for January, 12 for December): '))
	custom_day=int(raw_input('Day: '))
	span=int(raw_input('Span (in days): '))
	custom_select=(custom_year,custom_month,custom_day,0,0,0)
	custom_next=(custom_year,custom_month,custom_day+span,0,0,0)


###########################
#####    MAIN CODE    #####
###########################

# loop over species
for var in DATA:
	print 'Statistics for', var,':\n'
	CAS = DATA[var]['CAS']
	smooth = DATA[var]['smooth']
	TCCON = DATA[var]['TCCON']
	priori = DATA[var]['priori']

	Station_range=range(len(CAS))

	if var == 'co2':
		gas = '$XCO_2$ (ppm)'
		units = '(ppm)'
	if var == 'ch4':
		gas = '$XCH_4$ (ppb)'
		units = '(ppb)'
	# save data in a csv file.
	csvpath=os.path.join(save_path,var+'_compare.csv')
	csv=open(csvpath,'w')
	for i in range(0,len(CAS)):
		for j in range(0,len(CAS[i])):
			csv.write(str(CAS[i][j])+','+str(smooth[i][j])+','+str(TCCON[i][j])+'\n')
	csv.close()

	csvpath=os.path.join(save_path,var+'_stats.csv')
	csv=open(csvpath,'w')

	N_station=[[] for i in range(len(CAS))]

	if isinstance(all_check,basestring)==True:
		ALL_TCCON=[]
		ALL_COMPA=[]
		Bdev=[]

	############################
	#####    BOKEH PLOT    #####
	############################

	stations = Label_station[1:-1]


	iterable = [(stations[k],{'smooth':smooth[k],'TCCON':TCCON[k],'DMF_time':[datetime(*time.gmtime(i)[:6]) for i in DMF_time[k]],'TCCON_time':[datetime(*time.gmtime(i)[:6]) for i in TCCON_time[k]]}) for k in range(len(stations))]

	STRUCT = {key:value for key,value in iterable if len(value['smooth'])!=0}

	temp = [0 for site in STRUCT]

	latlist = [28.30,45.94,67.37,49.10,36.60,47.48,53.23,-34.41,47.97,-12.43,-45.06,-45.05,80.05]

	lat_ordered_ID = [latlist.index(lat) for lat in sorted(latlist)[::-1]]

	lat_ordered_sites = [[site for site in STRUCT][i] for i in lat_ordered_ID]

	names = [Tnames[stations.index(site)] for site in lat_ordered_sites]
	for i in range(len(names)):
		if names[i]=='Lauder':
			names[i]='Lauder 120'
		if names[i]=='LLauder':
			names[i]='Lauder 125'

	table_source = ColumnDataSource( data = {'site':names,'latitude':sorted(latlist)[::-1],'N':temp,'RMS':temp,'Bias':temp,'Scatter':temp,'R':temp} )

	for site in lat_ordered_sites:
		ID = STRUCT.keys().index(site)
		color = kelly_colors.keys()[ID]
		STRUCT[site]['color'] = kelly_colors[color]

	sources = {}
	sources2 = {}
	cor_sources = {}
	count = 0
	for site in lat_ordered_sites:
		sources[site] = ColumnDataSource(data = {'x':STRUCT[site]['DMF_time'],'y':STRUCT[site]['smooth']})
		sources2[site] = ColumnDataSource(data = {'x':STRUCT[site]['TCCON_time'],'y':STRUCT[site]['TCCON']})
		cor_sources[site] = ColumnDataSource(data={'x':[],'y':[]})
		count += 1

	columns = [ 
				TableColumn(field='site',title='Site'),
				TableColumn(field='latitude',title='Latitude'),
				TableColumn(field='N',title='N'),
				TableColumn(field='RMS',title='RMS '+units),
				TableColumn(field='Bias',title='Bias '+units),
				TableColumn(field='Scatter',title='Scatter '+units),
				TableColumn(field='R',title='R'),
				]

	### TCCON.html
	min_x = min([STRUCT[site]['DMF_time'][0] for site in lat_ordered_sites])
	max_x = max([STRUCT[site]['DMF_time'][-1] for site in lat_ordered_sites])

	min_y = min([min(STRUCT[site]['TCCON']) for site in lat_ordered_sites])
	max_y = max([max(STRUCT[site]['TCCON']) for site in lat_ordered_sites])

	ampl_y = max_y-min_y

	min_y = min_y - ampl_y*10/100
	max_y = max_y + ampl_y*10/100

	milestone = time.time()

	TOOLS = ["pan,wheel_zoom,box_zoom,undo,redo,reset,box_select,save"]

	bokeh_fig = figure(webgl = True, title = var+' time series', y_range=[min_y,max_y], plot_width = 750, plot_height = 250+20*len(lat_ordered_sites), tools = TOOLS, toolbar_location = 'right', x_axis_type='datetime', x_range = Range1d(min_x,max_x) )

	bokeh_fig.tools[-2].dimensions = 'width'

	plots=[]
	for site in lat_ordered_sites:
		plots.append( bokeh_fig.scatter(x='x',y='y',color=STRUCT[site]['color'],alpha=0.5,source=sources[site]) )
		plots.append( bokeh_fig.scatter(x='x',y='y',color='black',alpha=0.5,source=sources2[site]) )


	data_table = DataTable(source=table_source, columns=columns, width= 700, height=400)

	count = 0
	for site in lat_ordered_sites:
		sources[site].callback = CustomJS(args = dict(s2=sources2[site],dt=data_table,scor=cor_sources[site]), code="""
		var inds = cb_obj.get('selected')['1d'].indices;
		var d1 = cb_obj.get('data');
		var d2 = s2.get('data');
		var tab = dt.get('source').get('data');
		var dcor = scor.get('data');

		var difm = 0;
		var difm2 = 0;
		var scat = 0;

		var ym1 = 0;
		var ym2 = 0;

		var T1 = 0;
		var T2 = 0;
		var T3 = 0;

		tab['N']["""+str(count)+"""] = inds.length;

		if (inds.length == 0) {
			tab['RMS']["""+str(count)+"""] = 0;
			tab['Bias']["""+str(count)+"""] = 0;
			tab['Scatter']["""+str(count)+"""] = 0;
			tab['R']["""+str(count)+"""] = 0;
			dt.trigger('change');
			return;
		}

		for (i=0; i < inds.length; i++){
			difm += d1['y'][inds[i]] - d2['y'][inds[i]];
			difm2 += Math.pow(d1['y'][inds[i]] - d2['y'][inds[i]],2);
			ym1 += d1['y'][inds[i]];
			ym2 += d2['y'][inds[i]];

			dcor['x'].push(d2['y'][inds[i]]);
			dcor['y'].push(d1['y'][inds[i]]);
		}

		difm /= inds.length;
		difm2 /= inds.length;
		ym1 /= inds.length;
		ym2 /= inds.length;

		for (i=0; i < inds.length; i++){
			scat += Math.pow(d1['y'][inds[i]] - d2['y'][inds[i]] - difm,2);
		}

		for (i=0; i < inds.length; i++){
			T1 += (d1['y'][inds[i]] - ym1)*(d2['y'][inds[i]] - ym2);
			T2 += Math.pow(d1['y'][inds[i]] - ym1,2);
			T3 += Math.pow(d2['y'][inds[i]] - ym2,2);
		}

		tab['RMS']["""+str(count)+"""] = Math.sqrt(difm2).toFixed(2);
		tab['Bias']["""+str(count)+"""] = difm.toFixed(2);
		tab['Scatter']["""+str(count)+"""] = Math.sqrt(scat/(inds.length -1)).toFixed(2);
		tab['R']["""+str(count)+"""] = (T1/Math.sqrt(T2*T3)).toFixed(2);

		dt.trigger('change');
		scor.trigger('change');
		""")
		count += 1

	N_plots = range(len(plots))

	N_plots2 = range(len(plots)/2)

	even_plots = [i for i in N_plots if i%2 == 0]

	#legend=Legend(legends=[('GEM-MACH-GHG',[plots[0]]),('TCCON',[plots[1]])],location=(0,0))
	legend=Legend(legends=[('TCCON',[plots[1]])]+[([site for site in lat_ordered_sites][i],[plots[even_plots[i]]]) for i in range(len([site for site in lat_ordered_sites]))],location=(0,0))
	bokeh_fig.add_layout(legend,'right')
	bokeh_fig.yaxis.axis_label = 'x'+var+' '+units
	bokeh_fig.xaxis.axis_label = 'Time'

	# correlation plot
	bokeh_corfig = figure(webgl = True, title = 'Correlations', plot_width = 400, plot_height = 400, x_range = [min_y,max_y], y_range = [min_y,max_y]) 
	bokeh_corfig.toolbar.logo = None
	bokeh_corfig.toolbar_location = None
	bokeh_corfig.xaxis.axis_label = ' '.join(['TCCON',var])
	bokeh_corfig.yaxis.axis_label = ' '.join(['GEM-MACH-GHG',var])

	corplots = []
	linerange = list(np.arange(0,int(max_y),int(max_y)/10.0))+[max_y]
	for site in lat_ordered_sites:
		corplots.append( bokeh_corfig.scatter(x='x',y='y',color=STRUCT[site]['color'],alpha=0.5,source=cor_sources[site]) )
	
	bokeh_corfig.line(x=linerange,y=linerange,color='black')

	N_corplots = range(len(corplots))

	checkbox = CheckboxGroup(labels=[Tnames[stations.index(sitename)] for sitename in lat_ordered_sites],active=range(len(lat_ordered_sites)),width=200)

	####
	checkbox_iterable = [('p'+str(i),plots[i]) for i in N_plots]+[('pcor'+str(i),corplots[i]) for i in N_corplots]+[('checkbox',checkbox)]
	checkbox_code = ''.join(['p'+str(i)+'.visible = checkbox.active.includes('+str(i/2)+');p'+str(i+1)+'.visible = checkbox.active.includes('+str(i/2)+');pcor'+str(i/2)+'.visible = checkbox.active.includes('+str(i/2)+');' for i in range(0,len(N_plots),2)])
	checkbox.callback = CustomJS(args={key: value for key,value in checkbox_iterable}, code=checkbox_code)

	# button to uncheck all checkboxes
	clear_button = Button(label='Clear all',width=120)
	clear_button_code = """checkbox.active=[];"""+checkbox_code
	clear_button.callback = CustomJS(args={key: value for key,value in checkbox_iterable}, code=clear_button_code)

	# button to check all checkboxes
	check_button = Button(label='Check all',width=120)
	check_button_code = """checkbox.active="""+str(N_plots)+""";"""+checkbox_code
	check_button.callback = CustomJS(args={key: value for key,value in checkbox_iterable}, code=check_button_code)


	download_button = Button(label='Save Table to CSV', width = 200)
	download_button.callback = CustomJS(args=dict(dt=data_table),code="""
	var tab = dt.get('source').get('data');
	var filetext = 'site,latitude,N,RMS,Bias,Scatter,R'+String.fromCharCode(10);
	for (i=0; i < tab['site'].length; i++) {
	    var currRow = [tab['site'][i].toString(),
	                   tab['latitude'][i].toString(),
	                   tab['N'][i].toString(),
	                   tab['RMS'][i].toString(),
	                   tab['Bias'][i].toString(),
	                   tab['Scatter'][i].toString(),
	                   tab['R'][i].toString()+String.fromCharCode(10)];

	    var joined = currRow.join();
	    filetext = filetext.concat(joined);
	}

	var filename = 'data_result.csv';
	var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });

    var link = document.createElement("a");
    link = document.createElement('a')
    link.href = URL.createObjectURL(blob);
    link.download = filename
    link.target = "_blank";
    link.style.visibility = 'hidden';
    link.dispatchEvent(new MouseEvent('click'))

	""")

	group = widgetbox(checkbox,clear_button,check_button)

	grid = gridplot([[bokeh_fig,group],[download_button],[data_table,bokeh_corfig]],toolbar_location='left')

	tabs.append(Panel(child=grid,title=var))


	#################################
	#####    Loop over sites    #####
	#################################
	for k in range(0,len(CAS)):
		m_path=os.path.join(save_path,Label_station[k+1])
		if not os.path.exists(m_path):
			os.mkdir(m_path)

		print Label_station[k+1]

		#Number of matches.
		N_station[k]=len(CAS[k])
		print '  N is', N_station[k]

		#limits:
		if B=='y':
			COMP=smooth
		else:
			COMP=CAS
		
		sc=[0]*len(CAS[k])
		for i in range(0,len(CAS[k])):
			sc[i]=CAS[k][i]-smooth[k][i]

		di=[0]*len(CAS[k])
		for i in range(0,len(CAS[k])):
			di[i]=COMP[k][i]-TCCON[k][i]
		try:
			scmin,scmax=round(min(sc)-0.5),round(max(sc)+0.5)
			dimin,dimax=round(min(di)-0.5),round(max(di)+0.5)
			gasmin,gasmax=round(min([min(COMP[k]),min(TCCON[k])])-0.5),round(max([max(COMP[k]),max(TCCON[k])])+0.5)
		except ValueError:
			print '  '+Label_station[k+1]+' has no measurements'

		csv.write(Label_station[k+1]+','+str(N_station[k])+'\n'+''+','+'RMS'+','+'Bias'+','+'Scatter'+','+'R'+'\n')

		###################################################################################################################################################################
		###################################################################################################################################################################
		###########################################################				Custom				###################################################################
		###################################################################################################################################################################
		###################################################################################################################################################################

		if isinstance(all_check,basestring)==True and sc!=[]:

			custom_station=k
			custom_tags=[]
			for i in TCCON_time[custom_station]:
				if i>=calendar.timegm(custom_select) and i<calendar.timegm(custom_next):
					custom_tags.append(TCCON_time[custom_station].index(i))

			custom_priori=[priori[custom_station][i] for i in custom_tags]
			custom_TCCON=[TCCON[custom_station][i] for i in custom_tags]
			custom_CAS=[CAS[custom_station][i] for i in custom_tags]
			custom_smooth=[smooth[custom_station][i] for i in custom_tags]
			custom_TCCON_time=[TCCON_time[custom_station][i] for i in custom_tags]
			custom_DMF_time=[DMF_time[custom_station][i] for i in custom_tags]

			if B=="y":
				COMPA=custom_smooth
			else:
				COMPA=custom_CAS

			#######################
			#####    STATS    #####
			#######################

			#Bias and Root Mean Square
			Bdif=[0]*len(custom_CAS)
			for i in range(0,len(custom_CAS)):
				Bdif[i]=COMPA[i] - custom_TCCON[i]
			Bdif2=[i**2 for i in Bdif]
			Bias=sum(Bdif)/len(custom_CAS)
			RMS=sqrt(sum(Bdif2)/len(custom_CAS))
			
			#Scatter
			difScatter=[0]*len(custom_CAS)
			for i in range(0,len(custom_CAS)):
				difScatter[i]=((COMPA[i] - custom_TCCON[i])-Bias)**2
			Scatter=sqrt((1.0/(len(custom_CAS)-1))*sum(difScatter))

			#Correlation
			meanCAS=sum(COMPA)/len(custom_CAS)
			meanTCCON=sum(custom_TCCON)/len(custom_TCCON)
			T1=[0]*len(custom_CAS)
			T2=[0]*len(custom_CAS)
			T3=[0]*len(custom_CAS)
			for i in range(0,len(custom_CAS)):
				T1[i]=(COMPA[i]-meanCAS)*(custom_TCCON[i]-meanTCCON)
				T2[i]=(COMPA[i]-meanCAS)**2
				T3[i]=(custom_TCCON[i]-meanTCCON)**2
			R=sum(T1)/sqrt(sum(T2)*sum(T3))

			csv.write(str(RMS)+','+str(Bias)+','+str(Scatter)+','+str(R)+'\n')


			###################################
			#####    Correlation plots    #####
			###################################
			max1=int(max(COMPA)+0.5)
			min1=int(min(COMPA)-0.5)
			max2=int(max(custom_TCCON)+0.5)
			min2=int(min(custom_TCCON)-0.5)
			if max(COMPA)>max(custom_TCCON):
				maxi=int(max(COMPA)+1)
			else:
				maxi=int(max(custom_TCCON)+1)

			f=range(maxi)
			linfit=np.polyfit(COMPA,custom_TCCON,1)
			f1=[linfit[0]*j+linfit[1] for j in f]

			fig, ax=plt.subplots()
			fig.set_size_inches(6,4)
			plt.plot(f,f,c='red')
			plt.plot(f,f1,c='green')
			plt.hexbin(COMPA,custom_TCCON, bins='log',gridsize=(int(max1/hexsize),int(max2/hexsize)),cmap=plt.cm.Blues)

			cb=plt.colorbar(extend='max',spacing='uniform')
			cb.set_label('log$_{10}$(Density per hex)')

			plt.legend(['y=x','y='+str(round(linfit[0],2))+'x +'+str(round(linfit[1],2))+'\n $R^2$='+str(round(pearsonr(COMPA,custom_TCCON)[0]*pearsonr(custom_smooth,custom_TCCON)[0],3))],fontsize=12, loc='upper left')
			plt.ylabel('TCCON '+gas,fontsize=12)
			plt.xlabel('GEM-MACH-GHG '+gas,fontsize=12)	
			plt.tick_params(axis='both',which='major',labelsize=12)
			plt.xlim(min1,max1)
			plt.ylim(min2,max2)
			plt.title(Label_station[custom_station+1]+', N='+str(len(custom_TCCON))+',  RMS='+str(round(RMS,2))+', b='+str(round(Bias,2))+', s='+str(round(Scatter,2))+', R='+str(round(R,2)), fontsize=12)
			#plt.subplots_adjust(bottom=0.15,left=0.18)
			saveName=os.path.join(save_path,Label_station[custom_station+1]+var+'_Correlation.png')
			fig.savefig(saveName,bbox_inches='tight')
			plt.clf()
			plt.close()

			if Label_station[custom_station+1]!='Lau':
				print '  ALL +', Label_station[custom_station+1]
				ALL_TCCON+=custom_TCCON
				ALL_COMPA+=COMPA
				Bdev+=[Bias]

			print "  RMS=",RMS
			print "  Bias=",Bias
			print "  Scatter=",Scatter
			print "  R=",R


	##########################
	#####   ALL STATS    #####
	##########################
	if isinstance(all_check,basestring)==True:
			#Bias and Root Mean Square
			Bdif=[0]*len(ALL_TCCON)
			for i in range(0,len(ALL_TCCON)):
				Bdif[i]=ALL_COMPA[i] - ALL_TCCON[i]
			Bdif2=[i**2 for i in Bdif]
			Bias=sum(Bdif)/len(ALL_TCCON)
			RMS=sqrt(sum(Bdif2)/len(ALL_TCCON))
			
			#Scatter
			difScatter=[0]*len(ALL_TCCON)
			for i in range(0,len(ALL_TCCON)):
				difScatter[i]=((ALL_COMPA[i] - ALL_TCCON[i])-Bias)**2
			Scatter=sqrt((1.0/(len(ALL_TCCON)-1))*sum(difScatter))

			#Correlation
			meanCAS=sum(ALL_COMPA)/len(ALL_TCCON)
			meanTCCON=sum(ALL_TCCON)/len(ALL_TCCON)
			T1=[0]*len(ALL_TCCON)
			T2=[0]*len(ALL_TCCON)
			T3=[0]*len(ALL_TCCON)
			for i in range(0,len(ALL_TCCON)):
				T1[i]=(ALL_COMPA[i]-meanCAS)*(ALL_TCCON[i]-meanTCCON)
				T2[i]=(ALL_COMPA[i]-meanCAS)**2
				T3[i]=(ALL_TCCON[i]-meanTCCON)**2
			R=sum(T1)/sqrt(sum(T2)*sum(T3))

			meanBdev=sum(Bdev)/len(Bdev)
			Bscatter=0
			for i in range(0,len(Bdev)):
				Bscatter+=(Bdev[i]-meanBdev)**2
			Bscatter=sqrt(Bscatter/(len(Bdev)-1))


			print 'Overall RMS is:',RMS
			print 'Overall Bias is:',Bias
			print 'Overall Scatter is:',Scatter
			print 'Overall R is:',R
			print 'Bias deviation is:',Bscatter,'\n'
	csv.close()

final=Tabs(tabs=tabs)

milestone = time.time()
print 'Plotting series.html ...'
outfile=open(os.path.join(save_path,'series.html'),'w')
outfile.write(file_html(final,CDN,'GEM-MACH-GHG vs TCCON'))
outfile.close()
print 'series.html DONE in', time.time()-milestone,'seconds'


		


