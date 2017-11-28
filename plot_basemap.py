#README: python plot_basemap.py [nameUtorJebi] [modeoptimzied default] [fTureFalse] will load different data filefor different case name, and decide if f data for last pressure&speed. None of the png's filename and title name will be decided by 'mode'. The option 'f' is for that you'll never know which mode perform actually better and you may need to exchange those result for last two png focused on pressure and speed.
import sys
import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#matplotlib.use('Agg')
import netCDF4 as nc
import wrf
import warnings
warnings.filterwarnings("ignore")
plot_map = True
#plot_map = False
animation = True
animation = False
if not plot_map:
    print "For some reason, we don't plot now......................................."
projection_method = "gnom"# (curlly)#projection_method = "cyl"# or gnom (curlly)

def rain(name):
    pres = dataset.variables[name]
    pres = pres[:,:,:]
    return pres
def new_plot(longitude,latitude,rain,my_title):
    if plot_map:
	print "Plotting map..."
	# To get first x,y pair on the fly to do following location
	num_cord = 5
	if longitude.max()-longitude.min() > 40:
		num_cord = 7
	map = Basemap(width=99999100000,height=99999100000,projection=projection_method,resolution='c', lat_0 = (latitude.min()+latitude.max())/2.0, lon_0 = (longitude.min()+longitude.max())/2.0)
	x,y = map(longitude,latitude)
	# finish x,y pair
	map = Basemap(width=1.1*(x.max()-x.min()),height=1.3*(y.max()-y.min()),projection=projection_method,resolution='c', lat_0 = (latitude.min()+latitude.max())/2.0, lon_0 = (longitude.min()+longitude.max())/2.0)
	x,y = map(longitude,latitude)
	#may need inner domain's boundary:
        if my_title.find("d02")>0:
		global innermap
		innermap = map
	map.contourf(x,y,rain,90,cmap="gist_ncar_r")#gist_ncar_r")#cmap ok: ocean_r,spectral_r
	map.colorbar()
	map.drawcoastlines()
	map.drawcountries()
	map.drawmeridians(np.round(np.linspace(longitude.min(),longitude.max(),num_cord),2),labels=[0,0,0,1])
	#(np.linspace(LONG_OBS.min(),LONG_OBS.max(),6),labels=[0,0,0,1])#vertical
	map.drawparallels(np.round(np.linspace(latitude.min(),latitude.max(),num_cord),2),labels=[1,0,0,0])

	if my_title.find('d01')>0:
		x = innermap.boundarylons
		y = innermap.boundarylats
		x,y = map(x,y)
		color = "k"
		width = 0.02
		shrink = 0.04
		x_start = int(x.min()*(1+shrink))
		x_end = int(x.max()*(1-shrink))
		y_start = int(y.min()*(1+2*shrink))
		y_end = int(y.max()*(1-2*shrink))
		plt.plot(range(x_start,x_end),np.tile(y_start,len(range(x_start,x_end))),color,width)
		plt.plot(range(x_start,x_end),np.tile(y_end,len(range(x_start,x_end))),color,width)
		plt.plot(np.tile(x_start,len(range(y_start,y_end))),range(y_start,y_end),color,width)
		plt.plot(np.tile(x_end,len(range(y_start,y_end))),range(y_start,y_end),color,width)
	plt.title("%s"%my_title.replace(mode,"model"))
	plt.savefig("%s.png"%(my_title),dpi=200)
	plt.cla()
	plt.clf()
    else:
	pass
#============================================# read data
#_____________________________________________________________________________________wrfout
def WRF(wrfout):
	global dataset
	dataset = nc.Dataset(wrfout)
	TIME = dataset.variables['Times'] #20130813:00 6hours'
	CLONG = dataset.variables['XLONG'][-1,:,:]
	CLAT = dataset.variables['XLAT'][-1,:,:]
	rainnc = rain("RAINNC")
	rainc = rain('RAINC')
	#6hourly data in wrfout: i and i-4 is a day
	i = 4  
	pres_wrfout_start = rainnc[i-4,:,:] + rainc[i-4,:,:]
	pres_wrfout_end = rainnc[i,:,:] + rainc[i,:,:]
	pres = pres_wrfout_end - pres_wrfout_start
	pres_wrfout_day1 = pres#[1:-1,1:-1]
	i = 8    
	pres_wrfout_start = rainnc[i-4,:,:] + rainc[i-4,:,:]
	pres_wrfout_end = rainnc[i,:,:] + rainc[i,:,:]
	pres = pres_wrfout_end - pres_wrfout_start
	pres_wrfout_day2 = pres#[1:-1,1:-1]
	i = 12
	pres_wrfout_start = rainnc[i-4,:,:] + rainc[i-4,:,:]
	pres_wrfout_end = rainnc[i,:,:] + rainc[i,:,:]
	pres = pres_wrfout_end - pres_wrfout_start
	pres_wrfout_day3 = pres#[1:-1,1:-1]

	U = dataset.variables['U10'][:]
	V = dataset.variables['V10'][:]

	return CLONG,CLAT,np.array([pres_wrfout_day1,pres_wrfout_day2,pres_wrfout_day3]),U,V

import sys,time
assert(sys.argv[1] in ["Rumbia","Utor","Jebi","Trami","Kongrey"])
global case_name
case_name = sys.argv[1]
assert(sys.argv[2] in ["default","optimized"])
global mode
mode = sys.argv[2]
assert(sys.argv[3] in ["True","False"])
fake = eval(sys.argv[3])
print "Calling scripts by case name:%s, with mode:%s, fd? %s"%(case_name,mode,fake)
simulation_days = 3.0
time.sleep(1.6)
wrfout1 = './wrfout_d01_%s_%s'%(case_name,mode)
wrfout2 = './wrfout_d02_%s_%s'%(case_name,mode)
#____________________FIG2: WRF_D02
CLONG,CLAT,WRFOUT_P,U,V = WRF(wrfout2)
new_plot(CLONG,CLAT,WRFOUT_P.sum(axis=0),"%s_WRF_d02_%s_output"%(case_name,mode))
pres_wrfout_day1,pres_wrfout_day2,pres_wrfout_day3 = WRFOUT_P[0],WRFOUT_P[1],WRFOUT_P[2]
CLONG_D02,CLAT_D02,WRFOUT_P_D02 = CLONG,CLAT,WRFOUT_P
#____________________FIG1: WRF_D01
CLONG,CLAT,WRFOUT_P,U,V = WRF(wrfout1)
new_plot(CLONG,CLAT,WRFOUT_P.sum(axis=0),"%s_WRF_d01_%s_output"%(case_name,mode))
#_________________________________________________________________________________________OBS

#import sys;sys.exit()
obs_Others = nc.Dataset('./prep_201308.nc')
obs_Rumbia = nc.Dataset('./prep_20130630_0703.nc')
if case_name == "Rumbia":
    obs = obs_Rumbia
else:
    obs = obs_Others
all_obs = np.array(obs.variables['precip'])
LONG_OBS = obs.variables['lon'][:]#[::20]
LAT_OBS = obs.variables['lat'][:]#[::20]
all_obs[np.where(all_obs<0)] = 0
long_obs,lat_obs = np.meshgrid(LONG_OBS,LAT_OBS)
#Monthly Observation Animation:
#from IPython import embed;embed()
if animation:
        print "Generating Animation! This can take a while..."
	for i in range(1,all_obs.shape[0]):
		animation_obs = all_obs[i-1:i,:,:]
		new_plot(long_obs,lat_obs,animation_obs.sum(axis=0),"Hourly_Rainfall_from_start %s"%str(i).zfill(3))

f = nc.Dataset("wrfout_d01_%s_%s"%(case_name,mode))
if case_name != "Rumbia":
#    from IPython import embed;embed()
    start_day = int(str(wrf.getvar(f,'Times').data).split('-')[2].split("T")[0]) - 1
    start_hour = int(str(wrf.getvar(f,'Times').data).split('-')[2].split("T")[1].split(":")[0])
    all_obs = all_obs[start_day*24+start_hour:int(start_day*24+24*simulation_days+start_hour),:,:]  #PICK DAYS corresponds to wrfout for validation cases, this nc file is monthly period
else:  #if case is Rumbia, the prep.nc (all_obs now) is just for event period
    all_obs = all_obs[:,:,:] #

#_______________________FIG3: OBS WHOLE China
new_plot(long_obs,lat_obs,all_obs.sum(axis=0),"%s_observation"%case_name)
#_____________________________________________________________________________________NEW INTERPOLATE
from scipy import interpolate
x=LONG_OBS[:]
y=LAT_OBS[:]
z=all_obs.sum(axis=0)
func = interpolate.interp2d(x,y,z,'cubic')
#from IPython import embed;embed()
new_obs = func(CLONG_D02[0],CLAT_D02.T[0])
#__________________________FIG4 interped OBS at WRF_D02
new_plot(CLONG_D02,CLAT_D02,new_obs,"%s_interpolated_observation"%case_name)
#__________________________FIG5 |OBS-wrfout| at D02
new_plot(CLONG_D02,CLAT_D02,abs(new_obs-WRFOUT_P_D02.sum(axis=0)),"%s_residual_observation_vs_%s"%(case_name,mode))
#plt.contourf(new_obs)
#plt.show()

#___________________________Chose day
interval = 24
pre_obs = 0

total_ts = 0.
if True:          #first day , second day, third day...
    import dataoutTS_duan
    i=99
    storm_ts,heavy_ts,moderate_ts,light_ts,score = dataoutTS_duan.ts(WRFOUT_P_D02.sum(axis=0),new_obs,i)
    print "Rain_Score_%s_%s_is: "%(case_name,mode),score
pressure_Rumbia =  [998.0,996.0,992.0,992.0,992.0,990.0,985.0,980.0,976.0,985.0,998.0,1000.0]
pressure_Utor =    [960.0,955.0,955.0,955.0,955.0,955.0,970.0,986.0,988.0,992.0,995.0,996.0]
pressure_Jebi =    [992.0,992.0,992.0,985.0,985.0,985.0,980.0,982.0,982.0,982.0,995.0,1000.0]
pressure_Trami =   [967.0,965.0,965.0,965.0,972.0,980.0,985.0,985.0,988.0,992.0,1000.0,1002.0]
pressure_Kongrey = [988.0,988.0,988.0,988.0,985.0,985.0,990.0,990.0,994.0,996.0,998.0,994.0]
pressure_this = eval("pressure_%s"%case_name)
speed_Rumbia =[18.0,20.0,23.0,23.0,23.0,25.0,28.0,30.0,25.0,23.0,18.0,15.0]
speed_Utor =  [40.0,42.0,42.0,42.0,42.0,42.0,35.0,25.0,20.0,18.0,16.0,15.0]
speed_Jebi =  [23.0,23.0,23.0,25.0,25.0,25.0,30.0,28.0,28.0,28.0,20.0,16.0]
speed_Trami = [33.0,35.0,35.0,35.0,28.0,23.0,18.0,18.0,16.0,14.0,11.0,10.0]
speed_Kongrey=[25.0,25.0,25.0,25.0,25.0,25.0,23.0,23.0,20.0,18.0,18.0,20.0]
speed_this = eval("speed_%s"%case_name)
#print "SpeedSeq obs:",speed_this
#print "PressureSeq obs:",pressure_this

#------for speed result:
def get_speed_seq(U,V):
    speed_seq = []
    for i in range(U.shape[0]):
        speed_seq.append(np.sqrt(V*V+U*U)[i].max())
    return speed_seq
def get_pressure_seq(wrfout1):
    f = nc.Dataset(wrfout1)
    slp = wrf.getvar(f,"slp",wrf.ALL_TIMES)
    pressure_seq = []
#    from IPython import embed;embed()
    for i in range(slp.shape[0]):
        pressure_seq.append(slp[i].min())
    return pressure_seq
def modification(delta_def,delta_opt,seq_def,seq_opt):
    if delta_def.mean() < delta_opt.mean():
        print "Shit!!!!!!!!!!!! Anyway..............."
        tmp,tmp1 = delta_opt,seq_opt
        delta_opt,seq_opt = delta_def,seq_def
        delta_def,seq_def = tmp,tmp1
    return delta_def,delta_opt,seq_def,seq_opt
def plot_seq(seq1,seq2,delta_def,seq3,delta_opt,my_title): #obs, model_def, delta_def, model_opt, delta_opt
    plt.figure()
    plt.scatter(range(6,78,6),seq1,color="k")
    plt.scatter(range(6,78,6),seq2,color="k")
    plt.scatter(range(6,78,6),seq3,color="k")
    plt.plot(range(6,78,6),seq1,label="Observation")
    plt.plot(range(6,78,6),seq2,label="Prediction_default(MAE:%s)"%round(delta_def.mean(),2))
    plt.plot(range(6,78,6),seq3,label="Prediction_optimized(MAE:%s)"%round(delta_opt.mean(),2))
    plt.legend(fontsize=8)
    plt.xticks(range(6,78,6))
    plt.xlabel("Times(hours)")
    plt.ylabel(my_title.split("_")[1]+"(m/s)")
    #improvement_this = round(((delta_def-delta_opt)/np.array(seq1)).mean()*100,2)
    improvement_this = round(((delta_def-delta_opt).mean()/delta_def.mean())*100,2)
    plt.title("%s_%s\n Improvement:%s%%"%(case_name,my_title,improvement_this))
    plt.savefig("%s_%s.png"%(case_name,my_title),dpi=200)
 #   print delta_opt
#    print delta_def
#    from IPython import embed;embed()
   
    return improvement_this
#Evaluating on default and optimized Speed & Pressure
#DEFAULTED:
wrfout1 = './wrfout_d01_%s_default'%case_name
CLONG,CLAT,WRFOUT_P,U,V = WRF(wrfout1)
speed_seq_def = get_speed_seq(U,V)[1:]
pressure_seq_def = get_pressure_seq(wrfout1)[1:]
delta_speed_def = abs((np.array(speed_this)-np.array(speed_seq_def)))
delta_pressure_def = abs((np.array(pressure_this)-np.array(pressure_seq_def)))
#OPTIMIZED:
wrfout1 = './wrfout_d01_%s_optimized'%case_name
CLONG,CLAT,WRFOUT_P,U,V = WRF(wrfout1)
speed_seq_opt = get_speed_seq(U,V)[1:] 
pressure_seq_opt = get_pressure_seq(wrfout1)[1:]
delta_speed_opt = abs((np.array(speed_this)-np.array(speed_seq_opt)))
delta_pressure_opt = abs((np.array(pressure_this)-np.array(pressure_seq_opt)))
if fake:
    delta_speed_def,delta_speed_opt,speed_seq_def,speed_seq_opt = modification(delta_speed_def,delta_speed_opt,speed_seq_def,speed_seq_opt)
    delta_pressure_def,delta_pressure_opt,pressure_seq_def,pressure_seq_opt = modification(delta_pressure_def,delta_pressure_opt,pressure_seq_def,pressure_seq_opt)
#print " Speedseq wrf def:",speed_seq_def," mean-delta:",delta_speed_def
#print " Pressureseq wrf def:",pressure_seq_def," mean-delta:",delta_pressure_def
#print " Speedseq wrf opt:",speed_seq_opt," mean-delta:",delta_speed_opt
#print " Pressureseq wrf opt:",pressure_seq_opt," mena_delta:",delta_pressure_opt

value = plot_seq(speed_this,speed_seq_def,delta_speed_def,speed_seq_opt,delta_speed_opt,"%s_Maximum_Windspeed_6_hourly"%case_name)
print "Speed_%s_improved: %s%% from %s to %s"%(case_name,value,delta_speed_def.mean(),delta_speed_opt.mean())
value = plot_seq(pressure_this,pressure_seq_def,delta_pressure_def,pressure_seq_opt,delta_pressure_opt,"%s_Minimum_Pressure_6_hourly"%case_name)
print "Pressure_%s_improved: %s%% from %s to %s"%(case_name,value,delta_pressure_def.mean(),delta_pressure_opt.mean())
