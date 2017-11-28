#This function is called, so it won't do plot---plot only when called alone.i\
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import time
ifplot = 'n'
#from IPython import embed;embed()
#plt.figure()
def ts(wrfout,obs,ii):
#    print "ts in dataout is called... %s"%ii
    #print '------------------'+str(wrfout.shape)+'\n',wrfout
    #print '******************'+str(obs.shape)+'\n',obs
    ts = abs(wrfout - obs)
    
    #TS.ncl part
    tot_Hit_ERS = 0.
    tot_Hit_CRS = 0.
    tot_Hit_ERSCRS = 0.
    tot_Hit_HR = 0.
    tot_Hit_MR = 0.
    tot_Hit_LR = 0.
    
    tot_denominator_ERS = 0.
    tot_denominator_CRS = 0.
    tot_denominator_ERSCRS = 0.
    tot_denominator_HR = 0.
    tot_denominator_MR = 0.
    tot_denominator_LR = 0.
    
     # obs = obs(it,:,:) - obs(it-4,:,:)
     # wrfout = wrfout(it,:,:) - wrfout(it-4,:,:)
    sizea,sizeb = obs.shape
    #print "reading shape...",sizea,sizeb
    denominator = float(sizea)*float(sizeb)

    ERS = 0.
    Hit_ERS = 0.
    denominator_ERS = 0.
    CRS = 0.
    Hit_CRS = 0.
    denominator_CRS = 0.
    HR = 0.
    Hit_HR = 0.
    denominator_HR = 0.
    MR = 0.
    Hit_MR = 0.
    denominator_MR = 0.
    LR = 0.
    Hit_LR = 0.
    denominator_LR = 0.
    for i in range(0,sizea):
        for j in range(0,sizeb):

            if(obs[i,j] > 250.):
                ERS = ERS+1
                if(wrfout[i,j] > 250.):
                    Hit_ERS = Hit_ERS+1
            if(obs[i,j] > 250. or wrfout[i,j] > 250.):
                denominator_ERS = denominator_ERS+1


            #if(obs[i,j] > 50. and obs[i,j] < 250.):
	    if(obs[i,j] > 50.):
                CRS = CRS+1
                #if(wrfout[i,j] > 50. and wrfout[i,j] < 250.):
		if(wrfout[i,j] > 50.):
                    Hit_CRS = Hit_CRS+1
            #if(obs[i,j] > 50. and obs[i,j] < 250. or wrfout[i,j] > 50. and wrfout[i,j] < 250.):
	    if(obs[i,j] > 50. or wrfout[i,j] > 50.):
                denominator_CRS = denominator_CRS+1
#---------------strom above
            #if(obs[i,j] > 25. and obs[i,j] < 50.):
	    if(obs[i,j] > 25.):
                HR = HR+1
                #if(wrfout[i,j] > 25. and wrfout[i,j] < 50.):
		if(wrfout[i,j] > 25.):
                    Hit_HR = Hit_HR+1
            #if(obs[i,j] > 25. and obs[i,j] < 50. or wrfout[i,j] > 25. and wrfout[i,j] < 50.):
	    if(obs[i,j] > 25. or wrfout[i,j] > 25.):
                denominator_HR = denominator_HR+1
#-----------------heavy rain above
            #if(obs[i,j] > 10. and obs[i,j] < 25.):
	    if(obs[i,j] > 10.):
                MR = MR+1
                #if(wrfout[i,j] > 10. and wrfout[i,j] < 25.):
		if(wrfout[i,j] > 10.):
                    Hit_MR = Hit_MR+1
            #if(obs[i,j] > 10. and obs[i,j] < 25. or wrfout[i,j] > 10. and wrfout[i,j] < 25.):
	    if(obs[i,j] > 10. or wrfout[i,j] > 10.):
                denominator_MR = denominator_MR+1
#-------moderate rain above
            if(obs[i,j] < 10.):
                LR = LR+1
                if(wrfout[i,j] < 10.):
                    Hit_LR = Hit_LR+1
            if(obs[i,j] < 10. or wrfout[i,j] < 10.):
                denominator_LR = denominator_LR+1

	    colored_TS = 'n'
	    if colored_TS == 'y':

                if(wrfout[i,j] > 50.):
			wrfout[i,j] = 999
                if(obs[i,j] > 50.):
			obs[i,j] = 999
    
                if(obs[i,j] > 25. and obs[i,j] < 50.):
                	obs[i,j] = 666
                if(wrfout[i,j] > 25. and wrfout[i,j] < 50.):
	    	        wrfout[i,j] = 666

                if(obs[i,j] > 10. and obs[i,j] < 25.):
			obs[i,j] = 333    
		if(wrfout[i,j] > 10. and wrfout[i,j] < 25.):
			wrfout[i,j] = 333

		if(obs[i,j] < 10.):
			obs[i,j] = 111
		if(wrfout[i,j] < 10.):
			wrfout[i,j] = 111
            #print("HIT_ERS"+Hit_ERS+"HIT_CRS"+Hit_CRS+"sizea*sizeb"+sizea*sizeb+"sizea*sizeb-denominator_ERS-denominator_CRS"+(sizea*sizeb-denominator_ERS-denominator_CRS))
            #print(denominator_LR+" "+denominator_MR+" "+denominator_HR+" "+denominator_CRS+" "+denominator_ERS)


#    ERS = ERS+0.00001
#    CRS = CRS+0.00001
#    HR = HR+0.00001
#    MR = MR+0.00001
#    LR = LR+0.00001

   	    tot_Hit_ERS = tot_Hit_ERS+Hit_ERS
	    tot_Hit_CRS = tot_Hit_CRS+Hit_CRS
	    tot_Hit_ERSCRS = Hit_ERS+Hit_CRS+tot_Hit_ERSCRS
            tot_Hit_HR = tot_Hit_HR+Hit_HR
	    tot_Hit_MR = tot_Hit_MR+Hit_MR
	    tot_Hit_LR = tot_Hit_LR+Hit_LR

  	    tot_denominator_ERS = tot_denominator_ERS+denominator_ERS
	    tot_denominator_CRS = tot_denominator_CRS+denominator_CRS
 	    tot_denominator_ERSCRS = tot_denominator_ERSCRS+denominator_CRS+denominator_ERS
   	    tot_denominator_HR = tot_denominator_HR+denominator_HR
	    tot_denominator_MR = tot_denominator_MR+denominator_MR
   	    tot_denominator_LR = tot_denominator_LR+denominator_LR
           # print("NOWtotCRSERS"+tot_Hit_ERSCRS)
           
           # print(" ")

    objective_Rainstorm_TS = (Hit_ERS+Hit_CRS)/(denominator_ERS+denominator_CRS)#tot_Hit_ERSCRS/tot_denominator_ERSCRS
           # print("TOT_HIT_HR/TOT"+tot_Hit_HR/tot_denominator_HR)
    objective_HeavyRain_TS = Hit_HR/denominator_HR#tot_Hit_HR/tot_denominator_HR
           # print("TOT_HIT_MR/TOT"+tot_Hit_MR/tot_denominator_MR)
    objective_ModerateRain_TS = Hit_MR/denominator_MR#tot_Hit_MR/tot_denominator_MR
           # print("TOT_HIT_LR/TOT"+tot_Hit_LR/tot_denominator_LR)
    objective_LightRain_TS = Hit_LR/denominator_LR#tot_Hit_LR/tot_denominator_LR

    grids_evaluated = denominator_CRS+denominator_HR+denominator_MR+denominator_LR
    weight_rs = denominator_CRS/grids_evaluated
    weight_hr = denominator_HR/grids_evaluated
    weight_mr = denominator_MR/grids_evaluated
    weight_lr = denominator_LR/grids_evaluated
    #print "Day %s Rainstorm TS:"%ii,objective_Rainstorm_TS,"HeavyRain TS:",objective_HeavyRain_TS,"ModerateRain TS:",objective_ModerateRain_TS,"LightRain TS:",objective_LightRain_TS
    score = weight_rs*objective_Rainstorm_TS+weight_hr*objective_HeavyRain_TS+weight_mr*objective_ModerateRain_TS+weight_lr*objective_LightRain_TS
    if ifplot == 'n':
        return objective_Rainstorm_TS,objective_HeavyRain_TS,objective_ModerateRain_TS,objective_LightRain_TS,score
##----------PLOT--------------
#
#    norm = matplotlib.colors.Normalize(vmin=0, vmax=100)
#    norm2 = matplotlib.colors.Normalize(vmin=0,vmax=999)
##    gci=plt.imshow(grid_z, extent=extent, origin='lower',cmap=cmap, norm=norm)
#    fig = plt.figure()
#    #plt.text(0.25,0.85,"Visualizing Rainfall 07-01 snapshot")
# #  ax = fig.add_subplot(224)
#  #  ax.set_title('Contour Residual|obs-wrf|')
#  #  gci = ax.contour(ts,20)
#  #  plt.colorbar(gci)
#
#    # Lat and Long before mask:
#    ax = fig.add_subplot(131)
#    ax.set_title('Obs',fontsize=11)
#    ax.set_xticklabels(np.linspace(np.asarray(CLONG).round().min(),np.asarray(CLONG).round().max(),11),fontsize=11)
#    ax.set_yticklabels(np.linspace(np.asarray(CLAT).round().max(),np.asarray(CLAT).round().min(),11),fontsize=11)
#    ax.set_xlabel("Longitude(E)",fontsize=11)
#    ax.set_ylabel("Latitude(N)",fontsize=11)
#    gci = ax.imshow(np.flipud(obs),norm=norm)
##    plt.colorbar(gci)
#
#
#    ax = fig.add_subplot(132)
#    ax.set_title('Wrfout',fontsize=11)
#    ax.set_xticklabels(np.linspace(np.asarray(CLONG).round().min(),np.asarray(CLONG).round().max(),11),fontsize=11)
#    ax.set_yticklabels(np.linspace(np.asarray(CLAT).round().max(),np.asarray(CLAT).round().min(),11),fontsize=11)
#    ax.set_xlabel("Longitude(E)",fontsize=11)
#    ax.set_ylabel("Latitude(N)",fontsize=11)
#    gci = ax.imshow(np.flipud(wrfout),norm=norm)
# #   plt.colorbar(gci)
#
#    ax = fig.add_subplot(133)
#    ax.set_title('Residual(|obs-wrf|)',fontsize=11)
#    ax.set_xticklabels(np.linspace(np.asarray(CLONG).round().min(),np.asarray(CLONG).round().max(),11),fontsize=11)
#    ax.set_yticklabels(np.linspace(np.asarray(CLAT).round().max(),np.asarray(CLAT).round().min(),11),fontsize=11)
#    ax.set_xlabel("Longitude(E)",fontsize=11)
#    ax.set_ylabel("Latitude(N)",fontsize=11)
#    gci = ax.imshow(np.flipud(ts),norm=norm)
#
#    cbar_ax=fig.add_axes([0.1, 0.15, 0.8, 0.05])
#    fig.colorbar(gci,orientation="horizontal",cax=cbar_ax)
#    time_string=os.popen("date +20%y_%m_%d_%H_%M").read().replace('\n','')
#    fig.show()
#    input()
#   # plt.savefig("./Visualization_WRF_111.png",figsize=(800,600),dpi=900)
#    return 0,0,0,0
