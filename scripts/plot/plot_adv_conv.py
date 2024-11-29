import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors, colorbar
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import colorsys
import os.path
import math

graphdir='/scratch/cimes/ls9640/graphs_solo_sw/'
datadir='/scratch/cimes/ls9640/solo_sw/'
graphdir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/graphs_solo_sw/'
datadir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/solo_sw/'
figformat='eps'

#--------------------------------------------------------------------------------------------------------
# test case
tc = -5

# 1d advection scheme
hords = (0,8)
#hords = (0, )

#grid type 0-equiedge; 2-equiangular
#gtypes = (2,)
#gtypes = (0,)
gtypes = (0,2)

#duogrid
dgs = (2,2,2,2,2,2,2,2)

# 2d advection scheme
advs = (1,2)
#advs = (1,1,2,2)

# mass fixer
mfs = (1,1,1,1)
#mfs = (0,1,0,1)

# values of N
Ns = (48, 96, 192, 384, 768)
#Ns = (48, 96, 192, 384)
#Ns = (48, 96, 192, )
#Ns = (96, 192, 384, 768)
#Ns = (48, 96, )
ngrids = len(Ns)
#--------------------------------------------------------------------------------------------------------


# basename used in outputs
if tc==1:
    basename='cosine-zonal'
    alpha = 45  # rotation angle
    Tf = 12

elif tc==-3:
    basename='gaussian-zonal'
    alpha = 45  # rotation angle
    Tf = 12

elif tc==-4:
    basename='geobalance'
    alpha = 45  # rotation angle
    Tf = 12

elif tc==-5 :
    basename='twogaussians-ndiv'
    alpha = 0 # rotation angle
    Tf = 12
    hmin, hmax = 0.0, 1.0
    
elif tc==-6 :
    basename='twogaussians-div'
    alpha = 0 # rotation angle
    Tf = 12
    hmin, hmax = 0.0, 5.5

else:
   print('ERROR: invalid initial condition')

# initialize error arrays and counter
M = len(advs)
error_linf, error_l1, error_l2 = np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))) 
den_linf, den_l1, den_l2 = np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))),np.zeros((len(Ns),M,len(hords),len(gtypes)))

#--------------------------------------------------------------------------------------------------------
# Loop over all schemes - to get errors
for g in range(0, len(gtypes)):
  gtype = gtypes[g]
  for k in range(0, len(hords)):
    hord = hords[k]
    if hord==0:
       hord_name='UNLIM'
    elif hord==8:
       hord_name='MONO'
   
    for m in range(0, M):
        dg = dgs[m]
        adv = advs[m]
        mf = mfs[m]
        
        # grid name
        if gtype==0:
            gname = 'equiedge'
        elif gtype==2:
            gname = 'equiangular'

        # duogrid name
        if dg==1:
            dg = 'dg1'
        elif dg==2:
            dg = 'dg2'
        else:
            dg = 'kinked'

        # adv name
        if adv==1:
           advname = 'PL'
        elif adv==2:
           advname = 'LT'
        n = 0
        # Now let us loop over all N values and compute the errors
        #------------------------------------------------------------------------------------------------
        for N in Ns:
            # Directory where the netcdf files are
            filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
            +".g"+str(gtype)+"."+dg+".adv"+str(adv)+".hord"+str(hord)+".mf"+str(mf)+".tf"+str(Tf)+"/rundir/"
            #print(filepath)
            # Loop over tiles
            for tile in range(1,7):
                # Files to be opened
                atmos_file = filepath+"atmos_daily.tile"+str(tile)+".nc"
                grid_file  = filepath+"grid_spec.tile"+str(tile)+".nc"

                # Check if they exist
                files_exist = os.path.exists(grid_file) and os.path.exists(atmos_file)

                if files_exist:
                    # Load the files
                    data = xr.open_dataset(atmos_file, decode_times=False)
                    grid = xr.open_dataset(grid_file , decode_times=False)
                    areas = grid.area[:,:].values

                    # Final timestep
                    t = int(len(data.time.values)-1)
             
                    # Variable to be plotted (ps = fluid depth in the SW model)
                    h_fv3 = data['ps'][t,:,:].values

                    # Exact fluid depth (initial condition, in this case)
                    h_ref = data['ps_ic'][:,:].values

                    # compute the errors
                    error = h_ref - h_fv3

                    #linf
                    den_linf_temp = np.amax(abs(h_ref))
                    den_linf[n,m,k,g] = max(den_linf[n,m,k,g], den_linf_temp)  
                    error_linf_temp = np.amax(abs(error))
                    error_linf[n,m,k,g] = max(error_linf[n,m,k,g], error_linf_temp)

                    #l1
                    den_l1  [n,m,k,g] = den_l1  [n,m,k,g] + np.sum(abs(h_ref)*areas)
                    error_l1[n,m,k,g] = error_l1[n,m,k,g] + np.sum(abs(error)*areas)

                    #l2
                    den_l2  [n,m,k,g] = den_l1  [n,m,k,g] + np.sum(abs(h_ref**2)*areas)
                    error_l2[n,m,k,g] = error_l2[n,m,k,g] + np.sum(abs(error)**2*areas)

                else: # if the file does not exist, we assum it to be unstable
                     error_linf[n,m,k,g] = -10000.0
                     error_l1[n,m,k,g] = -10000.0
                     error_l2[n,m,k,g] = -10000.0


            error_linf[n,m,k,g] = error_linf[n,m,k,g]/den_linf[n,m,k,g]
            error_l1  [n,m,k,g] = error_l1  [n,m,k,g]/den_l1  [n,m,k,g]
            # take the square root for L_2 norm
            if error_l2[n,m,k,g] > 0.0:
                error_l2[n,m,k,g] = np.sqrt(error_l2[n,m,k,g])/np.sqrt(den_l2[n,m,k,g])


            #print(gname, 'hord'+str(hord), advname, N, error_linf[n,m,k,g], error_l1[n,m,k,g], error_l2[n,m,k,g])
            print(gname, 'hord '+str(hord_name), advname, 'mf', mf, N, error_linf[n,m,k,g])

            # update counter
            n = n + 1
        print()
#------------------------------------------------------------------------------------------------

#exit()


#------------------------------------------------------------------------------------------------
# get max errors - for using same plotting scale
emax = max(np.amax(error_linf),np.amax(error_l1), np.amax(error_l2))
emax = 1.5*emax

# make unstable methods have large error
error_linf[error_linf<0] = float('NaN')
error_l1[error_l1<0] = float('NaN')
error_l2[error_l2<0] = float('NaN')

# get min errors
emin = min(np.amin(error_linf),np.amin(error_l1), np.amin(error_l2))
emin =  0.7*emin
#print(emin, emax)


#------------------------------------------------------------------------------------------------
# Now let us plot the error graph
errors = [error_linf, error_l1, error_l2]
names = [r'$L_{\infty}$',r'$L_1$',r'$L_2$']
enames = ['linf','l1','l2']
colors = ('lightgreen', 'darkgreen', 'lightblue', 'darkblue', 'orange', 'black', 'brown', 'gray')
colors = ('orange', 'blue', 'blue', 'red', 'green', 'red', 'blue')
#colors = ('green', 'red', 'blue', 'red', 'blue')
markers = ('*','+','x','*', '+', 'x', '*', '+')
lines_style = ('-','--')

#-------------------------------------------------------------------------------------------------
# Plot ONLY error graph
dpi=100

#for k in range(0, M):
for l in range(0, len(errors)):
  #fig, axs = plt.subplots(1, 2, figsize=(12, 6))  # Creating subplots, 1 row, 2 columns
  emin, emax = np.min(errors[l][:,:,:,:]), np.amax(errors[l][:,:,:,:])
  emin, emax = 0.5*emin, 1.5*emax
  for g in range(0, len(gtypes)):
    # if gtypes[g]==0:
    gtype = gtypes[g]

    # grid name
    if gtype==0:
      gname = 'equiedge'
    elif gtype==2:
      gname = 'equiangular'
    #title = names[l] + " error - TC" + str(tc)
    #plt.suptitle(title)
    ax = plt.gca()
    #ax = axs[g]

    for k in range(0, len(hords)):
      hord = str(hords[k])
      if hord=='0':
         hord_name='UNLIM'
      elif hord=='8':
         hord_name='MONO'
  
      # Loop over all schemes
      for m in range(0,M):
         dg = dgs[m]
         adv = advs[m]
         mf = str(mfs[m])

         # duogrid name
         if dg==1:
           dg = 'dg1'
         elif dg==2:
           dg = 'dg2'
         else:
           dg = 'kinked'

         # adv name
         if adv==1:
           advname = 'PL07'
           #advname = 'FV3'
         elif adv==2:
           advname = 'LT2'
           #advname = 'MOD'
         
         # subtitle
         subtitle = gname+' - '+dg+'-'+advname+'.hord'+hord+'.mf'+mf
         subtitle = gname+' - '+dg+'-'+advname+'.hord'+hord
         # Get the maximum error and plot
         #plt.ylim(emin, emax)
         error = errors[l][:,m,k,g]

         # convergence rate
         n = len(Ns)-1
         CR = (np.log(error[n-1])-np.log(error[n]))/np.log(2.0)
         CR = str("{:2.1f}".format(CR))

         #axs[g].loglog(Ns, error, lines_style[k], color=colors[m], marker=markers[m], \
         ax.loglog(Ns, error, lines_style[k], color=colors[m], marker=markers[m], \
         label = advname+'.'+hord_name+" - order "+CR)
         #label = advname+'.'+hord_name+".mf"+mf+" - order "+CR)

    # plot reference lines
    eref = 10*np.amin(error)
    order1 = [eref, eref/2.0]
    order2 = [eref, eref/4.0]
    order3 = [eref, eref/8.0]
    if gtype==0:
      ax.loglog(Ns[ngrids - 2:ngrids], order1, '-',  color='black', label='1st order')
      ax.loglog(Ns[ngrids - 2:ngrids], order2, '--', color='black', label='2nd order')
      ax.loglog(Ns[ngrids - 2:ngrids], order3, '-.', color='black', label='3rd order')
    else:
      ax.loglog(Ns[ngrids - 2:ngrids], order1, '-',  color='black')
      ax.loglog(Ns[ngrids - 2:ngrids], order2, '--', color='black')
      ax.loglog(Ns[ngrids - 2:ngrids], order3, '-.', color='black')
 
    ax.set_xlabel('$N$', fontsize=14)
    ax.set_ylabel('Error',fontsize=14)
    ax.set_xlim(0, 1000)
    ax.set_ylim(emin,emax)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=14)
    #ax.set_title(names[l] + " error - "+gname, fontsize=14)
    ax.set_title(names[l] + r" relative error for $\phi$", fontsize=14)

    ax.legend(fontsize=12)
    ax.grid(True, which="both")
    plt.tight_layout()
    filename = graphdir+enames[l]+"error_tc"+str(tc)+'_alpha'+str(alpha)
    #plt.savefig(filename+'.'+figformat, format=figformat)
    plt.savefig(filename+'.'+gname+'.'+figformat, format=figformat)
    plt.close()
#------------------------------------------------------------------------------------------------
