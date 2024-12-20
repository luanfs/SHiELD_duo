import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors, colorbar
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import colorsys
#import dask
import os.path
import math
graphdir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/graphs_solo_sw/'
datadir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/solo_sw/'
#graphdir='/scratch/cimes/ls9640/graphs_solo_sw/'
#datadir='/scratch/cimes/ls9640/solo_sw/'
figformat='png'

#--------------------------------------------------------------------------------------------------------
# test case
tc = 122

# 1d advection scheme
hords = (5, 8)
#hords = (8, )

# 2d advection scheme
advs = (1,)
advs = (1,2)
#advs = (1,1,1,1)
#advs = (2,2,2,2)

#grid type 0-equiedge; 2-equiangular
#gtypes = (2,)
#gtypes = (0, )
#gtypes = (2,2,2,2)
#gtypes = (0,0,0,0)
gtypes = (0,2)
#gtypes = (0,2,0,2)
#gtypes = (0,2)

#duogrid
dgs = (2, 2, 2, 2)
#dgs = (1, 1, 1, 1)
#dgs = (0, 0, 0, 0)

#mass fixers
mfs = (1, 1, 1, 1)

#divergence damping (only for sw)
if tc == 2:
   dds = (0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12)
else:
   #dds = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
   dds = (0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12)
#dds = (0, 0.12, 0, 0.12)
#dds = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# vorticity damping (same size as hords)
vds = (0.06,0)

# values of N
Ns = (48, 96, 192, 384, 768)
#Ns = (48, 96, 192, 384)
#Ns = (48, 96, 192, )
#Ns = (48, 96, )
ngrids = len(Ns)
#--------------------------------------------------------------------------------------------------------


# basename used in outputs
if tc==122 :
    basename='thinlayer'
    alpha = 45
    Tf = 1
 
else:
   print('ERROR: invalid initial condition')

# initialize error arrays and counter
M = len(advs)
error_linf_h, error_l1_h, error_l2_h = np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes)))
error_linf_u, error_l1_u, error_l2_u = np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))) 
error_linf_v, error_l1_v, error_l2_v = np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))) 

den_linf_h, den_l1_h, den_l2_h = np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))),np.zeros((len(Ns),M,len(hords),len(gtypes)))
den_linf_u, den_l1_u, den_l2_u = np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))),np.zeros((len(Ns),M,len(hords),len(gtypes)))
den_linf_v, den_l1_v, den_l2_v = np.zeros((len(Ns),M,len(hords),len(gtypes))), np.zeros((len(Ns),M,len(hords),len(gtypes))),np.zeros((len(Ns),M,len(hords),len(gtypes)))

#--------------------------------------------------------------------------------------------------------
# Loop over all schemes - to get errors
for g in range(0, len(gtypes)):
  gtype = gtypes[g]
  for k in range(0, len(hords)):
    hord = hords[k]
    vd = str(vds[k])
    for m in range(0, M):
        dg = dgs[m]
        dd = str(dds[m])
        adv = advs[m]
        mf = str(mfs[m])

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
           advname = 'PL07'
        elif adv==2:
           advname = 'LT2'
        n = 0
        # Now let us loop over all N values and compute the errors
        #------------------------------------------------------------------------------------------------
        for N in Ns:
            # Directory where the netcdf files are
            filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
            +".g"+str(gtype)+"."+dg+".adv"+str(adv)+".hord"+str(hord)+'.dd'+dd+'.vd'+vd\
            +".mf"+mf+".tf"+str(Tf)+"/rundir/"
            print(filepath)
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
                    h_fv3 = data['ps'][:,:].values
                    u_fv3 = data['ucomp'][:,:].values
                    v_fv3 = data['vcomp'][:,:].values

                    # Exact fluid depth (initial condition, in this case)
                    h_ref = data['ps_ic'][:,:].values
                    u_ref = data['ua_ic'][:,:].values
                    v_ref = data['va_ic'][:,:].values
 
                    # compute the errors
                    error_h = h_ref - h_fv3
                    error_u = u_ref - u_fv3
                    error_v = v_ref - v_fv3

                    #--------------------------------------------------------------------------------------
                    # h errors
                    #linf
                    den_linf_temp = np.amax(abs(h_ref))
                    den_linf_h[n,m,k,g] = max(den_linf_h[n,m,k,g], den_linf_temp)  
                    error_linf_temp = np.amax(abs(error_h))
                    error_linf_h[n,m,k,g] = max(error_linf_h[n,m,k,g], error_linf_temp)

                    #l1
                    den_l1_h  [n,m,k,g] = den_l1_h  [n,m,k,g] + np.sum(abs(h_ref)*areas)
                    error_l1_h[n,m,k,g] = error_l1_h[n,m,k,g] + np.sum(abs(error_h)*areas)

                    #l2
                    den_l2_h  [n,m,k,g] = den_l1_h  [n,m,k,g] + np.sum(abs(h_ref**2)*areas)
                    error_l2_h[n,m,k,g] = error_l2_h[n,m,k,g] + np.sum(abs(error_h)**2*areas)
                    #--------------------------------------------------------------------------------------


                    #--------------------------------------------------------------------------------------
                    # u errors
                    #linf
                    den_linf_temp = np.amax(abs(u_ref))
                    den_linf_u[n,m,k,g] = max(den_linf_u[n,m,k,g], den_linf_temp)  
                    error_linf_temp = np.amax(abs(error_u))
                    error_linf_u[n,m,k,g] = max(error_linf_u[n,m,k,g], error_linf_temp)

                    #l1
                    den_l1_u  [n,m,k,g] = den_l1_u  [n,m,k,g] + np.sum(abs(u_ref)*areas)
                    error_l1_u[n,m,k,g] = error_l1_u[n,m,k,g] + np.sum(abs(error_u)*areas)

                    #l2
                    den_l2_u  [n,m,k,g] = den_l1_u  [n,m,k,g] + np.sum(abs(u_ref**2)*areas)
                    error_l2_u[n,m,k,g] = error_l2_u[n,m,k,g] + np.sum(abs(error_u)**2*areas)
                    #--------------------------------------------------------------------------------------

                    #--------------------------------------------------------------------------------------
                    # v errors
                    #linf
                    den_linf_temp = np.amax(abs(v_ref))
                    den_linf_v[n,m,k,g] = max(den_linf_v[n,m,k,g], den_linf_temp)  
                    error_linf_temp = np.amax(abs(error_v))
                    error_linf_v[n,m,k,g] = max(error_linf_v[n,m,k,g], error_linf_temp)

                    #l1
                    den_l1_v  [n,m,k,g] = den_l1_v  [n,m,k,g] + np.sum(abs(v_ref)*areas)
                    error_l1_v[n,m,k,g] = error_l1_v[n,m,k,g] + np.sum(abs(error_v)*areas)

                    #l2
                    den_l2_v  [n,m,k,g] = den_l1_v  [n,m,k,g] + np.sum(abs(v_ref**2)*areas)
                    error_l2_v[n,m,k,g] = error_l2_v[n,m,k,g] + np.sum(abs(error_v)**2*areas)

                    #--------------------------------------------------------------------------------------  


                else: # if the file does not exist, we assum it to be unstable
                     error_linf_h[n,m,k,g] = -10000.0
                     error_l1_h[n,m,k,g] = -10000.0
                     error_l2_h[n,m,k,g] = -10000.0
                     den_linf_h[n,m,k,g] = 1.0
                     den_l1_h[n,m,k,g] = 1.0
                     den_l2_h[n,m,k,g] = 1.0

                     error_linf_u[n,m,k,g] = -10000.0
                     error_l1_u[n,m,k,g] = -10000.0
                     error_l2_u[n,m,k,g] = -10000.0
                     den_linf_u[n,m,k,g] = 1.0
                     den_l1_u[n,m,k,g] = 1.0
                     den_l2_u[n,m,k,g] = 1.0
 
                     error_linf_v[n,m,k,g] = -10000.0
                     error_l1_v[n,m,k,g] = -10000.0
                     error_l2_v[n,m,k,g] = -10000.0
                     den_linf_v[n,m,k,g] = 1.0
                     den_l1_v[n,m,k,g] = 1.0
                     den_l2_v[n,m,k,g] = 1.0

            error_linf_h[n,m,k,g] = error_linf_h[n,m,k,g]/den_linf_h[n,m,k,g]
            error_linf_u[n,m,k,g] = error_linf_u[n,m,k,g]/den_linf_u[n,m,k,g]
            error_linf_v[n,m,k,g] = error_linf_v[n,m,k,g]/den_linf_v[n,m,k,g]

            error_l1_h[n,m,k,g] = error_l1_h[n,m,k,g]/den_l1_h[n,m,k,g]
            error_l1_u[n,m,k,g] = error_l1_u[n,m,k,g]/den_l1_u[n,m,k,g]
            error_l1_v[n,m,k,g] = error_l1_v[n,m,k,g]/den_l1_v[n,m,k,g]

            # take the square root for L_2 norm
            if error_l2_h[n,m,k,g] > 0.0:
                error_l2_h[n,m,k,g] = np.sqrt(error_l2_h[n,m,k,g])/np.sqrt(den_l2_h[n,m,k,g])

            if error_l2_u[n,m,k,g] > 0.0:
                error_l2_u[n,m,k,g] = np.sqrt(error_l2_u[n,m,k,g])/np.sqrt(den_l2_u[n,m,k,g])

            if error_l2_v[n,m,k,g] > 0.0:
                error_l2_v[n,m,k,g] = np.sqrt(error_l2_v[n,m,k,g])/np.sqrt(den_l2_v[n,m,k,g])

            #print(gname, 'hord'+str(hord), advname, N, error_linf[n,m,k,g], error_l1[n,m,k,g], error_l2[n,m,k,g])
            print(gname, 'hord'+str(hord), advname, 'mf', mf, N, error_linf_h[n,m,k,g], error_linf_u[n,m,k,g], error_linf_v[n,m,k,g])

            # update counter
            n = n + 1
        print()
#------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------
# get max errors - for using same plotting scale
emax_h = max(np.amax(error_linf_h),np.amax(error_l1_h), np.amax(error_l2_h))
emax_u = max(np.amax(error_linf_u),np.amax(error_l1_u), np.amax(error_l2_u))
emax_v = max(np.amax(error_linf_v),np.amax(error_l1_v), np.amax(error_l2_v))

emax_h = 1.5*emax_h
emax_u = 1.5*emax_u
emax_v = 1.5*emax_v

# make unstable methods have large error
'''
error_linf_h[error_linf_h<0] = float('NaN')
error_l1_h[error_l1_h<0] = float('NaN')
error_l2_h[error_l2_h<0] = float('NaN')

error_linf_u[error_linf_u<0] = float('NaN')
error_l1_u[error_l1_u<0] = float('NaN')
error_l2_u[error_l2_u<0] = float('NaN')

error_linf_v[error_linf_v<0] = float('NaN')
error_l1_v[error_l1_v<0] = float('NaN')
error_l2_v[error_l2_v<0] = float('NaN')
'''

# get min errors
emin_h = min(np.amin(error_linf_h),np.amin(error_l1_h), np.amin(error_l2_h))
emin_u = min(np.amin(error_linf_u),np.amin(error_l1_u), np.amin(error_l2_u))
emin_v = min(np.amin(error_linf_v),np.amin(error_l1_v), np.amin(error_l2_v))

emin_h =  0.7*emin_h
emin_u =  0.7*emin_u
emin_v =  0.7*emin_v
#print(emin, emax)


#------------------------------------------------------------------------------------------------
names = [r'$L_{\infty}$',r'$L_1$',r'$L_2$']
enames = ['linf','l1','l2']
#colors = ('lightgreen', 'darkgreen', 'lightblue', 'darkblue', 'orange', 'black', 'brown', 'gray')
colors = ('orange', 'blue', 'orange', 'blue')
#colors = ('green', 'red', 'blue', 'red', 'blue')
markers = ('*','+','x','*', '+', 'x', '*', '+')
lines_style = ('-','--')

#-------------------------------------------------------------------------------------------------
# Plot ONLY error graph
dpi=100

# Now let us plot the error graph
errors_h = [error_linf_h, error_l1_h, error_l2_h]
errors_u = [error_linf_u, error_l1_u, error_l2_u]
errors_v = [error_linf_v, error_l1_v, error_l2_v]

errors_list = [errors_h, errors_v, errors_u]
fields_name = ('h','u','v')

for f in range(0,len(fields_name)):
 errors = errors_list[f]
 for l in range(0, len(errors)):
  fig, axs = plt.subplots(1, 2, figsize=(12, 6))  # Creating subplots, 1 row, 2 columns
  emin, emax = np.min(errors[l][:,:,:,:]), np.amax(errors[l][:,:,:,:])
  emin, emax = 0.5*emin, 1.5*emax
  for g in range(0, len(gtypes)):
    gtype = gtypes[g]

    # grid name
    if gtype==0:
      gname = 'equiedge'
    elif gtype==2:
      gname = 'equiangular'
    title = fields_name[f]+' - '+names[l] + " error - TC" + str(tc)
    fig.suptitle(title)

    for k in range(0, len(hords)):
        hord = str(hords[k])
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
             advname = 'PL'
           elif adv==2:
             advname = 'LT'

           # subtitle
           subtitle = gname+' - '+dg+'-'+advname+'.hord'+hord+'.mf'+mf

           # Get the maximum error and plot
           #plt.ylim(emin, emax)
           error = errors[l][:,m,k,g]

           # convergence rate
           n = len(Ns)-1
           CR = (np.log(error[n-1])-np.log(error[n]))/np.log(2.0)
           CR = str("{:2.1f}".format(CR))
           axs[g].loglog(Ns, error, lines_style[k], color=colors[m], marker=markers[m], \
           label = advname+'.hord'+hord+".mf"+mf+" - order "+CR)

    # plot reference lines
    eref = 10*np.amin(error)
    order1 = [eref, eref/2.0]
    order2 = [eref, eref/4.0]
    #order3 = [eref, eref/8.0]
    if g==1:
       axs[g].loglog(Ns[ngrids - 2:ngrids], order1, '-',  color='black', label='1st order')
       axs[g].loglog(Ns[ngrids - 2:ngrids], order2, '--', color='black', label='2nd order')
       #axs[g].loglog(Ns[ngrids - 2:ngrids], order3, '-.', color='black', label='3rd order')
    else:
       axs[g].loglog(Ns[ngrids - 2:ngrids], order1, '-',  color='black')
       axs[g].loglog(Ns[ngrids - 2:ngrids], order2, '--', color='black')
       #axs[g].loglog(Ns[ngrids - 2:ngrids], order3, '-.', color='black')
 

    # Set a common title
    #if tc==1 or tc==2:
    #    title = names[l]+" error - TC"+str(tc)+", $\\alpha$ = "+str(alpha)+', '+str(Tf)+' days'
    #else:
    #    title = names[l]+" error - TC"+str(tc)+', '+str(Tf)+' days'

    # Label
    axs[g].set_xlabel('$N$')
    axs[g].set_ylabel('Error')
    axs[g].set_xlim(0, 1000)
    axs[g].set_ylim(emin,emax)
    axs[g].set_title(gname)
    #if g==1:
    axs[g].legend()
    axs[g].grid(True, which="both")
  filename = graphdir+fields_name[f]+'_'+enames[l]+"error_tc"+str(tc)+'_alpha'+str(alpha)
  plt.savefig(filename+'.'+figformat, format=figformat)
  plt.close()
#------------------------------------------------------------------------------------------------
