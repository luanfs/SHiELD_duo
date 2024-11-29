import numpy as np
import matplotlib.pyplot as plt
import os.path

#graphdir='/scratch/cimes/ls9640/graphs_solo_sw/'
#datadir='/scratch/cimes/ls9640/solo_sw/'
graphdir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/graphs_solo_sw/'
datadir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/solo_sw/'
figformat='eps'

#--------------------------------------------------------------------------------------------------------
# test case
tc = -3

# value of N
#N = 192
#N = 384
N = 768
#N = 48
#N = 96

# advection scheme
hords = (0,8)
#hords = (8,)

#grid type 0-equiedge; 2-equiangular
gtypes = (0, 2)

#duogrid scheme
#dgs = (0, 1, 1,)
#dgs = (1, 2, 2)
dgs = (2, 2, 2, 2)
#dgs = (1, 1, 1, 1)

#2d adv scheme
#advs = (1, 1, 2, 2)
advs = (1, 2)

# mass fixers
#mfs = (0, 1, 0, 1)
mfs = (1, 1, 1, 1)

#divergence damping (only for sw)
if tc == 2:
   dds = (0.12,)
   #dds = (0, 0.12,)
   #dds = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
   #dds = (0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12)
else:
   dds = (0,)

#--------------------------------------------------------------------------------------------------------

# basename used in outputs
if tc==1 :
    basename='cosine-zonal'
    alpha = 45  # rotation angle
    Tf = 12

elif tc==2 :
    basename='geobalance'
    alpha = 45
    #alpha = 0
    Tf = 30

elif tc==-3:
    basename='gaussian-zonal'
    alpha = 45
    Tf = 12

elif tc==-4:
    basename='geobalance'
    alpha = 45
    Tf = 12

else:
   print('ERROR: invalid initial condition')

#--------------------------------------------------------------------------------------------------------
# Error lists
errors_linf = []
errors_l1   = []
errors_l2   = []

#--------------------------------------------------------------------------------------------------------
# Loop over all schemes - to get errors
M = len(advs)
schemes_label = []
for g in range(0, len(gtypes)):
 gtype = gtypes[g]
 for k in range(0, len(hords)):
   hord = hords[k]
   for d in range(0, len(dds)):
     dd = str(dds[d])
     for m in range(0, M):
       dg = dgs[m]
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
 
       if hord==0:
          hord_name='UNLIM'
       elif hord==8:
          hord_name='MONO'

       # adv name
       if adv==1:
         advname = 'PL07'
         #advname = 'FV3'
       elif adv==2:
         advname = 'LT2'
         #advname = 'MOD'
 
       #print(adv, advname)
       #------------------------------------------------------------------------------------------------
       # Directory where the netcdf files are
       if tc>1:
         filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
           +".g"+str(gtype)+"."+dg+".adv"+str(adv)+".hord"+str(hord)+'.dd'+str(dd)+".mf"+str(mf)+".tf"+str(Tf)+"/rundir/"
       else:
         filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
       +".g"+str(gtype)+"."+dg+".adv"+str(adv)+".hord"+str(hord)+".mf"+str(mf)+".tf"+str(Tf)+"/rundir/"

       schemes_label.append(advname+"."+hord_name+".mf"+str(mf))
       #schemes_label.append(advname+"."+hord_name)
       #schemes_label.append(advname+"."+hord_name+".dd"+str(dd))
       print("g"+str(gtype)+'.'+advname+"."+hord_name+".mf"+str(mf))
       #print("g"+str(gtype)+'.'+advname+"."+hord_name+".dd"+str(dd))
       file = filepath+"error_delp.txt"
       errors = np.loadtxt(file)

       # get errors
       errors_linf.append(errors[:,0])
       errors_l1.append(errors[:,1])
       errors_l2.append(errors[:,2])
#------------------------------------------------------------------------------------------------

errors = [errors_linf, errors_l1, errors_l2]
names = [r'$L_{\infty}$',r'$L_1$',r'$L_2$']
enames = [r'$L_{\infty}$',r'$L_1$',r'$L_2$']
etitle = ['linf','l1','l2']
colors = ('orange', 'blue', 'orange', 'blue')
#colors = ('lightgreen','lightblue','darkgreen', 'darkblue', )
#colors = ('darkgreen', 'darkblue', )
lines_style = ('-','--')

for l in range(0, len(errors)):
  error = errors[l]
  emin, emax = np.amin(error), np.amax(error)
  emax = 10**(np.floor(np.log10(emax)+1))
  emin = 10**(np.floor(np.log10(emin)-1))
  #emax = 1*10**(-4)
  #emin = 1*10**(-6)
  #fig, axs = plt.subplots(1, 2, figsize=(12, 6))  # Creating subplots, 1 row, 2 columns
  #title = names[l] + " error - TC" + str(tc) +", N = "+str(N)
  #fig.suptitle(title)

  #ax = axs[g]

  c = 0
  for g in range(0, len(gtypes)):
    ax = plt.gca()
    gtype = gtypes[g]
    # grid name
    if gtype==0:
        gname = 'equiedge'
    elif gtype==2:
        gname = 'equiangular'

    for k in range(0, len(hords)):
      hord = str(hords[k])
      color_count = 0
      for d in range(0, len(dds)):
        dd = str(dds[d])
        for m in range(0,M):
         e = error[c]
         Nsteps = np.shape(e)[0]
         time = np.linspace(0,Tf,Nsteps+1)[1:]
         ax.semilogy(time[0:Nsteps:4], e[0:Nsteps:4], lines_style[k], color=colors[color_count], label = schemes_label[c])
         c = c+1
         color_count = color_count+1

    ax.set_xlabel('Time (days)', fontsize=14)
    ax.set_ylabel('Error', fontsize=14)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    ax.set_ylim(emin, emax)
    ax.legend(fontsize=12)
    ax.grid(True, which="both")
    plt.tight_layout()
    #ax.set_title(names[l]+' error evolution for $\phi$ on the '+gname+' grid', fontsize=14)
    ax.set_title(names[l]+r' error evolution for $\phi$ - $N=768$ (13km)', fontsize=14)
    plt.savefig(graphdir+'tc'+str(tc)+'_C'+str(N)+'_'+etitle[l]+'_errors_'+gname+'.'+figformat, format=figformat)
    plt.close()

