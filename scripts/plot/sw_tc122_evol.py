import numpy as np
import matplotlib.pyplot as plt
import os.path

#graphdir='/scratch/cimes/ls9640/graphs_solo_sw/'
#datadir='/scratch/cimes/ls9640/solo_sw/'
graphdir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/graphs_solo_sw/'
datadir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/solo_sw/'
figformat='png'

#--------------------------------------------------------------------------------------------------------
# test case
tc = 122

# value of N
#N = 192
#N = 384
#N = 768
N = 48
N = 96

# advection scheme
hords = (5,8)
#hords = (8,)
hords = (5,)
hord = hords[0]
#hords = (5,)

#grid type 0-equiedge; 2-equiangular
gtypes = (0, 2)
gtypes = (0, )
gtype = gtypes[0]

#duogrid scheme
#dgs = (0, 1, 1,)
#dgs = (1, 2, 2)
#dgs = (2, 2, 2, 2)
dgs = (1, 1, 1, 1)

#2d adv scheme
#advs = (1, 1, 2, 2)
advs = (1, 2)
#dds = (0.12,0.12)
dds = (0.15,0.15)
vds = (0.06,0.06)
#vds = (0,0.06)
#vds = (0,0)


# mean depths
mds = (1000, 100, 25)
mds = (100, 25, 10, 5)
mds = (25, 10, 5, 1, 0.5, 0.1, 0.05, 0.01)

#--------------------------------------------------------------------------------------------------------

# basename used in outputs
basename='thinlayer'
alpha = 45
alpha = 0
Tf = 100
 
#--------------------------------------------------------------------------------------------------------
# Error lists
errors_linf = []
errors_l2   = []
times = []

#--------------------------------------------------------------------------------------------------------
# Loop over all schemes - to get errors
M = len(advs)
schemes_label = []
for md in mds:
  for m in range(0, M):
     dg = dgs[m]
     adv = advs[m]
     dd = dds[m]
     vd = vds[m]

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
  
     if hord==5:
        hord_name='UNLIM'
     elif hord==8:
        hord_name='MONO'

     # adv name
     if adv==1:
       advname = 'PL07'
     elif adv==2:
       advname = 'LT2'
       #advname = 'MOD'
  
     #print(adv, advname)
     #------------------------------------------------------------------------------------------------
     # Directory where the netcdf files are
     filepath = datadir+"C"+str(N)+".sw."+basename+".alpha"+str(alpha)\
     +".g"+str(gtype)+"."+dg+".adv"+str(adv)+".hord"+str(hord)+'.dd'+str(dd)+'.vd'+str(vd)+".tf"+str(Tf)+".meandepth"+str(md)+"/rundir/"
     print(filepath)
     #filepath = datadir+"C"+str(N)+".sw."+basename+".tc"+str(tc)+".alpha"+str(alpha)\
     #+".g"+str(gtype)+"."+dg+".adv"+str(adv)+".hord"+str(hord)+'.dd'+dd+'.vd'+vd\
  

     #schemes_label.append(advname+"."+hord_name)
     schemes_label.append(''+str(md)+'m')
     #schemes_label.append(advname+"."+hord_name+".dd"+str(dd))
     print(filepath)
     #print("g"+str(gtype)+'.'+advname+"."+hord_name+".dd"+str(dd))
     file = filepath+"error_delp.txt"
     errors = np.loadtxt(file)

     # Replace NaN values with 100000
     errors[np.isnan(errors)] = 100000

     # Replace -Inf values with 100000
     errors[np.isinf(errors) & (errors < 0)] = 100000

     # get errors
     times.append(errors[:,0])
     errors_linf.append(errors[:,1])
     errors_l2.append(errors[:,2])
#exit()
#------------------------------------------------------------------------------------------------
errors = [errors_linf, errors_l2]
errors = [errors_linf, ]
names  = [r'$L_{\infty}$',r'$L_2$']
enames = [r'$L_{\infty}$',r'$L_2$']
etitle = ['linf', 'l2']
colors = ('orange', 'blue', 'orange', 'blue')
colors = ('orange', 'blue', 'red', 'magenta', 'green', 'brown')

# Create a colormap by interpolating the colors from your array
cmap = plt.cm.colors.LinearSegmentedColormap.from_list("custom_cmap", colors, N=15)

# Generate the 15 colors
color_variation = [cmap(i / 14) for i in range(15)]
colors = color_variation

#colors = ('lightgreen','lightblue','darkgreen', 'darkblue', )
#colors = ('darkgreen', 'darkblue', )
lines_style = ('-','--')
lines_style = ('-','-')

for l in range(0, len(errors)):
  error = errors[l]
  #emin, emax = np.amin(error), np.amax(error)
  #emax = 10**(np.floor(np.log10(emax)+1))
  #emin = 10**(np.floor(np.log10(emin)-1))
  emax = 1*10**(1)
  emin = 1*10**(-3)
  fig, axs = plt.subplots(1, 2, figsize=(12, 6))  # Creating subplots, 1 row, 2 columns
  title = names[l] + " error - Thin layer shallow water test - C"+str(N) +' - hord='+str(hord) 
  fig.suptitle(title, fontsize=16)

  c = 0
  #ax = plt.gca()
  # grid name
  if gtype==0:
      gname = 'equiedge'
  elif gtype==2:
      gname = 'equiangular'

  for k, md in enumerate(mds):
      hord = str(hords[0])
      for m in range(0,M):
        e = error[c]
        Nsteps = np.shape(e)[0]
        time = times[c]
        axs[m].semilogy(time[0:Nsteps:4], e[0:Nsteps:4], lines_style[m], color=colors[k], label = schemes_label[c])
        c = c+1
        print(etitle[l], schemes_label)

  axs[0].set_xlabel('Time (days)', fontsize=14)
  axs[1].set_xlabel('Time (days)', fontsize=14)
  axs[0].set_ylabel(enames[l]+' error', fontsize=14)

  axs[0].set_title('PL07 - div_damp='+str(dds[0])+', vort_damp='+str(vds[0]), fontsize=14)
  axs[1].set_title('LT2  - div_damp='+str(dds[1])+', vort_damp='+str(vds[1]), fontsize=14)

  axs[0].tick_params(axis='x', labelsize=14)
  axs[1].tick_params(axis='x', labelsize=14)
  axs[0].tick_params(axis='y', labelsize=14)

  axs[0].set_ylim(emin, emax)
  axs[1].set_ylim(emin, emax)

  axs[1].set_yticks([])


  axs[1].legend(fontsize=12)

  axs[0].grid(True, which="both")
  axs[1].grid(True, which="both")

  plt.tight_layout()
  #ax.set_title(names[l]+' error evolution for $\phi$ on the '+gname+' grid', fontsize=14)
  #ax.set_title(names[l]+r' error evolution for $\phi$ - $N=768$ (13km)', fontsize=14)
  plt.savefig(graphdir+'tc'+str(tc)+'_C'+str(N)+'_'+etitle[l]+'_hord'+str(hord)+\
  '_dd'+str(dds[0])+'_'+str(dds[1])+'_vd'+str(vds[0])+'_'+str(vds[1])+'_errors_'+gname+'.'+figformat, format=figformat)
  plt.close()
