import numpy as np
import matplotlib.pyplot as plt
import os.path

#graphdir='/scratch/cimes/ls9640/graphs_solo_sw/'
#datadir='/scratch/cimes/ls9640/solo_sw/'
graphdir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/graphs_solo_sw/'
datadir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/solo_sw/'
figformat='png'

#--------------------------------------------------------------------------------------------------------
# value of N
#N = 192
#N = 384
N = 768
#N = 48
#N = 96

# advection scheme
#hords = (5,5,8,8,5)
hords = (5,5,8,8)
#hords = (5,5,5)

#2d adv scheme
#advs = (1,2,1,2,2)
advs = (1,2,1,2)
#advs = (2,1,2)

#grid type 0-equiedge; 2-equiangular
gtypes = (0, 2)
#gtypes = (0,)
#gtypes = (2,)


#vort damp
#vds = (0,0.04,0,0,0)
#vds = (0,0,0.04)
vds = (0,0,0,0,0)

#duogrid scheme
dg = 1

#--------------------------------------------------------------------------------------------------------
basename='geobalance.tc2'
alpha = 45
Tf = 1

#--------------------------------------------------------------------------------------------------------
# Error lists
errors_linf = []
errors_l1   = []
errors_l2   = []

#--------------------------------------------------------------------------------------------------------
# Loop over all schemes - to get errors
schemes_label = []
for gtype in gtypes:
   for adv, hord, vd in zip(advs, hords, vds):
      # grid name
      if gtype==0:
          gname = 'equiedge'
      elif gtype==2:
          gname = 'equiangular'

      # duogrid name
      if dg==1:
         dgname = 'dg1'
      elif dg==2:
         dgname = 'dg2'
      else:
         dgname = 'kinked'
 
      if hord==5:
         hord_name='UNLIM'
         hord_name='hord5'
      elif hord==8:
         hord_name='MONO'
         hord_name='hord8'

      # adv name
      if adv==1:
        advname = 'PL07'
      elif adv==2:
        advname = 'LT2'
 
      #------------------------------------------------------------------------------------------------
      # Directory where the netcdf files are
      if vd ==0:
         filepath = datadir+"C"+str(N)+".sw."+basename+".alpha"+str(alpha)\
         +".g"+str(gtype)+"."+dgname+".adv"+str(adv)+".hord"+str(hord)+".tf"+str(Tf)+"/rundir/"
      else:
         filepath = datadir+"C"+str(N)+".sw."+basename+".alpha"+str(alpha)\
         +".g"+str(gtype)+"."+dgname+".adv"+str(adv)+".hord"+str(hord)+".vd"+str(vd)+".tf"+str(Tf)+"/rundir/"
 
      schemes_label.append(advname+"."+hord_name+".vtdm."+str(vd))
      print(filepath)

      print("g"+str(gtype)+'.'+advname+"."+hord_name)
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
colors = ('orange', 'blue', 'orange', 'blue', 'blue')
#colors = ('lightgreen','lightblue','darkgreen', 'darkblue', )
#colors = ('darkgreen', 'darkblue', )
lines_style = ('solid','solid','dashed','dashed','dotted')
lines_style = ('dashed','solid','solid','dotted')

for l in range(0, len(errors)):
  error = errors[l]
  emin, emax = np.amin(error), np.amax(error)
  #emax = 10**(np.floor(np.log10(emax))+1)
  emax = 1.2*emax
  emin = 10**(np.floor(np.log10(emin)))
  #emax = 1*10**(-4)
  #emin = 1*10**(-6)
  #fig, axs = plt.subplots(1, 2, figsize=(12, 6))  # Creating subplots, 1 row, 2 columns
  #title = names[l] + " error - TC" + str(tc) +", N = "+str(N)
  #fig.suptitle(title)

  #ax = axs[g]

  c = 0
  for gtype in gtypes:
     for adv, hord, Color, scheme_label, e, ls in zip(advs, hords, colors, schemes_label, error, lines_style):
        ax = plt.gca()
        # grid name
        if gtype==0:
           gname = 'equiedge'
        elif gtype==2:
           gname = 'equiangular'

        Nsteps = np.shape(e)[0]
        time = np.linspace(0,Tf,Nsteps+1)[1:]
        gap=10
        ax.semilogy(time[0:Nsteps:gap], e[0:Nsteps:gap], linestyle=ls, color=Color, label=scheme_label)

     ax.set_xlabel('Time (days)', fontsize=14)
     ax.set_ylabel('Error', fontsize=14)
     ax.tick_params(axis='x', labelsize=14)
     ax.tick_params(axis='y', labelsize=14)
     ax.set_ylim(emin, emax)
     ax.legend(fontsize=12)
     ax.grid(True, which="both")
     ax.set_title('Steady geostrophic flow - N='+str(N)+'\n'\
     +names[l]+' error evolution for $h$ on the '+gname+' grid', fontsize=14)
     plt.savefig(graphdir+basename+'_C'+str(N)+'_'+etitle[l]+'_errors_'+gname+'.png', format='png')
     plt.close()
