import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors, colorbar
from matplotlib.colors import ListedColormap
#from plot_cs3 import *

#----------------------------------------------------------------------------------------------
#source https://unidata.github.io/python-gallery/examples/Precipitation_Map.html
cmap_data = [(1.0, 1.0, 1.0),
             (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
             (0.0, 1.0, 1.0),
             (0.0, 0.8784313797950745, 0.501960813999176),
             (0.0, 0.7529411911964417, 0.0),
             (0.501960813999176, 0.8784313797950745, 0.0),
             (1.0, 1.0, 0.0),
             (1.0, 0.6274510025978088, 0.0),
             (1.0, 0.0, 0.0),
             (1.0, 0.125490203499794, 0.501960813999176),
             (0.9411764740943909, 0.250980406999588, 1.0),
             (0.501960813999176, 0.125490203499794, 1.0),
             (0.250980406999588, 0.250980406999588, 1.0),
             (0.125490203499794, 0.125490203499794, 0.501960813999176),
             (0.125490203499794, 0.125490203499794, 0.125490203499794)]
 
# Create the colormap using ListedColormap
cmap_precp = ListedColormap(cmap_data)
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
def plot_scalarfield(qs, subtitles, filename, title, filepath, colormap, qmin, qmax):
    # Set figure DPI and format
    dpi = 100
    figformat = 'png'

    # Create a figure with 1 row and 2 columns (side by side plots)
    fig, ax = plt.subplots(1, 2, figsize=(24, 9), subplot_kw={'projection': ccrs.PlateCarree()})  # Define projection here

    # Loop through each scalar field
    for k, (q, subtitle) in enumerate(zip(qs, subtitles)):
        # Add gridlines with labels
        gl = ax[k].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                             linewidth=1, color='gray', alpha=0.5, linestyle='--')
        gl.xlabel_style = {'size': 19, 'color': 'black'}
        gl.ylabel_style = {'size': 19, 'color': 'black'}

        # Define color for each cube (6 tiles)
        colors = ('black', 'black', 'black', 'black', 'black', 'black')
        N = np.shape(q)[0]

        # Plot each tile (cube section)
        for tile in range(6):
            # Grid file to be loaded for the specific tile
            grid_file = filepath + "grid_spec.tile" + str(tile + 1) + ".nc"

            # Load the grid dataset
            grid = xr.open_dataset(grid_file, decode_times=False)

            # Get longitude and latitude arrays
            lon = grid['grid_lon'] + 180  # Adjust longitude to be in [0, 360]
            lat = grid['grid_lat']

            # Plot the cube edges
            A_lon, A_lat = lon[0, 0], lat[0, 0]
            B_lon, B_lat = lon[N, 0], lat[N, 0]
            C_lon, C_lat = lon[N, N], lat[N, N]
            D_lon, D_lat = lon[0, N], lat[0, N]
            lw = 0.2
            ax[k].plot([A_lon, B_lon], [A_lat, B_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
            ax[k].plot([B_lon, C_lon], [B_lat, C_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
            ax[k].plot([C_lon, D_lon], [C_lat, D_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
            ax[k].plot([D_lon, A_lon], [D_lat, A_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)

            # Plot the scalar field data for each tile
            im = ax[k].pcolormesh(lon, lat, q[:, :, tile], alpha=1, transform=ccrs.PlateCarree(),
                                   zorder=10, vmin=qmin, vmax=qmax, cmap=colormap)

        # Set plot limits
        ax[k].set_xlim(-20, 20)
        ax[k].set_ylim(-10, 30)

        # Set the title for each plot
        ax[k].set_title(subtitle, fontsize=19)

        # Add colorbar
        cax, kw = colorbar.make_axes(ax[k], orientation='horizontal', fraction=0.046, pad=0.04, format='%.1e')
        cb = plt.colorbar(im, cax=cax, extend='both', **kw)

        ticks = np.linspace(qmin, qmax, num=5)
        cb.set_ticks(ticks)
        cb.ax.tick_params(labelsize=22)

    # Add common title for the entire figure
    fig.suptitle(title, fontsize=22)

    # Save the figure only once, after both plots are generated
    plt.savefig(filename + '.' + figformat, format=figformat)
    print('Plotted ' + filename)

    plt.close()
#----------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------
# directories
#hydro_option="nh"
hydro_option="hydro"
#inputdir='/scratch/cimes/ls9640/solo_'+hydro_option+'/'
#outputdir = '/scratch/cimes/ls9640/graphs_solo_'+hydro_option+'/'
inputdir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/solo_'+hydro_option+'/'
outputdir = '/gpfs/f5/scratch/Luan.Santos/gfdl_w/graphs_solo_'+hydro_option+'/'
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
#simulation parameters
N=96
N=192
#N=384
gtype = 0
hords = (5,)
advs  = (1,2) # plot for PL07 and LT2 advection schemes
testname='tc'
dg='dg1'
basename= "C"+str(N)+"."+hydro_option+"."+testname+".g"+str(gtype)+"."+str(dg)
#----------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------
# grid name
if gtype==0:
    gname = 'equiedge'
elif gtype==2:
    gname = 'equiangular'
#-----------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------
# Create filepaths list
filepaths = []
for hord in hords: 
    datas = []
    for adv in advs:
        filepath = inputdir+basename+'.adv'+str(adv)+'.hord'+str(hord)+"/rundir/"
        filepaths.append(filepath)
#-----------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------
# Get the number of plots
atmos_file = filepaths[0]+"atmos_4xdaily.tile1.nc"
data = xr.open_dataset(atmos_file, decode_times=False)
times = data.time.values
nplots = len(times)
dtplot = times[1]-times[0]

#-----------------------------------------------------------------------------------------
# Arrays to store the data that will be plotted
# FV3 data
pwat = np.zeros((N,N,6,nplots,len(advs)))
vort = np.zeros((N,N,6,nplots,len(advs)))
u = np.zeros((N,N,6,nplots,len(advs)))
v = np.zeros((N,N,6,nplots,len(advs)))
delp = np.zeros((N,N,6,nplots,len(advs)))

# open & read files
tgap=1
 
for k, filepath in enumerate(filepaths):
    for tile in range(0,6):
        # Files to be opened
        atmos_file = filepath+"atmos_4xdaily.tile"+str(tile+1)+".nc"
        grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

        # Load the files
        print('loading ' + atmos_file)
        data = xr.open_dataset(atmos_file, decode_times=False)
        grid = xr.open_dataset(grid_file , decode_times=False)
        #for d in data.keys():
        #   print(d)
        #exit()

        for t in range(0,nplots,tgap):
            # Variables to be plotted (ps = fluid depth in the SW model)
            u[:,:,tile,t,k] = data['u850'][t,:,:].values
            v[:,:,tile,t,k] = data['v850'][t,:,:].values
            vort[:,:,tile,t,k] = data['vort850'][t,:,:].values
            pwat[:,:,tile,t,k] = data['PWAT'][t,:,:].values

# plot
time = times[0]
dif_u, dif_v = np.zeros(nplots), np.zeros(nplots)
for t in range(0,nplots,tgap):
    print(t)
    ##############################################################################################
    fref = 0
    field_adv1 = vort[:,:,:,t,0]-fref
    field_adv2 = vort[:,:,:,t,1]-fref

    fmin = min(np.amin(field_adv1), np.amin(field_adv2))
    fmax = max(np.amax(field_adv1), np.amax(field_adv2))
    fabs = max(abs(fmin),abs(fmax))
    fmin, fmax = -fabs, fabs

    field='vort'
    title = field+'_'+testname+'_t'+str(t)+'_'+basename+'.hord'+str(hord)
    filename = outputdir+title
    subtitle1 = 'PL07\n'
    subtitle2 = 'LT2\n'

    plot_scalarfield([field_adv1, field_adv2], [subtitle1, subtitle2], filename, title, \
                     filepaths[0], 'seismic', fmin, fmax)
    ##############################################################################################

    ##############################################################################################
    field_adv1 = pwat[:,:,:,t,0]
    field_adv2 = pwat[:,:,:,t,1]
    fmin, fmax = np.amin(pwat), np.amax(pwat)
    fmin = min(np.amin(field_adv1), np.amin(field_adv2))
    fmax = max(np.amax(field_adv1), np.amax(field_adv2))
    #fmin, fmax = -fabs, fabs

    field='pwat'
    title = field+'_'+testname+'_t'+str(t)+'_'+basename+'.hord'+str(hord)
    filename = outputdir+title
    subtitle1 = 'PL07\n'
    subtitle2 = 'LT2\n'

    plot_scalarfield([field_adv1, field_adv2], [subtitle1, subtitle2], filename, title, \
                     filepaths[0], cmap_precp, fmin, fmax)
    ##############################################################################################
    # time update
    time = time + dtplot 
print('bye')
