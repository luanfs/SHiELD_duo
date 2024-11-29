#----------------------------------------------------------------------------------------------
# Module for plotting routines
#----------------------------------------------------------------------------------------------
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors, colorbar
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import math
import matplotlib as mpl
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
graphdir = '../graphs/'
#----------------------------------------------------------------------------------------------
# Plot the scalar field q given in a A grid.
#----------------------------------------------------------------------------------------------
def plot_scalarfield(q, map_projection, title, filename, filepath, colormap, qmin, qmax,\
    dpi, figformat):
    # map projection
    if map_projection == "mercator":
        #plt.figure(figsize=(1832/dpi, 977/dpi), dpi=dpi)
        # gaussian at corner
        #plt.figure(figsize=(1000/dpi, 1000/dpi), dpi=dpi)
        plt.figure(figsize=(1000/dpi, 1050/dpi), dpi=dpi)
        #plt.figure(figsize=(1600/dpi, 1000/dpi), dpi=dpi)
        plateCr = ccrs.PlateCarree()
    elif map_projection == "sphere":
        #plt.figure(figsize=(800/dpi, 800/dpi), dpi=dpi) 
        plt.figure(figsize=(1000/dpi,1000/dpi), dpi=dpi) 
        plateCr = ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0)
        #plateCr = ccrs.Orthographic(central_longitude=0.25*180, central_latitude=180/6.0)
        #plateCr = ccrs.Orthographic(central_longitude=-45, central_latitude=30)

    projection=ccrs.PlateCarree(central_longitude=0)
    plateCr._threshold = plateCr._threshold/10.

    ax = plt.axes(projection=plateCr)
    #ax.set_global()
    #ax.stock_img()
    #ax.gridlines(draw_labels=True)
    #ax.xlabel_style = {'size': 18}
    #ax.ylabel_style = {'size': 18}
    #gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    #                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
    #gl.xlabel_style = {'size': 19, 'color': 'black'}
    #gl.ylabel_style = {'size': 19, 'color': 'black'}

 
    # Color of each cubed panel
    colors = ('black','black','black','black','black','black')
    N = np.shape(q)[0]
    # plot for each tile
    for tile in range(0,6):
        # Grid file to be opened
        grid_file  = filepath+"/rundir/grid_spec.tile"+str(tile+1)+".nc"
        print('plotting tile ', tile)

        # Load the file
        grid = xr.open_dataset(grid_file , decode_times=False)

        # Get grid
        lon = grid['grid_lon']
        lat = grid['grid_lat']

        # Plot cube edges
        A_lon, A_lat = lon[0,0], lat[0,0]
        B_lon, B_lat = lon[N, 0], lat[N, 0]
        C_lon, C_lat = lon[N, N], lat[N, N]
        D_lon, D_lat = lon[0, N], lat[0, N]
        lw = 0.2
        plt.rcParams["axes.axisbelow"] = True

        ax.plot([A_lon, B_lon], [A_lat, B_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([B_lon, C_lon], [B_lat, C_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([C_lon, D_lon], [C_lat, D_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([D_lon, A_lon], [D_lat, A_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)

        # Plot scalar field
        im = ax.pcolormesh(lon, lat, q[:,:,tile], alpha=1, transform=ccrs.PlateCarree(), \
        zorder=10, vmin = qmin, vmax=qmax,  cmap=colormap)
    # gaussian at corner
    #plt.xlim(45-25,45+25)
    #plt.ylim(10,45+15)
    # cylinder
    #plt.xlim(10,80)
    #plt.ylim(5,65)
    #plt.tick_params(axis='both', labelsize=18)
    plt.title(title,  fontsize=19)
    #ax.coastlines()
    # Plot colorbar


    r0 = np.pi/9.
    deg2rad = np.pi/180
    rad2deg = 1.0/deg2rad
    lon0 = -np.pi/2 + np.pi/4.0
    lat0 = np.pi/6.
    lat0, lon0, r0 = lat0*rad2deg, lon0*rad2deg, r0*rad2deg
    circle1 = plt.Circle((lon0, lat0), r0, color='black', fill=False, zorder=11, linestyle='--',linewidth=0.3)
    #ax.add_patch(circle1)

    #cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.046, pad=0.04, shrink=0.9, format='%.0e')
    # gaussian
    #cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.046, pad=0.04, format='%.0e')
    # cilinder
    #cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.046, pad=0.04, format='%.1e')
    cax,kw = colorbar.make_axes(ax,orientation='vertical' , fraction=0.09, pad=0.04, shrink=0.9, format='%.0e')

    cb=plt.colorbar(im, cax=cax, extend='both',**kw)

    ticks = np.linspace(qmin, qmax, num=5)
    cb.set_ticks(ticks)
    cb.ax.tick_params(labelsize=22)
    plt.savefig(filename+'.'+figformat, format=figformat)
    plt.close()
#-----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
# Plot the scalar field q given in a A grid.
#----------------------------------------------------------------------------------------------
def plot_scalarfield_wind(q, ua, va, map_projection, title, filename, filepath, colormap, qmin, qmax,\
    dpi, figformat):
    # map projection
    if map_projection == "mercator":
        #plt.figure(figsize=(1832/dpi, 977/dpi), dpi=dpi)
        #plt.figure(figsize=(1000/dpi, 1000/dpi), dpi=dpi)
        plt.figure(figsize=(1600/dpi, 1000/dpi), dpi=dpi)
        plateCr = ccrs.PlateCarree()
    elif map_projection == "sphere":
        plt.figure(figsize=(800/dpi, 800/dpi), dpi=dpi) 
        plateCr = ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0)

    projection=ccrs.PlateCarree(central_longitude=0)
    plateCr._threshold = plateCr._threshold/10.

    ax = plt.axes(projection=plateCr)
    ax.set_global()
    #ax.stock_img()
    ax.gridlines(draw_labels=True)

    # The third example illustrates the use of custom length colorbar
    # extensions, used on a colorbar with discrete intervals.
    #colormap = mpl.colors.ListedColormap([[0., .4, 1.], [0., .8, 1.],[1., .8, 0.], [1., .4, 0.]])
    #colormap.set_over((1., 0., 0.))
    #colormap.set_under((0., 0., 1.))

    #bounds = [-1., -.5, 0., .5, 1.]
    #norm = mpl.colors.BoundaryNorm(bounds, colormap.N)

    # Color of each cubed panel
    colors = ('black','black','black','black','black','black')
    #colors = ('gray','gray','gray','gray','gray','gray')
    #colors = ('white','white','white','white','white','white')
    N = np.shape(q)[0]
    # plot for each tile
    for tile in range(0,6):
        # Grid file to be opened
        grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

        # Load the file
        grid = xr.open_dataset(grid_file , decode_times=False, chunks={'time': 5}).squeeze()

        # Get grid
        lon = grid['grid_lon']
        lat = grid['grid_lat']

        # Plot cube edges
        A_lon, A_lat = lon[0,0], lat[0,0]
        B_lon, B_lat = lon[N, 0], lat[N, 0]
        C_lon, C_lat = lon[N, N], lat[N, N]
        D_lon, D_lat = lon[0, N], lat[0, N]
        lw = 0.2
        plt.rcParams["axes.axisbelow"] = True

        ax.plot([A_lon, B_lon], [A_lat, B_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([B_lon, C_lon], [B_lat, C_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([C_lon, D_lon], [C_lat, D_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([D_lon, A_lon], [D_lat, A_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)

        # Plot scalar field
        im = ax.pcolormesh(lon, lat, q[:,:,tile], alpha=1, transform=ccrs.PlateCarree(), \
        zorder=10, vmin = qmin, vmax=qmax,  cmap=colormap)
        #lata = np.zeros((N,N))
        #lat1, lon1 = np.zeros((N+1,N+1)),np.zeros((N+1,N+1))
        #lat1[:,:], lon1[:,:] = lat[:,:], lon[:,:]
        #A_lat, B_lat, C_lat, D_lat = lat1[0:N,0:N], lat1[1:N+1,0:N], lat1[0:N,1:N+1], lat1[1:N+1,1:N+1]
        #lat0 = (A_lat+B_lat+C_lat+D_lat)*0.5
        #A_lon, B_lon, C_lon, D_lon = lon1[0:N,0:N], lon1[1:N+1,0:N], lon1[0:N,1:N+1], lon1[1:N+1,1:N+1]
        #lon0 = (A_lon+B_lon+C_lon+D_lon)*0.5
        #print(np.shape(lat0), np.shape(A_lat), np.shape(B_lat), np.shape(C_lat), np.shape(D_lat))
        #lona = (lon[0:N,0:N] + lon[0:N,1:N+1] + lon[1:N+1,0:N] + lon[1:N+1,1:N+1])*0.25
        #print(N,np.shape(lata), np.shape(lat[0:N+2,0:N+1] + lat[0:N+2,1:N+2]), np.shape(lat[0:N,0:N]), np.shape(lat[0:N,1:N+1]))#, np.shape(lat[1:N+1,0:N]), np.shape(lat[1:N+1,1:N+1]), np.shape(q[:,:,tile]))
       # print(np.shape(lon),np.shape(lat[1:N+1,1:N+1]),np.shape(lata),np.shape(q[:,:,tile]))
        #exit()
        #plt.contour(lon0, lat0, q[:,:,tile], 4, colors='k')
        gap = N//8
        plt.quiver(lon[0:N:gap, 0:N:gap], lat[0:N:gap, 0:N:gap], ua[0:N:gap, 0:N:gap, tile], va[0:N:gap, 0:N:gap, tile], zorder=12, width=0.0005, color='white')
    """
    r0 = np.pi/9.
    deg2rad = np.pi/180
    rad2deg = 1.0/deg2rad
    lon0 = -np.pi/2 + np.pi/4.0
    lat0 = np.pi/6.
    lat0, lon0, r0 = lat0*rad2deg, lon0*rad2deg, r0*rad2deg
    circle1 = plt.Circle((lon0, lat0), r0, color='black', fill=False, zorder=13, linestyle='--',linewidth=0.3)
    ax.add_patch(circle1)
    """
    #plt.xlim(45-33,45+33)
    #plt.ylim(5,45+20)
 
    plt.title(title,  fontsize=17)
    #cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.046, pad=0.04, shrink=0.9, format='%.1e')
    #cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.046, pad=0.04, format='%.1e')
    cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.09, pad=0.04, shrink=0.9, format='%.0e')
    cb=plt.colorbar(im, cax=cax, extend='both',**kw)
    cb.ax.tick_params(labelsize=17)
    plt.savefig(filename+'.'+figformat, format=figformat)
    plt.close()
#-----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
# Plot the scalar field q given in a A grid in a panel for all time steps
#----------------------------------------------------------------------------------------------
def plot_scalarfield_panel(q, map_projection, filename, filepath, title, colormap, qmin, qmax,\
    nplots, dt, dpi, figformat):

    # panel dimensions
    ncols = 2
    nlines = math.ceil((nplots+1) / ncols)
    fig = plt.figure(figsize=(nlines*500/dpi, ncols*4000/dpi),dpi=dpi)
    gs = fig.add_gridspec(nlines, ncols, width_ratios=[1] * ncols)
    k = 0
    time = 0
    for i in range(0, nlines):
        for j in range(0, ncols):
            if k<=nplots:
                ax = fig.add_subplot(gs[i, j], projection=ccrs.PlateCarree())  # Use PlateCarree projectioni
                #ax.set_global()
                plateCr = ccrs.PlateCarree()        
                projection=ccrs.PlateCarree(central_longitude=0)
                plateCr._threshold = plateCr._threshold/10.
                #ax.stock_img()
                ax.gridlines(draw_labels=True)
                # plot for each tile
                for tile in range(0,6):
                    # Grid file to be opened
                    grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

                    # Load the file
                    grid = xr.open_dataset(grid_file , decode_times=False, chunks={'time': 5}).squeeze()

                    # Get grid
                    lon = grid['grid_lon']
                    lat = grid['grid_lat']

                    # Plot scalar field
                    im = ax.pcolormesh(lon, lat, q[:,:,tile,k],alpha=1,transform=ccrs.PlateCarree(),\
                    zorder=10, vmin = qmin, vmax=qmax,  cmap=colormap)

                dmin = str("{:.1e}".format(np.amin(abs(q[:,:,:,k]))))
                dmax = str("{:.1e}".format(np.amax(abs(q[:,:,:,k]))))
                Time = str("{:.1e}".format(time))
                if np.amin(abs(q[:,:,:,k]))>0.0000000000001:
                    subtitle =str(Time)+" days, min = "+ dmin+", max = "+dmax
                else:
                    subtitle =str(Time)+" days, max = "+dmax

                # Add coastlines to the subfigure
                #ax.add_feature(ccrs.COASTLINE, edgecolor='k', linewidth=0.8)
                #ax.add_feature(cfeature.COASTLINE, edgecolor='k', linewidth=0.8)
     
                # Add a subtitle to the subfigure
                ax.set_title(subtitle)

                # time update
                time = time + dt

            #else:
            #    print('bye!')
            k = k+1

    #plt.tight_layout()  # Ensure proper spacing between subplotsp

    # Add a single colorbar
    #cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.5])  # Define the position and size of the vertical colorbar
    #cbar = fig.colorbar(im, cax=cbar_ax, fraction=0.046, pad=0.04, shrink=0.8, format='%.1e')
    #cbar.ax.tick_params(labelsize=8)

    #fig.suptitle(title+'\n\n', fontsize=30)

    plt.savefig(filename+'.'+figformat, format=figformat)
    #plt.show()
    plt.close()
    #exit()


#-----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------
# Plot the scalar field q given in a A grid.
#----------------------------------------------------------------------------------------------
def plot_scalarfield2(q, lons, lats, map_projection, title, filename, colormap, dpi, figformat, qmin, qmax):
    # map projection
    if map_projection == "mercator":
        plt.figure(figsize=(1600/dpi, 1000/dpi), dpi=dpi)
        plateCr = ccrs.PlateCarree()        
    elif map_projection == "sphere":
        plt.figure(figsize=(800/dpi, 800/dpi), dpi=dpi) 
        plateCr = ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0)
       # plateCr = ccrs.Orthographic(central_longitude=0.25*180, central_latitude=180/6.0)

    projection=ccrs.PlateCarree(central_longitude=0)
    plateCr._threshold = plateCr._threshold/10.

    ax = plt.axes(projection=plateCr)
    #ax.coastlines()
    ax.set_global()
    ax.stock_img()
    #ax.gridlines(draw_labels=True)
    # Color of each cubed panel
    colors = ('black','black','black','black','black','black')


    # The third example illustrates the use of custom length colorbar
    # extensions, used on a colorbar with discrete intervals.
    #colormap = mpl.colors.ListedColormap([[0., .4, 1.], [0., .8, 1.], [1., .8, 0.], [1., .6, 0.],  [1., .3, 0.]])
    #colormap.set_over((1., 0., 0.))
    #colormap.set_under((0., 0., 1.))

    #bounds = [-1., -.5, 0., .5, 1.]
    #norm = mpl.colors.BoundaryNorm(bounds, colormap.N)
    N = np.shape(q)[0]
    # plot for each tile
    for tile in range(0,6):
        # Get grid
        lon = lons[:,:,tile]
        lat = lats[:,:,tile]

        # Plot cube edges
        A_lon, A_lat = lon[0,0], lat[0,0]
        B_lon, B_lat = lon[N, 0], lat[N, 0]
        C_lon, C_lat = lon[N, N], lat[N, N]
        D_lon, D_lat = lon[0, N], lat[0, N]
        lw = 0.2
        plt.rcParams["axes.axisbelow"] = True

        #ax.plot([A_lon, B_lon], [A_lat, B_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        #ax.plot([B_lon, C_lon], [B_lat, C_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        #ax.plot([C_lon, D_lon], [C_lat, D_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        #ax.plot([D_lon, A_lon], [D_lat, A_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)

        # Plot scalar field
        im = ax.pcolormesh(lon, lat, q[:,:,tile], alpha=1, transform=ccrs.PlateCarree(), \
        zorder=10, vmin = qmin, vmax=qmax,  cmap=colormap)
 

    plt.title(title,  fontsize=17)
    # Plot colorbar
    #cax,kw = colorbar.make_axes(ax,orientation='vertical' , fraction=0.046, pad=0.04, shrink=0.8, format='%.1e')
    cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.09, pad=0.04, shrink=0.9, format='%.0e')
    cb=plt.colorbar(im, cax=cax, extend='both',**kw)
    #cb.ax.tick_params(labelsize=17)
    plt.savefig(graphdir+filename+'.'+figformat, format=figformat)
    plt.close()
#-----------------------------------------------------------------------------------------
