#----------------------------------------------------------------------------------------------
# Module for reference solution computation
#----------------------------------------------------------------------------------------------
import xarray as xr
import numpy as np


####################################################################################
# Convert from spherical (lat,lon) to cartesian coordinates (x,y,z)
# on the unit sphere.
# Inputs: latitude (lat), longitude (lon)
####################################################################################
def sph2cart(lon, lat):
    x = np.cos(lat)*np.cos(lon)
    y = np.cos(lat)*np.sin(lon)
    z = np.sin(lat)
    return x, y, z

####################################################################################
# Convert from cartesian (x,y,z) to spherical coordinates (lat,lon)
# on the unit sphere.
# Outputs: latitude (lat), longitude (lon)
####################################################################################
def cart2sph(X, Y, Z):
    hypotxy = np.hypot(X, Y)
    lat = np.arctan2(Z, hypotxy)
    lon = np.arctan2(Y, X)
    return lon, lat

####################################################################################
# Reference resolution for h
####################################################################################
def hexact(lon, lat, t, tc, alpha):
    alpha_rad = alpha*np.pi/180.0 #convert from degree to radians
    #alpha_rad = alpha

    # Gaussian hill or cosine bell
    if tc==-3 or tc==1:
        # Wind speed
        u0 =  2.0*np.pi/(12.0*86400.0) # Wind speed
        ws = -u0
        wt = ws*t

        #Rotation parameters
        cosa  = np.cos(alpha_rad)
        cos2a = cosa*cosa
        sina  = np.sin(alpha_rad)
        sin2a = sina*sina
        coswt = np.cos(wt)
        sinwt = np.sin(wt)

        X, Y, Z = sph2cart(lon, lat)

        rotX = (coswt*cos2a+sin2a)*X -sinwt*cosa*Y + (coswt*cosa*sina-cosa*sina)*Z
        rotY =  sinwt*cosa*X + coswt*Y + sina*sinwt*Z
        rotZ = (coswt*sina*cosa-sina*cosa)*X -sinwt*sina*Y + (coswt*sin2a+cos2a)*Z

        # Gaussian center
        if tc==-3:
            c = 1.0/np.sqrt(3.0)
            X0, Y0, Z0 = c, c, c
            #lon0, lat0 = np.pi/4.0, np.pi/6.0
            #X0, Y0, Z0 = sph2cart(lon0, lat0)
            h = 0.1 + 0.9*np.exp(-10*((rotX-X0)**2+ (rotY-Y0)**2 + (rotZ-Z0)**2))

        # Cosine center
        elif tc==1:
            lon0, lat0 =  np.pi*0.25, np.pi*0.5 - np.arccos(1.0/np.sqrt(3.0))
            lon , lat  = cart2sph(rotX, rotY, rotZ)
            N = np.shape(lon)[0]
            erad = 6371200
            r0 = erad/3.0
            r = 2.0*erad*np.arcsin(np.sqrt(np.sin((lat-lat0)/2)**2 + np.cos(lat)*np.cos(lat0)*np.sin((lon-lon0)/2.)**2 ))
            h = 0.5*np.ones((N,N))
            mask = r<=r0
            h[mask] = 0.5 + 0.5*(1.0+np.cos(np.pi*r[mask]/r0))
    return h

####################################################################################
# Reference resolution for h on A grid using B grid points
####################################################################################
def href_Agrid(grid, t, dt, tc, alpha):
    # dt is given in days
    day2sec = 86400
    time = t*dt*day2sec

    # Get B grid
    rad2deg = 180.0/np.pi              # Radians to degrees conversion
    deg2rad = 1.0/rad2deg  # Degrees to radians conversion
    lonB = (grid['grid_lon'].values)*deg2rad
    latB = grid['grid_lat'].values*deg2rad
    XB, YB, ZB = sph2cart(lonB, latB)

    # Get A grid points
    N = np.shape(lonB)[0]
    XA = XB[0:N-1,0:N-1] + XB[1:N,0:N-1] + XB[0:N-1,1:N] + XB[1:N,1:N]
    YA = YB[0:N-1,0:N-1] + YB[1:N,0:N-1] + YB[0:N-1,1:N] + YB[1:N,1:N]
    ZA = ZB[0:N-1,0:N-1] + ZB[1:N,0:N-1] + ZB[0:N-1,1:N] + ZB[1:N,1:N]
    # project onto sphere
    normA = np.sqrt(XA*XA + YA*YA + ZA*ZA)
    XA, YA, ZA = XA/normA, YA/normA, ZA/normA
    # get latlon
    lonA, latA = cart2sph(XA, YA, ZA)

    # Gaussian hill
    if tc==-3 or tc==1:
        href = hexact(lonA, latA, time, tc, alpha)
    return href

def topography(grid, tc):
    if tc == 5:
        r0 = np.pi/9.
        deg2rad = np.pi/180
        lon0 = -np.pi/2 + np.pi/4.0
        lat0 = np.pi/6.

        lonB = grid['grid_lon'].values*deg2rad
        latB = grid['grid_lat'].values*deg2rad
        XB, YB, ZB = sph2cart(lonB, latB)

        # Get A grid points
        N = np.shape(lonB)[0]
        XA = XB[0:N-1,0:N-1] + XB[1:N,0:N-1] + XB[0:N-1,1:N] + XB[1:N,1:N]
        YA = YB[0:N-1,0:N-1] + YB[1:N,0:N-1] + YB[0:N-1,1:N] + YB[1:N,1:N]
        ZA = ZB[0:N-1,0:N-1] + ZB[1:N,0:N-1] + ZB[0:N-1,1:N] + ZB[1:N,1:N]

        # project onto sphere
        normA = np.sqrt(XA*XA + YA*YA + ZA*ZA)
        XA, YA, ZA = XA/normA, YA/normA, ZA/normA

        # get latlon
        lonA, latA = cart2sph(XA, YA, ZA)

        r = (lonA-lon0)*(lonA-lon0) + (latA-lat0)*(latA-lat0)
        r = np.minimum(r0*r0,r)
        r = np.sqrt(r)
        b = 2000.0*(1.0-(r/r0))
        
    return b

