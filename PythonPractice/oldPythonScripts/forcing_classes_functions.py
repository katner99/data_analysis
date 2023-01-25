import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

# ------------------ 00100 CLASS TO CONTAIN ALL VARIABLE DATA -------------------------------------------
    # fileName: Name of NCEP file for processing
    # nc_id: the ncetCDF id for the input file, initialised as -1
    # varId: variable name for the output file
    # limits: applying these limits to the interpolated data (avoids cubic fits exceeding these)
    # factor: a factor to apply. NOT used, superceded by function
    # inBulk: flag if output data is to go into the bulk file or a different one (i.e. vector one)
    # directConversion: whether the values are directly interpolated or converted using a function, or pre-processing step
    #       (i.e the conversionFunc is set to None if a simple post-interpolation conversion is applied at the same 
    #        stage as applying limits. It True and no function exists a pre-processing step specific to the variable is applied)
    # conversionFunc: the conversion function imported from classes and functions file

class forceVariable:
    def __init__(self, fileName, nc_id, varId, inVarId, limits, factor, inBulk, directConversion, units, conversionFunc, 
                type_tag, type_var, type_units, type_factor, type_time_res,type_steps_from_midnight):
        self.fileName = fileName
        self.nc_id = nc_id
        self.varId = varId
        self.inVarId = inVarId
        self.limits = limits
        self.inBulk = inBulk
        self.factor = factor
        self.directConversion = directConversion
        self.units = units
        self.conversionFunc = conversionFunc
        self.type_tag = type_tag
        self.type_var = type_var
        self.type_units = type_units
        self.type_factor = type_factor
        self.type_time_res = type_time_res
        self.type_steps_from_midnight = type_steps_from_midnight

# ------------------ 00200 CONVERSION FOR KELVIN TO DEG C -------------------------------------------
def KelvinToDegrees(valueInKelvin):
    valueInCelcius = valueInKelvin
    #  apply a secondary check value is not already in degrees C
    if valueInKelvin > 100:
        valueInCelcius = valueInKelvin - 273.15
    return valueInCelcius

# ------------------ 00300 CONVERSION FOR DEG C TO KELVIN -------------------------------------------
def DegreesToKelvin(valueInDegrees):
    #  apply a secondary check value is not already in degrees K
    valueInKelvin = valueInDegrees + 273.15
    return valueInKelvin

# ------------------ 00400 CONVERSION FROM PERCENTAGE TO FACTOR (0-1) -------------------------------------------
def PCtoFraction(valueInPC):
    return valueInPC/100.0

# ------------------ 00500 CONVERT A GEOGRAPHIC POLAR REGION TO A FLAT CARTESIAN PLANE -------------------------------------------
def ConvertPolarToCartesian(lat, lng, u, v):
    x = np.tan(lat) * np.sin(lng)
    y = -np.tan(lat) * np.cos(lng)
    pos = [x,y]

    u_new = u * np.cos(lng) - v * np.sin(lng)
    v_new = v * np.cos(lng) + u * np.sin(lng)

    dir = [u_new, v_new]

    return pos,dir

# ------------------ 00600 CONVERT A ROW ON A WOA GRID TO ARRAY OF VECTORS (U,V) -------------------------------------------
def ConvertRowOfVectorsToPolar(vect_u, vect_v, lngs):
    new_u = vect_u * np.cos(np.radians(lngs)) - vect_v * np.sin(np.radians(lngs))
    new_v = vect_v * np.cos(np.radians(lngs)) + vect_u * np.sin(np.radians(lngs))
    return(new_u, new_v)

# ------------------ 00700 ARRAY OF VECTORS (U,V) TO A ROW ON A WOA GRID -------------------------------------------
def MakeARowFromAvgFromPolar(u_val, v_val, lngs):
    new_u = u_val * np.cos(np.radians(lngs)) + v_val * np.sin(np.radians(lngs))
    new_v = v_val * np.cos(np.radians(lngs)) - u_val * np.sin(np.radians(lngs))
    return(new_u, new_v)

# ------------------ 00800 THE ANGLE NORTH WRT ORCA GRID -------------------------------------------
def fx(glam, gphi):
    # z = double( 2. * cos( rad*plam ) * tan( rpi/4. - rad*pphi/2. ))
    res = 2 * np.cos( np.radians(glam) ) * np.tan( np.radians(180)/4. - np.radians(gphi)/2)
    return res
def fy(glam, gphi):
    # z = double( 2. * sin( rad*plam ) * tan( rpi/4. - rad*pphi/2. ))
    res = 2 * np.sin( np.radians(glam) ) * np.tan( np.radians(180)/4. - np.radians(gphi)/2)
    return res

# ------------------ 00900 INTERPOLATION -------------------------------------------
def Interpolate(out_dims, xLim, yLim, tLim, data_in, lons, lats, nc_lat, nc_lon, directConversion, conversionFunc, limits, type, point):

    # out_dims, dimensions of array to output
    # xLim, x dimension of input array
    # yLim, y dimension of input array
    # tLim, time dimension
    # data_in, input data
    # lons, longitude array of input data
    # lats, latitude array of input data
    # nc_lat, grid of latitudes for output grid
    # nc_lon, grid of longitudes for output grid
    # directConversion, flag for whether the number is interpolated directly or req pre-processing
    # conversionFunc, conversion funcion if required
    # limits, limits of outputvalues that are applied post-interpolation

    # Set output array
    data_out = np.zeros(out_dims,dtype=float)

    print(type,tLim, data_out.shape,data_in.shape)

    if type == 'ncep':
        for t in range(0,tLim): # t is lon
            print("Progress: day ",t,"/",tLim,end='\r') # progress update to screeen

            data_to_use = data_in[t,:,:] # split out time point array of data

            # extend the grid longitudinally by 4 points each end to allow wrapping of data for fitting splines
            data_to_add_before = data_to_use[:,-4:]
            data_to_add_after = data_to_use[:,:4]
            data_in_extended = np.concatenate((data_to_add_before,data_to_use,data_to_add_after),axis = 1)
            # add nonsensical longtitudes to allow for the fit if needed
            lons_extended = np.append(lons[-4:]-360,np.append(lons,lons[:4]+360))

            # extend the grid in latitude and flip data as it moves over the pole
            data_to_add_before = np.flip(data_in_extended[:4,:],axis=0)
            data_to_add_after = np.flip(data_in_extended[-4:,:],axis=0)
            data_in_extended = np.concatenate((data_to_add_before,data_in_extended,data_to_add_after),axis = 0)
            # add nonsensical latitudes to allow for contraining fit if needed
            lats_before = 180. - np.flip(lats[:4],axis=0)
            lats_after = -1*np.flip(lats[-4:],axis=0) - 180.
            lats_extended = np.concatenate((lats_before,lats,lats_after),axis = 0)

            # create a grid of latitudes and longitudes that match the input data
            grid_lat,grid_lon = np.meshgrid(lats_extended,lons_extended)
            # then convert these to a list of points
            points = np.column_stack([grid_lat.ravel(), grid_lon.ravel()])

            # create a list of the values that match the points fron the input data
            vals = []
            for ix in range(0,len(lons_extended)):
                for iy in range(0,len(lats_extended)):
                    vals.append(data_in_extended[iy,ix])
            vals = np.array(vals) # and convert to array

            # interpolate the established grid of points and values from the input data
            # defined above. The points extracted from the grid of data is defined by
            # the nc_lat and nc_lon defined by the mask file. Method is type of interpolation
            # here we use cubid to maintain the peaks/troughs as much as possible in the data
            # linear, would overly remove these.
            data_out[t,:,:] = griddata(points,vals,(nc_lat, nc_lon), method='cubic')

            # tidy up the data for cubic fits allowing out of range data
            # and apply a conversion function if it exists fo this variable
            for ix in range(0,xLim): # x is lon
                for iy in range(0,yLim): # y is lat
                    valToUse = data_out[t,iy,ix]
                    if directConversion == False and conversionFunc != None:
                        valToUse = conversionFunc(valToUse)

                    if valToUse < limits[0]:
                        valToUse = limits[0]

                    if valToUse > limits[1]:
                        valToUse = limits[1]

                    data_out[t,iy,ix] =  valToUse



    else:
        # JRA55 or ERA5, read in conversion file and use this
        for t in range(0,tLim): # t is lon
            print("Progress: day ",t,"/",tLim,end='\r') # progress update to screeen
    
            data_to_use = data_in[t,:,:] # split out time point array of data
            # print(data_to_use.shape, out_dims)
            if point == "T":
                conversionFile = "t_jra_to_mesh.dat"
            if point == "U":
                conversionFile = "u_jra_to_mesh.dat"
            if point == "V":
                conversionFile = "v_jra_to_mesh.dat"
            conv_file = open(conversionFile, "r")
            conv_lines = conv_file.readlines()
            for line in conv_lines:
                words = line.split(',')
                ty = int(words[0]) 
                tx = int(words[1]) 
                count = int(words[2])
                contrib = []
                # print(words)
                for c in range(3,len(words)-1,2):
                    contrib.append( ( int(words[c]),int(words[c+1]) ) )

                valToUse = 0
                for p in contrib:
                    valToUse += data_to_use[p[0],p[1]]/count
                    # print(p, valToUse)
                    # print(contrib)
                    # input("--")
                # print(ty,tx,data_out.shape)

                val = valToUse #data_out[t,iy,ix]
                if directConversion == False and conversionFunc != None:
                    val = conversionFunc(val)

                if val < limits[0]:
                    val = limits[0]

                if val > limits[1]:
                    val = limits[1]

                data_out[t,ty,tx] =  val

                # data_out[t,ty,tx] = valToUse

    return data_out

def Interpolate_perDay(out_dims, xLim, yLim, tLim, data_in, lons, lats, nc_lat, nc_lon, directConversion, conversionFunc, limits, type, point):
    # JRA55, read in conversion file and use this
        # for t in range(0,2):#tLim): # t is lon
            # print("Progress: day ",t,"/",tLim,end='\r') # progress update to screeen
    
    data_out = np.zeros((out_dims[1],out_dims[2]),dtype=float)
    data_to_use = data_in[:,:] # split out time point array of data

    if point == "T":
        conversionFile = "t_jra_to_mesh.dat"
    if point == "U":
        conversionFile = "u_jra_to_mesh.dat"
    if point == "V":
        conversionFile = "v_jra_to_mesh.dat"
    conv_file = open(conversionFile, "r")
    conv_lines = conv_file.readlines()
    for line in conv_lines:
        words = line.split(',')
        ty = int(words[0]) 
        tx = int(words[1]) 
        count = int(words[2])
        contrib = []
        # print(words)
        for c in range(3,len(words)-1,2):
            contrib.append( ( int(words[c]),int(words[c+1]) ) )

        valToUse = 0
        for p in contrib:
            valToUse += data_to_use[p[0],p[1]]/count
            # print(p, valToUse)
            # print(contrib)
            # input("--")
        # print(ty,tx,data_out.shape)

        val = valToUse #data_out[t,iy,ix]
        # print(valToUse.shape, valToUse)
        if directConversion == False and conversionFunc != None:
            val = conversionFunc(val)

        if val < limits[0]:
            val = limits[0]

        if val > limits[1]:
            val = limits[1]

        data_out[ty,tx] =  val

        # data_out[t,ty,tx] = valToUse

    return data_out


# ------------------ 01000 NORTH POLE COHERANCE -------------------------------------------
def NorthPoleWindVectorCoherance(u_data_in, v_data_in, data_in_lon, data_in_lat, show):

    # Adjust north-pole row for cohereance
    for t in range(0,u_data_in.shape[0]): # 192 is a big one for the sake of an example (2018 data)
        # take the top row at the north pole
        u_nth = u_data_in[t,-1,:]
        v_nth = v_data_in[t,-1,:]
        # respective co-ords
        lon_nth = np.array(data_in_lon)[:]
        lat_nth = np.array(data_in_lat)[-1:]
        u_nth_pol, v_nth_pol = ConvertRowOfVectorsToPolar(u_nth, v_nth, lon_nth)
        u_mean = np.mean(u_nth_pol)
        v_mean = np.mean(v_nth_pol)

        u_nth_corr, v_nth_corr = MakeARowFromAvgFromPolar(u_mean, v_mean, lon_nth)

        if show == True: # set to false but code left if needed, will plot vectors for each time point
            # -------------- IF YOU NEED TO SEE IT IN POLAR REGION -----------
            u_nth_pol = []
            v_nth_pol = []
            # project these two to a 2D array centred on the pole 
            pol_pnts = []
            pol_v = []
            meanPosList = []
            for i in range(len(u_nth)):
                new_pos, new_dir = ConvertPolarToCartesian(np.radians(lat_nth[0]),np.radians(lon_nth[i]), u_nth[i], v_nth[i])
                pol_pnts.append(new_pos)
                pol_v.append(new_dir)
                meanPosList.append(new_dir)

            meanPosList = np.array(meanPosList)
            meanPos = [ np.mean(meanPosList[:,0]), np.mean(meanPosList[:,1])] 
            meanPos1 = [u_mean,v_mean]
            pol_pnts = np.array(pol_pnts)
            pol_v = np.array(pol_v)

            plt.quiver(pol_pnts[:,0],pol_pnts[:,1],pol_v[:,0],pol_v[:,1], scale=3, width=0.001, pivot='tail', headwidth=10)
            plt.quiver([0],[0], meanPos[0], meanPos[1], scale=3, width=0.001, pivot='tail', headwidth=10)
            plt.quiver([0],[0], meanPos1[0], meanPos1[1], scale=3, width=0.001, pivot='tail', headwidth=10, color='red')
            plt.show()
            # -------------- END: IF YOU NEED TO SEE IT IN POLAR REGION -----------

        # these two rows are then used to over-write the last row (at the north pole)
        u_data_in[t,-1,:] = u_nth_corr
        v_data_in[t,-1,:] = v_nth_corr

    # returns the adjusted input data
    return u_data_in, v_data_in
    