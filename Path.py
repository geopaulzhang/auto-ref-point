#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#### The main program for estimating the optimal reference point for each wetland unit
#### and output indices of interferograms to remove for SBAS processing
###  Including 7 steps in the script

import os
import h5py
import sys
import numpy as np
from mintpy.objects import ifgramStack
from mintpy.utils import readfile
from mintpy.utils import utils as ut
from mintpy.utils import network as pnet, ptime, readfile, utils0 as ut0
import inspect
import csv
from statistics import median, mean
from skimage import measure, morphology as morph, segmentation as seg
import shapefile as shp
from shapely.geometry import Point
from shapely.geometry import shape 
import pandas as pd
from math import sin, cos, sqrt, atan2, radians
import matplotlib.pyplot as plt
import datetime
from collections import Counter
from osgeo import gdal, gdal_array
import heapq

def label_conn_comp(mask, min_area, erosion_size=5, print_msg=False):
    """Label / clean up the conn comp (mask)
    Parameters: mask         - 2D np.ndarray of bool/int
                min_area     - float, minimum region/area size
                erosion_size - int (odd number), size of erosion structure
                               set to 0 to turn it off.
    Returns:    label_img    - 2d np.ndarray of int, labeled array where all
                               connected regions are assigned the same value
                num_label    - int, number of labeled regions
    """

    # label
    label_img, num_label = measure.label(mask, connectivity=1, return_num=True)

    ## remove small regions
    min_area = min(min_area, label_img.size * 3e-3)
    if print_msg:
        print(f'remove regions with area < {int(min_area)}')
    mask = morph.remove_small_objects(label_img, min_size=min_area, connectivity=1)
    label_img[mask == 0] = 0
    # update label
    label_img, num_label = measure.label(label_img, connectivity=1, return_num=True) # re-label

    ## remove regions that would disappear after erosion
    # to ensure the consistency between label_img and label_bound
    if erosion_size > 0:
        erosion_structure = np.ones((erosion_size, erosion_size))
        label_erosion_img = morph.erosion(label_img, erosion_structure).astype(np.uint8)

        erosion_regions = measure.regionprops(label_erosion_img)
        if len(erosion_regions) < num_label:
            if print_msg:
                print('regions lost during morphological erosion operation:')

            label_erosion = [reg.label for reg in erosion_regions]
            for orig_reg in measure.regionprops(label_img):
                if orig_reg.label not in label_erosion:
                    label_img[label_img == orig_reg.label] = 0
                    if print_msg:
                        print('label: {}, area: {}, bbox: {}'.format(orig_reg.label,
                                                                     orig_reg.area,
                                                                     orig_reg.bbox))

            # update label
            label_img, num_label = measure.label(label_img, connectivity=1, return_num=True)

    return label_img, num_label

def eulid_dist(xt, yt, xo, yo):
    return np.sqrt((xt - xo) * (xt - xo) + (yt - yo) * (yt - yo))

def close_point(latbt, lonbt):
    dist = 1000
    for i in range(xrg * yrg):   ## calculate the closest sar pixel 
        temp = eulid_dist(latbt, lonbt, lat_1d[i], lon_1d[i])
        if (temp < dist):
            dist = temp
            indexbt = i
    ybt = indexbt % yrg
    xbt = (indexbt - ybt) / yrg
    ybt = int(ybt)
    xbt = int(xbt)
    return indexbt, ybt, xbt
    
def close_point_wet(latbt, lonbt, wu):
    dist = 1000
    for i in range(wetnum):
        temp = eulid_dist(latbt, lonbt, lalo[i][0], lalo[i][1])
        if (temp < dist):
            dist = temp
            wbest = i
    return wbest ## index in the wetland

def wgs2utm(lat, lon):
    lat_rad = lat * (np.pi / 180)
    lon_rad = lon * (np.pi / 180)
    a = 6378137  # Semi-major axis of the WGS84 ellipsoid
    f = 1 / 298.257223563  # Flattening of the WGS84 ellipsoid
    k0 = 0.9996
    e_squared = 2 * f - f * f
    N = a / sqrt(1 - e_squared * sin(lat_rad)**2)
    T = np.tan(lat_rad)**2
    C = e_squared / (1 - e_squared) * cos(lat_rad)**2
    lon0 = -81 * (np.pi / 180)  # Central meridian for UTM Zone 17 North in radians
    A = (lon_rad - lon0) * cos(lat_rad)
    M = a * (
       (1 - e_squared / 4 - 3 * e_squared**2 / 64 - 5 * e_squared**3 / 256) * lat_rad -
       (3 * e_squared / 8 + 3 * e_squared**2 / 32 + 45 * e_squared**3 / 1024) * sin(2 * lat_rad) +
       (15 * e_squared**2 / 256 + 45 * e_squared**3 / 1024) * sin(4 * lat_rad) -
       (35 * e_squared**3 / 3072) * sin(6 * lat_rad))
    x = k0 * N * (A + (1 - T + C) * A**3 / 6 + (5 - 18 * T + T**2 + 72 * C - 58 * e_squared) * A**5 / 120)
    y = k0 * (M + N * np.tan(lat_rad) * (A**2 / 2 + (5 - T + 9 * C + 4 * C**2) * A**4 / 24 + (61 - 58 * T + T**2 + 600 * C - 330 * e_squared) * A**6 / 720))
    x_utm = x + 500000  # False easting for UTM
    y_utm = y   # False northing for Northern hemisphere
    return x_utm/30, y_utm/30 ## convert to pixel 

def binary_image(input_image, threshold):
    # Apply the threshold
    bi = (input_image > threshold).astype(np.uint8)

    return bi
    
def dilate_nb(x, y, size, img):
    startx = int(x - size)
    starty = int(y + size)
    for i in range(size * 2):
        for j in range(size * 2):
            img[startx + i, starty - i] = 1
            
    return img      
    
def dijkstra_with_max_avg_coherence(start, candidates, image):
    priority_queue = [(0, start, [start])]
    distances = {node: float('inf') for node in candidates}
    avg_coherences = {node: 0 for node in candidates}
    
    while priority_queue:
        dist, node, path = heapq.heappop(priority_queue)

        if node in candidates:
            if dist > distances[node]:
                continue

            for neighbor in candidates:
                if neighbor != node:
                    new_dist = dist + euclidean_distance(node, neighbor)
                    if new_dist < distances[neighbor]:
                        distances[neighbor] = new_dist
                        avg_coherences[neighbor] = avg_coherences[node] + image[neighbor[0]][neighbor[1]]
                        new_path = path + [neighbor]
                        heapq.heappush(priority_queue, (new_dist, neighbor, new_path))

    best_candidate = max(candidates, key=lambda c: avg_coherences[c] / distances[c])

    return best_candidate, avg_coherences[best_candidate] / len(path)
    
site_select = 1
usecc = 1 ## (0 for using coherence) using connectedComponent as label instead of coherence-based label, 
reducedt = 1 ## (0 for using the entire InSAR time series) trucated period of time based on nd and ndi 
refthd = 0.90
#refthd = 0.9 ## threshold for REFERENCE coherence

thdpt = 0.95
#thdpt = 0.9
#radius = 100 ### only find pixels within 100 pixels distance- a square, not the optimal
radius = 1000
tops = radius * radius  ## Select the top 30 reference points for each wetland pixel 
normct_thd = 0.8
wetthd = 30 ## minimum percentage of connection from one wetland pixel to the radius*radius window 
pxpercomp = 5 ## how many pixel selected from each component
paththd = 0.9

#wetthd = 5
cc20 = 'exclude'
cc200 = '200cc'
cclevel = 0
conn = 6
min_area = 30
#min_area = 150 ## not the optimal
bicoh_min_area = 2
dilation = 1
################# Other thresholds #########################
uthd = 0.35 ## usability threshold of coherence
hthd = 0.7 ## high quality threshold

connthd = 80 ## the threshold for percentage of connection between the selected reference point and wetland pixels.

np.set_printoptions(threshold=sys.maxsize)

site1 = 'SacFull'
site2 = 'DelFull'
site3 = 'ClsFull'
ws = ''

for wu in [0]: ## no finding in the initial high thresholds.
#for wu in range(0, 10, 1):
#for wu in [2]:
#for wu in [4]:
    if wu == 0:
        #trgoi = range(3, 19 + 1) ## SacT34c1
        trgoi = range(4, 6 + 1) ## SacT34c1
    elif wu == 31:
        trgoi = range(23, 35 + 1) ## SacT4c1 corresponding to 20171102 - 20180326
    elif wu == 3:
        trgoi = range(15, 21) ## SacT4c1 summer, corresponding to 20170729 - 20170927
    elif wu == 2:
        trgoi = range(6, 20 + 1) ## SacT4c2 corresponding 20170729 - 20170927 dry for SacT4c1 
    elif wu == 4:
        trgoi = range(9, 19 + 1) ## SacP1Ac1
    elif wu == 6:
        trgoi = range(8, 16 + 1) ## SacT9c2
    elif wu == 8:
        trgoi = range(15, 19 + 1) ## SacT41c1
    elif wu == 81:
        trgoi = range(22, 36 + 1) ## SacT41c1_w
    elif wu == 1:
        trgoi = range(13, 20) ## SacT29c2
    elif wu == 7:
        trgoi = range(8, 21) ## SacP11c2
    elif wu == 9:
        trgoi = range(8, 22)
    elif wu == 5:
        continue
    #    trgoi = range(21, 33 + 1) ## SacT41c1_w
    
    ##### Threshold combination for three specials #####
    if wu == 1: # bad condition
        normct_thd = 0.3
    elif wu == 2:
        thdpt = 0.9
    elif wu == 6 or wu == 8:
        normct_thd = 0.75
    elif wu == 7:
        thdpt = 0.9
        normct_thd = 0.75
        paththd = 0.65 ## even 0.7 with min_area 30 does not work
    elif wu == 9:
        #paththd = 0.7
        thdpt = 0.8
    # elif wu == 8:
        # min_area = 150
    #elif wu == 0:
    #    thdpt = 0.90
        
        

    #trgoi = range(8, 18 + 1) ## ClsT15c3 May 5 to Sep 3
    nd = trgoi ## the index for SAR acquisitions
    ndi = range(trgoi[0] * conn, trgoi[-1] * conn + 1) ## the index for interferograms
    print("Original index for interferograms of interest were: {ndi}".format(ndi = ndi))
    ndi = [x for x in ndi if x % conn == 0 or x % conn == 1]
    
    if wu == 2:
    #if wu == 2: ## remove the unit 0 - the problematic one
        ndi = [x for x in ndi if x % conn == 0] ## only 12 days baseline 
        
    print("Filtered index for interferograms of interest were: {ndi}".format(ndi = ndi))

    print("The index for dates of interest are:{nd}".format(nd = nd))

    print("The length of ndi is {x}".format(x = len(ndi)))
    print("ndi is {x}".format(x = ndi))

    if cclevel == 0:
        Path = '/beanstore/Personal/Paul/InSAR/timeseries/ts_115_exclude/crops/'
        Pathout = '/beanstore/Personal/Paul/InSAR/timeseries/ts_115_exclude/Analysis/'
        Path_cls = '/beanstore/Personal/Paul/InSAR/timeseries/ts_115_part/crops/'
    else:
        Path = '/beanstore/Personal/Paul/InSAR/timeseries/ts_115_200cc/crops/'
        Pathout = '/beanstore/Personal/Paul/InSAR/timeseries/ts_115_200cc/Analysis'
        Path_cls = '/beanstore/Personal/Paul/InSAR/timeseries/ts_115_200cc/crops/'

    #################################### Preparation #######################################
    print("site index is :")
    print(site_select)

    if site_select == 1:
        ifgram_stack = Path + site1 + '/inputs/ifgramStack.h5'
        geom_file = Path + site1 + '/inputs/geometryGeo.h5' 
        inc = Path + site1 + '/SacFull.csv'
        shpvec = Path + site1 + '/SacInc.npy'
        site = site1
        
    elif site_select == 2:
        ifgram_stack = Path + site2 + '/inputs/ifgramStack.h5'
        geom_file = Path + site2 + '/inputs/geometryGeo.h5' 
        inc = Path + site2 + '/DelFull.csv'
        shpvec = Path + site2 + '/DelInc.npy'
        site = site2
    else:
        ifgram_stack = Path_cls + site3 + '/inputs/ifgramStack.h5'
        geom_file = Path_cls + site3 + '/inputs/geometryGeo.h5' 
        inc = Path_cls + site3 + '/ClsFull.csv'
        shpvec = Path_cls + site3 + '/ClsInc.npy'
        site = site3
        Path = Path_cls
        
     
    Shapath = '/beanstore/Personal/Paul/InSAR/timeseries/ts_115_exclude/Shpfile/Project_sac.shp'
    shpfile = Shapath
    temp = readfile.read(ifgram_stack, datasetName='coherence', box=(0, 0, 1, 1))[0]
    stack_obj = ifgramStack(ifgram_stack)                    
    coh = stack_obj.read(datasetName='coherence') ## connections, y, x

    conncomp = stack_obj.read(datasetName='connectComponent')
    print("This is the shape of the connected Component {cc}".format(cc = conncomp.shape))

    xrg = coh.shape[2]
    yrg = coh.shape[1]
    k   = coh.shape[0]
    cohadj = np.zeros((k, xrg, yrg))
    connadj = np.zeros((k, xrg, yrg))
    connadj_2d = np.zeros((k, xrg * yrg))
    for i in range(xrg):
        for j in range(yrg):
            cohadj[:, i, j] = coh[:, j, i]
            connadj[:, i, j] = conncomp[:, j, i]

    meta = readfile.read_attribute(ifgram_stack)
    coord_obj = ut.coordinate(meta, geom_file)
    quality = np.zeros((xrg, yrg))
        
    #usda_veg = np.transpose(usda_veg)
    lat_stack_all = np.zeros((xrg, yrg))
    lon_stack_all = np.zeros((xrg, yrg))
    lat_1d = np.zeros((xrg * yrg))
    lon_1d = np.zeros((xrg * yrg))
    coh_1d = np.zeros((xrg * yrg))
    coh_2d = np.zeros((xrg, yrg))
    sf = shp.Reader(shpfile) 
     
    sfRec = sf.records()
    all_shapes = sf.shapes()
    all_records = sf.records() 
    feature = sf.shapeRecords()[0]
    first = feature.shape.__geo_interface__

    inclusion = np.zeros((len(all_shapes), xrg * yrg))


    ct = 0
    for i in range(xrg):
        for j in range(yrg):
            lat, lon = coord_obj.radar2geo(j, i)[0:2]
            lat_stack_all[i, j] = lat
            lon_stack_all[i, j] = lon
            lat_1d[ct] = lat
            lon_1d[ct] = lon
            connadj_2d[:, ct] = connadj[:, i, j] 
            ct = ct + 1
    
    
    print("x = {x}".format(x = xrg))
    print('In total we have this many number of polygons:')

    #top = np.zeros((len(all_shapes) + 1, tops)) ## initialize the top list
    top = np.zeros(tops) ## initialize the top list - ONLY FOR ONE WETLAND UNIT. 
    print("The total number of units is {n}".format(n = len(all_shapes)))


     
    shape_ind = np.load(shpvec)
    label = np.zeros((k, xrg,yrg))
    num_ary = np.zeros(k)
    mark = np.zeros((k, xrg, yrg))


    stack_obj = ifgramStack(ifgram_stack)
    stack_obj.open()
    date12List = stack_obj.get_date12_list(dropIfgram=False) ## This is the reference and secondary intferograms dates
    dateList = stack_obj.get_date_list(dropIfgram=False) ## This is the acquisiton dates

    output = np.vstack((dateList))
    #output = output.transpose()

    with open(Path + site + '/datelist.csv',"w") as f:
        np.savetxt(f, output, delimiter=",", header="ID,disp", 
               fmt="% s")

    #### start the iteration!!!!
  
    print("wetland unit is:")
    print(wu)     
    print("Saved the date list!!")
    print("List of dates is{l}".format(l = dateList))
    print("the type of the variable is{v}".format(v = type(dateList)))
    print("start of the water change {s}".format(s = dateList[23]))
    print("start of the water change inf {s}".format(s = date12List[23 * 6]))

    print("end of the water change {s}".format(s = dateList[35]))
    print("end of the water change inf {s}".format(s = date12List[35 * 6]))


    ct = 0
    inclusion = np.zeros((xrg * yrg))
    points = pd.read_csv(inc)

    #ID = points.iloc[:,0].values
    inclusion = pd.read_csv(inc)
    lat_within = inclusion.iloc[:, 1].values
    lon_within = inclusion.iloc[:, 2].values
    include = inclusion.iloc[:, 3].values
    include = np.where(include == 6, 0, include) ## remove SacP2 as background

    for i in range(k):
            temp = cohadj[i, :, :]
            mask = np.where(temp >= uthd, 1, 0) ## connection using 0.5 instead of 0.7
            mask = mask.astype(int)
            label_img, num = label_conn_comp(mask, min_area, 0, False)
            label[i, :, :] = label_img
            num_ary[i] = num

    ind = np.arange(0, xrg * yrg, 1)      
    out = np.vstack((ind, lat_1d, lon_1d, connadj_2d[64, :]))
    out = out.transpose()


    with open(Path + site + '/connection.csv',"w") as f:
            np.savetxt(f, out, delimiter=",", header="ID, lat,lon, conn_comp", 
                   fmt="%f") 

    if usecc == 1:
        label = connadj
        print("Using CC label instead of Coherence label!")
        
    if reducedt == 1:
        label = label[ndi, :, :]
        cohadj = cohadj[ndi, :, :]
        krt = len(ndi) ## important - change the temporal length
        print("Using REDECED PERIOD OF TIME!")
    
    coh_2d = np.mean(cohadj, axis=0) ## average coherence value
    
    
    ct = 0
    for i in range(xrg):
        for j in range(yrg):
            coh_1d[ct] = np.mean(cohadj[:, i, j]) ## recording the mean value of coh in time. 
            ct = ct + 1
    
    ############################# Filtering in the beginning to accelerate ########################
    ###############################################################################################
    np.save(Path + site + '/coh_2d.npy', coh_2d)
    print("the size of the label is {s}".format(s = label.shape))
    print("the size of the cohadj is {s}".format(s = cohadj.shape))

    ######################### Compare the best point and northwest point - each wetland pixel connectivity in time ###
    # latnw = 39.4382
    # lonnw = -122.1720 ## the northwest point
    # latbt = 39.435565
    # lonbt = -122.155205

    # latnw = 39.4382
    # lonnw = -122.172 ## on the top
    # #latnw = 39.4304
    # #lonnw = -122.155 ## on the right
    # latbt = 39.434755
    # lonbt = -122.157905
    # latne = 39.438
    # lonne = -122.16 ## northeastern point outside of P1Ac1 
    # lats = 39.423
    # lons = -122.162 ## southern point within P1Ac1
    
    lat_st = 39.409
    lon_st = -122.197 ## the single point on the street.
    
    lat_ref = 39.389665
    lon_ref = -122.182205

    # indexnw, ynw, xnw = close_point(latnw, lonnw)

    # indexbt, ybt, xbt = close_point(latbt, lonbt)
    ##indexne, yne, xne = close_point(latne, lonne)
    
    indexst, ynw, xnw = close_point(lat_st, lon_st)
    indexref, ybt, xbt = close_point(lat_ref, lon_ref)
    
    #print('indexne is {f}'.format(f = indexne))

    ct = 0
    latout = np.zeros(len(np.where(shape_ind == wu + 1)[0]))
    lonout = np.zeros(len(np.where(shape_ind == wu + 1)[0]))
    conn_nw = np.zeros(len(np.where(shape_ind == wu + 1)[0]))
    conn_bt = np.zeros(len(np.where(shape_ind == wu + 1)[0]))

    for x in range(xrg):
        for y in range(yrg):
            
            ind = yrg * x + y
            if shape_ind[ind] == wu + 1: ## just for SacP1Ac1 
                nonzero = np.where(label[:, x, y] != 0)[0] ## select the non-zero connected component
                latout[ct] = lat_1d[ind]
                lonout[ct] = lon_1d[ind]
                conn_nw[ct] = len(np.where(label[nonzero, xnw, ynw] == label[nonzero, x, y])[0])
                conn_bt[ct] = len(np.where(label[nonzero, xbt, ybt] == label[nonzero, x, y])[0])
                ct = ct + 1

    ind = np.arange(0, len(np.where(shape_ind == wu + 1)[0]), 1)      
    out = np.vstack((ind, latout, lonout, conn_nw, conn_bt))
    out = out.transpose()


    with open(Path + site + '/final_pk.csv',"w") as f:
            np.savetxt(f, out, delimiter=",", header="ID, lat,lon, conn_nw, conn_bt", 
                   fmt="%f")          
    ################################## Select the optimal reference point #######################
    ## Step 1: Value the pixels based on wetland units (Done in the script "Wetland_unit_inclusion.py")
    ## Step 2: Get the reference pixel candidates 
    print("STEP 1 DONE")
    print("STEP 2: start QUALITY analysis!!!")
    qua1d = np.zeros(xrg * yrg)
    pt1d = np.zeros(xrg * yrg)
    ct = 0
    for i in range(xrg):
        for j in range(yrg):
            
            temp = cohadj[:, i, j]    
            pt = len(np.where(temp > refthd)[0]) / krt ## do not use the percentage
            pt1d[ct] = pt
            #pt = coh_2d[i, j] ## instead, using the average coherence
            index = yrg * i + j ## the index for ct in inclusion, based on xrg first and then yrg for loops 
            if ((pt > thdpt) & (include[index] == 0)):
                quality[i][j] = 1
                qua1d[ct] = 1
            if(include[index] > 0): # wetland pixel
                pt1d[ct] = 0
            ct = ct + 1
            
    ind = np.arange(0, xrg * yrg, 1)      
    out = np.vstack((ind, lat_1d, lon_1d, pt1d))
    out = out.transpose()


    with open(Path + site + '/quality_rt.csv',"w") as f:
            np.savetxt(f, out, delimiter=",", header="ID, lat,lon, quality", 
                   fmt="%f")

                       
    print("end reference candidate analysis!!!")              
    #############################################################################################
    ## Step 3: using connected component as label (Done)
    ## Step 4: Searching for radius * radius neighbor pixels and calculate the number of connections
    print("STEP 3 and 4")

    wetnum = len(np.where(shape_ind == wu + 1)[0])
    #count = np.zeros((wetnum, radius * radius)) ## the frequency of the connected in time domain
    targ = np.zeros((wetnum, radius * radius)) ## the index of the target pixel
    normct = np.zeros((wetnum, radius * radius)) ## normalized the top list
    lalo = np.zeros((wetnum, 2))
    print('Starting calculating the frequency and reference point')
    wetind = 0 ## used a index for recording the wetland pixels. 
    for i in range(xrg):
        print(i)
        for j in range(yrg):
            
            ind = yrg * i + j
            
            if(shape_ind[ind] == wu + 1): ## a wetland pixel IN THE UNIT OF INTEREST, this is using the wetland unit as inclusion pixel value
                label_ts = label[:, i, j]
                nonzero = np.where(label_ts != 0)[0] ## select the non-zero connected component
                label_filt = label_ts[nonzero]
                lalo[wetind][0] = lat_1d[ind]
                lalo[wetind][1] = lon_1d[ind]
                ulx = int(i - radius/2)
                uly = int(j + radius/2) ## upper left corner pixel
                ct = 0
                for m in range(radius):
                    for n in range(radius):
                        xt = ulx + m
                        yt = uly - n
                       
                        if xt >= 0 and yt >= 0 and xt < xrg and yt < yrg:
                            indt = yrg * xt + yt ## the index of the target pixel
                            if quality[xt, yt] == 1: ## the target pixel has to be a qualified pixel
                                
                                label_tg = label[:, xt, yt] ## target pixel label ts
                                label_tg = label_tg[nonzero] ## filter as the same as the original pixel 
                                label_dif = label_filt - label_tg
                                #count[ind][ct] = len(np.where(label_dif == 0)[0]) ## dif == 0 indicates the target pixel is connected to the original pixel 
                                count = len(np.where(label_dif == 0)[0]) ## dif == 0 indicates the target pixel is connected to the original pixel 
                                targ[wetind][ct] = indt
                                #normct[ind][ct] = count[ind][ct]/ krt 
                                normct[wetind][ct] = count / krt
                                
                                ct = ct + 1
                            else:  ## the pixel is not qualified as a ref point. 
                                #count[ind][ct] = -1
                                targ[wetind][ct] = -1
                                normct[wetind][ct] = -1
                                ct = ct + 1
                        else: ## pixel is outside from the image!
                            #count[ind][ct] = -1
                            targ[wetind][ct] = -1
                            normct[wetind][ct] = -1
                            ct = ct + 1
                wetind = wetind + 1
    print('wetind is {x}'.format(x = wetind))
    #np.save(Path + site + '/normct.npy', normct)
    #np.save(Path + site + '/targ.npy', targ)
            # else:
                # #count[ind][:] = -1  ## the original pixel is not a wetland pixel 
                # targ[wetind][:] = -1
                # normct[wetind][:] = -1

    #####################################################################################
    # Step 5: Sort the radius * radius pixels based on connctions
    # Step 6: Intersect with all pixels within a wetland

    print('STEP 5 and 6. start selecting the optimal reference points!')
    # for t in range(1, len(all_shapes) + 1): ### because it is t + 1, so the output is FID + 1 !!!!!!
        # print(t)
    t = wu + 1 ## just for one wetland of interest
    ct = 0

    for x in range(xrg):
        for y in range(yrg):
            
            ind = yrg * x + y
            if shape_ind[ind] == t: ## only continue the analysis if the point belongs to this shapefile
                comb = zip(normct[ct], targ[ct])
                srt = sorted(comb, key=lambda k: (-k[0], k[1]))
                temp = [tple[1] for tple in srt] ## only for the reference pixels, second column
                temp = np.array(temp) ## very important step!
                #temp = temp[0:tops] ## select the top 100
                
                normtmp = [tple[0] for tple in srt] ## the first column
                normtmp = np.array(normtmp)
                goodref = np.where(normtmp > normct_thd)[0] ## thresholding based on normalized connections between the two points
                temp = temp[goodref] ## the filter based on normalized count. 
                #print("temp value is : {x}".format(x = temp))
                
                if (len(goodref) < wetthd):
                   continue
                   
                if (ct == 0): 
                   top[0:len(temp)] = temp      
                else:
                   tmp = np.intersect1d(top, temp) ## gradually reduce the candiates
                   top[0:len(tmp)] = tmp
                   top[len(tmp):] = 0 ## remove the rest of items # do not use -1, use :]
                   
                #contain = np.isin(top, indexne)
                
                #print(ind)
                #print("the top value is: {x}".format(x = top))
                ct = ct + 1

    np.save(Path + site + '/top.npy', top)
    #output = np.zeros(top.shape[0], 10) ## top 10.
    output = np.zeros(10) ## top 10.

    #for i in range(1, top.shape[0]):
    i = wu + 1
    tmp = top
    tmp = [x for x in tmp if x != 0]
    refcount = len(tmp)
    print('refcount is {x}'.format(x = refcount))
    latwet = np.where(shape_ind == i)[0]
    lonwet = np.where(shape_ind == i)[0]
    cohref = np.zeros(refcount) ## calculate mean value for coherence
    closedist = np.zeros(refcount)
    np.save(Path + site + '/latwet.npy', latwet)
    np.save(Path + site + '/tmp.npy', tmp)
    np.save(Path + site + '/top.npy', top)
    
    ############################## Create a new connected component ##################################
    ycc = [element % yrg for element in top]
    xcc = [(element - ytarg) // yrg for element, ytarg in zip(top, ycc)]
    ycc = np.array(ycc).astype(int)
    xcc = np.array(xcc).astype(int)
    bi_image = np.zeros((xrg, yrg))
    bi_image[xcc, ycc] = 1
    bi_image = bi_image.astype(int)
    label_cc, num_features = label_conn_comp(bi_image, min_area, 0, False)
    label_1d = np.zeros(xrg * yrg)
    ct = 0
    for i in range(xrg):
        for j in range(yrg):
              
            label_1d[ct] = label_cc[i, j]
            ct = ct + 1
    
    ind = np.arange(0, xrg * yrg, 1) 
    out = np.vstack((ind, lat_1d, lon_1d, label_1d))
    out = out.transpose()
    with open(Path + site + '/' + 'newlabel_' + str(wu) + '_' + str(min_area) + '_dil.csv',"w") as f:
            np.savetxt(f, out, delimiter=",", fmt="%f")
    print("********************************************************************")
    print("The number of connected component is:")
    print(num_features)
    
    ################################## Selecting 5 points for each component #######################
    print("Selection of points from each component nearest to the wetland central point")
    
    cwetla = np.mean(latout)
    cwetlo = np.mean(lonout)
    cwet = (cwetla, cwetlo)
    top_px_lat = np.zeros((pxpercomp, num_features))
    top_px_lon = np.zeros((pxpercomp, num_features))
    index_px = np.zeros((pxpercomp, num_features))
    dist_px = np.zeros((pxpercomp, num_features)) ## recording distance.
    coh_px = np.zeros((pxpercomp, num_features)) ## recording distance.
    
        
    for i in range(1, num_features + 1): ## +1 because in total there are num_features + 1 features (including the background)
        rg = np.where(label_1d == i)[0]
        rglen = len(np.where(label_1d == i)[0])
        coordinates = list(enumerate(zip(lat_1d[rg], lon_1d[rg]))) 
        cohlist = coh_1d[rg] ## coh_1d for those points, pt1d cannot distinguish the candidate points
        #print(cohlist)
        #distances = [(coord, eulid_dist(cwet, coord)) for coord in coordinates]
        distances = [(index, coord, cohlist, eulid_dist(cwetla, cwetlo, coord[0], coord[1])) for index, coord in coordinates]
        # Sort distances in ascending order
        sorted_distances = sorted(distances, key=lambda x: x[3])
        sorted_indices = sorted(range(len(distances)), key=lambda k: distances[k][3])
        
        coord_res = [coord for index, coord, cohlist, distance  in sorted_distances]
        disttemp = [distance for index, coord, cohlist, distance  in sorted_distances]
        #cohtemp = [cohlist for index, coord, cohlist, distance in sorted_distances]
        cohtemp = [cohlist[i] for i in sorted_indices]
        
        #index_res = [index for index, coord, distance in sorted_distances]
        index_res = [rg[i] for i in sorted_indices] ## sorted rg. 
        index_res = np.array(rg).astype(int)
        
        
        # np.save(Path + site + '/srt' + str(i) + '.npy', sorted_distances)
        # np.save(Path + site + '/cohlist' + str(i) + '.npy', cohlist)
        # np.save(Path + site + '/cohtemp' + str(i) + '.npy', cohtemp)
        #print(cohtemp[:pxpercomp])
        #print(cohtemp)
        
        # Get the top 5 coordinates with the shortest distances
        top_coords = coord_res[:pxpercomp]
        top_px_lat[:, i - 1] = [coord[0] for coord in top_coords]
        top_px_lon[:, i - 1] = [coord[1] for coord in top_coords] ## minus 1 because 0 is the background.
        index_px[:, i - 1] = index_res[:pxpercomp]
        dist_px[:, i - 1] = disttemp[:pxpercomp]
        coh_px[:, i - 1] = cohtemp[:pxpercomp]
        
    np.save(Path + site + '/top_px_lat.npy', top_px_lat)
    np.save(Path + site + '/top_px_lon.npy', top_px_lon)
    np.save(Path + site + '/coh_px.npy', coh_px)
    ######################### Dijkstra ! ###################################################################
    print("")
    index, ywet, xwet = close_point(cwetla, cwetlo)
    start = (xwet, ywet)
   
    index_flat = index_px.flatten()
    index_flat = np.array(index_flat).astype(int)
    dist_flat = dist_px.flatten()
    dist_flat = np.array(dist_flat)
    coh_flat = coh_px.flatten()
    coh_flat = np.array(coh_flat)
    print("dist_px is {x}".format(x = dist_px))
    print("dist_flat is {x}".format(x = dist_flat))
    
    y_cand = [element % yrg for element in index_flat]
    x_cand = [(element - y) // yrg for element, y in zip(index_flat, y_cand)]
    y_cand = np.array(y_cand).astype(int)
    x_cand = np.array(x_cand).astype(int)
    coherence_cand = coh_1d[index_flat]

    candidates = {}

    for x, y, coherence in zip(x_cand, y_cand, coherence_cand):
        candidates[(x, y)] = coherence
    
    print("candidates are {x}".format(x = candidates))
    #best_candidate, max_avg_coherence = dijkstra_with_max_avg_coherence(start, candidates.keys(), coh_2d) ## Dijkstra algorithm! start - center wetland pixel, candiates from each component

    #print(best_candidate)
    #print(max_avg_coherence)
    
    ######################## connectivity based on binary map ##############
    print("calculate connectivity based on binary map")
    ### Step 1 - build binary map #### which is bicod_2D
    ### Step 2 - Find the closest point for each candidate 
    print("Find the closest point for each candidate")
    
    clpt = np.zeros(index_flat.shape[0]) ## the closes point to each of the top candidate
    lalo = np.zeros((wetnum, 2))
    wetind = 0 ## used a index for recording the wetland pixels. 
    wetloc = np.zeros(wetnum) ## record the index for each wetland pixel
    toplat_flat = top_px_lat.flatten()
    toplon_flat = top_px_lon.flatten()
    for i in range(xrg):
            for j in range(yrg):
                ind = yrg * i + j
                
                if(shape_ind[ind] == wu + 1): ## a wetland pixel IN THE UNIT OF INTEREST, this is using the wetland unit as inclusion pixel value
                    lalo[wetind][0] = lat_1d[ind]
                    lalo[wetind][1] = lon_1d[ind]
                    wetloc[wetind] = ind
                    wetind = wetind + 1
                    
    for i in range(clpt.shape[0]):
        clpt[i] = close_point_wet(toplat_flat[i], toplon_flat[i], wu)

    clpt_ind = wetloc[np.array(clpt).astype(int)]
    clpt_ind = np.array(clpt_ind).astype(int)
    print("the closest points are {x}".format(x = clpt_ind))          

    y_cl = [element % yrg for element in clpt_ind]
    x_cl = [(element - y) // yrg for element, y in zip(clpt_ind, y_cand)]
    y_cl = np.array(y_cl).astype(int)
    x_cl = np.array(x_cl).astype(int)
    print("y and x of the closest points are:")
    print(y_cl)
    print(x_cl)
    
    ## Step 3 Check Conncectivity ##
    
    bicoh_2d= binary_image(coh_2d, paththd)
    print("paththd is {x}".format(x = paththd))
    ywet = [element % yrg for element in wetloc]
    xwet = [(element - y) // yrg for element, y in zip(wetloc, ywet)]

    for i in range(wetnum):
        bicoh_2d[int(xwet[i]), int(ywet[i])] = 1 ## give wetland connectivity
        

    for i in range(wetnum):
        bicoh_2d = dilate_nb(xwet[i], ywet[i], dilation, bicoh_2d)
                    
    bicoh_2d = bicoh_2d.astype(int)
    label_bicoh, num = label_conn_comp(bicoh_2d, bicoh_min_area, 0, False) 
    np.save(Path + site + '/label_bicoh.npy', label_bicoh)
    np.save(Path + site + '/bicoh_2d.npy', bicoh_2d)
    connectivity = np.zeros(index_flat.shape[0])
    
    
    for i in range(index_flat.shape[0]):
        if (label_bicoh[x_cand[i], y_cand[i]] == label_bicoh[x_cl[i], y_cl[i]]):
            connectivity[i] = 1
        else:
            connectivity[i] = 0
        
    print("connectivity of all candidiates are:")
    print(connectivity)
    
    # Find the indices where connectivity is equal to 1
    indices_connectivity_1 = np.where(connectivity == 1)[0]

    # Subset dist_flat and connectivity based on the indices with connectivity equal to 1
    dist_px_subset = dist_flat[indices_connectivity_1]
    connectivity_subset = connectivity[indices_connectivity_1]
    lat_sub = toplat_flat[indices_connectivity_1]
    lon_sub = toplon_flat[indices_connectivity_1]
    coh_sub = coh_flat[indices_connectivity_1]
    comb = zip(dist_px_subset, connectivity_subset, lat_sub, lon_sub, coh_sub)
    srt = sorted(comb, key=lambda k: (k[0])) ## use distance ascending order
    latsrt = [lat_sub for dist_px_subset, connectivity_subset, lat_sub, lon_sub, coh_sub in srt]
    lonsrt = [lon_sub for dist_px_subset, connectivity_subset, lat_sub, lon_sub, coh_sub  in srt]
    dissrt = [dist_px_subset for dist_px_subset, connectivity_subset, lat_sub, lon_sub, coh_sub  in srt]
    cohsrt = [coh_sub for dist_px_subset, connectivity_subset, lat_sub, lon_sub, coh_sub  in srt]
    #### Step 4. Export results with connectivity, lat and lon 
    bicoh_1d = np.zeros(xrg * yrg)
    label_1d = np.zeros(xrg * yrg)
    ct = 0
    for i in range(xrg):
            for j in range(yrg):
                bicoh_1d[ct] = bicoh_2d[i, j]
                label_1d[ct] = label_bicoh[i, j]
                ct = ct + 1
    ind = np.arange(0, index_flat.shape[0], 1)
    out = np.vstack((ind, toplat_flat, toplon_flat, connectivity, coh_flat))
    out = out.transpose()
    with open(Path + site + '/' + 'final_selection' + str(wu) + '_' + str(min_area) + 'dil.csv',"w") as f:
            np.savetxt(f, out, delimiter=",", fmt="%f")
            
    ind = np.arange(0, xrg * yrg, 1)
    out = np.vstack((ind, lat_1d, lon_1d, bicoh_1d, label_1d))
    out = out.transpose()
    with open(Path + site + '/' + 'check' + str(wu) + '_' + str(min_area) + 'dil.csv',"w") as f:
            np.savetxt(f, out, delimiter=",", fmt="%f")
            
    ind = np.arange(0, index_flat.shape[0], 1)
    out = np.vstack((ind, lat_1d[clpt_ind], lon_1d[clpt_ind]))
    out = out.transpose()
    with open(Path + site + '/' + 'clpt' + str(wu) + '_' + str(min_area) + 'dil.csv',"w") as f:
            np.savetxt(f, out, delimiter=",", fmt="%f")
            
    ind = np.arange(0, len(latsrt), 1)
    out = np.vstack((ind, latsrt, lonsrt, dissrt, cohsrt))
    out = out.transpose()
    with open(Path + site + '/' + 'thefinal' + str(wu) + '_' + str(min_area) + 'dil.csv',"w") as f:
            np.savetxt(f, out, delimiter=",", fmt="%f")
    #conn = is_path_connected(bicoh_2d, (x_cl[0], y_cl[0]), (x_cand[0], y_cand[0]), tolerance=2)
    # if refcount > 0:
        # for j in range(len(tmp)):
            # ytarg = tmp[j] % yrg
            # xtarg = (tmp[j] - ytarg) / yrg
            # closedist[j] = 100000
            # cohref[j] = np.mean(cohadj[:, int(xtarg), int(ytarg)])
            # #print('ref point is {x}'.format(x = lat_1d[int(tmp[j])]))
            # #print('target point is {x}'.format(x = latwet[m]))
            # for n in range(latwet.shape[0]): ## calculat the nearest distance between the ref and entire wetland
                # #print('target point is {x}'.format(x = lat_1d[latwet[n]]))
                # #print('target point is {x}'.format(x = lon_1d[lonwet[n]]))
                # x1, y1 = wgs2utm(lat_1d[int(tmp[j])], lon_1d[int(tmp[j])])
                # xw, yw = wgs2utm(lat_1d[latwet[n]], lon_1d[lonwet[n]])
                # #distance = eulid_dist(lat_1d[int(tmp[j])], lon_1d[int(tmp[j])], lat_1d[latwet[n]], lon_1d[lonwet[n]])
                # distance = eulid_dist(x1, y1, xw, yw)
                # #print('distance is {x}'.format(x = distance))
                # #print('closedj is {x}'.format(x = closedist[j]))
                
                # if distance < closedist[j]:
                    # closedist[j] = distance
        
        # comb = zip(closedist, tmp)
        # srt = sorted(comb, key=lambda k: (k[0])) ## use distance ascending order
        # # comb = zip(cohref, tmp)
        # # srt = ssorted(comb, key=lambda k: (-k[0])) ## use coh descending. 
        # print(type(srt))
        # srt = np.array(srt) 
        # tmptop = srt[:, 1]
       
        # output = tmptop[:10]
         
            
    # print('top 10 is {x}'.format(x = output))
    # latout = lat_1d[output.astype(int)]
    # lonout = lon_1d[output.astype(int)]
    # cohpt = pt1d[output.astype(int)]
    # out = np.vstack((latout, lonout, output, cohpt))
    # out = out.transpose()

    # #np.save(Path + site + '/srt.npy', srt)
    # with open(Path + site + '/' + 'output_allref' + str(wu) + '_ref.csv',"w") as f:
            # np.savetxt(f, out, delimiter=",", fmt="%f")

    # ############### all qualified RP pixels ####################
    # latout = lat_1d[tmptop.astype(int)]
    # lonout = lon_1d[tmptop.astype(int)]
    # cohpt = pt1d[tmptop.astype(int)]
    # out = np.vstack((latout, lonout, tmptop, cohpt, closedist))
    # out = out.transpose()

    # #np.save(Path + site + '/srt.npy', srt)
    # with open(Path + site + '/' + 'allref' + str(wu) + '_ref.csv',"w") as f:
            # np.savetxt(f, out, delimiter=",", fmt="%f")
    
    
    print('sub figure 2 - qualified pixels')
    ind = np.arange(0, len(lat_1d), 1) 
    out = np.vstack((ind, lat_1d, lon_1d, pt1d))
    out = out.transpose()
    with open(Path + site + '/inter_step2' + str(wu) + '.csv',"w") as f:
        np.savetxt(f, out, delimiter=",", header="step1", 
           fmt="%f")

    #lat_hw = 39.436
    #lon_hw = -122.163 ## the top wetland point in SacP1Ac1
    
    lat_hw = 39.386
    lon_hw = -122.19 ## the top wetland point in unit 0 or 1 later. 
    wind = close_point_wet(lat_hw, lon_hw, wu)
    windex, y_hw, x_hw = close_point(lat_hw, lon_hw)
    print('wind is:')
    print(wind)
            
    lat_st = 39.409
    lon_st = -122.197 ## the single point on the street. 
    sindex, y_st, x_st = close_point(lat_st, lon_st)

    # indexhw, yhw, xhw = close_point(lat_hw, lon_hw)
    # indexst, yst, xst = close_point(lat_st, lon_st)

    print('check the inclusivity of this point!!!')
    print(include[int(windex)])
    
    print('check the quality of street point!!!')
    print(quality[x_st, y_st])
    #print(targ[wind])
    filtered_array = [x for x in targ[wind] if x != -1]

    length = len(filtered_array)

    print(length)
    rg = np.where(targ[wind] != -1 ) ## only left the ones are not -1
    hw_outlat = lat_1d[np.array(targ[wind]).astype(int)][rg]
    hw_outlon = lon_1d[np.array(targ[wind]).astype(int)][rg] ## randomly pick pixel #0. 
    hw_normct = normct[int(wind)][rg]
    hw_coh = pt1d[np.array(targ[wind]).astype(int)][rg]
    
    ind = np.arange(0, len(hw_normct), 1) 
    out = np.vstack((ind, hw_outlat, hw_outlon, hw_normct, hw_coh))
    out = out.transpose()
    with open(Path + site + '/inter_step3_1_new' + '.csv',"w") as f:
        np.savetxt(f, out, delimiter=",", header="ID, lat, lon, normct", 
               fmt="%f")
    
    print('the original weltand pixel lat lon are:')
   
    print(lalo[0][0])
    print(lalo[0][1])
    
    print('the new weltand pixel lat lon are:')
    
    print(lalo[wind][0])
    print(lalo[wind][1])
    
    # ct = 0
    # numc = np.zeros(wetnum)
    # for i in range(wetnum):
        # numc[i] = len(np.where(targ[i] != -1))
    
    counts = []

    # Iterate through the first dimension (1-5)
    for i in range(wetnum):
        # Count the number of values not equal to -1 in the second dimension (1-100)
        count_not_minus_one = sum(1 for value in targ[i] if value != -1)
        counts.append(count_not_minus_one)
    
    
    with open(Path + site + '/conn_num.csv',"w") as f:
            np.savetxt(f, counts, delimiter=",", header="numc", 
                   fmt="%f")
                   
    print("figure step4_1")
    
    rg = np.where(top > 0) ## only need those top greater than 0!
    outlat = lat_1d[np.array(top).astype(int)][rg]
    outlon = lon_1d[np.array(top).astype(int)][rg]
    outcoh = pt1d[np.array(top).astype(int)][rg]
    ind = np.arange(0, len(outlat), 1) 
    out = np.vstack((ind, outlat, outlon, outcoh))
    out = out.transpose()
    with open(Path + site + '/step4_1_' + str(wu) + '.csv',"w") as f:
            np.savetxt(f, out, delimiter=",", header="ID, lat, lon, coh", 
                   fmt="%f")
    print('THE END.')
