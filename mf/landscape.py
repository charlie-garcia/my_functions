import numpy as np
from skimage.segmentation import watershed

def GetWatershed(static_deformation, rho, t, P):
    
    # Watershed (label==0 & regions)
    static_deformation = np.nan_to_num(static_deformation)  
    labels = watershed(-static_deformation, connectivity = 1, watershed_line=True)
    
    # remove out of region watershed.
    waterlines = np.int16(labels==0)
    wl_index = np.argwhere(labels == 0)
    
    # Remove false positive region
    ref = np.sum(labels==0)
    for ii in range(labels.max()):
        wtotal = np.sum( labels==ii )    
        if wtotal < ref:
            l_idx = np.nonzero(labels == ii)
            labels[l_idx] = ii+1
    
    # count number of regions
    npk = np.unique(labels)      
    npk = npk[npk != 0]
    
    # separate sections given watershed in each region.
    N = len(npk)
    segment, locs= ([0]*N for ii in range(2))
    peaks = np.zeros((N,1))
    jj=0

    # Compute landscape function
    landscape = static_deformation*rho*t/P

    for ii in npk:
        segment[jj] = labels == ii
        zone        = landscape*segment[jj]
        peaks[jj]   = np.max(zone)
        locs[jj]    = [np.where(zone == np.amax(zone))[0][0], np.where(zone == np.amax(zone))[1][0]]
        jj+=1
    return landscape, peaks, waterlines, wl_index, locs, segment
