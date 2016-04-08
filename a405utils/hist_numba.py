import numba
import numpy as np

# nopython=True means an error will be raised
# if fast compilation is not possible.
@numba.jit(nopython=True)
def fill_counts(y_centers,x_centers,y_indices,x_indices):
    """
    Use numba to to count datapoints in a 2d histogram grid
 
    parameters
    ----------

    y_centers: vector (float)
        bin centers for row variable

    x_centers: vector (float)
       bin centers for column variable

    y_indices: vector (float)
       vector of bin number for each data point's y value


    x_indices: vector (float)
       vector of bin number for each data point's x value
    
    Returns
    -------

    hist_array: 2d array of shape (len(y_centers),len(x_centers))
       counts of data falling into each row(y), col(x)  bin
    """
    num_xbins=int(x_centers.shape[0])
    num_ybins=int(y_centers.shape[0])
    num_y=y_indices.shape[0]   #number of x and y data points
    hist_array=np.zeros((num_ybins,num_xbins),dtype=np.float32)
    for n in range(num_y): #y_indices and x_indices both size of raw data
        if x_indices[n] > 0 and y_indices[n] > 0 and \
            x_indices[n] <= num_xbins and y_indices[n] <= num_ybins:
            bin_row=y_indices[n]-1 # '-1' to get the index of the bin center
            bin_col=x_indices[n]-1
            hist_array[bin_row, bin_col] += 1
    rows,cols=hist_array.shape
    for row in range(rows):
        for col in range(cols):
            if hist_array[row,col] < 1.:
                hist_array[row,col]=np.nan
    return hist_array

            
def numba_hist2d(x_raw,y_raw,x_edges,y_edges):
    """
    produce a 2-d historgram of x (columns) and y (rows)
    """
    print('in numba')
    x_centers=(x_edges[:-1] + x_edges[1:])/2.
    y_centers=(y_edges[:-1] + y_edges[1:])/2.
    x_indices=np.asarray(np.searchsorted(x_edges, x_raw.flat, 'right'),dtype=np.int64)
    y_indices=np.asarray(np.searchsorted(y_edges, y_raw.flat, 'right'),dtype=np.int64)
    hist_array=fill_counts(y_centers,x_centers,y_indices,x_indices)
    return hist_array,x_centers,y_centers

