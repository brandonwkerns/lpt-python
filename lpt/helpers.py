import numpy as np
from scipy.signal import convolve2d

## These functions are used for LPT.

############################################################################
######################  Gaussian Smoothing Functions  ######################
############################################################################

def gauss_smooth_kernel(nx,ny=-999,stdx=-999,stdy=-999):

    """
    ## Originally written in Matlab.
    %
    % F = gaussSmoothKernel(nx,ny,stdx,stdy)
    %    -- or --
    % F = gaussSmoothKernel(stdev)
    %
    % Returns 2-D Gaussian smoothing kernel
    % With x and y lengths specified.
    % along with standard deviation on in the x and y direction.
    % (standard deviations are in terms of gridpoints.)
    %
    % !!!! nx and ny must be odd, or else 1 will be added to
    % them. !!!!
    %
    % New simplified option: only specify stdx (which is used in x
    % and y directions).  It is also in terms of gridpoints.
    %
    % ----------------------------------------------
    % 2010.11.11  BK  bkerns@rsmas.miami.edu
    % 2012.11.12  BK  added simplified option to only specify the stdx
    ## Python port: 1.28.2019
    """

    ## Process the input, if needed.
    if ny < -900:
        stdx=nx
        stdy=nx
        nx = 6*stdx + 1
        ny = 6*stdy + 1

    if int(nx) / 2 == 0:
        print('Warning: nx should be odd! Adding one to it.')
        nx += 1

    if int(ny) / 2 == 0:
        print('Warning: ny should be odd! Adding one to it.')
        ny += 1

    if stdx < -900 or stdy < -900:
        stdx = (nx - 1) / 6
        stdy = (ny - 1) / 6

    ## Now, do the calculation.
    F = np.ones([nx,ny])

    halfX = (nx - 1) / 2
    x = np.arange(-1*halfX, halfX + 1)
    halfY = (ny - 1) / 2
    y = np.arange(-1*halfY, halfY + 1)
    X, Y = np.meshgrid(x,y)

    F = np.exp(-1.0*(np.power(X,2))/(2*stdx**2) + -1.0*(np.power(Y,2))/(2*stdy**2))

    ## Noramlize
    F /= np.sum(F)

    print((stdx,stdx))
    print(np.power(X,2))
    print(X)
    print(Y)
    return F


def gauss_smooth(data,nx,ny=-999,stdx=-999,stdy=-999):

    """
    ## Originally written in Matlab.
    % dataSmooth = gaussSmooth(data,nx,ny,stdx,stdy)
    % -- OR --
    % dataSmooth = gaussSmooth(data,stdev)
    %
    % Data is a 2-D array
    % The Gaussian kernel is defined by: nx,ny,stdx,stdy
    %  using the function gaussSmoothKernel
    %    nx = full width in x direction
    %    ny = full width in y direction
    %    stdx = standard deviation in x direction (gridpoints)
    %    stdy = standard deviation in y direction (gridpoints)
    %
    %  NOTE:  Function gaussSmooth requires that nx and ny be odd.
    %         If you give even numbers, 1 is added to make them odd.
    %
    %  ALSO: dataSmooth has the same dimensions as data. The way this works is
    %        Matlab's built-in 'conv2' is called with the SHAPE option set to 'same'
    %
    % UPDATE 2012.11.12 you can call this in the simplified way:
    %   dataSmooth = gaussSmooth(data,stdev)
    %  and it will use stdev in both x and y direction, extending the
    %  kernel out to 3 stdev's.
    ## Python port: 1.28.2019
    """
    
    kernel = gauss_smooth_kernel(nx,ny,stdx,stdy)
    data_smooth = convolve2d(data, kernel, 'same')
    return data_smooth
