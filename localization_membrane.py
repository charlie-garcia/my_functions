from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
from skimage.segmentation import watershed
import scipy.interpolate as inter
from mf.plots import PlotSettings

def test_fenics(V, bc, fmesh):
	# Define variational problem
	u = TrialFunction(V)
	v = TestFunction(V)
	c0 = 1
	a = inner(grad(u), grad(v))*dx

	# --> Solve landscape here
	f = Constant(1)
	L = f*v*dx

	# Compute solution
	u = Function(V)
	solve(c0*a == L, u, bc)
	# =============================================================================

	npoints = 300
	offset = 0
	nodal_values_u = u.vector()
	array_u = np.array(nodal_values_u)

	# Acces to coordinates
	n = V.dim()                                                                      
	d = fmesh.geometry().dim()                                                        
	dof_coordinates = V.tabulate_dof_coordinates().reshape(n,d)

	xx = dof_coordinates[:,0]
	yy = dof_coordinates[:,1]

	[X, Y] = np.meshgrid(np.linspace(np.min(xx)-offset, np.max(xx) + offset, npoints),np.linspace(np.min(yy)-offset, np.max(yy)+offset,npoints));
	landscape0 = inter.griddata(dof_coordinates,array_u,(X,Y), method = 'cubic')

	landscape = np.nan_to_num(landscape0)  
	# =============================================================================
	labels = watershed(-landscape, connectivity = 1, watershed_line=True)

	# remove out of region watershed.
	# It has to be a better way !!!!!!!!!
	Nx = np.shape(X)[0]
	Ny = np.shape(X)[1]
	ROI = np.zeros([Nx,Ny])
	tol = 1e-16
	for ix in range(0, Nx):
	    for iy in range(0, Ny):
	        if landscape[ix,iy] >=0-tol:
	            ROI[ix,iy] = 1                # this could work also with the wtaershed algorithm

	labels = labels*ROI
	# # reduce order if multiples regios are found
	min_value_loc = np.where(np.unique(labels>0))[0]
	min_value = np.unique(labels)[min_value_loc][0]
	max_value = np.max(labels)
    
	labels = labels - min_value + 1
	labels = np.where(labels<0, 0, labels)

	waterlines = labels==0
	waterlines = waterlines*ROI
	wl_index = np.argwhere(labels == 0)
	
	# separate sections given watershed in each region.
	npk = np.unique(labels)
	npk = npk[npk != 0]                  #npk(npk==0) = []
	N = int(np.max(npk))
	zone= [0]*N
	peaks = np.zeros([1,N], 'double')
	locs = [0]*N
	for ii in range(0,N):
	    segment = labels == ii+1
	    zone[ii] = landscape*segment;
	    peaks[0,ii] = np.max(zone[ii])
	    locs[ii]  = [np.where(zone[ii] == np.amax(zone[ii]))[0][0], np.where(zone[ii] == np.amax(zone[ii]))[1][0]]

	fig, ax = plt.subplots(1)
	plt.pcolormesh(X, Y, landscape0,shading='auto')
	#plot(fmesh)

	# Plot the surface
	PlotSettings(fig,ax)
	ax.set_aspect('equal', 'box')
	plt.scatter(X[ wl_index[:,0],wl_index[:,1] ], Y[ wl_index[:,0],wl_index[:,1] ], waterlines[ wl_index[:,0],wl_index[:,1] ], color='white')
	plt.title('Landscape and watershed lines', color ='w')

	for ii in range(0,N):
	    plt.scatter(X[locs[ii][0], locs[ii][1]], Y[locs[ii][0], locs[ii][1]], s=20, color='red', marker ='x')

	max_landscape= np.max(landscape)

	f = ((4/np.pi)/np.sqrt(peaks))/(2*np.pi)
	f = np.sort(f)

	for ii in range(0,N):
# 		print('f0 from landscape: {:.2f} '.format(f[0,ii]))
		print('f0 from landscape %.2f Hz' %f[0,ii])
