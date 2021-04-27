import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inter
# from fenics import *

def set_fontsize(N):
	plt.rcParams.update({'font.size': N})

def PlotSettings(fig,ax):
    from mf.math import check_array

    fig.set_facecolor('#19232d')
        
    if check_array(ax) is True:
        cc = get_axes_coord(ax)
        idx=0
        nlen = ax.size
        for ii in range(0, nlen):
            ax[cc[idx]].set_facecolor('#19232d')
            ax[cc[idx]].tick_params(color='#C0C0C0', labelcolor='w')
            # We change the fontsize of minor ticks label 
            ax[cc[idx]].tick_params(axis='both', which='major', labelsize=7)
            ax[cc[idx]].tick_params(axis='both', which='minor', labelsize=6)
            ax[cc[idx]].xaxis.label.set_color('w')
            ax[cc[idx]].yaxis.label.set_color('w')
            if hasattr(ax[cc[idx]], 'get_zlim'): 
                ax[cc[idx]].zaxis.label.set_color('w')

            ax[cc[idx]].title.set_color('w')
            # set imshow outline
            for spine in ax[cc[idx]].spines.values():
                spine.set_edgecolor('w') 

            idx=idx+1

    elif isinstance(ax, list) is True:
        nlen = len(ax)
        for ii in range(0, nlen):
            ax[ii].set_facecolor('#19232d')
            ax[ii].tick_params(color='#C0C0C0', labelcolor='w')
            # We change the fontsize of minor ticks label 
            ax[ii].tick_params(axis='both', which='major', labelsize=7)
            ax[ii].tick_params(axis='both', which='minor', labelsize=6)
            ax[ii].xaxis.label.set_color('w')
            ax[ii].yaxis.label.set_color('w')
            if hasattr(ax[ii], 'get_zlim'): 
                ax[ii].zaxis.label.set_color('w')
            ax[ii].title.set_color('w')
            for spine in ax[ii].spines.values():
                spine.set_edgecolor('w') 
    else:        
        ax.set_facecolor('#19232d')
        ax.tick_params(color='#C0C0C0', labelcolor='w')
        # We change the fontsize of minor ticks label 
        set_axis_parameters(ax)

def set_axis_parameters(ax):
    ax.tick_params(axis='both', which='major', labelsize=7)
    ax.tick_params(axis='both', which='minor', labelsize=6)
    ax.xaxis.label.set_color('w')
    ax.yaxis.label.set_color('w')
    if hasattr(ax, 'get_zlim'): 
            ax.zaxis.label.set_color('w')
    ax.title.set_color('w')
    for spine in ax.spines.values():
        spine.set_edgecolor('w') 
        
def PlotSettingsSmall(fig,ax):
    from mf.math import check_array
        
    if check_array(ax) is True:
        cc = get_axes_coord(ax)
        idx=0
        nlen = ax.size
        for ii in range(0, nlen):
            ax[cc[idx]].tick_params(color='#262626')
            ax[cc[idx]].tick_params(axis='both', which='major', labelsize=7)
            ax[cc[idx]].tick_params(axis='both', which='minor', labelsize=7)
            idx=idx+1

    elif isinstance(ax, list) is True:
        nlen = len(ax)
        for ii in range(0, nlen):
            ax[ii].tick_params(color='#262626')
            # We change the fontsize of minor ticks label 
            ax[ii].tick_params(axis='both', which='major', labelsize=7)
            ax[ii].tick_params(axis='both', which='minor', labelsize=7)
            
    else:        
        ax.tick_params(color='#262626')
        # We change the fontsize of minor ticks label 
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.tick_params(axis='both', which='minor', labelsize=7)

def PlotSettingsSize(fig,ax, set_size):
    from mf.math import check_array
        
    if check_array(ax) is True:
        cc = get_axes_coord(ax)
        idx=0
        nlen = ax.size
        for ii in range(0, nlen):
            ax[cc[idx]].tick_params(color='#262626')
            ax[cc[idx]].tick_params(axis='both', which='major', labelsize=set_size)
            ax[cc[idx]].tick_params(axis='both', which='minor', labelsize=set_size)
            ax[cc[idx]].xaxis.label.set_size(set_size)
            ax[cc[idx]].yaxis.label.set_size(set_size)
            idx=idx+1

    elif isinstance(ax, list) is True:
        nlen = len(ax)
        for ii in range(0, nlen):
            ax[ii].tick_params(color='#262626')
            # We change the fontsize of minor ticks label 
            ax[ii].tick_params(axis='both', which='major', labelsize=set_size)
            ax[ii].tick_params(axis='both', which='minor', labelsize=set_size)
            ax[ii].xaxis.label.set_size(set_size)
            ax[ii].yaxis.label.set_size(set_size)
            
    else:        
        ax.tick_params(color='#262626')
        # We change the fontsize of minor ticks label 
        ax.tick_params(axis='both', which='major', labelsize=set_size)
        ax.tick_params(axis='both', which='minor', labelsize=set_size)
        ax.xaxis.label.set_size(set_size)
        ax.yaxis.label.set_size(set_size)
        

def PolarPlotZoom(fig, loc, theta, Pfield, k, set_scale, label):    
    import numpy as np
    import mpl_toolkits.axisartist.angle_helper as angle_helper
    from matplotlib.projections import PolarAxes
    from matplotlib.transforms import Affine2D
    from mpl_toolkits.axisartist import SubplotHost, ParasiteAxesAuxTrans, GridHelperCurveLinear
    
    # see demo_curvelinear_grid.py for details
    # tr = Affine2D().scale(np.pi / 180., 1.) + PolarAxes.PolarTransform()
    tr = Affine2D().scale(np.pi / 180., 1)  + PolarAxes.PolarTransform()
    
    extreme_finder = angle_helper.ExtremeFinderCycle(5,                         # On x axis
                                                     5,                         # On y axis
                                                     lon_cycle=360,             # Arc 
                                                     lat_cycle=None,
                                                     lon_minmax=None,
                                                     lat_minmax=(0,np.inf),
                                                     )
    
    grid_locator1 = angle_helper.LocatorDMS(12)
    tick_formatter1 = angle_helper.FormatterDMS()
    grid_helper = GridHelperCurveLinear(tr,
                                        extreme_finder=extreme_finder,
                                        grid_locator1=grid_locator1,
                                        tick_formatter1=tick_formatter1
                                        )
    
    n_row, n_col, axes_id = loc[0], loc[1], loc[2]
    ax1 = SubplotHost(fig, n_row, n_col, axes_id, grid_helper=grid_helper)
    
    fig.add_subplot(ax1)
    
    # A parasite axes with given transform
    ax2 = ParasiteAxesAuxTrans(ax1, tr)#, "auto")
    ax1.parasites.append(ax2)
    
    ax2.plot(theta*90/(np.pi/2), Pfield, linewidth=1.0, label=label)
    if set_scale == 'fixed':
        ax1.set_aspect(1)
    ax1.set_xlim(0, 1)
    if k<0.1:
        my_ylim = 1
    elif k>=0.1 and k<=1:
        my_ylim = -0.83*k + 1.08
    elif k>1:
        my_ylim = 0.25    
    
    ax1.set_ylim(-my_ylim, my_ylim)

    ax1.grid(True)
    return ax1, ax2
    
def ColorbarSettings(ax, cb, colorbar_label):
	cb.set_label(colorbar_label, color='w')            			   # set colorbar label plus label color
	cb.ax.yaxis.set_tick_params(color='w')                         # set colorbar tick color
	cb.outline.set_edgecolor('w')                                  # set colorbar edgecolor 
	plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color='w')       # set colorbar ticklabels


def get_axes_coord(ax):
    ii=0
    Ni, Nj = ax.shape
    # cc = np.zeros((Ni*Nj, 1))
    cc = [0]*Ni*Nj
    for ix in range(0,Ni):
        for iy in range(0,Nj):
            cc[ii] = (ix,iy)
            ii=ii+1
    return cc

def connect_squares2D_gui(x,y,N,idx,connect, plotInfo, ax):
    cx, cy = ([0]*N**2 for i in range(2))
    for ii in range(0,connect.shape[0]):
        cx[ii] = np.mean(x[0,idx[connect[ii,:]]])
        cy[ii] = np.mean(y[0,idx[connect[ii,:]]])
        
        if plotInfo in 'yes':
            ax.plot( [x[0, idx[connect[ii,0]] ], x[0, idx[connect[ii,1]] ] ], [y[0, idx[connect[ii,0]] ], y[0, idx[connect[ii,1]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            ax.plot( [x[0, idx[connect[ii,1]] ], x[0, idx[connect[ii,2]] ] ], [y[0, idx[connect[ii,1]] ], y[0, idx[connect[ii,2]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            ax.plot( [x[0, idx[connect[ii,2]] ], x[0, idx[connect[ii,3]] ] ], [y[0, idx[connect[ii,2]] ], y[0, idx[connect[ii,3]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            ax.plot( [x[0, idx[connect[ii,3]] ], x[0, idx[connect[ii,0]] ] ], [y[0, idx[connect[ii,3]] ], y[0, idx[connect[ii,0]] ] ], ':', color=[.6,.6,.6], lw=0.5)
        
        cs = np.c_[np.array(cx),np.array(cy)]
    return cs

def create_plate2D_gui(Lx,Ly,N,plot_info, ax):
    x = np.linspace(0,Lx,N+1)    # 0 : dx : Lx ;              
    y = np.linspace(0,Ly,N+1)    # 0 : dy : Ly ;
    
    X,Y = np.meshgrid(x,y)
    xx = X.reshape(1, (N+1)**2)
    yy = Y.reshape(1, (N+1)**2)
    
    idx = np.array(np.r_[:xx.shape[1]],dtype='int')
    connect =np.zeros([4,N**2], dtype='int')
    
    for ii in range(N):
        for jj in range(4):
            connect[0,N*ii:N*(ii+1)] =  np.r_[(ii+0) * N+1 + ii: (ii+1) * N+1 + ii]
            connect[1,N*ii:N*(ii+1)] =  np.r_[(ii+0) * N+2 + ii: (ii+1) * N+2 + ii]
            connect[2,N*ii:N*(ii+1)] =  np.r_[(ii+1) * N+3 + ii: (ii+2) * N+3 + ii]
            connect[3,N*ii:N*(ii+1)] =  np.r_[(ii+1) * N+2 + ii: (ii+2) * N+2 + ii]
    
    connect = connect.T - 1
    cs = connect_squares2D_gui(xx, yy, N,idx, connect, plot_info, ax)
    
    return cs, x, y

def connect_squares2D(x,y,N,idx,connect, plotInfo):
    cx, cy = ([0]*N**2 for i in range(2))
    for ii in range(0,connect.shape[0]):
        cx[ii] = np.mean(x[0,idx[connect[ii,:]]])
        cy[ii] = np.mean(y[0,idx[connect[ii,:]]])
        
        if plotInfo in 'yes':
            plt.plot( [x[0, idx[connect[ii,0]] ], x[0, idx[connect[ii,1]] ] ], [y[0, idx[connect[ii,0]] ], y[0, idx[connect[ii,1]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            plt.plot( [x[0, idx[connect[ii,1]] ], x[0, idx[connect[ii,2]] ] ], [y[0, idx[connect[ii,1]] ], y[0, idx[connect[ii,2]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            plt.plot( [x[0, idx[connect[ii,2]] ], x[0, idx[connect[ii,3]] ] ], [y[0, idx[connect[ii,2]] ], y[0, idx[connect[ii,3]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            plt.plot( [x[0, idx[connect[ii,3]] ], x[0, idx[connect[ii,0]] ] ], [y[0, idx[connect[ii,3]] ], y[0, idx[connect[ii,0]] ] ], ':', color=[.6,.6,.6], lw=0.5)
        
        cs = np.c_[np.array(cx),np.array(cy)]
    return cs

def create_plate2D(Lx,Ly,N,plot_info):
    x = np.linspace(0,Lx,N+1)    # 0 : dx : Lx ;              
    y = np.linspace(0,Ly,N+1)    # 0 : dy : Ly ;
    
    X,Y = np.meshgrid(x,y)
    xx = X.reshape(1, (N+1)**2)
    yy = Y.reshape(1, (N+1)**2)
    
    idx = np.array(np.r_[:xx.shape[1]],dtype='int')
    connect =np.zeros([4,N**2], dtype='int')
    
    for ii in range(N):
        for jj in range(4):
            connect[0,N*ii:N*(ii+1)] =  np.r_[(ii+0) * N+1 + ii: (ii+1) * N+1 + ii]
            connect[1,N*ii:N*(ii+1)] =  np.r_[(ii+0) * N+2 + ii: (ii+1) * N+2 + ii]
            connect[2,N*ii:N*(ii+1)] =  np.r_[(ii+1) * N+3 + ii: (ii+2) * N+3 + ii]
            connect[3,N*ii:N*(ii+1)] =  np.r_[(ii+1) * N+2 + ii: (ii+2) * N+2 + ii]
    
    connect = connect.T - 1
    cs = connect_squares2D(xx, yy, N,idx, connect, plot_info)
    
    return cs, x, y

def areaTriangle(x,y):
    area = 1/2 * np.abs( x[0] * ( y[1]-y[2] ) + x[1] * ( y[2]-y[0] ) + x[2] * ( y[0]-y[1] ) )
    return area


def getCenters(x,y,z):
    lv = x.shape[0]
    cx, cy, cz,area = (np.zeros((lv,1)) for i in range(4))
    for ii in range(0,lv):
        cx[ii] = np.mean(x[ii,:])
        cy[ii] = np.mean(y[ii,:])
        cz[ii] = np.mean(z[ii,:])
        area[ii] = areaTriangle(x[ii,:], y[ii,:])
        plt.plot(cx[ii], cy[ii], 'x')
    cs = np.c_[cx,cy,cz]
    
    return cs, area

def CreatePiston(radius, N):
    phi, cx, cy = (np.zeros((4*N**2,N)) for i in range(3))
    for ii in range(1,N+1):
        r = radius*(2*ii-1)/(2*N)
        for jj in range(1, 4*(2*ii-1)+1):
            phi[jj-1,ii-1] = np.pi*(2*jj-1)/(4*(2*ii-1))
            cx[jj-1,ii-1] = r*np.cos(phi[jj-1,ii-1])
            cy[jj-1,ii-1] = r*np.sin(phi[jj-1,ii-1])
            
    idx = phi!=0
    
    cx = cx.T[idx.T]
    cy = cy.T[idx.T]
    cs = np.vstack((cx, cy, 0*cx))  
    cs = cs.T
    area = np.ones((4*N**2,1))*(radius/N)**2*np.pi/4
    return cs, area    

def plotFig(filename,fignr=1):
    from scipy.io import loadmat
    from numpy import size
    from matplotlib.pyplot import plot,figure,xlabel,ylabel,show,clf,xlim,legend
    d = loadmat(filename,squeeze_me=True, struct_as_record=False)
    matfig = d['hgS_070000']
    childs = matfig.children
    ax1 = [c for c in childs if c.type == 'axes']
    if(len(ax1) > 0):
        ax1 = ax1[0]
    legs = [c for c in childs if c.type == 'scribe.legend']
    if(len(legs) > 0):
        legs = legs[0]
    else:
        legs=0
    pos = matfig.properties.Position
    size = np.array([pos[2]-pos[0],pos[3]-pos[1]])/96
    plt.figure(fignr,figsize=size)
    plt.clf()
    # plt.hold(True)
    counter = 0    
    for line in ax1.children:
        if line.type == 'graph2d.lineseries':
            if hasattr(line.properties,'Marker'):
                mark = "%s" % line.properties.Marker
                if(mark != "none"):
                    mark = mark[0]
            else:
                mark = '.'
            if hasattr(line.properties,'LineStyle'):
                linestyle = "%s" % line.properties.LineStyle
            else:
                linestyle = '-'
            if hasattr(line.properties,'Color'):
                r,g,b =  line.properties.Color
            else:
                r = 0
                g = 0
                b = 1
            if hasattr(line.properties,'MarkerSize'):
                marker_size = line.properties.MarkerSize
            else:
                marker_size = -1                
            x = line.properties.XData
            y = line.properties.YData
            if(mark == "none"):
                plt.plot(x,y,linestyle=linestyle,color=[r,g,b])
            elif(marker_size==-1):
                plt.plot(x,y,marker=mark,linestyle=linestyle,color=[r,g,b])
            else:
                plt.plot(x,y,marker=mark,linestyle=linestyle,color=[r,g,b],ms=marker_size)
        elif line.type == 'text':
            if counter == 0:
                plt.xlabel("$%s$" % line.properties.String,fontsize =16)
            elif counter == 1:
                plt.ylabel("$%s$" % line.properties.String,fontsize = 16)
            elif counter == 3:
                plt.title("$%s$" % line.properties.String,fontsize = 16)
            counter += 1        
    # plt.grid(ax1.properties.XGrid)

    if(hasattr(ax1.properties,'XTick')):
        if(hasattr(ax1.properties,'XTickLabelRotation')):
            plt.xticks(ax1.properties.XTick,ax1.properties.XTickLabel,rotation=ax1.properties.XTickLabelRotation)
        else:
            plt.xticks(ax1.properties.XTick,ax1.properties.XTickLabel)
    if(hasattr(ax1.properties,'YTick')):
        if(hasattr(ax1.properties,'YTickLabelRotation')):
            plt.yticks(ax1.properties.YTick,ax1.properties.YTickLabel,rotation=ax1.properties.YTickLabelRotation)
        else:
            plt.yticks(ax1.properties.YTick,ax1.properties.YTickLabel)
    plt.xlim(ax1.properties.XLim)
    plt.ylim(ax1.properties.YLim)
    if legs:        
        leg_entries = tuple(['$' + l + '$' for l in legs.properties.String])
        py_locs = ['upper center','lower center','right','left','upper right','upper left','lower right','lower left','best','best']
        MAT_locs=['North','South','East','West','NorthEast', 'NorthWest', 'SouthEast', 'SouthWest','Best','none']
        Mat2py = dict(zip(MAT_locs,py_locs))
        location = legs.properties.Location
        plt.legend(leg_entries)#,loc=Mat2py[location])
    # plt.hold(False)
    plt.show()