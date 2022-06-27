import numpy as np
import matplotlib.pyplot as plt


### The following are outward facing functions which may be called

### The main function of the package, returns the RJ-values of a 2D image (2d-array) as well as the classification of the structure.
def get_rj(image):

    J1,J2,com,theta = moments(image)

    RJ1 = (J1 - J2)/np.sqrt(2)
    RJ2 = (J1 + J2)/np.sqrt(2)
    RJ2 = RJ2/(-RJ1 + np.sqrt(2))

    cla = classification(RJ1,RJ2)

    return RJ1,RJ2,cla

### Uses the classification scheme described in the accompanying paper. 
### 1 - Quasi-circular, centrally over-dense
### 2 - Quasi-circular, centrally under-dense
### 3 - Elongated, centrally over-dense
### 4 - Elongated, centrally under-dense
def classification(RJ1,RJ2):

    z = np.poly1d([-0.31717289,  0.00112143,  0.47870963])
    cla = -2

    if(RJ1>0.47870963):
        if(RJ2>0):
            cla = 3
        else:
            cla = 4

    else:
        if(RJ1 < z(RJ2) and RJ1>-z(RJ2)):
            if(RJ2>0):
                cla = 1
            else:
                cla = 2

        else:
            if(RJ2>0):
                cla = 3
            else:
                cla = 4

    if(cla==-2):
        print("Something went wrong with the classification",RJ1,RJ2)

    return cla


### Creates and saves an RJ-plot given an array of RJ-values. The section which produces the axes, sets the limits and displays the colour-coded quadrants in highlighted so that it may copied to produce bespoke RJ-plots.
def make_rj_plot(RJ1,RJ2,plotname,dpi=100,markersize=1,markercolor="C0"):

    fig = plt.figure(1,figsize=(10,10))
    plt.plot(RJ1,RJ2,marker="o",color=markercolor,ls="None",ms=markersize)

    ### This following section of the code can be taken and copied to produce the axis lines and the four colour-coded quadrants
    ### From here:
    z = np.poly1d([-0.31717289,  0.00112143,  0.47870963])
    ans = 1.2303055886801646
    yy = np.linspace(-ans,ans,100)
    plt.plot(z(yy),yy,"r",lw=1.5)

    xx = np.linspace(-0.05,np.sqrt(2),100)
    plt.plot(xx,np.zeros_like(xx),"k",lw=1.5)
    xx = np.linspace(-np.sqrt(2),np.sqrt(2),100)
    plt.plot(np.zeros_like(xx),xx,"k",lw=1.5)

    plt.xlim(-0.05,np.sqrt(2))
    plt.ylim(-1,1)

    y2 = np.linspace(0,ans,100)
    plt.fill_between(z(y2),y2,alpha=0.15, color="y")
    plt.fill_between(z(y2),-y2,0,alpha=0.15, color="b")

    xx = np.linspace(np.amax(z(y2)),np.sqrt(2),100)
    q = np.zeros_like(xx)

    x2 = z(y2)[::-1]
    y2 = y2[::-1]

    x1 = np.concatenate((x2,xx))
    y1 = np.concatenate((y2,q))

    plt.fill_between(x1,y1,np.sqrt(2),alpha=0.15,color="g")
    plt.fill_between(x1,-y1,-np.sqrt(2),alpha=0.15,color="m")

    plt.xlabel("R$_1$")
    plt.ylabel("R$_2$")
    ### To here.

    plt.savefig(plotname,dpi=dpi)



### This function is from J-plots and will return the two J-values, the centre of weight of a structure, as well as its position angle.
def moments(g):
    com=find_com(g)
    Atot=np.count_nonzero(g)
    Mtot=np.sum(g)
    I1,I2,t1=i12(g,com)
    first=min([I1,I2])
    second=max([I1,I2])
    ma=Atot*Mtot
    II1=(ma - 4*np.pi*first)/(ma + 4*np.pi*first)
    II2=(ma - 4*np.pi*second)/(ma + 4*np.pi*second)

    return II1,II2,com,t1
    






### The following are functions used internally and are not necessarily meant to be called outside of RJ-plots. As such they are commented.

def find_com(grid):
    com = np.array([0,0],dtype=float)
    nx=grid.shape[0]
    ny=grid.shape[1]

    ax = np.linspace(0,nx-1,nx)
    ay = np.linspace(0,ny-1,ny)

    com[0] = np.sum(np.sum(grid,axis=1)*ax)
    com[1] = np.sum(np.sum(grid,axis=0)*ay)

    com=com/np.sum(grid)

    return com
    
def m11(grid,com):
    nx=grid.shape[0]
    ny=grid.shape[1]
    dr11=0

    ay = np.linspace(0,ny-1,ny)
    dr11 = np.sum(np.sum(grid,axis=0)*(ay-com[1])**2)

    return dr11
    
def m22(grid,com):
    nx=grid.shape[0]
    ny=grid.shape[1]
    dr22=0

    ax = np.linspace(0,nx-1,nx)
    dr22 = np.sum(np.sum(grid,axis=1)*(ax-com[0])**2)

    return dr22
    
def m12(grid,com):
    nx=grid.shape[0]
    ny=grid.shape[1]
    dr12=0

    ax = np.linspace(0,nx-1,nx) - com[0]
    ay = np.linspace(0,ny-1,ny) - com[1]

    ax2 = np.array([ax,]*ny).T
    ay2 = np.array([ay,]*nx)

    dr12 = np.sum(grid*ay2*ax2)

    return dr12
    
def theta(M11,M22,M12):

    th=0.5*np.arctan2((2*M12),(M11-M22))

    return th
    
def i12(grid,com):
    M11=m11(grid,com)
    M22=m22(grid,com)
    M12=m12(grid,com)
    tt=theta(M11,M22,M12)

    i1=(M11*np.cos(tt)*np.cos(tt)) + \
        (M22*np.sin(tt)*np.sin(tt)) + \
        (M12*np.sin(2*tt))

    i2=(M11*np.sin(tt)*np.sin(tt)) + \
        (M22*np.cos(tt)*np.cos(tt)) - \
        (M12*np.sin(2*tt))

    return i1,i2,tt
 