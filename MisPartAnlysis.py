import numpy as np
import matplotlib.pyplot as plt
import CalcMSD as MSD

def split(data,cutoff,wanted = [1,1,1,1],showloss = True):
#wanted: 0: List of Trajectories as lists of particle positions; 1:trajectory length info, minimum, maximum, average;
#   2: List of Frame numbers for each particle in each trajectory; 3:min/max x and y position throughout all frames

    for i in xrange(4-len(wanted)):
        wanted.append(0)
    Traj = []
    ntraj = int(max(data[:,3]))
    mintrajlen = 100
    maxtrajlen = 0
    ave = 0
    count = 0
    Len = []
    Frame = []
    minx, miny = 0 , 0
    maxx, maxy = 0, 0
    for i in range(ntraj):
        ind = np.where(data[:,3] == i+1)
        traj = data[ind,0:2][0]    
        Len.append(len(traj))
        if(len(traj)<=cutoff): count = count+1
        if(len(traj)>cutoff):
            if(wanted[0] == 1 or wanted[1] == 1):
                Traj.append(traj)
            if(wanted[2] == 1):
                Frame.append(data[ind,2][0])
            if(wanted[1] == 1):
                mintrajlen = min(mintrajlen,len(traj))
                maxtrajlen = max(maxtrajlen,len(traj))
                ave = ave + len(traj)
            if(wanted[3] == 1):
                mintraj = np.min(traj,axis=0)
                if(minx > mintraj[0]): minx = mintraj[0]
                if(miny > mintraj[1]): miny = mintraj[1]
                maxtraj = np.max(traj,axis=0)
                if(maxx < maxtraj[0]): maxx = maxtraj[0]
                if(maxy < maxtraj[1]): maxy = maxtraj[1]
    if(showloss): print( count/ntraj,'lost')

    output = []
    if(wanted[0] == 1): output.append(Traj)
    if(wanted[1] == 1):
        L = len(Traj)
        ave = ave/L
        ave = mintrajlen
        output = output + [mintrajlen,maxtrajlen,ave]
    if(wanted[2] == 1): output.append(Frame)
    if(wanted[3] == 1): output = output + [[minx,maxx],[miny,maxy]]
    return output

def Dense(data): #returns number of particles per frame
    frames = data[:,2]
    Npart = []
    for i in range(int(max(frames))):
        Npart.append(len(np.where(frames == i)[0]))
    return Npart
        
def Crowd(Frame,Npart,Traj): #returns average number of particles of trajectory
    Npart = np.array(Npart)
    crowd = []
    for i in range(len(Traj)):
        crowd.append(np.ave(Frame.pop()))
    return crowd

def vor_density(Loc): #density given by voronoi tesselations
    from scipy.spatial import Voronoi, voronoi_plot_2d
    vor = Voronoi(Loc)
    edges = np.array(vor.regions)[vor.point_region]
    vert = vor.vertices
    Dens = []
    ind_reg = -1
    for edg in edges:
     if(-1 in edg):
        Dens.append(-1)
     elif(len(edg) > 1):
        ind_reg = ind_reg+1
        pt = []
        if(-1 in edg):
            edg.remove(-1)
        l = len(edg)
        for cnt in xrange(l):
            v = edg[cnt]
            pt.append(vert[v])
        pt.sort(key=lambda pos: pos[0])
        pt = np.array(pt)
        y = np.array(pt[:,1])
        x = np.array(pt[:,0])
        m = (y[-1]-y[0])/(x[-1]-x[0])
        b = y[0]-m*x[0]
        upind = np.where(y+0.0001>=x*m+b)[0]
        downind = np.where(y-0.0001<= x*m+b)[0]
        x2 = x[upind]
        y2 = y[upind]
        y3 = x2*m+b
        area = np.sum((x2[1:]-x2[:-1])*((y2[:-1]-y3[:-1])-(y2[:-1]-y2[1:])/2-(y3[1:]-y3[:-1])/2))
        x2 = x[downind]
        y2 = y[downind]
        y3 = x2*m+b
        area = area+np.sum((x2[:-1]-x2[1:])*((y2[:-1]-y3[:-1])-(y2[:-1]-y2[1:])/2-(y3[1:]-y3[:-1])/2))
        Dens.append(1/area)
     else: Dens.append(0)
    return Dens
    
def find_poly(Loc): #identify polygon that encloses datapoints, Loc=array([[x_1,y_1],[x_2,y_2]...])
    def slope_int(p1,p2):
        if(p1[0] == p2[0]):
            p = p1[0]+0.0000001
            m = (p1[1]-p2[1])/(p-p2[0])
            b = p1[1]-m*p
            print "warning: p1=p2"
        else:
            m = (p1[1]-p2[1])/(p1[0]-p2[0])
            b = p1[1]-m*p1[0]
        return [m,b]
    ind1 = np.argmax(Loc[:,0])
    ind2 = np.argmax(Loc[:,1])
    ind3 = np.argmin(Loc[:,0])
    ind4 = np.argmin(Loc[:,1])
    topx = Loc[ind2,0]
    botx = Loc[ind4,0]
    Upedge = [slope_int(Loc[ind1],Loc[ind2]),slope_int(Loc[ind3],Loc[ind2])]
    Uppoint = [Loc[ind1],Loc[ind2],Loc[ind3]]
    Downedge = [slope_int(list(Loc[ind1]),Loc[ind4]),slope_int(Loc[ind3],Loc[ind4])]
    Downpoint = [Loc[ind1],Loc[ind4],Loc[ind3]]
    if(ind2 == ind3):
        Upedge.pop(1)
        Uppoint.pop(1)
    if(ind1 == ind2):
        Upedge.pop(0)
        Uppoint.pop(0)
    if(ind3 == ind4):
        Downedge.pop(1)
        Downpoint.pop(1)
    if(ind1 == ind4):
        Downedge.pop(0)
        Downpoint.pop(0)
    X = Loc[:,0]
    Y = Loc[:,1]
    lenedg = len(Upedge)
    i = 0
    while(i<lenedg):
        edg = Upedge[i]
        ind = np.argmax(Y - X*edg[0]+edg[1])
        x = X[ind]
        y = Y[ind]
        pt0,pt1 = Uppoint[i],Uppoint[i+1]
        x0, x1 = x,x
        if(pt1[0] == x):
            if(x<topx):
                pt1[0] = pt1[0]-0.00001 
            else:
                x1 = x1+0.00001
        if(pt0[0] == x): 
            if(x<topx): x0=x0-0.00001
            else: pt0[0] = pt0[0]+0.00001
        if(y-0.001 > x*edg[0]+edg[1]):
            Upedge[i] = slope_int([x1,y],pt1)
            Uppoint.insert(i+1,[x,y])
            Upedge.insert(i,slope_int(pt0,[x0,y]))
            lenedg = lenedg+1
            i = i-1
        i = i+1
    lenedg = len(Downedge)
    i = 0
    while(i<lenedg):
        edg = Downedge[i]
        ind = np.argmax(X*edg[0]+edg[1]-Y)
        x = X[ind]
        y = Y[ind]
        pt0 = list(Downpoint[i])
        pt1 = list(Downpoint[i+1])
        x_1,x_0 = x, x
        ind = 0
        if(pt0[0] == x):
            if(pt0[1] == y): ind = 1
            if(x<botx):
                x_0 = x-0.00001
            else:
                pt0[0]=x+0.00001
        if(pt1[0] == x):
            if(pt1[1] == y): ind = 1
            if(x<botx):
                pt1[0] = x-0.0001
            else:
                x_1 = x+0.0001
        if(x*edg[0]+edg[1] > y-0.0001):
         if(ind == 0):
            Downpoint.insert(i+1,[x,y])
            Downedge[i] = slope_int([x_1,y],pt1)
            Downedge.insert(i,slope_int(pt0,[x_0,y]))
            lenedg = lenedg+1
            i = i-1
        i = i+1
    return [[Upedge,Downedge],[Uppoint,Downpoint]]

def dist_from_edge(Loc,Upedge,Downedge): #identify shortest distance from point array to polygon 
    X = Loc[:,0]
    Y = Loc[:,1]
    Dist = []
    for edg in Upedge:
        [m,b] = edg
        if(m == 0):
            Dist.append((Y-b)**2)
        else:
            Dist.append(np.linalg.norm([(Y-m*X+b)/2/m, (m*X-Y+3*b)/2],axis=0))
    for edg in Downedge:
        [m,b] = edg
        if(m == 0):
            Dist.append((Y-b)**2)
        else:
            Dist.append(np.linalg.norm([(Y-m*X+b)/2/m, (m*X-Y+3*b)/2],axis=0))
    return np.min(Dist,axis=0)

def find_all_poly(data,t=0): #data = array([[x_1,y_1,t_1,...]...])
#identify polygon encompassing datapoints in each frame and the t frames before/after it. Takes awhile
    maxframes = int(max(data[:,2]))
    Output = []
    for i in xrange(maxframes):
        if(i%20 == 0): print i
        ind = []
        for j in range(-t,t+1):
            ind = ind+list(np.where(data[:,2] == i+j)[0])
        if(ind != []):
            Output.append(find_poly(data[ind,:2]))
    return Output

def find_all_dist(data,t=0):#data = array([[x_1,y_1,t_1,...]...])
#Identify shortest distance of each datapoint from the Polygon for the frame its in
    maxframes = int(max(data[:,2]))
    Output = data[:,1]
    Poly = find_all_poly(data,t)
    for i in xrange(maxframes):
        ind = np.where(data[:,2] == i)[0]
        if(ind != []):
            [Upedg,Dnedg] = Poly[i][0]
            Output[ind] = dist_from_edge(data[:,:2],Upedg,Dnedg)
    return Output

def showPoly(Uppoint,Downpoint,Loc,specialedge=[],specialpt = []):
    Uppoint = np.array(Uppoint)
    plt.plot(Uppoint[:,0],Uppoint[:,1],'k')
    Downpoint = np.array(Downpoint)
    plt.plot(Downpoint[:,0],Downpoint[:,1],'k')
    for edg in specialedge:
        edg = np.array(edg)
        plt.plot(edg[:,0],edg[:,1],'r')
    plt.plot(Loc[:,0],Loc[:,1],'bo')
    for pt in specialpt:
        plt.plot(pt[0],pt[1],'ro')
    plt.show()

def showPoly2(Uppoint,Downpoint,Loc,maxx,minx,specialedge=[],specialpt = []):
    x=np.arange(21)*(maxx-minx)/20+minx
    for edg in Uppoint:
        edg = np.array(edg)
        plt.plot(edg[:,0],edg[:,1],'k')
    for edg in Downpoint:
        edg = np.array(edg)
        plt.plot(edg[:,0],edg[:,1],'k')
    for edg in specialedge:
        plt.fill_between(x,edg[0]*x+edg[1])
        print edg[0]
    plt.plot(Loc[:,0],Loc[:,1],'bo')
    for pt in specialpt:
        plt.plot(pt[0],pt[1],'ro')
    plt.ylim(0,500)
    plt.show()


if(__name__ == '__main__'): #For testing functions
#    Loc = np.transpose([30*np.random.random(30),20*np.random.random(30)])

    f = 'Bot_large.txt' #scrum data
    data_l = np.genfromtxt(f,delimiter = ',')
    f = 'Bot_small.txt' #particle data
    data_s = np.genfromtxt(f,delimiter = ',')
    #generating necesary info from data
    maxframe = int(np.max(data_s[:,2]))
#    [Traj_s,mintrajlen,maxtrajlen,ave,Frame_s,[minx,maxx],[miny,maxy]] = split(data_s,5)
#    [Traj_l,Frame_l] = split(data_l,5,[1,0,1,0])
    data = np.concatenate([data_s,data_l])
##    Celledg = find_all_poly(data,3)
    loc = data[:50,:2]
    print vor_density(loc)
