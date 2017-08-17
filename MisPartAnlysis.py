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

def Dense(data):
    frames = data[:,2]
    Npart = []
    for i in range(int(max(frames))):
        Npart.append(len(np.where(frames == i)[0]))
    return Npart
        
def Crowd(Frame,Npart):
    Npart = np.array(Npart)
    crowd = []
    for i in range(len(Traj)):
        crowd.append(np.var(Frame.pop()))
    return crowd

def find_poly(Loc): #identify polygon that encloses datapoints
    def slope_int(p1,p2):
        if(p1[0] == p2[0]):
            p = p1[0]+0.0000001
            m = (p1[1]-p2[1])/(p-p2[0])
            b = p1[1]-m*p
  #          print "warning: p1=p2"
        else:
            m = (p1[1]-p2[1])/(p1[0]-p2[0])
            b = p1[1]-m*p1[0]
        return [m,b]
    maxx = max(Loc[:,0])
    maxy = max(Loc[:,1])
    minx = min(Loc[:,0])
    miny = min(Loc[:,1])
    ind1 = np.where(Loc[:,0]==maxx)[0][0]
    ind2 = np.where(Loc[:,1]==maxy)[0][0]
    ind3 = np.where(Loc[:,0]==minx)[0][0]
    ind4 = np.where(Loc[:,1]==miny)[0][0]
    Upedge = [slope_int(Loc[ind1],Loc[ind2]),slope_int(Loc[ind3],Loc[ind2])]
    Uppoint = [[Loc[ind1],Loc[ind2]],[Loc[ind3],Loc[ind2]]]
    Downedge = [slope_int(list(Loc[ind1]),Loc[ind4]),slope_int(Loc[ind3],Loc[ind4])]
    Downpoint = [[list(Loc[ind1]),Loc[ind4]],[Loc[ind3],Loc[ind4]]]
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
    for j in xrange(1):
        X = Loc[:,0]
        Y = Loc[:,1]
        lenedg = len(Upedge)
        i = 0
        while(i < lenedg):
            edg = Upedge[i]
            ind = np.argmax(Y - X*edg[0]+edg[1])
            x = X[ind]
            y = Y[ind]
            if(y > x*edg[0]+edg[1]+.000001):
                pt = Uppoint[i]
                if(x*edg[0]+edg[1]>y+0.00001):
                    Uppoint[i] = [pt[1],[x,y]]
                    Upedge[i] = slope_int(pt[1],[x,y])
                    Uppoint.insert(i,[pt[0],[x,y]])
                    Upedge.insert(i,slope_int(pt[0],[x,y]))
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
            if(x*edg[0]+edg[1]>y+0.00001):
                pt = Downpoint[i]
                if(pt[1][0] != x and pt[0][0] != x):
                    Downpoint.pop(i)
                    Downpoint.insert(i,[pt[1],[x,y]])
                    Downedge[i] = slope_int(pt[1],[x,y])
                    Downpoint.insert(i,[pt[0],[x,y]])
                    Downedge.insert(i,slope_int(pt[0],[x,y]))
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

if(__name__ == '__main__'): #For testing functions
#    Loc = np.transpose([30*np.random.random(30),20*np.random.random(30)])
    f = 'Bot_large.txt'
    data = np.genfromtxt(f,delimiter = ',')
    Loc = data[np.where(data[:,2] == 1)[0],:2]
    maxx = max(Loc[:,0])
    minx = min(Loc[:,0])
    Output = find_all_poly(data)
    [Upedge,Downedge] = Output[1][1]
    for edg in Upedge:
        edg = np.array(edg)
        plt.plot(edg[:,0],edg[:,1],'k')
    for edg in Downedge:
        edg = np.array(edg)
        plt.plot(edg[:,0],edg[:,1],'k')
    plt.plot(Loc[:,0],Loc[:,1],'o')
    plt.show()
