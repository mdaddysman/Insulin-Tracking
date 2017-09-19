import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d,Delaunay, delaunay_plot_2d
import matplotlib.pyplot as plt
import time
def split(rdgvert,rdgpt):
    outerboundry = set([])
    innerboundry = set([])
    inside = set([])
    i=0
    for edg in rdgvert:
        if(-1 in edg):
            outerboundry.add(rdgpt[i][0])
            outerboundry.add(rdgpt[i][1])
        i=i+1
    i=0
    for pt in rdgpt:
        if(not(outerboundry.issuperset(pt))):
            if(pt[0] in outerboundry):
               innerboundry.add(pt[1])
            elif(pt[1] in outerboundry):
               innerboundry.add(pt[0])
            else:
               inside=inside.union(pt)
    return (outerboundry,innerboundry,inside)

def rmEdg(Regions,Vert):
    k=0
    for reg in Regions:
        k=k+1
        v = np.round(Vert[reg],3)
        for i in range(1,len(v)):
            ind = np.where(np.sum(v[i:] == v[:-i],axis=1)==2)[0]
            for j in ind:
                reg.pop(j)
                v[j] = np.array([-5,-5])
    return Regions

def chkEdgs(v):
    extraEdg = []
    extras = set([])
    l = len(v)
    for i in xrange(1,l):
        ind = np.where(np.linalg.norm(v[i:] -v[:-i],axis=1)<0.005)[0]
        for j in ind:
            extraEdg.append([j,i+j])
            extras.add(j)
            extras.add(i+j)
    return (extraEdg,extras)

def neigh(Regions): #number of neighbors
    Neigh = []
    for reg in Regions:
        Neigh.append(len(reg))
    return Neigh

def area(poly): #surveyer's formula
    return np.absolute(np.sum(poly[:-1,0]*poly[1:,1]-poly[1:,0]*poly[-1,1]))

def area_per_reg(Regions,Vert):
    Area = []
    for reg in Regions:
        v = Vert[reg]
        Area.append(area(v))
    return Area

def degree(ridge_vert,extraEdg=None,extras=None):
    Degree = []
    if(isinstance(extras,set)):
        k=0
        for v in ridge_vert:
            if(v in extraEdg):
                ridge_vert.pop(k)
            else:
                k=k+1
    ridge_vert = np.concatenate(ridge_vert)
    extraEdg=np.array(extraEdg)
    for i in xrange(max(ridge_vert)+1):
        d = sum(ridge_vert==i)
        if(isinstance(extras,set) and i in extras):
            ind = np.where(extraEdg == i)
            if(np.sum(ind[0])==0):
                for j in ind[1]:
                    d=d+sum(ridge_vert == j)
                Degree.append(d)
        else:
            Degree.append(d)
    return Degree

def dist_from_neigh(Loc,ridge_pt,ptset=None):
    dist = []
    for ptind in ridge_pt:
        if(ptset is None or ptset.issuperset(ptind)):
            pts = Loc[ptind]
            dist.append(np.linalg.norm(pts[0]-pts[1]))
    return dist

def perimeter(Regions,Vert):
    Perim = []
    for v in Regions:
        Perim.append(np.linalg.norm(Vert[v]-Vert[np.roll(v,1)]))
    return Perim

def plot_vor(Loc,subsets):
    for s in subsets:
        s = np.array ( list(s))
        pts = Loc[s]
        plt.plot(pts[:,0],pts[:,1],'o')
        
def vor_stats(data):
    Stats = np.empty([len(data),7])
    for i in xrange(int(max(data[:,2]))):
        ind = np.where(data[:,2] == i)[0]
        Loc = data[ind,:2]
        vor = Voronoi(Loc)
        Regions = np.array(vor.regions)[vor.point_region]
        (outerboundry,innerboundry,inside) = split(vor.ridge_vertices,vor.ridge_points)
        Regions = rmEdg(Regions,vor.vertices)
        (extraEdg,extras) = chkEdgs(vor.vertices)
        Dist = dist_from_neigh(Loc,vor.ridge_points)
        Degree = degree(vor.ridge_vertices,extraEdg,extras)
        DistAve = []
        DistVar = []
        DegAve = []
        for j in xrange(len(ind)):
            dist = np.array(Dist)[np.where(vor.ridge_points == j)[0]]
            DistAve.append(np.mean(dist))
            DistVar.append(np.var(dist))
#            DegAve.append(np.mean(np.array(Degree)[np.where(vor.ridge_vertices == j)[0]]))
        cnt = 0
        Stats[ind,cnt] = neigh(Regions)
        cnt = cnt+1
        Stats[ind,cnt] = area_per_reg(Regions,vor.vertices)
        cnt = cnt+1
        Stats[ind,cnt] = perimeter(Regions,vor.vertices)
        cnt = cnt+1
        Stats[ind,cnt] = DistAve
        cnt = cnt+1
        Stats[ind,cnt] = DistVar
        cnt = cnt+1
        Stats[ind,cnt] = Degree
        cnt = cnt+1
        Stats[ind[list(outerboundry)],cnt] = 2
        Stats[ind[list(innerboundry)],cnt] = 1
        Stats[ind[list(inside)],cnt] = 0
        if(i/10*10==i): print i
    return Stats

def vor_stats2(data):
    Stats = []
    for i in xrange(1,int(max(data[:,2]))):
        ind = np.where(data[:,2] == i)[0]
        Loc = data[ind,:2]
        vor = Voronoi(Loc)
        vor.close()
        ridgevert = vor.ridge_vertices
        Regions = np.array(vor.regions)[vor.point_region]
        (outerboundry,innerboundry,inside) = split(vor.ridge_vertices,vor.ridge_points)
        Regions = rmEdg(Regions,vor.vertices)
        (extraEdg,extras) = chkEdgs(vor.vertices)
        Dist = dist_from_neigh(Loc,vor.ridge_points)
        DistAve = []
        DistVar = []
        DegAve = []
        vor.regions.remove([])
        for j in xrange(len(ind)):
            ind2 = np.where(vor.ridge_points == j)[0]
            dist = np.array(Dist)[ind2]
            DistAve.append(np.mean(dist))
            DistVar.append(np.var(dist))
        stats = []
        stats.append(neigh(Regions))
        stats.append(area_per_reg(Regions,vor.vertices))
        stats.append(perimeter(Regions,vor.vertices))
        stats.append(DistAve)
#        stats.append(DistVar)
#        stats = np.hstack([np.zeros([len(ind),1]),np.transpose(stats)])
#        stats[list(outerboundry),0] = 1
        Stats.append(np.transpose(stats))
        if(i/10*10==i): print i
    return Stats

if (__name__ == '__main__'):
    f = '../../Downloads/Bot_small.txt'
 #   f = '../../Downloads/modeltraj_1.txt'
    data = np.genfromtxt(f,delimiter=',')
    A = []
    for i in xrange(int(max(data[:,3]))):
        A.append(len(np.where(data[:,3]==i)[0]))
    plt.hist(A,bins)
    plt.show()
    sdf
    stats = vor_stats2(data)
    A = np.mean(np.concatenate(stats),axis=0)
    for e in A:
        digits = int(np.log10(e))
        if(digits>3):
            print round(e/10**digits,3),'$\\times 10^{',str(digits-3), '}$ &',
        elif(digits>0):
            print round(e/10**(digits-1),4-digits), '&',
        else:
            print round(e/10**(digits-3),4),'&',
    print '\\\\'
    A = []
    for stat in stats:
        A.append(np.var(stat,axis=0))
    A = np.mean(A,axis=0)
    for e in A:
        digits = int(np.log10(e))
        if(digits>3):
            print round(e/10**digits,3),'$\\times 10^{',str(digits-3), '}$ &',
        elif(digits>0):
            print round(e/10**(digits-1),4-digits), '&',
        else:
            print round(e/10**(digits-3),4),'&',
    print '\\\\'
    A = []
    for stat in stats:
        A.append(np.mean(stat,axis=0))
    A = np.var(A,axis=0)
    for e in A:
        digits = int(np.log10(e))
        if(digits>3):
            print round(e/10**digits,3),'$\\times 10^{',str(digits-3), '}$ &',
        elif(digits>0):
            print round(e/10**(digits-1),4-digits), '&',
        else:
            print round(e/10**(digits-3),3),'&',
    print '\\\\'
 #   print np.mean(stats,axis=[0,1])
    
##    Traj = [np]*(np.max(data[:,3])-1)
        
##        j=0
##        for traj in Traj:
##            j=j+1
##            ind2 = np.where(data[ind,3] == j)[0]
##            if(traj==[]): traj = [Loc[ind2][0],Loc[ind2][1],i,Neigh[ind2],Area[ind2],dist[ind2]]
##            else:
##                if(data[ind[ind2]+1,3]==j):
##                    traj = np.array(traj)

##    if(np.mean(degree(vor.ridge_vertices,extraEdg,extras)) != 3):
##        print i,np.mean(degree(vor.ridge_vertices,extraEdg,extras))

#print np.average(A), np.var(A)


##for i in xrange(3):
##x=list(np.random.random(100))
##y = list(np.random.random(100))
##for i in xrange(100):
##    if(x[99-i]>0.8 and y[99-i]>0.8):
##        x.pop(99-i)
##        y.pop(99-i)
##Loc=np.transpose([x,y])[:60]
#    Regions = np.array(vor.regions)[vor.point_region]
#    Regions = rmEdg(Regions,vor.vertices)
