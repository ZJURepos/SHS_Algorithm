#coding=utf-8
from sklearn.cluster import DBSCAN
from sklearn.cluster import OPTICS
from sklearn.cluster import cluster_optics_dbscan

from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering

from open3d import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import DistanceMetric
from mpl_toolkits.mplot3d import Axes3D
import math
import random
from itertools import chain


def creatcolor():
    r=random.random()
    g=random.random()
    b=random.random()
    # print(r,g,b)
    return r,g,b

def avedis(listin):
    dis=np.mean(listin)
    return dis

def transpose(listin):
    return list(zip(*listin))

def ComputeConvariance(pts):
    xall=0
    yall=0
    zall=0
    psize=np.asarray(pts).size / 3
    # for i in range (0,int(psize)):
    #     xall += pts[i][0]
    #     yall += pts[i][1]
    #     zall += pts[i][2]
    # print([xall/psize,yall/psize,zall/psize])
    x = pts[:,0]
    y = pts[:,1]
    z = pts[:,2]
    mx=x.mean()
    my=y.mean()
    mz=z.mean()
    covxx = np.mean((x-mx)**2)
    covyy = np.mean((y-my)**2)
    covzz = np.mean((z-mz)**2)
    covxy = np.mean((x-mx)*(y-my))
    covxz = np.mean((x-mx)*(z-mz))
    covyz = np.mean((y-my)*(z-mz))

    out=[[covxx,covxy,covxz],[covxy, covyy, covyz],[covxz, covyz, covzz]]
    out2=np.asarray(out)
    return out

def ComputeMahalanobisDistance(p1,p2,con):
    x = p1[0] - p2[0]
    y = p1[1] - p2[1]
    z = p1[2] - p2[2]
    vec=np.asarray([x,y,z])
    a11 = x * con[0][0] + y * con[1][0] + z * con[2][0]
    a12 = x * con[0][1] + y * con[1][1] + z * con[2][1]
    a13 = x * con[0][2] + y * con[1][2] + z * con[2][2]
    A = np.matrix(con)
    s_vec = np.linalg.inv(A)
    in1=np.inner(vec.transpose(), s_vec)
    out=np.inner(in1,vec)
    # print(out)
    return math.sqrt(out)

def MahDisvec(incloud):
    pts1=np.asarray(incloud.points)
    psize1 = np.asarray(pts1).size / 3
    pts = pts1[0:int(psize1 - 1)]
    disvec = []
    con = ComputeConvariance(pts)
    for i in  range (0,int(psize1-1)):
        distemp=ComputeMahalanobisDistance(pts[i], pts1[int(psize1-1)], con)
        # print("disance is {}".format(distemp))
        disvec.append(distemp)
    return disvec

def getDisList(incloud):
    points = incloud.points
    listoutfinal = []
    listoutsorted = []
    ptssize = int(np.asarray(points).size / 3)
    for i in range(0, ptssize - 1):
        listtemp = []
        sortedtemp = []
        for j in range(0, ptssize - 1):
            xt = points[i][0] - points[j][0]
            yt = points[i][1] - points[j][1]
            zt = points[i][2] - points[j][2]
            dis = math.sqrt(xt * xt + yt * yt + zt * zt)

            listtemp.append(dis)
            sortedtemp.append(dis)
        listoutfinal.append(listtemp)
        sortedtemp.sort()
        listoutsorted.append(sortedtemp)
    return listoutfinal, listoutsorted

    # 马氏距离
    # pts1 = np.asarray(incloud.points)
    # con = ComputeConvariance(pts1)
    # listoutfinal = []
    # listoutsorted = []
    # ptssize = int(pts1.size / 3)
    # for i in range(0, ptssize - 1):
    #     listtemp = []
    #     sortedtemp = []
    #     for j in range(0, ptssize - 1):
    #         # 欧氏距离
    #         dis = ComputeMahalanobisDistance(pts1[i], pts1[j], con)
    #         listtemp.append(dis)
    #         sortedtemp.append(dis)
    #     listoutfinal.append(listtemp)
    #     sortedtemp.sort()
    #     listoutsorted.append(sortedtemp)
    # return listoutfinal, listoutsorted


def AdpatedDBscan(incloud):

    [dislist, dissortlist] = getDisList(incloud)
    distrans = transpose(dissortlist)

    for i in range (0,len(dissortlist[0])):
        print(dissortlist[0][i])

    disvec = []
    for i in range(0, len(dissortlist)):
        dis = avedis(distrans[i])
        disvec.append(dis)

    disall = list(chain.from_iterable(dissortlist))
    disall.sort()

    currentpos = 0
    for i in range(0, len(disall)):
        if (disall[i] == 0):
            currentpos += 1
            continue;

    disall2 = []
    for i in range(currentpos, len(disall)):
        disall2.append(disall[i])

    minptsvec = [0]
    minpts_current = 0
    currentpos2 = 0
    for i in range(1, len(disvec)):
        minptstemp = 0
        minpts = 0
        for j in range(currentpos2, len(disall2)):
            if (disall2[j] <= disvec[i]):
                minptstemp += 1
            else:
                currentpos2 = j
                minpts = minptstemp + minpts_current
                minpts_current = minpts
                break
        minptsvec.append(minpts / len(dissortlist))

    # Density
    Density = [0.0]
    densityx = [0.0]
    for i in range(1, len(minptsvec)):
        densitytemp = minptsvec[i] / (math.pi * disvec[i] * disvec[i])
        Density.append(densitytemp)
        densityx.append(i)
    return disvec,minptsvec,Density

def DBscan_ad1(incloud,data,psize1,eps,minpts,density):
    X = data[["x", "y", "z"]]
    X_train=np.asarray(incloud.points)
    pts1=np.asarray(incloud.points)
    color=np.asarray(incloud.colors)
    X_train2 = pts1[0:int(psize1 - 1)]
    db = DBSCAN(eps=eps, min_samples=minpts).fit(X)


    labels = db.labels_
    data['id_ini'] = 0
    data['cluster_label'] = labels

    samelabelcounts =  data['cluster_label'].value_counts()
    print(samelabelcounts)
    print(samelabelcounts.max())
    if(samelabelcounts.max()==psize1):
        return 1

    # classification with color
    datacolor = pd.DataFrame(color, columns=['r', 'g', 'b'])
    datacolor['id_ini'] = 0
    for i in range(0, int(data.size / 5)):
        datacolor.loc[i, 'id_ini'] = i

    datacolor['cluster_label'] = labels
    cenlabel = datacolor.iloc[int(psize1 - 1), 4]
    d2 = data.sort_values('cluster_label')

    # label
    print("------------")
    dc2 = datacolor.sort_values('cluster_label')
    labels2 = dc2['cluster_label']
    # print(labels2)
    # print(dc2)
    label_begin = dc2['cluster_label'][0]
    label_size = 1
    [r, g, b] = creatcolor()
    dc2.loc[0, 'r'] = r
    dc2.loc[0, 'g'] = g
    dc2.loc[0, 'b'] = b
    pointsize = int(data.size / 5)
    for i in range(1, pointsize, 1):
        label_temp = dc2.iloc[i, 4]
        # print(label_temp)
        if (label_temp == cenlabel):
            dc2.iloc[i, 0] = 1
            dc2.iloc[i, 1] = 0
            dc2.iloc[i, 2] = 0
            continue
        if (label_temp == label_begin):
            dc2.iloc[i, 0] = r
            dc2.iloc[i, 1] = g
            dc2.iloc[i, 2] = b
            label_begin = label_temp
        else:
            [r, g, b] = creatcolor()
            dc2.iloc[i, 0] = r
            dc2.iloc[i, 1] = g
            dc2.iloc[i, 2] = b
            label_begin = label_temp
            label_size = label_size + 1
    # print(label_size)
    print("------------")

    x = pts1[:, 0]
    y = pts1[:, 1]
    z = pts1[:, 2]

    for i in range(1, int(dc2.size / 5)):
        id = dc2['id_ini'][i]
        incloud.colors[id][0] = dc2['r'][i]
        incloud.colors[id][1] = dc2['g'][i]
        incloud.colors[id][2] = dc2['b'][i]

    ax4 = plt.axes(projection='3d')
    ax4.scatter3D(x, y, z)
    plt.figure()

    return 0

def Optis_ad1(incloud,data,psize1,eps,minpts,density):
    X = data[["x", "y", "z"]]
    pts1=np.asarray(incloud.points)
    color=np.asarray(incloud.colors)
    db = OPTICS(min_samples=5, xi=.05, min_cluster_size=.05).fit(X)

    labels_050 = cluster_optics_dbscan(reachability=db.reachability_,
                                       core_distances=db.core_distances_,
                                       ordering=db.ordering_, eps=eps)

    labels = labels_050
    data['id_ini'] = 0
    data['cluster_label'] = labels

    samelabelcounts =  data['cluster_label'].value_counts()
    print(samelabelcounts)
    print(samelabelcounts.max())
    if(samelabelcounts.max()==psize1):
        return 1

    datacolor = pd.DataFrame(color, columns=['r', 'g', 'b'])
    datacolor['id_ini'] = 0
    for i in range(0, int(data.size / 5)):
        datacolor.loc[i, 'id_ini'] = i

    datacolor['cluster_label'] = labels
    cenlabel = datacolor.iloc[int(psize1 - 1), 4]
    d2 = data.sort_values('cluster_label')

    print("------------")
    dc2 = datacolor.sort_values('cluster_label')
    labels2 = dc2['cluster_label']
    # print(labels2)
    # print(dc2)
    label_begin = dc2['cluster_label'][0]
    label_size = 1
    [r, g, b] = creatcolor()
    dc2.loc[0, 'r'] = r
    dc2.loc[0, 'g'] = g
    dc2.loc[0, 'b'] = b
    pointsize = int(data.size / 5)
    for i in range(1, pointsize, 1):
        label_temp = dc2.iloc[i, 4]
        # print(label_temp)
        if (label_temp == cenlabel):
            dc2.iloc[i, 0] = 1
            dc2.iloc[i, 1] = 0
            dc2.iloc[i, 2] = 0
            continue
        if (label_temp == label_begin):
            dc2.iloc[i, 0] = r
            dc2.iloc[i, 1] = g
            dc2.iloc[i, 2] = b
            label_begin = label_temp
        else:
            [r, g, b] = creatcolor()
            dc2.iloc[i, 0] = r
            dc2.iloc[i, 1] = g
            dc2.iloc[i, 2] = b
            label_begin = label_temp
            label_size = label_size + 1
    # print(label_size)
    print("------------")

    x = pts1[:, 0]
    y = pts1[:, 1]
    z = pts1[:, 2]

    for i in range(1, int(dc2.size / 5)):
        id = dc2['id_ini'][i]
        incloud.colors[id][0] = dc2['r'][i]
        incloud.colors[id][1] = dc2['g'][i]
        incloud.colors[id][2] = dc2['b'][i]

    ax4 = plt.axes(projection='3d')
    ax4.scatter3D(x, y, z, c=incloud.colors)
    plt.figure()

    return 0


def Agglomerative_ad1(incloud,data,psize1,eps,minpts,density):
    X = data[["x", "y", "z"]]
    X_train=np.asarray(incloud.points)
    pts1=np.asarray(incloud.points)
    color=np.asarray(incloud.colors)
    row_clusters = linkage(pdist(data, metric='euclidean'), method='complete')

    db = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='complete')
    labels = db.fit_predict(X)

    labels = labels
    data['id_ini'] = 0
    data['cluster_label'] = labels

    samelabelcounts =  data['cluster_label'].value_counts()
    print(samelabelcounts)
    print(samelabelcounts.max())
    if(samelabelcounts.max()==psize1):
        return 1

    datacolor = pd.DataFrame(color, columns=['r', 'g', 'b'])
    datacolor['id_ini'] = 0
    for i in range(0, int(data.size / 5)):
        datacolor.loc[i, 'id_ini'] = i

    datacolor['cluster_label'] = labels
    cenlabel = datacolor.iloc[int(psize1 - 1), 4]
    d2 = data.sort_values('cluster_label')

    print("------------")
    dc2 = datacolor.sort_values('cluster_label')
    labels2 = dc2['cluster_label']
    # print(labels2)
    # print(dc2)
    label_begin = dc2['cluster_label'][0]
    label_size = 1
    [r, g, b] = creatcolor()
    dc2.loc[0, 'r'] = r
    dc2.loc[0, 'g'] = g
    dc2.loc[0, 'b'] = b
    pointsize = int(data.size / 5)
    for i in range(1, pointsize, 1):
        label_temp = dc2.iloc[i, 4]
        # print(label_temp)
        if (label_temp == cenlabel):
            dc2.iloc[i, 0] = 1
            dc2.iloc[i, 1] = 0
            dc2.iloc[i, 2] = 0
            continue
        if (label_temp == label_begin):
            dc2.iloc[i, 0] = r
            dc2.iloc[i, 1] = g
            dc2.iloc[i, 2] = b
            label_begin = label_temp
        else:
            [r, g, b] = creatcolor()
            dc2.iloc[i, 0] = r
            dc2.iloc[i, 1] = g
            dc2.iloc[i, 2] = b
            label_begin = label_temp
            label_size = label_size + 1
    # print(label_size)
    print("------------")

    x = pts1[:, 0]
    y = pts1[:, 1]
    z = pts1[:, 2]

    for i in range(1, int(dc2.size / 5)):
        id = dc2['id_ini'][i]
        incloud.colors[id][0] = dc2['r'][i]
        incloud.colors[id][1] = dc2['g'][i]
        incloud.colors[id][2] = dc2['b'][i]

    ax4 = plt.axes(projection='3d')
    ax4.scatter3D(x, y, z, c=incloud.colors)
    plt.figure()

    return 0

