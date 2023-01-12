# coding=utf-8
from sklearn.cluster import DBSCAN
from sklearn.cluster import OPTICS
from sklearn.cluster import cluster_optics_dbscan

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
import sys
import tools
sys.path.append("..")

from open3d import *
import numpy as np
import copy
import random
import MyRansac as rs
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from sklearn.neighbors import DistanceMetric
from mpl_toolkits.mplot3d import Axes3D
import math
# from ransac import *
from itertools import chain

def getFloatData(data):
    # step=len(data[0])
    psize = len(data)
    dataout = []
    for i in range(psize):
        dataout.append([float(data[i][0]), float(data[i][1]), float(data[i][2]), int(data[i][3]), int(data[i][4])]);
    return dataout

def loadDatadet(infile):
    f = open(infile, 'r')
    sourceInLine = f.readlines()
    dataset = []
    psize = 0
    for line in sourceInLine:
        temp1 = line.strip('\n')
        temp1 = temp1.strip('')
        temp2 = temp1.split('\t')
        dataset.append(temp2)
        psize += 1
    dataout = getFloatData(dataset)
    return dataout, psize

def creatcolor():
    r = random.random()
    g = random.random()
    b = random.random()
    return r, g, b

def computeCon(pts1):
    cov_mat = np.cov(pts1.T)
    eigen_vals, eigen_vecs = np.linalg.eig(cov_mat)  #
    return cov_mat, eigen_vals, eigen_vecs

def computeLamda(eigen_val):
    l1 = eigen_val[0]
    l2 = eigen_val[1]
    l3 = eigen_val[2]
    l = l1 / (l1 + l2 + l3)
    return l, l1

def line_pcd(p1, p2):

    linepoints = np.array([[p1[0], p1[1], p1[2]], [p2[0], p2[1], p2[2]]])
    lines = [[0, 1]]  # Right leg
    colors = [[0, 0, 1] for i in range(len(lines))]
    point_pcd = open3d.geometry.PointCloud()
    point_pcd.points = open3d.utility.Vector3dVector(linepoints)

    line_pcd = open3d.geometry.LineSet()
    line_pcd.lines = open3d.utility.Vector2iVector(lines)
    line_pcd.colors = open3d.utility.Vector3dVector(colors)
    line_pcd.points = open3d.utility.Vector3dVector(linepoints)
    return line_pcd

def line_pcds2(ptsin, tempsize, cloudin):
    pts = []
    lines = []
    for i in range(0, tempsize):
        pts.append(cloudin.points[i])
    linepoints = np.array(pts)
    lines = ptsin  # Right leg
    c = [1, 0, 0]
    colors = [c for i in range(len(lines))]
    point_pcd = open3d.geometry.PointCloud()
    point_pcd.points = utility.Vector3dVector(linepoints)

    line_pcd = open3d.geometry.LineSet()
    line_pcd.lines = utility.Vector2iVector(lines)
    line_pcd.colors = utility.Vector3dVector(colors)
    line_pcd.points = utility.Vector3dVector(linepoints)
    return line_pcd

def computeAngle(v1, v2):
    modv1 = math.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2])
    modv2 = math.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2])
    cosA = math.fabs((v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (modv1 * modv2))
    return cosA

def giveColor(SortedPts, psize):
    colors = np.zeros((int(psize), 5))
    label_begin = SortedPts[0][3]
    [r, g, b] = creatcolor()
    colors[0][0] = r
    colors[0][1] = g
    colors[0][2] = b
    # colors[0][3] = int(label_begin)
    for i in range(0, int(psize)):
        label_temp = SortedPts[i][3]
        if (label_temp == label_begin):
            colors[i][0] = r
            colors[i][1] = g
            colors[i][2] = b
            colors[i][3] = SortedPts[i][3]  # label
            colors[i][4] = SortedPts[i][4]  # id
            # colors[i][3]=int(label_temp)
        else:
            label_begin = label_temp
            [r, g, b] = creatcolor()
            colors[i][0] = r
            colors[i][1] = g
            colors[i][2] = b
            colors[i][3] = SortedPts[i][3]  # label
            colors[i][4] = SortedPts[i][4]  # id
    return colors

def computeDis(p1, p2):
    p1p2=p1-p2
    dis = math.sqrt(p1p2[0] * p1p2[0] + p1p2[1] * p1p2[1] + p1p2[2] * p1p2[2])
    return dis

def computeDireDis(p1, p2, p1dire):
    p1p2 = p2-p1
    p1p2_p1dire = (p1p2[0] * p1dire[0] + p1p2[1] * p1dire[1] + p1p2[2] * p1dire[2])
    p1p2mod = math.sqrt(p1p2[0] * p1p2[0] + p1p2[1] * p1p2[1] + p1p2[2] * p1p2[2])
    p1diremod = math.sqrt(p1dire[0] * p1dire[0] + p1dire[1] * p1dire[1] + p1dire[2] * p1dire[2])
    # p1p2mod=math.sqrt(p1p2[0] * p1p2[0] + p1p2[1] * p1p2[1] + p1p2[2] * p1p2[2])
    cosA = p1p2_p1dire / (p1p2mod * p1diremod)
    if(cosA>1):
        cosA=1
    if(cosA<-1):
        cosA=-1
    sinA = math.sqrt(1 - cosA * cosA)
    disEu = math.sqrt(p1p2[0] * p1p2[0] + p1p2[1] * p1p2[1] + p1p2[2] * p1p2[2])
    angle =(math.acos(cosA)) * 180 / math.pi
    dis = float(p1p2_p1dire / p1diremod)
    #
    if (math.fabs(float(sinA / cosA)) > 5):
        dis = float("inf")
    else:
        dis = float(p1p2_p1dire / p1diremod)
    if (disEu > 20):
        dis = float("inf")
    if (dis > 20 or dis <-20):
        dis = float("inf")
    disEuDir = float(p1p2_p1dire / math.fabs(p1p2_p1dire)) * disEu
    # dis=math.sqrt(x*x+y*y+z*z)
    return dis, disEuDir, angle

def computeDisPN(p1, p2, p1dire):
    p1p2 = p2 - p1
    p1p2_p1dire = (p1p2[0] * p1dire[0] + p1p2[1] * p1dire[1] + p1p2[2] * p1dire[2])
    p1p2mod = math.sqrt(p1p2[0] * p1p2[0] + p1p2[1] * p1p2[1] + p1p2[2] * p1p2[2])
    p1diremod = math.sqrt(p1dire[0] * p1dire[0] + p1dire[1] * p1dire[1] + p1dire[2] * p1dire[2])
    dis = float(p1p2_p1dire / (p1diremod))
    return dis

def computeDireDis_pure(p1, p2, p1dire):
    p1p2 = []
    p1p2.append(p2[0] - p1[0])
    p1p2.append(p2[1] - p1[1])
    p1p2.append(p2[2] - p1[2])
    p1p2_p1dire = (p1p2[0] * p1dire[0] + p1p2[1] * p1dire[1] + p1p2[2] * p1dire[2])
    p1p2mod = math.sqrt(p1p2[0] * p1p2[0] + p1p2[1] * p1p2[1] + p1p2[2] * p1p2[2])
    p1diremod = math.sqrt(p1dire[0] * p1dire[0] + p1dire[1] * p1dire[1] + p1dire[2] * p1dire[2])
    cosA = p1p2_p1dire / (p1p2mod * p1diremod)
    sinA = math.sqrt(1 - cosA * cosA)
    disEu = math.sqrt(p1p2[0] * p1p2[0] + p1p2[1] * p1p2[1] + p1p2[2] * p1p2[2])
    # dis=0
    angle = (math.acos(math.fabs(cosA))) * 180 / math.pi
    dis = float(p1p2_p1dire / p1diremod)
    #
    if (angle > 15):
        dis = float("inf")
    disEuDir = float(p1p2_p1dire / math.fabs(p1p2_p1dire)) * disEu
    # dis=math.sqrt(x*x+y*y+z*z)
    return dis, disEuDir, angle

def bulid_graph_dire_singleP(currnentP, neiborsize, ske, cloudin):
    neibos = currnentP.Neibors
    graph = np.zeros(neiborsize)
    Eugraph = np.zeros(neiborsize)
    angleSet = np.zeros(neiborsize)
    p1dire = currnentP.direciton
    for i in range(0, neiborsize):
        id = neibos[i]
        p2 = cloudin.points[id]
        if (i == 0):
            dis = 0
            disEudir = 0
            angle = 0
        else:
            [dis, disEudir, angle] = computeDireDis(currnentP.point, p2, p1dire)
        graph[i] = dis
        Eugraph[i] = disEudir
        angleSet[i] = angle
        # disvectemp=sortWithDis(disvectemp)
    return graph, Eugraph, angleSet

def computeAngleJiaodu(v1, v2, isfbs=False):
    v1v2 = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])
    v1mod = math.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2])
    v2mod = math.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2])
    cosA = v1v2  / (v1mod*v2mod)
    if (cosA > 1):
        cosA = 1
    if (cosA<-1):
        cosA=-1
    if(isfbs==True):
        cosA=math.fabs(cosA)
    angle = (math.acos(cosA)) * 180 / math.pi
    return angle

def dij_dire_singleP_withClusterLabel(currnentp, neibosize, ske, graph, Eugraph, angleSet):
    anglethreshold=120
    vec1 = []
    vec2 = []
    if (currnentp.lines == 2):
        return vec1,vec2
    id1 = 0
    id2 = 0
    min = float("inf")
    if (currnentp.hasLine1 != "true"):
        for i in range(0, neibosize):
            idtemp = currnentp.Neibors[i]
            if(idtemp==currnentp.id):
                continue
            temp = ske[idtemp]
            labeltemp = temp.label
            if (labeltemp != currnentp.label):
                continue
            if (graph[i] == 0):
                continue
            if (graph[i] == float("inf")):
                continue
            if (graph[i] < 0):
                continue
            dis2Points = graph[i]
            if (dis2Points < min):
                min = dis2Points
                id1 = idtemp

    max = -float("inf")
    if (currnentp.hasLine2 != "true"):
        for i in range(0, neibosize):
            idtemp = currnentp.Neibors[i]
            temp = ske[idtemp]
            if (temp.lines == 2):
                continue
            if (graph[i] == 0):
                continue
            if (graph[i] == -float("inf")):
                continue
            if (graph[i] >= 0):
                continue

            dis2Points = graph[i]
            if (dis2Points > max):
                max = dis2Points
                id2 = idtemp

    if (min==float("inf") and max==-float("inf")):
        vec1=vec1
        vec2=vec1
    elif (min!=float("inf") and max!=-float("inf")):
        vec1.append([currnentp.id, id1])
        currnentp.hasLine1 = "true"
        currnentp.lineset.append([currnentp.id, id1])
        currnentp.lines = currnentp.lines + 1
        vec2.append([currnentp.id, id2])
        currnentp.hasLine2 = "true"
        currnentp.lineset.append([currnentp.id, id2])
        currnentp.lines = currnentp.lines + 1
    elif (min!=float("inf") and max==-float("inf")):
        vec1.append([currnentp.id, id1])
        currnentp.hasLine1 = "true"
        currnentp.lineset.append([currnentp.id, id1])
        currnentp.lines = currnentp.lines + 1
        vec2=vec2
    elif (min==float("inf") and max!=-float("inf")):
        vec1 = vec1
        vec2.append([currnentp.id, id2])
        if (currnentp.hasLine1 != "true"):
            currnentp.hasLine1 = "true"
        else:
            currnentp.hasLine2 = "true"
        currnentp.lineset.append([currnentp.id, id2])
        currnentp.lines = currnentp.lines + 1
    return vec1, vec2

def countTimes(connectlist):
    idtimes = np.zeros(pts_size)
    for j in range(0, len(connectlist)):
        a = connectlist[j]
        id1 = a[0]
        id2 = a[1]
        idtimes[id1] = idtimes[id1] + 1
        idtimes[id2] = idtimes[id2] + 1

def computeVectorAngle(v1, v2):
    v1v2 = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])
    v1mod = math.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2])
    v2mod = math.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2])

    cosA = v1v2  # / (v1mod*v2mod)
    if (math.fabs(cosA) > 1):
        cosA = 1

    angle = (math.acos(math.fabs(cosA))) * 180 / math.pi
    return angle

def computeVecDisAngle(p1, p2, e, p1dire):
    p = p2 - p1
    p1p2_p1dire = (p[0] * p1dire[0] + p[1] * p1dire[1] + p[2] * p1dire[2])
    pmod = math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
    p1diremod = math.sqrt(p1dire[0] * p1dire[0] + p1dire[1] * p1dire[1] + p1dire[2] * p1dire[2])
    cosA = p1p2_p1dire / (pmod * p1diremod)
    if (cosA > 1):
        cosA = 1
    if (cosA < -1):
        cosA = -1
    angleA = math.acos(math.fabs(cosA)) * 180 / math.pi
    if(angleA>20):
        dis=float("inf")
    else:
        dis = math.fabs(float(p1p2_p1dire / p1diremod))

    e1e2 = (p[0] * e[0] + p[1] * e[1] + p[2] * e[2])
    emod = math.sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2])
    cosB = e1e2 / (pmod * emod)
    if (cosB > 1):
        cosB = 1
    if (cosB < -1):
        cosB = -1
    angleB = math.acos(cosB) * 180 / math.pi
    if (angleB <= 120 and angleB >= 0):
        dis = float("inf")
    # dis=math.sqrt(x*x+y*y+z*z)
    return dis, angleA, angleB

def connectEnds_withLabel(cp, ske, graphtemp, Eugraph, anglegraph):
    vec1 = []
    if(cp.lines==2):
        return vec1
    neiborsize = len(cp.Neibors)
    min = float("inf")
    dis2Points = 0
    id1 = 0
    # 得到边1
    if (cp.lineset[0][0] != cp.id):
        p1id1 = cp.lineset[0][0]
    elif (cp.lineset[0][1] != cp.id):
        p1id1 = cp.lineset[0][1]
    v1 = ske[p1id1].point - cp.point

    for i in range(1, neiborsize):
        idtemp = cp.Neibors[i]
        if(idtemp==p1id1):
            continue
        ptemp = ske[idtemp]
        if(ptemp.isend=="false"):
            continue
        labeltemp = ptemp.label
        if (cp.label != labeltemp):
            continue
        if (ptemp.lines == 2):
            continue
        if (ptemp.lines == 1):
            v2=ptemp.point-cp.point
            id3=ptemp.lineset[0][1]
            v3=ptemp.point-ske[id3].point
            angle = computeAngleJiaodu(v3, v2)
            if (angle < 135):
                continue
        [dis, angleA, angleB] = computeVecDisAngle(cp.point, ptemp.point, v1, cp.direciton)
        dis=math.fabs(dis)

        if (dis == float("inf")):
            continue
        if (dis <= 0):
            continue
        dis2Points = dis
        if (dis2Points < min):
            min = dis2Points
            id1 = idtemp

    if (min != float("inf")):
        if (id1 == 1833 or id1==1455):
        if(ske[id1].lines<2):
            vec1.append([cp.id, id1])
            cp.hasLine2 = "true"
            cp.lines = cp.lines + 1
            cp.lineset.append([cp.id, id1])
            if (ske[id1].hasLine1 == "true"):
                ske[id1].hasLine2 = "true"
                ske[id1].lineset.append([id1, cp.id])
                ske[id1].lines = ske[id1].lines + 1
            else:
                ske[id1].hasLine1 = "true"
                ske[id1].lineset.append([id1, cp.id])
                ske[id1].lines = ske[id1].lines + 1

    return vec1

def connectEnds_withLabelV2(cp, ske, graphtemp, Eugraph, anglegraph):
    vec1 = []
    if(cp.lines==2):
        return vec1
    neiborsize = len(cp.Neibors)
    min = float("inf")

    dis2Points = 0
    id1 = 0
    if(cp.lines==0):
        return vec1
    if (cp.lineset[0][0] != cp.id):
        p1id1 = cp.lineset[0][0]
    elif (cp.lineset[0][1] != cp.id):
        p1id1 = cp.lineset[0][1]
    v1 = ske[p1id1].point - cp.point

    for i in range(1, neiborsize):
        idtemp = cp.Neibors[i]

        if(idtemp==p1id1):
            continue
        ptemp = ske[idtemp]
        if(ptemp.isend=="false"):
            continue
        labeltemp = ptemp.label
        if (cp.label != labeltemp):
            continue
        if (ptemp.lines == 2):
            continue
        v2=ptemp.point-cp.point
        # a=ptemp.lineset
        # aa=ptemp.lineset[0]
        id3=ptemp.lineset[0][1]
        v3=ptemp.point-ske[id3].point
        [dis, angleA, angleB] = computeVecDisAngle(cp.point, ptemp.point, v1, cp.direciton)
        # dis = computeDis(cp.point, ptemp.point)
        dis=math.fabs(dis)
        angle=computeAngleJiaodu(v3,v2)
        if(angle<135):
            continue
        if (dis == float("inf")):
            continue
        if (dis <= 0):
            continue
        dis2Points = dis
        if (dis2Points < min):
            min = dis2Points
            id1 = idtemp
    if (min != float("inf")):
        if(ske[id1].lines<2):
            vec1.append([cp.id, id1])
            cp.hasLine2 = "true"
            cp.lines = cp.lines + 1
            cp.lineset.append([cp.id, id1])
            if (ske[id1].hasLine1 == "true"):
                ske[id1].hasLine2 = "true"
                ske[id1].lineset.append([id1, cp.id])
                ske[id1].lines = ske[id1].lines + 1
            else:
                ske[id1].hasLine1 = "true"
                ske[id1].lineset.append([id1, cp.id])
                ske[id1].lines = ske[id1].lines + 1

    return vec1

def checkAllPts(skeTemp):
    for c in skeTemp:
        temp = []
        if (len(c.lineset) == 0):
            c.lines = 0
            c.lineset = []
            continue
        if (len(c.lineset) == 1):
            c.lines = 1
            continue
        if (len(c.lineset) == 2):
            for item in c.lineset:
                item2 = [item[1], item[0]]
                if not item in temp and not item2 in temp:
                    temp.append(item)
            c.lineset = temp
            c.lines = len(temp)

    return skeTemp

def breakedges(item,ske):
    breaklist=[]
    if(item.lines==2):
        id0=item.id
        l1=item.lineset[0]
        id1=l1[1]
        l2 = item.lineset[1]
        id2 = l2[1]
        v1 = ske[id1].point - ske[id0].point
        v2 = ske[id2].point - ske[id0].point
        angle=computeAngleJiaodu(v1,v2)
        if(angle<135):
            item.lines=0
            item.lineset.clear()
            breaklist.append([id0,id1])
            breaklist.append([id0,id2])
            ske[id1].lines-=1
            ske[id1].lineset.remove([id1, id0])
            ske[id2].lines -= 1
            ske[id2].lineset.remove([id2, id0])
    return breaklist

class lineCounts:
    def __init__(self,id):
        self.id=id
        self.linenums=0
        self.lineset=[]
        self.isSet=False
        self.disset=[]

def computeEdges(skeTemp):
    pts = []
    endSet = []
    for i in range(0, len(skeTemp)):
        tp = lineCounts(i)
        pts.append(tp)
    for i in range(0,len(skeTemp)):
        tp = skeTemp[i]
        if (tp.lines==0):
            continue
        for item in tp.lineset:
            id1=item[0]
            id2=item[1]
            if not [id1,id2] in pts[id1].lineset:
                pts[id1].linenums += 1
                pts[id1].lineset.append([id1,id2])
            if not [id2, id1] in pts[id2].lineset:
                pts[id2].linenums += 1
                pts[id2].lineset.append([id2,id1])
    return pts

def optiEdges(pts,skeTemp):
    ptsout=copy.deepcopy(pts)
    for item in ptsout:
        max=-float("inf")
        min=float("inf")
        saveminid = -1
        savemaxid = -1
        id = item.id
        if(item.linenums<=2):
            continue

        disset=[]
        for item2 in item.lineset:
            id1 = item2[0]
            id2 = item2[1]
            p1=skeTemp[id1].point
            p2=skeTemp[id2].point
            dire=skeTemp[id1].direciton
            dis=computeDisPN(p1,p2,dire)
            disset.append([dis,id2])
        for item2 in disset:
            if(item2[0]>0):
                if(item2[0]<min):
                    min=item2[0]
                    saveminid=item2[1]
        for item2 in disset:
            if (item2[0] < 0):
                if (item2[0] > max):
                    max = item2[0]
                    savemaxid = item2[1]
        if (min!=float("inf") and max !=-float("inf")):
            l1=[id,saveminid]
            l2=[id,savemaxid]
            item.linenums=2
            for item2 in item.lineset:
                if(item2[1]==saveminid or item2[1]==savemaxid):
                    continue
                ptsout[item2[1]].lineset.remove([item2[1],id])
            item.lineset.clear()
            item.lineset.append(l1)
            item.lineset.append(l2)
            item.disset.append(min)
            item.disset.append(max)
        elif(min==float("inf") and max !=-float("inf")):
            l2 = [id, savemaxid]
            for item2 in item.lineset:
                if(item2[1]==savemaxid):
                    continue
                ptsout[item2[1]].lineset.remove([item2[1],id])
            item.linenums = 1
            item.lineset.clear()
            item.lineset.append(l2)
            item.disset.append(max)
        elif (min != float("inf") and max == -float("inf")):
            l1 = [id, saveminid]
            for item2 in item.lineset:
                if(item2[1]==saveminid):
                    continue
                ptsout[item2[1]].lineset.remove([item2[1],id])
            item.linenums = 1
            item.lineset.clear()
            item.lineset.append(l1)
            item.disset.append(min)

    for item in ptsout:
        item.linenums=len(item.lineset)

    return ptsout

def updateOpti(ptsopti,skeTemp):
    for item in ptsopti:
        id=item.id
        skeTemp[id].lineset.clear()
        skeTemp[id].lineset=item.lineset
        skeTemp[id].lines=item.linenums

def connectSinglePfromBranch(branch,skeTemp):
    id1=branch[0][0]
    id2=branch[-1][1]
    endset=[id1,id2]
    addlist=[]
    anthre1=15
    anthre2=165
    for i in range(0,len(endset)):
        for idp in skeTemp[endset[i]].Neibors:
            id = endset[i]
            item = skeTemp[id]
            Vcdir = item.direciton
            ptemp=skeTemp[idp]
            if(id==idp):
                continue
            if (item.lines == 2):
                continue
            if(ptemp.lines==2):
                continue
            if(item.label!=ptemp.label):
                continue
            if(ptemp.lines==0):
                Vpdir = ptemp.direciton
                angle1 = computeAngleJiaodu(Vcdir, Vpdir, isfbs=True)
                if (angle1 > anthre1):
                    continue
                V0_1 = skeTemp[id].point - skeTemp[item.lineset[0][1]].point
                V0_2 = item.point - ptemp.point
                angle2 = computeAngleJiaodu(V0_1, V0_2, isfbs=False)
                if (angle2 < anthre2):
                    continue
                p1 = skeTemp[idp]
                skeTemp[id].lineset.append([id, idp])
                skeTemp[id].lines += 1
                skeTemp[id].isend = "false"
                skeTemp[idp].lineset.append([idp, id])
                skeTemp[idp].lines += 1
                skeTemp[idp].branchid = skeTemp[id].branchid
                skeTemp[idp].isend = "true"
                skeTemp[idp].issingle2end = True
                addlist.append([id, idp])
                p2=skeTemp[idp]
                if (id == endset[0]):
                    endset[0] = idp
                else:
                    endset[1] = idp
            elif(ptemp.lines==1):
                V0_1 = skeTemp[id].point - skeTemp[item.lineset[0][1]].point
                V0_2 = item.point - ptemp.point
                angle1 = computeAngleJiaodu(V0_1, V0_2, isfbs=False)
                if (angle1 < anthre2):
                    continue
                V3_2 = skeTemp[ptemp.lineset[0][1]].point - ptemp.point
                angle2 = computeAngleJiaodu(V3_2, V0_2, isfbs=False)
                if (angle2 < anthre2):
                    continue
                skeTemp[id].lineset.append([id, idp])
                skeTemp[id].lines += 1
                skeTemp[id].isend = "false"
                skeTemp[idp].lineset.append([idp, id])
                skeTemp[idp].lines += 1
                skeTemp[idp].branchid = skeTemp[id].branchid
                skeTemp[idp].isend = "true"
                addlist.append([id, idp])
                if (id == endset[0]):
                    endset[0] = idp
                else:
                    endset[1] = idp
    return addlist

def findBranches(input, pts_size):
    branchlist = []
    for item in input:
        item2 = [item[1], item[0]]
        if not item in branchlist and not item2 in branchlist:
            branchlist.append(item)
    pts = []
    endSet = []
    for i in range(0, pts_size):
        tp = BranchNode(i)
        pts.append(tp)
    for i in range(len(branchlist)):
        tp = branchlist[i]
        id1 = tp[0]
        id2 = tp[1]
        pts[id1].connectPnum += 1
        pts[id1].connectP.append(id2)
        pts[id2].connectPnum += 1
        pts[id2].connectP.append(id1)

    for i in range(0, pts_size):
        if (pts[i].connectPnum == 0):
            continue
        if (pts[i].connectPnum == 1):
            pts[i].connectP1 = pts[i].connectP[0]
            pts[i].isEnd = True
            endSet.append(pts[i])
        if (pts[i].connectPnum == 2):
            pts[i].connectP1 = pts[i].connectP[0]
            pts[i].connectP2 = pts[i].connectP[1]

    ptscopy = pts
    ptsout = []
    for item in endSet:
        id = item.id
        if (ptscopy[id].isInBranch == True):
            continue
        btemp = Branch()
        stemp = Stack()
        stemp.push(ptscopy[id])
        while (stemp.items != []):
            cnode = stemp.pop()
            idtemp = cnode.id
            ptscopy[idtemp].isInBranch = True
            btemp.pts.append(ptscopy[idtemp])
            if (ptscopy[idtemp].connectP1 != -1):
                idt1 = ptscopy[idtemp].connectP1
                if (ptscopy[idt1].isInBranch == False):
                    stemp.push(ptscopy[idt1])
            if (ptscopy[idtemp].connectP2 != -1):
                idt2 = ptscopy[idtemp].connectP2
                if (ptscopy[idt2].isInBranch == False):
                    stemp.push(ptscopy[idt2])
        if (btemp != []):
            ptsout.append(btemp)
    return ptsout

def updateBranch(branchesSet,skeTemp):
    for i in range(0,len(branchesSet)):
        cbran=branchesSet[i]
        for item in cbran:
            skeTemp[item[0]].branchid=i
            skeTemp[item[0]].isend = "false"
        skeTemp[cbran[-1][1]].branchid = i
        skeTemp[cbran[0][0]].isend = "true"
        skeTemp[cbran[-1][1]].isend = "true"

def resetSketempFromBranch(branchesSet,skeTemp):
    for item in skeTemp:
        item.lines=0
        item.lineset=[]
        item.branchid=-1
    for i in range(0,len(branchesSet)):
        cbran=branchesSet[i]
        for item in cbran:
            id1=item[0]
            id2=item[1]
            skeTemp[id1].branchid = i
            skeTemp[id2].branchid = i
            skeTemp[id1].lines += 1
            skeTemp[id2].lines += 1
            skeTemp[id1].lineset.append([id1,id2])
            skeTemp[id2].lineset.append([id2,id1])
        skeTemp[cbran[0][0]].isend = "true"
        skeTemp[cbran[-1][1]].isend = "true"
    for i in range (0,len(skeTemp)):
        if (skeTemp[i].lines==1):
            skeTemp[i].hasLine1="true"
            skeTemp[i].hasLine2="false"
        elif (skeTemp[i].lines==2):
            skeTemp[i].hasLine1 = "true"
            skeTemp[i].hasLine2 = "true"
        elif (skeTemp[i].lines==0):
            skeTemp[i].hasLine1 = "false"
            skeTemp[i].hasLine2 = "false"

def getEdgesFromBranch(connectlist,pts_size,skeTemp):

    branchlistout2 = findBranches(connectlist, pts_size)
    branchesSet2 = []
    ConnectLines2 = []

    for item in branchlistout2:
        edgestemp = item.getBranchedges()
        branchesSet2.append(edgestemp)
        for item2 in edgestemp:
            ConnectLines2.append(item2)
    updateBranch(branchesSet2, skeTemp)
    return branchesSet2,ConnectLines2,branchlistout2

class endConnect:
    def __init__(self):
        self.id1 = -1
        self.id2 = -1
        self.angle =float("inf")
        self.lineOk=False
    def __init__(self,id1,id2,angle,lineok):
        self.id1 = id1
        self.id2 = id2
        self.angle = angle
        self.lineOk = lineok

def connectBranchSameLabel(branch,skeTemp):
    conlist=[]
    id1 = branch[0][0]
    id1_spare = branch[0][1]
    p1 = skeTemp[id1]
    p1_spare=skeTemp[id1_spare]
    for idp1 in p1.Neibors:
        if(p1.lines==2 or p1.isend!="true"):
            break
        if (ptemp.label != p1.label):
            continue
        if (ptemp.lines == 2):
            continue
        if (ptemp.isend == "false"):
            continue
        ptemp=skeTemp[idp1]
        ptemp_spare=skeTemp[ptemp.lineset[0][1]]

def connectBigDis(endBigDis,endPts,skeTemp,anglethreshold):
    outlist=[]
    for item in endBigDis:
        cp=skeTemp[item[0]]
        if(cp.bigdisconnect==True):
            continue
        v1=item[-1]
        v2=cp.direciton

        endPts_Dis=[]
        for id2 in endPts:
            tp=skeTemp[id2]
            distemp=computeDis(cp.point,tp.point)
            endPts_Dis.append([id2,distemp])
        endPts_Dis2=sorted(endPts_Dis,key=lambda x: x[1])

        for item in endPts_Dis2:
            id2=item[0]
            tp=skeTemp[id2]
            if(tp.bigdisconnect==True):
                continue
            if(tp.branchid==cp.id):
                continue
            vt=cp.point-tp.point
            vt2=tp.point-cp.point
            angle3=180
            if(tp.lineset!=[]):
                id3=tp.lineset[0][1]
                vt3=tp.point-skeTemp[id3].point
                angle3 = computeAngleJiaodu(vt2, vt3)
            angle1=computeAngleJiaodu(v1,vt)
            angle2=computeAngleJiaodu(v2,vt)
            angle4=computeAngleJiaodu(cp.direciton,tp.direciton,True)
            if(angle1<anglethreshold or angle3<anglethreshold):
                continue
            if (angle4 > 15):
                continue
            dis=computeDis(cp.point,tp.point)
            dissin1 = dis * math.sin(angle4)
            dissin = dis * math.sin(angle2)
            discos1 = dis * math.cos(angle4)
            discos = dis * math.cos(angle2)

            if(dis>200):
                continue

            if(dissin<0.01):
                skeTemp[cp.id].bigdisconnect=True
                skeTemp[tp.id].bigdisconnect=True
                outlist.append([cp.id,tp.id])
                break
    return outlist

def connectSmallBranchBigDis(endPts,skeTemp,anglethreshold):
    outlist=[]
    for item in endPts:
        cp=skeTemp[item]
        if(cp.bigdisconnect==True):
            continue
        id12=cp.lineset[0][1]
        v1=cp.point-skeTemp[id12].point
        v2=cp.direciton
        endPts_Dis=[]
        for id2 in endPts:
            tp=skeTemp[id2]
            distemp=computeDis(cp.point,tp.point)
            endPts_Dis.append([id2,distemp])
        endPts_Dis2=sorted(endPts_Dis,key=lambda x: x[1])

        for item in endPts_Dis2:
            id2=item[0]
            tp=skeTemp[id2]
            if(tp.bigdisconnect==True):
                continue
            if(tp.branchid==cp.id):
                continue
            vt=cp.point-tp.point
            vt2=tp.point-cp.point
            angle3=180
            if(tp.lineset!=[]):
                id3=tp.lineset[0][1]
                vt3=tp.point-skeTemp[id3].point
                angle3 = computeAngleJiaodu(vt2, vt3)
            angle1=computeAngleJiaodu(v1,vt)
            angle2=computeAngleJiaodu(v2,vt)
            angle4=computeAngleJiaodu(cp.direciton,tp.direciton,True)
            if(angle1<anglethreshold or angle3<anglethreshold):
                continue
            if (angle4 > 30):
                continue
            dis=computeDis(cp.point,tp.point)
            dissin1 = dis * math.sin(angle4*math.pi/180)
            dissin = dis * math.sin(angle2*math.pi/180)
            discos1 = dis * math.cos(angle4*math.pi/180)
            discos = dis * math.cos(angle2*math.pi/180)

            if(dis>30):
                continue
            if(dissin<10):
                skeTemp[cp.id].bigdisconnect=True
                skeTemp[tp.id].bigdisconnect=True
                outlist.append([cp.id,tp.id])
                break
    return outlist

def connectBigDisBranch(endBigDis,endPts,skeTemp,anglethreshold):
    outlist=[]
    for item in endBigDis:
        cp=skeTemp[item[0]]
        v1=item[-1]
        v2=cp.direciton
        for id2 in endPts:
            tp=skeTemp[id2]
            if(tp.branchid==cp.id):
                continue
            vt=cp.point-tp.point
            vt2=tp.point-cp.point
            angle3=180
            if(tp.lineset!=[]):
                id3=tp.lineset[0][1]
                vt3=tp.point-skeTemp[id3].point
                angle3 = computeAngleJiaodu(vt2, vt3)
            angle1=computeAngleJiaodu(v1,vt)
            angle2=computeAngleJiaodu(v2,vt)
            if(angle1<anglethreshold or angle3<anglethreshold):
                continue
            dis=computeDis(cp.point,tp.point)
            dissin=dis*math.sin(angle2)
            if(dissin<0.01):
                outlist.append([cp.id,tp.id])
                break
    return outlist

def getfittingline(branch,skeTemp):
    id00 = branch[0][0]
    p00=skeTemp[id00]
    pts_set=[]
    pts_set.append(p00.point)
    for item in branch:
        idt=item[1]
        pts_set.append(skeTemp[idt].point)
    temp = np.asarray(pts_set)
    con, eigen_val, eigen_vec = computeCon(temp)
    dir_vec = [eigen_vec[0][0], eigen_vec[1][0], eigen_vec[2][0]]
    if (dir_vec[0] < 0):
        dir_vec[0] = -dir_vec[0]
        dir_vec[1] = -dir_vec[1]
        dir_vec[2] = -dir_vec[2]
    return dir_vec

def getfittinglineFromPts(pts,skeTemp):

    pts_sorted=sorted(pts,key=lambda x:x.posinbranch)
    dir_vec_set=[]
    for i in range(0,len(pts_sorted)-1):
        p1t = pts_sorted[i].point
        p2t = pts_sorted[i+1].point
        v1 = p1t[0] - p2t[0]
        v2 = p1t[1] - p2t[1]
        v3 = p1t[2] - p2t[2]
        dir_vec_set.append([v1,v2,v3])
    dir_vec=[0,0,0]
    for item in dir_vec_set:
        dir_vec[0] += item[0]
        dir_vec[1] += item[1]
        dir_vec[2] += item[2]
    dir_vec[0] = dir_vec[0]/(len(pts_sorted)-1)
    dir_vec[1] = dir_vec[1]/(len(pts_sorted)-1)
    dir_vec[2] = dir_vec[2]/(len(pts_sorted)-1)

    if (dir_vec[0] < 0):
        dir_vec[0] = -dir_vec[0]
        dir_vec[1] = -dir_vec[1]
        dir_vec[2] = -dir_vec[2]
    return dir_vec

def getNeiborBranchesPoints(branch,skeTemp,cloudin,pcd_tree):

    id00 = branch[0][0]
    p00 = skeTemp[id00]
    pts_set = []
    pts_set.append(p00)
    branchid00 = p00.branchid
    for item in branch:
        idt = item[1]
        pts_set.append(skeTemp[idt])
    neipts_id=[]
    for item in pts_set:
        id_in_src=item.id
        [k0, idx0, distance0] = pcd_tree.search_radius_vector_3d(cloudin.points[id_in_src], 10)
        for idj in idx0:
            if not idj in neipts_id:

                if (skeTemp[idj].branchid != -1 and skeTemp[idj].branchid != branchid00):
                    if(skeTemp[idj].isNegaCombined==True or skeTemp[idj].isCombined == False):
                        neipts_id.append(idj)
                    # if(skeTemp[idj].isCombined==False): #or skeTemp[idj].isNegaCombined==True):
    neipts=[]
    branchidtset=[]
    for id in neipts_id:
        branchidt=skeTemp[id].branchid
        if not branchidt in branchidtset:
            branchidtset.append(branchidt)
        neipts.append(skeTemp[id])
    branchidtset_sorted=sorted(branchidtset,reverse=True)
    neipts_sorted = sorted(neipts, key=lambda x: x.branchid,reverse=True)
    return branchidtset_sorted, neipts_sorted

def setCombined(neiBranchIDOut,branchset,skeTemp):
    lensets=[]
    for brid in neiBranchIDOut:
        lengtht=len(branchset[brid])
        lensets.append([brid,lengtht])
    len_sorted=sorted(lensets,key=lambda x: x[1],reverse=True)
    outbranid=len_sorted[0][0]
    for brid in neiBranchIDOut:
        item=branchset[brid]
        id0t=item[0][0]
        if(skeTemp[id0t].branchid==outbranid):
            skeTemp[id0t].isNegaCombined=True
        skeTemp[id0t].branchid = outbranid
        skeTemp[id0t].isCombined=True
        for item2 in item:
            idt=item2[1]
            if (skeTemp[idt].branchid == outbranid):
                skeTemp[idt].isNegaCombined = True
            skeTemp[idt].branchid = outbranid
            skeTemp[idt].isCombined = True

def combineNeighborBranch(branch1,branchset,skeTemp,cloudin,pcd_tree,id0):

    neiBranchIDOut = []
    fitting0=getfittingline(branch1, skeTemp)
    [branchsIdset,neipts_sorted]=getNeiborBranchesPoints(branch1, skeTemp, cloudin, pcd_tree)
    if(len(branchsIdset)==0):
        return neiBranchIDOut


    neiBranchIDOut.append(id0)
    for brid in branchsIdset:
        fitneipts=[]
        fittingt=[]
        for item in neipts_sorted:
            if(item.branchid==brid):
                fitneipts.append(item)
        if(len(fitneipts)<=1):
            continue
        elif(len(fitneipts)<10):
            fittingt = getfittinglineFromPts(fitneipts,skeTemp)
        else:
            branchneiTemp = branchset[brid]
            fittingt = getfittingline(branchneiTemp, skeTemp)
        angle=computeAngleJiaodu(fitting0,fittingt,True)
        if(angle<30):
            neiBranchIDOut.append(brid)
    if(len(neiBranchIDOut)<=1):
        return []
    setCombined(neiBranchIDOut,branchset,skeTemp)
    return neiBranchIDOut

def findClosestPointFromOtherBranch(p,skeTemp,cloudin,pcd_tree):
    outp=[]
    id_in_src=p.id
    [k0, idx0, distance0] = pcd_tree.search_knn_vector_3d(cloudin.points[id_in_src], 100)
    for id in idx0:
        pt=skeTemp[id]
        if(pt.id==id_in_src):
            continue
        if(pt.branchid==p.branchid):
            if(pt.isNegaCombined==True):
                outp=pt
                break
    return outp

def insertBranch(combinBranch,branchesSet9,skeTemp,cloudin,pcd_tree):
    brancht = branchesSet9[0]
    beinsertedID = skeTemp[brancht[0][0]].branchid
    for item in combinBranch:
        if(item==beinsertedID):
            continue
        p1id = branchesSet9[item][0][0]
        p2id = branchesSet9[item][-1][1]
        p1 = skeTemp[p1id]
        p2 = skeTemp[p2id]
        outp1 = findClosestPointFromOtherBranch(p1, skeTemp, cloudin, pcd_tree)
        outp2 = findClosestPointFromOtherBranch(p2, skeTemp, cloudin, pcd_tree)
        if(outp1==[] and outp2==[]):
            continue
        elif(outp1==[] and outp2!=[]):
            insertbranch = copy.deepcopy(branchesSet9[item])
            if(outp2.posinbranch==0):
                insertbranch.append(0, [p2.id, outp2.id])
                branchesSet9[beinsertedID][0:0] = insertbranch
            else:
                insertbranch.reverse()
                insertbranch.insert(0, [outp2.id, p2.id])
                for item2 in insertbranch:
                    branchesSet9[beinsertedID].append(item2)
        elif(outp1!=[] and outp2==[]):
            insertbranch = copy.deepcopy(branchesSet9[item])
            if (outp1.posinbranch == 0):
                insertbranch.reverse()
                insertbranch.append(0, [p1.id, outp1.id])
                branchesSet9[beinsertedID][0:0] = insertbranch
            else:
                insertbranch.insert(0, [outp1.id, p1.id])
                for item2 in insertbranch:
                    branchesSet9[beinsertedID].append(item2)
        else:
            if (outp1.posinbranch < outp2.posinbranch):
                insertbranch = copy.deepcopy(branchesSet9[item])
                insertbranch.insert(0, [outp1.id, p1.id])
                insertbranch.append([p2.id, outp2.id])
                branchesSet9[beinsertedID][outp1.posinbranch + 1:outp1.posinbranch + 1] = insertbranch
            else:
                insertbranch = copy.deepcopy(branchesSet9[item])
                insertbranch.reverse()
                insertbranch.insert(0, [outp2.id, p2.id])
                insertbranch.append([p1.id, outp1.id])
                branchesSet9[beinsertedID][outp2.posinbranch + 1:outp2.posinbranch + 1] = insertbranch

def insertBranches(combineBranchesSet,branchesSet9,skeTemp,cloudin,pcd_tree):
    for item in combineBranchesSet:
        insertBranch(item,branchesSet9,skeTemp,cloudin,pcd_tree)

def deleteRedundancyBranch(branchesSet9,skeTemp):
    out=[]
    for brid in range(0,len(branchesSet9)):
        item=branchesSet9[brid]
        id0t=item[0][0]
        pt=skeTemp[id0t]
        if(pt.branchid!=brid):
            continue
        out.append(item)
    return out

def comuteBranchLen(branch,skeTemp):
    dis=0
    for item in branch:
        p1=skeTemp[item[0]].point
        p2=skeTemp[item[1]].point
        dist=computeDis(p1,p2)
        dis+=dist
    return dis

def connectBranches_sameLabels(branchesset,skeTemp):
    endlist=[]
    for item in branchesset:
        id1=item[0][0]
        id2=item[-1][1]
        endlist.append(skeTemp[id1])
        endlist.append(skeTemp[id2])

def findClosestPts(id1_1,id2_1,it1_1,it2_1,skeTemp):
    p1=skeTemp[id1_1].point
    p2=skeTemp[id2_1].point
    pt1=skeTemp[it1_1].point
    pt2=skeTemp[it2_1].point
    d1 = computeDis(p1, pt1)
    d2 = computeDis(p1, pt2)
    d3 = computeDis(p2, pt1)
    d4 = computeDis(p2, pt2)
    dmin=min(d1,d2,d3,d4)
    if(dmin==d1):
        outdisid = [id1_1,it1_1]
    elif(dmin==d2):
        outdisid = [id1_1,it2_1]
    elif(dmin==d3):
        outdisid = [id2_1,it1_1]
    else:
        outdisid = [id2_1,it2_1]
    return outdisid[0],outdisid[1],dmin

def connectBranches(endBlist,branchesset,skeTemp):
    athre=135
    endDelete=[]
    endConnect=[]

    for item in endBlist:
        [id1_1,id1_2,id1_3,id2_1,id2_2,id2_3]=item
        dislist = []
        for item2 in endBlist:
            if (item2==item):
                continue
            [it1_1, it1_2,it1_3,it2_1, it2_2,it2_3] = item2
            endid1_3 = id1_3
            endid2_3 = it1_3
            [endid1,endid2,dmin]=findClosestPts(id1_1,id2_1,it1_1,it2_1,skeTemp)
            if(dmin>=10):
                continue
            dislist.append(dmin)
            if(endid1==id2_1):
                endid1_3=id2_3
            if(endid2==it2_1):
                endid2_3=it2_3
            p1=skeTemp[endid1]
            p1_2=skeTemp[p1.lineset[0][1]]
            p2 = skeTemp[endid2]
            p2_2 = skeTemp[p2.lineset[0][1]]
            vp1p2=p1.point-p2.point
            vp2p1=p2.point-p1.point
            vp1p12=p1.point-p1_2.point
            vp2p22=p2.point-p2_2.point
            vp12p22 = p1_2.point - p2_2.point
            vp22p12 = p2_2.point - p1_2.point
            vp12p13= p1_2.point-skeTemp[endid1_3].point
            vp22p23=p2_2.point-skeTemp[endid2_3].point
            a1 = computeAngleJiaodu(vp1p2, vp1p12)
            a2 = computeAngleJiaodu(vp2p1, vp2p22)
            a3 = computeAngleJiaodu(vp12p22,vp12p13)
            a4 = computeAngleJiaodu(vp22p12,vp22p23)
            if(a3<athre or a4 <athre):
                continue
            if(a1<athre or a2 <athre):
                endDelete.append([endid1, p1.lineset[0][1]])
                endDelete.append([endid2, p2.lineset[0][1]])
                endConnect.append([p1.lineset[0][1],p2.lineset[0][1]])
            else:
                endConnect.append([endid1,endid2])
    return endConnect,endDelete

def connectBranches_DiffLabels(branches,skeTemp):
    conlist=[]
    angleThreshold=165
    for i in range(0,len(branches)):
        id1=branches[i][0][0]
        id1_spare = branches[i][0][1]
        p1 = skeTemp[id1]
        if(p1.lines==2):
            continue
        p1_spare = skeTemp[id1_spare]
        for idp1 in p1.Neibors:
            ptemp=skeTemp[idp1]
            if(ptemp.branchid==p1.branchid):
                continue
            if(ptemp.lines==2):
                continue
            if(ptemp.isend=="false"):
                continue
            pp_spare=skeTemp[ptemp.lineset[0][1]]
            vp1_spare = p1.point - p1_spare.point
            vpp_spare = ptemp.point - pp_spare.point
            vp1_pp=p1.point-ptemp.point
            vpp_p1=ptemp.point-p1.point
            angle1 = computeAngleJiaodu(vp1_spare, vp1_pp, False)
            angle2 = computeAngleJiaodu(vpp_spare, vpp_p1, False)
            if(angle1>angleThreshold and angle2 >angleThreshold):
                skeTemp[id1].lineset.append([id1,idp1])
                skeTemp[id1].lines+=1;
                skeTemp[id1].isend="false"
                skeTemp[idp1].lineset.append([idp1,id1])
                skeTemp[idp1].lines += 1;
                skeTemp[idp1].isend = "false"
                conlist.append([id1,idp1])

    for i in range(0, len(branches)):

        id2 = branches[i][-1][1]
        id2_spare = branches[i][-1][0]
        p2 = skeTemp[id2]
        if (p2.lines == 2):
            continue
        p2_spare = skeTemp[id2_spare]
        for idp2 in p2.Neibors:
            ptemp = skeTemp[idp2]
            if (ptemp.branchid == p2.branchid):
                continue
            if (ptemp.lines == 2):
                continue
            if (ptemp.isend == "false"):
                continue
            pp_spare = skeTemp[ptemp.lineset[0][1]]
            vp2_spare = p2.point - p2_spare.point
            vpp_spare = ptemp.point - pp_spare.point
            vp2_pp = p2.point - ptemp.point
            vpp_p2 = ptemp.point - p2.point
            angle1 = computeAngleJiaodu(vp2_spare, vp2_pp, False)
            angle2 = computeAngleJiaodu(vpp_spare, vpp_p2, False)

            if (angle1 > angleThreshold and angle2 > angleThreshold):
                skeTemp[id2].lineset.append([id2, idp2])
                skeTemp[id2].lines += 1;
                skeTemp[id2].isend = "false"
                skeTemp[idp2].lineset.append([idp2, id2])
                skeTemp[idp2].lines += 1;
                skeTemp[idp2].isend = "false"
                conlist.append([id2, idp2])
        return conlist

def updateBranchDireVec(branch,skeTemp):
    branchptsnum=len(branch)
    for i in range(0, branchptsnum):
        if(i==0 or i==branchptsnum-1):
            id1 = branch[i][0]
            id2 = branch[i][1]
            x = skeTemp[id1].point[0] - skeTemp[id2].point[0]
            y = skeTemp[id1].point[1] - skeTemp[id2].point[1]
            z = skeTemp[id1].point[2] - skeTemp[id2].point[2]
            if(x<0):
                x=-x
                y=-y
                z=-z
            pmod=math.sqrt(x*x+y*y+z*z)
            if(i==0):
                skeTemp[id1].direciton = [x/pmod,y/pmod,z/pmod]
            else:
                skeTemp[id2].direciton = [x/pmod,y/pmod,z/pmod]
            continue
        id1=branch[i-1][0]
        id2=branch[i][1]
        idc=branch[i][0]
        x = skeTemp[id1].point[0] - skeTemp[id2].point[0]
        y = skeTemp[id1].point[1] - skeTemp[id2].point[1]
        z = skeTemp[id1].point[2] - skeTemp[id2].point[2]
        if (x < 0):
            x = -x
            y = -y
            z = -z
        skeTemp[idc].direciton=[x/pmod,y/pmod,z/pmod]

def resetBranchID(branchesSet9,skeTemp):
    for item in skeTemp:
        item.branchid=-1
    for brid in range(0,len(branchesSet9)):
        item=branchesSet9[brid]
        p0id=item[0][0]
        skeTemp[p0id].branchid=brid
        for item2 in item:
            ptid=item2[1]
            skeTemp[ptid].branchid=brid

def setPointPosInBranch(branch,skeTemp):
    id00=branch[0][0]
    skeTemp[id00].posinbranch=0
    for i in range(0,len(branch)):
        idt=branch[i][1]
        skeTemp[idt].posinbranch=i+1


class skePoint:
    def __init__(self, point, id, id_in_cluster, color, dir_vec, visited, neibors):
        self.point = point
        self.id = id
        self.direciton = dir_vec
        self.isvisit = visited
        self.id_in_cluster = id_in_cluster
        self.color = color

        self.Neibors = neibors

        self.hasLine1 = "false"
        self.hasLine2 = "false"
        self.issame = "false"
        self.lines = 0
        self.isend = "false"
        self.endconnect = "false"
        self.line1 = []
        self.line2 = []
        self.lineset = []
        self.label = -1
        self.branchid = -1
        self.isnoise=False
        self.issingle2end=False
        self.bigdisconnect=False
        self.isCombined=False
        self.isNegaCombined=False
        self.posinbranch=-1

    def setVisit(self, visited):
        self.isvisit = visited

    def setNeibors(self, neibors):
        self.Neibors = neibors

    def setLineSet(self):
        if (self.lines == 1):
            self.line1.append(self.lineset[0])
        elif (self.lines == 2):
            self.line1.append(self.lineset[0])
            self.line2.append(self.lineset[1])

    def setLabel(self, label):
        self.label = label

class Stack(object):

    def __init__(self):
        self.items = []

    def is_empty(self):
        return self.items == []

    def push(self, item):
        self.items.append(item)

    def pop(self):
        return self.items.pop()

    def peek(self):
        return self.items[len(self.items) - 1]

    def size(self):
        return len(self.items)

class BranchNode:
    def __init__(self, id):
        self.id = id
        self.connectP = []
        self.connectP1 = -1
        self.connectP2 = -1
        self.isInBranch = False
        self.connectPnum = 0
        self.isEnd = False
        self.edges = []
        self.branchId = -1

    def getEdges(self):
        if (self.connectPnum == 0):
            self.edges = []
        if (self.connectPnum == 1):
            self.edges = [[self.id, self.connectP1]]
        if (self.connectPnum == 2):
            self.edges.append([self.id, self.connectP1])
            self.edges.append([self.id, self.connectP2])
        return self.edges

class Branch:
    def __init__(self, nodes):
        self.pts = nodes
        self.edges = []
        self.id = -1
        self.isvisited=False

    def __init__(self):
        self.pts = []
        self.edges = []
        self.id = -1
        self.isvisited = False

    def sort(self):
        self.pts.sort(key=lambda x: x.id)

    def getBranchedges(self):
        list = []
        for item in self.pts:
            edges = item.getEdges()
            for e in edges:
                list.append(e)
        listout = []
        for item1 in list:
            item2 = [item1[1], item1[0]]
            if not item1 in listout and not item2 in listout:
                listout.append(item1)
        self.edges = listout
        return listout

    def setBranchID(self, id):
        self.id = id
        for item in self.pts:
            item.branchId = id

def file_name(file_dir,filetype):
    names=[]
    for root, dirs, files in os.walk(file_dir):
        for i in files:
            if os.path.splitext(i)[1] == filetype:
                names.append(i)
    purenames = []
    for i in range(len(names)):
        if not "-DB" in names[i]:
            names.append(names[i])
            t = names[i].split('.')[0]
            purenames.append(t)
    return purenames

if __name__ == "__main__":

    filepath = "F:/1-DesktopFile_Pointget/01Silique/SHScodetest/ske/"
    names = file_name("F:/1-DesktopFile_Pointget/01Silique/SHScodetest/ske", ".txt")


    eachsiliquenums=[]
    for ijj in range(0,10):#len(names)):
        itemname=names[ijj]
        print(itemname)
        path = filepath + itemname + ".ply"
        intxtpath = filepath  + itemname + "-DB.txt"
        path2 = filepath[:-5] + itemname[:-4] + ".ply"
        outtxtpath = filepath[:-5] + itemname+ "_classfied.txt"

        cloudsrc = open3d.io.read_point_cloud(path2)
        cloudin = open3d.io.read_point_cloud(path)
        points = np.asarray(cloudin.points)
        pts_size = int(points.size / 3)
        sortedtxtpath = intxtpath
        [SortedPts, psize] = loadDatadet(sortedtxtpath)
        colors = giveColor(SortedPts, psize)

        for i in range(0, pts_size):
            idtemp = SortedPts[i][4]
            labeltemp = SortedPts[i][3]
            if (labeltemp == -1):
                c=idtemp
                cloudin.points[idtemp][0] = 0
                cloudin.points[idtemp][1] = 0
                cloudin.points[idtemp][2] = 0

        # give color
        for i in range(0, psize):
            id = int(colors[i][4])
            cloudin.colors[id][0] = colors[i][0]
            cloudin.colors[id][1] = colors[i][1]
            cloudin.colors[id][2] = colors[i][2]

        pts_src = np.array(points)
        pts_src = (pts_src.T)[0:3].T
        pts_test = np.copy(np.array(pts_src))
        pcd_tree = open3d.geometry.KDTreeFlann(cloudin)

        # open3d
        [h1, m1, s1] = tools.getClockTime()
        ms1 = tools.getMiroTimestick()

        direction = []
        skeSetPts = []
        skesize = 0
        labelini = int(colors[0][3])
        skeTemp = []
        distvec = []
        for i in range(0, pts_size):
            idinsrc = int(colors[i][4])
            label = int(colors[i][3])
            [k, idx, distance] = pcd_tree.search_knn_vector_3d(cloudin.points[idinsrc], 5)

            temppoints = []
            tempdis = []
            for j in range(0, k):
                temppoints.append(cloudin.points[idx[j]])
                tempdis.append(distance[j])
            distvec.append(tempdis)
            temp = np.asarray(temppoints)
            con, eigen_val, eigen_vec = computeCon(temp)
            lprecent, maxlamda = computeLamda(eigen_val)
            dir_vec = [eigen_vec[0][0], eigen_vec[1][0], eigen_vec[2][0]]
            if (dir_vec[0] < 0):
                dir_vec[0] = -dir_vec[0]
                dir_vec[1] = -dir_vec[1]
                dir_vec[2] = -dir_vec[2]

            direction.append(dir_vec)
            ptemp = skePoint(cloudin.points[idinsrc], idinsrc, i, cloudin.colors[idinsrc], dir_vec, "false", [])
            ptemp.setLabel(label)
            skeTemp.append(ptemp)
            skesize += 1
        skeTemp.sort(key=lambda x: x.id)

        for i in range(0, pts_size):
            a = distvec[i]
            distvec[i].sort()

        for i in range(0, pts_size):
            if (skeTemp[i].direciton[0] < 0):
                skeTemp[i].direciton[0] = -skeTemp[i].direciton[0]
                skeTemp[i].direciton[1] = -skeTemp[i].direciton[1]
                skeTemp[i].direciton[2] = -skeTemp[i].direciton[2]

        # Dij``````````````
        connectlistBefore = []
        for i in range(0, pts_size):
            id = skeTemp[i].id
            [k, idx, distance] = pcd_tree.search_knn_vector_3d(cloudin.points[id], 25)
            skeTemp[i].setNeibors(idx)
            neiborszie = len(skeTemp[i].Neibors)
            [graphtemp, Eugraph, anglegraph] = bulid_graph_dire_singleP(skeTemp[i], neiborszie, skeTemp, cloudin)
            [vec1, vec2] = dij_dire_singleP_withClusterLabel(skeTemp[i], neiborszie, skeTemp, graphtemp, Eugraph,
                                                             anglegraph)
            if (vec1 != []):
                connectlistBefore.append(vec1[0])
            if (vec2 != []):
                connectlistBefore.append(vec2[0])
        lines = line_pcds2(connectlistBefore, pts_size, cloudin)
        # vis.add_geometry(lines)
        pedges = computeEdges(skeTemp)
        ptsopti = optiEdges(pedges, skeTemp)
        updateOpti(ptsopti, skeTemp)

        checkAllPts(skeTemp)
        print("Check initial")

        # 剔除异常点，导致多出了孤点
        breaklist = []
        nums = 0
        for item in skeTemp:
            btemp = breakedges(item, skeTemp)
            if (btemp != []):
                breaklist.append(btemp[0])
                breaklist.append(btemp[1])
                nums += 1
        for item in skeTemp:
            if (item.lines == 1):
                item.hasLine1 = "true"
                item.hasLine2 = "false"
            elif (item.lines == 2):
                item.hasLine1 = "true"
                item.hasLine2 = "true"
            elif (item.lines == 0):
                item.hasLine1 = "false"
                item.hasLine2 = "false"
        connectlist = []
        for item in skeTemp:
            if (item.lineset == []):
                continue
            for item2 in item.lineset:
                connectlist.append(item2)
        lines3 = line_pcds2(connectlist, pts_size, cloudin)

        checkAllPts(skeTemp)
        print("Check optimization and break")

        branchlistout = findBranches(connectlist, pts_size)
        branchesSet = []
        ConnectLines = []
        for item in branchlistout:
            edgestemp = item.getBranchedges()
            branchesSet.append(edgestemp)
            for item2 in edgestemp:
                ConnectLines.append(item2)
        updateBranch(branchesSet, skeTemp)
        linesske = line_pcds2(ConnectLines, pts_size, cloudin)
        singleLinesList = []
        for item in branchesSet:
            addlist = connectSinglePfromBranch(item, skeTemp)
            if (addlist != []):
                for item2 in addlist:
                    connectlist.append(item2)
                linesadd = line_pcds2(addlist, pts_size, cloudin)

        checkAllPts(skeTemp)
        connectlist2 = []
        for item in skeTemp:
            if (item.lineset == []):
                continue
            for item2 in item.lineset:
                connectlist2.append(item2)
        [branchesSet2, ConnectLines2, branchlistout2] = getEdgesFromBranch(connectlist, pts_size, skeTemp)
        updateBranch(branchesSet2, skeTemp)
        linesske2 = line_pcds2(ConnectLines2, pts_size, cloudin)

        deletelines = []
        invalidnums = 0
        for item in branchesSet2:
            if (len(item) <= 2):
                invalidnums += 1
                for item2 in item:
                    deletelines.append(item2)
                    for item3 in item2:
                        skeTemp[item3].point = [0, 0, 0]
            else:
                [r, g, b] = creatcolor()
                for item2 in item:
                    for item3 in item2:
                        skeTemp[item3].color = [r, g, b]
        for item in skeTemp:
            if (item.lines == 0):
                skeTemp[item.id].point = [0, 0, 0]
        print("invalid branches is:", invalidnums)

        AA = cloudin.points[0]
        for i in range(0, pts_size):
            if (all(skeTemp[i].point == np.asarray([0, 0, 0]))):
                cloudin.points[i] = [0, 0, 0]
            if (cloudin.points[i][0] == 0 and cloudin.points[i][1] == 1 and cloudin.points[i][2] == 0):
                cloudin.colors[i] = [1, 1, 1]
            else:
                cloudin.colors[i] = skeTemp[i].color


        endBlist = []
        for item in branchesSet2:
            if (len(item) <= 2):
                continue
            id1_1 = item[0][0]
            id1_2 = item[0][1]
            id1_3 = item[1][1]
            id2_1 = item[-1][1]
            id2_2 = item[-1][0]
            id2_3 = item[-2][0]
            endBlist.append([id1_1, id1_2, id1_3, id2_1, id2_2, id2_3])
        [endsConnect, endsDelete] = connectBranches(endBlist, branchesSet2, skeTemp)
        for item in endsDelete:
            if not item in deletelines:
                deletelines.append(item)
        ConnectLines3 = []
        for item in ConnectLines2:
            item2 = [item[1], item[0]]
            if not item in deletelines and not item2 in deletelines:
                ConnectLines3.append(item)
        for item in endsConnect:
            ConnectLines3.append(item)
        lines3 = line_pcds2(ConnectLines3, pts_size, cloudin)

        [branchesSet3, ConnectLines4, branchlistout3] = getEdgesFromBranch(ConnectLines3, pts_size, skeTemp)
        lines4 = line_pcds2(ConnectLines4, pts_size, cloudin)
        updateBranch(branchesSet3, skeTemp)
        resetSketempFromBranch(branchesSet3, skeTemp)

        bSet3_2 = sorted(branchesSet3, key=lambda i: len(i), reverse=True)

        for i in range(0, len(skeTemp)):
            tp = skeTemp[i]
            if (skeTemp[i].branchid == -1):
                continue
            id = tp.id
            [k, idx, distance] = pcd_tree.search_knn_vector_3d(cloudin.points[id], 30)
            temppoints = []
            for j in range(0, k):
                if (skeTemp[idx[j]].branchid != tp.branchid):
                    continue
                temppoints.append(cloudin.points[idx[j]])
            temp = np.asarray(temppoints)
            if (len(temp) <= 3):
                continue
            con, eigen_val, eigen_vec = computeCon(temp)
            dir_vec = [eigen_vec[0][0], eigen_vec[1][0], eigen_vec[2][0]]
            if (dir_vec[0] < 0):
                dir_vec[0] = -dir_vec[0]
                dir_vec[1] = -dir_vec[1]
                dir_vec[2] = -dir_vec[2]
            skeTemp[i].direciton = dir_vec
            tempx3 = cloudin.points[id][0] + eigen_vec[0][0] * 10
            tempy3 = cloudin.points[id][1] + eigen_vec[1][0] * 10
            tempz3 = cloudin.points[id][2] + eigen_vec[2][0] * 10
            p3 = [tempx3, tempy3, tempz3]
            linesx2 = line_pcd(cloudin.points[id], p3)

        endBigDis = []
        for item in bSet3_2:
            if (len(item) >= 50):
                id11 = item[0][0]
                id12 = item[0][1]
                id13 = item[1][1]
                id21 = item[-1][1]
                id22 = item[-1][0]
                id23 = item[-2][0]
                v1 = skeTemp[id11].point - skeTemp[id12].point
                v2 = skeTemp[id21].point - skeTemp[id22].point
                endBigDis.append([id11, id12, id13, v1])
                endBigDis.append([id21, id22, id23, v2])
        endPts = []
        for item in bSet3_2:
            id1 = item[0][0]
            id2 = item[-1][1]
            endPts.append(id1)
            endPts.append(id2)

        anthre2 = 90
        eelist = connectBigDis(endBigDis, endPts, skeTemp, anthre2)

        ConnectLines5 = copy.deepcopy(ConnectLines4)
        for item in eelist:
            ConnectLines5.append(item)

        [branchesSet6, ConnectLines6, branchlistout6] = getEdgesFromBranch(ConnectLines5, pts_size, skeTemp)
        updateBranch(branchesSet6, skeTemp)
        resetSketempFromBranch(branchesSet6, skeTemp)
        lines6 = line_pcds2(ConnectLines6, pts_size, cloudin)

        bSet6 = sorted(branchesSet6, key=lambda i: len(i), reverse=True)
        endPts2 = []
        for item in bSet6:
            id1 = item[0][0]
            id2 = item[-1][1]
            endPts2.append(id1)
            endPts2.append(id2)
        anthre3 = 145
        eelist2 = connectSmallBranchBigDis(endPts2, skeTemp, anthre3)
        ConnectLines7 = copy.deepcopy(ConnectLines6)
        for item in eelist2:
            ConnectLines7.append(item)
        [branchesSet8, ConnectLines8, branchlistout8] = getEdgesFromBranch(ConnectLines7, pts_size, skeTemp)
        updateBranch(branchesSet8, skeTemp)
        resetSketempFromBranch(branchesSet8, skeTemp)
        lines8 = line_pcds2(ConnectLines8, pts_size, cloudin)


        branchesSet9 = sorted(branchesSet8, key=lambda x: len(x), reverse=False)
        resetBranchID(branchesSet9, skeTemp)
        for item in branchesSet9:
            setPointPosInBranch(item, skeTemp)


        combineBranchesSet = []
        for brid in range(0, len(branchesSet9) - 1):
            item = branchesSet9[brid]
            id00 = item[0][0]
            branchid00 = skeTemp[id00].branchid
            combineBranches = combineNeighborBranch(item, branchesSet9, skeTemp, cloudin, pcd_tree, brid)
            if (combineBranches == []):
                continue
            if (len(combineBranches) <= 1):
                continue
            combineBranchesSet.append(combineBranches)

        insertBranches(combineBranchesSet, branchesSet9, skeTemp, cloudin, pcd_tree)
        branchesSet10 = deleteRedundancyBranch(branchesSet9, skeTemp)

        updateBranch(branchesSet10, skeTemp)

        branchlenset = []
        for item in branchesSet10:
            distemp = comuteBranchLen(item, skeTemp)
            branchlenset.append(distemp)

        branchwithLen = []
        for i in range(0, len(branchesSet10)):
            bran = branchesSet10[i]
            dis = branchlenset[i]
            branchwithLen.append([bran, dis])

        branchwithLen.sort(key=lambda x: x[1])
        brlen=[bl[-1] for bl in branchwithLen]
        stdt=np.std(brlen,ddof=1)
        avelen=np.mean(brlen)
        mediapos = int(len(branchwithLen) / 2)
        mediaLen = branchwithLen[mediapos][1]
        brlen2=[]
        for itemt in brlen:
            if itemt > 2*avelen:
                continue
            else:
                brlen2.append(itemt)
        avelen2=np.mean(brlen2)

        stdtheta=np.std(brlen2,ddof = 1)

        uselessbranch = 0;
        dishist=[]
        # for i in range(0, len(branchesSet10)):
        for i in range (len(brlen2)):
            # dis = branchwithLen[i][1]
            dis = branchwithLen[i][1]
            # if (dis <= mediaLen - 20 or dis >= mediaLen + 40):
            minthreh=0
            if mediaLen-20<0 or mediaLen - 1.96 * stdtheta<0:
                minthreh=0.1*mediaLen
            else:
                minthreh=mediaLen-20
                minthreh = mediaLen - 1.96 * stdtheta
            maxthreh1=1.5*mediaLen + 40
            maxthreh2=mediaLen+1.96*stdtheta
            maxthreh=max(maxthreh1,maxthreh2)
            maxthreh=maxthreh2
            if (dis <= minthreh or dis <=15 or dis >= maxthreh):
            # dis = brlen[i]
            # if (dis <= 15 or dis >= 1.5*mediaLen+40):
                uselessbranch += 1
                # dishist.append(0)
            else:
                dishist.append(dis)


        print("mid dis: ", mediaLen, "ave dis: ", avelen2)
        print("All siliques and stems number is: ", len(branchesSet10), " and siliques number is: ",
              len(dishist))
        [h2, m2, s2] = tools.getClockTime()
        [ho, mo, so] = tools.getCostTime([h1, m1, s1], [h2, m2, s2])
        ms2 = tools.getMiroTimestick()
        print("Time Use：", ho, "h", mo, "m", so, "s")
        print("Time Use：", round(float((ms2 - ms1) * 1.0 / 1000), 2))
        outinfor = itemname + " All branches: " + str(len(branchesSet10)) + " Siliques : " + str(len(dishist)) + '\n '
        # outinfor=itemname+" All branches: "+ str(len(branchesSet10))+" Siliques : " +str(len(branchesSet10) - uselessbranch)+'\n '
        outinfor=outinfor+"mid dis: "+ str(mediaLen)+ " ave dis: "+ str(avelen2) + "\n "
        outinfor=outinfor+" Time used: "+ str(round(float((ms2 - ms1) * 1.0 / 1000), 2))
        eachsiliquenums.append([outinfor])

        for i in range(0, len(branchesSet10)):
            item = branchesSet10[i]
            if (len(item) <= 9 and branchlenset[i] <= 15):
                continue
            linest = line_pcds2(item, pts_size, cloudin)


        # set color
        for i in range(0, pts_size):
            id = skeTemp[i].id
            if (skeTemp[i].lines == 0):
                cloudin.points[id] = [0, 0, 0]
                cloudin.colors[id] = [1, 1, 1]

        branchcolors = []
        for i in range(0, len(branchesSet9)):
            [r, g, b] = creatcolor()
            branchcolors.append([r, g, b])

        for item in branchesSet10:
            id0 = skeTemp[item[0][0]].id
            branchid = skeTemp[item[0][0]].branchid
            [r, g, b] = branchcolors[branchid]
            cloudin.colors[id0] = [r, g, b]
            for item2 in item:
                idt = skeTemp[item2[1]].id
                cloudin.colors[idt] = [r, g, b]

        outSkePointCloud = []
        for item in branchesSet10:
            id0 = skeTemp[item[0][0]].id
            branchid = skeTemp[item[0][0]].branchid
            p = skeTemp[item[0][0]]
            pos = str(p.point[0]) + ',' + str(p.point[1]) + ',' + str(p.point[2])
            color = str(p.color[0]) + ',' + str(p.color[1]) + ',' + str(p.color[2])
            op = [pos, color, str(branchid)]
            outSkePointCloud.append(op)
            for item2 in item:
                tbranchid = skeTemp[item2[1]].branchid
                tp = skeTemp[item2[1]]
                tpos = str(tp.point[0]) + ',' + str(tp.point[1]) + ',' + str(tp.point[2])
                tcolor = str(tp.color[0]) + ',' + str(tp.color[1]) + ',' + str(tp.color[2])
                top = [tpos, tcolor, str(tbranchid)]
                outSkePointCloud.append(top)
        skeTempOut = sorted(skeTemp, key=lambda x: x.branchid)


        f = open(outtxtpath, 'w')
        for item in outSkePointCloud:
            f.write(item[0] + " " + item[1] + " " + item[2] + '\n')
        f.close
       

        n, bins, patches = plt.hist(dishist, bins=100, normed=1, facecolor='grey', alpha=0.75)
        # n, bins, patches = plt.hist(brlen2, bins=100, normed=1, facecolor='grey', alpha=0.75)
        # plt.show()

        print("==================================\n")
        print("\n")
    print("Finish")

    for item in eachsiliquenums:
        print(item)




