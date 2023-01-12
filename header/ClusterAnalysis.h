#pragma once
#include <iostream>
#include <cmath>
#include "DBscan.h"
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "skeleton.h"

using namespace Eigen;
using namespace std;
using namespace ske;

class ClusterAnalysis
{
public:
	vector<DataPoint> dataset; 
	vector<DataPoint> sortedDataset; 
	int dimNum; 
	float radius; 
	int dataNum;  
	int minpts;  
	int kernelClusterId; 
	vector<float> eigenvec;
	float eigenvalue;
	float GetDistance(DataPoint &p1, DataPoint &p2); 
	float GetMaDistance(DataPoint &p1, DataPoint &p2); 
	Eigen::Matrix3d ComputeGeneralConvariance(); 
	Eigen::Matrix3d Convariance; 
	Eigen::Matrix3d convariance_inverse; 
	float NormalizePara; 
	bool isNormalized; 
	void SetArrivalPoints(DataPoint &p); 
	void KeyPointCluster(int i, int clusterId);
public:
	ClusterAnalysis(){}
	bool Init(string fileName, float radius, int minpts, bool isnormalize);
	bool transPCL2DBpoint(ske::deltPoint& p, float radius, int minPTs, bool isnormalize);
	bool preProcess(bool isnormalize);
	bool DoDBscan(); 
	bool DataAddClusterID(); 
	bool WriteToFile(string filename, bool issorted);
	void Normalize(); 
	void meandeal(); 
	void SortedByCluster();
	vector<int> GetKernelCluster(); 

	bool MyRansac();

};