#include "ClusterAnalysis.h"
#include <fstream>
#include <iosfwd>
#include <math.h>


Eigen::Matrix3d ClusterAnalysis::ComputeGeneralConvariance()
{
	float ave_x = 0, ave_y = 0, ave_z = 0;
#pragma omp parallel  for
	for (int i = 0; i < dataset.size(); i++)
	{
		ave_x += dataset.at(i).GetDimension()[0];
		ave_y += dataset.at(i).GetDimension()[1];
		ave_z += dataset.at(i).GetDimension()[2];
	}
	ave_x = ave_x / dataset.size();
	ave_y = ave_y / dataset.size();
	ave_z = ave_z / dataset.size();

	double xx11 = 0, xy12 = 0, xz13 = 0;
	double yy22 = 0, yz23 = 0, zz33 = 0;
	int num = dataset.size();
	for (int i = 0; i < dataset.size(); i++)
	{
		double deltax = ave_x - dataset.at(i).GetDimension()[0];
		double deltay = ave_y - dataset.at(i).GetDimension()[1];
		double deltaz = ave_z - dataset.at(i).GetDimension()[2];
		xx11 += deltax * deltax;
		xy12 += deltax * deltay;
		xz13 += deltax * deltaz;
		yy22 += deltay * deltay;
		yz23 += deltay * deltaz;
		zz33 += deltaz * deltaz;
	}

	xx11 = xx11 * 1.0 / num;
	xy12 = xy12 * 1.0 / num;
	xz13 = xz13 * 1.0 / num;
	yy22 = yy22 * 1.0 / num;
	yz23 = yz23 * 1.0 / num;
	zz33 = zz33 * 1.0 / num;

	double yx21 = xy12, zx31 = xz13, zy32 = yz23;
	Eigen::Matrix3d tempMatrix;
	tempMatrix << xx11, xy12, xz13,
		yx21, yy22, yz23,
		zx31, zy32, zz33;
	this->Convariance = tempMatrix;
	this->convariance_inverse= this->Convariance.inverse();
	return tempMatrix;
}


Eigen::Vector3d VecMutipleMat(Eigen::Vector3d &vec, Eigen::Matrix3d &matric)
{
	float a1 = vec.data()[0] * matric.data()[0] + 
		vec.data()[1] * matric.data()[3] + 
		vec.data()[2] * matric.data()[6];
	float a2 = vec.data()[0] * matric.data()[1] + 
		vec.data()[1] * matric.data()[4] + 
		vec.data()[2] * matric.data()[7];
	float a3 = vec.data()[0] * matric.data()[2] + 
		vec.data()[1] * matric.data()[5] + 
		vec.data()[2] * matric.data()[8];
	Eigen::Vector3d out;
	out << a1, a2, a3;
	return out;
}

float ClusterAnalysis::GetMaDistance(DataPoint &p1, DataPoint &p2)
{
	
	float y = p1.GetDimension()[1] - p2.GetDimension()[1];
	float z = p1.GetDimension()[2] - p2.GetDimension()[2];
	Eigen::Vector3d diff_temp, out,out2;
	diff_temp << x, y, z;
	out2 = VecMutipleMat(diff_temp, convariance_inverse);
	float outnum = out2.dot(diff_temp);
	return sqrt(outnum);
}

float ClusterAnalysis::GetDistance(DataPoint &p1, DataPoint &p2)
{
	float sum = 0;
	sum += pow(p1.GetDimension()[0] - p2.GetDimension()[0], 2);
	sum += pow(p1.GetDimension()[1] - p2.GetDimension()[1], 2)*0.5;
	sum += pow(p1.GetDimension()[2] - p2.GetDimension()[2], 2);
	return pow(sum, 0.5);
}

void ClusterAnalysis::SetArrivalPoints(DataPoint &p)
{
	float tempradius = radius *1.0 / this->NormalizePara;

	for (int i = 0; i < dataNum; i++)
	{
		float distance = GetMaDistance(dataset[i], p);
		if (distance <= tempradius&& i != p.GetDpId())
		{
			p.GetArrivalPoints().push_back(i);
		}
	}

	if (p.GetArrivalPoints().size() >= minpts)
	{
		p.SetKey(true);
		return;
	}
	p.SetKey(false);
}

void ClusterAnalysis::KeyPointCluster(int id, int clusterId)
{
	DataPoint& srcp = dataset[id];
	if (!srcp.isVisited())
		return;
	vector<int>& arrvalPoints = srcp.GetArrivalPoints();
	for (int i = 0; i < arrvalPoints.size(); i++)
	{
		DataPoint &desDp = dataset[arrvalPoints[i]]; 
		if (!desDp.isVisited())
		{
			desDp.SetClusterId(clusterId);      
			desDp.SetVisited(true);               
			if (desDp.IsKey())                    
			{
				KeyPointCluster(desDp.GetDpId(), clusterId);
			}
		}
	}
}

void ClusterAnalysis::Normalize()
{
	vector<float> xvec;
	vector<float> yvec;
	vector<float> zvec;
	for (int i = 0; i < this->dataset.size(); i++)
	{
		xvec.push_back(dataset[i].GetDimension()[0]);
		yvec.push_back(dataset[i].GetDimension()[1]);
		zvec.push_back(dataset[i].GetDimension()[2]);
	}
	sort(xvec.begin(), xvec.end());
	sort(yvec.begin(), yvec.end());
	sort(zvec.begin(), zvec.end());
	float xrange = xvec.at(0) - xvec.at(xvec.size() - 1);
	float yrange = yvec.at(0) - yvec.at(yvec.size() - 1);
	float zrange = zvec.at(0) - zvec.at(zvec.size() - 1);
	this->NormalizePara = (pow(xrange*xrange + yrange*yrange + zrange*zrange, 0.5));
	//this->normalizeDataset = this->dataset;
	for (int i=0;i<this->dataset.size();i++)
	{
		float uppoint[3] = { (dataset[i].GetDimension()[0]-xvec.at(0)) * 1.0 / xrange,
			(dataset[i].GetDimension()[1] -yvec.at(0)) * 1.0 / yrange,
			(dataset[i].GetDimension()[2]-zvec.at(0)) * 1.0 / zrange };
		this->dataset.at(i).SetDimension(uppoint);
	}
}

void ClusterAnalysis::meandeal()
{
	float xmean=0;
	float ymean=0;
	float zmean=0;
	for (int i = 0; i < this->dataset.size(); i++)
	{
		xmean +=dataset[i].GetDimension()[0];
		ymean +=dataset[i].GetDimension()[1];
		zmean +=dataset[i].GetDimension()[2];
	}
	for (int i = 0; i < this->dataset.size(); i++)
	{
		float uppoint[3] = { (dataset[i].GetDimension()[0] - xmean * 1.0 / this->dataset.size()),
			(dataset[i].GetDimension()[1] - ymean * 1.0 / this->dataset.size()),
			(dataset[i].GetDimension()[2] - zmean * 1.0 / this->dataset.size()) };
		this->dataset.at(i).SetDimension(uppoint);
	}
}

bool ClusterAnalysis::transPCL2DBpoint(ske::deltPoint& pts,
	float radius, int minPTs, bool isnormalize)
{
	this->radius = radius;
	this->minpts = minPTs;
	this->dimNum = Dim_num;
	this->isNormalized = isnormalize;
	this->dataNum = pts.neibor_point_src->size() + 1;
	for (int i=0;i< pts.neibor_point_src->size();i++)
	{
		DataPoint tempDP;
		float tempDimData[Dim_num];
		tempDimData[0] = pts.neibor_point_src->at(i).x;
		tempDimData[1] = pts.neibor_point_src->at(i).y;
		tempDimData[2] = pts.neibor_point_src->at(i).z;
		tempDP.SetDimension(tempDimData);
		tempDP.SetDpId(i);
		tempDP.SetPidinSrc(pts.neibor_point_src->at(i).id_in_src);
		tempDP.SetVisited(false);
		tempDP.SetClusterId(-1);
		tempDP.SetKernel(false);
		dataset.push_back(tempDP);
	}
	DataPoint tempDP;
	float tempCenter[Dim_num];
	tempCenter[0] = pts.p.x;
	tempCenter[1] = pts.p.y;
	tempCenter[2] = pts.p.z;
	tempDP.SetDimension(tempCenter);
	tempDP.SetDpId(pts.neibor_point_src->size());
	tempDP.SetPidinSrc(pts.id_in_src);
	tempDP.SetVisited(false);
	tempDP.SetClusterId(-1);
	tempDP.SetKernel(true);
	dataset.push_back(tempDP);
	return true;
}

bool ClusterAnalysis::preProcess(bool isnormalize)
{
	if (isnormalize == true)
	{
		meandeal();
		Normalize();
	}
	else
	{
		meandeal();
		this->NormalizePara = 1;
	}
	ComputeGeneralConvariance();
	for (int i = 0; i < dataNum; i++)
	{
		SetArrivalPoints(dataset[i]);
	}
	return true;
}

bool ClusterAnalysis::Init(string fileName, float radius, int minPTs, bool isnormalize)
{
	this->radius = radius;
	this->minpts = minPTs;
	this->dataNum = Dim_num;
	this->isNormalized = isnormalize;
	ifstream ifs(fileName);
	if (!ifs.is_open())
	{
		cout << "Error opening file" << endl;
		getchar();
		exit(-1);
	}
	
	int datacouts = 0; 
	while (!ifs.eof()) 
	{
		DataPoint tempDP; 
		float tempDimData[Dim_num];
		for (int j = 0; j < Dim_num; j++)
		{
			ifs >> tempDimData[j];
		}
		tempDP.SetDimension(tempDimData);
		tempDP.SetDpId(datacouts);
		tempDP.SetVisited(false);
		tempDP.SetClusterId(-1);
		dataset.push_back(tempDP);
		datacouts++;
	}
	ifs.close();
	dataNum = datacouts;
	if (isnormalize == true)
	{
		meandeal();
		Normalize();
	}
	else
	{
		meandeal();
		this->NormalizePara = 1;	
	}
	ComputeGeneralConvariance();
	for (int i = 0; i < dataNum; i++)
	{
		SetArrivalPoints(dataset[i]);
	}
	return true;

}

bool ClusterAnalysis::DoDBscan()
{
	int clusterId = 0; 

	for (int i = 0; i < dataNum; i++)
	{
		DataPoint &dp = dataset[i];
		if (!dp.isVisited() && dp.IsKey())
		{
			dp.SetClusterId(clusterId);
			dp.SetVisited(true);
			KeyPointCluster(i, clusterId);
			clusterId++;
		}
	}
	
	return true;
}

bool ClusterAnalysis::DataAddClusterID()
{
	for (int i = 0; i < dataNum; i++)
	{
		this->dataset.at(i).DemisionAddClusterId();
	}
	return true;
}

vector<int> ClusterAnalysis::GetKernelCluster()
{
	DataPoint temp = this->dataset.at(dataset.size() - 1);
	vector<int> tempvec;
	try
	{
		if (temp.GetIsKernel() != true)
		{
			throw ("New Kernel point");
		}
		this->kernelClusterId = temp.GetClusterId();
	}
	catch (string v)
	{
		cout << v << endl;
		tempvec.push_back(-100);
		return tempvec;
	}
	for (int i=0;i<this->dataset.size()-1;i++)
	{
		int tempcluster = dataset.at(i).GetClusterId();
		if (tempcluster == this->kernelClusterId)
			tempvec.push_back(i);
		
	}
	return tempvec;
}

bool ClusterAnalysis::WriteToFile(string filename, bool issorted)
{
	int dimsize= Dim_num + 1;
	vector<DataPoint> tempvec;
	if (issorted == true)
	{
		SortedByCluster();
		tempvec = this->sortedDataset;
	}
	else
	{
		tempvec = this->dataset;

	}
	ofstream of1(filename);
	for (int i = 0; i < dataNum; i++)
	{
		for (int d = 0; d < dimsize; d++)
		{
			of1 << tempvec[i].GetDimensionWithClusterId()[d] << '\t';
		}
		of1 << endl;
	}
	of1.close();
	return true;
}

bool isSmaller(DataPoint&p1,DataPoint &p2)
{
	if (p1.GetClusterId() < p2.GetClusterId())
		return true;
	else
		return false;
}

void ClusterAnalysis::SortedByCluster()
{

	this->sortedDataset = this->dataset;
	sort(this->sortedDataset.begin(), this->sortedDataset.end(), isSmaller);
}

