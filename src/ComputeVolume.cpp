#pragma once
#include "PCLlibrary.h"
#include "getfile.h"
#include <iostream>
#include <string>
#include <map>
#include <omp.h>
#include "geometry.h"
#include "skeleton.h"
#include "tools.h"
#include "process.h"
#include <float.h>
#include <thread>
#include <vtkOutputWindow.h>  
#include <vtkSetGet.h>
#include "connect_ske.h"
#include "ClusterAnalysis.h"
#include <math.h>
#include <string.h>
#include <ctime>
#include <vtkPlaneSource.h> 
#include <Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h>
Geo::ProPlane ptt;
Geo::ProPlane ptt1;
Geo::ProPlane ptt2;
Geo::ProPlane ptt3;
pcl::PointCloud<pcl::PointXYZRGB>::Ptr neipts_proj_test(new pcl::PointCloud<pcl::PointXYZRGB>);
vector<pcl::PointXYZRGB> bod_l;
vector<pcl::PointXYZRGB> bod_w;

Files::ConfigFile parameters;
const int demension = 7;
void loadfromtxt(string txtpath, vector<string> &plist)
{
	//cout << "test" << endl;
	ifstream infile(txtpath);
	if (!infile.is_open())
	{
		cout << "Error opening file" << endl;
		getchar();
		exit(-1);
	}

	int datacouts = 0; //数据个数统计
	while (!infile.eof())
	{

		string s;
		while (getline(infile, s))
		{
			plist.push_back(s);
			//cout << s << endl;
		}
	}
	infile.close();
}

void getdatafromstring(string input,vector<float> &output)
{
	stringstream temp(input);
	//依次输出到result中，并存入res中 
	string allstring;
	vector<string> allset;
	while (getline(temp, allstring, ' '))
		allset.push_back(allstring);
	string pos;
	string color;
	string branchid;
	//pos
	stringstream posstream;
	posstream << allset.at(0);
	while (getline(posstream, pos, ',')) 
	{
		output.push_back(atof(pos.c_str()));
	}
	//color
	stringstream colorstream;
	colorstream << allset.at(1);
	while (getline(colorstream, color, ','))
	{
		output.push_back(atof(color.c_str()));
	}
	//branchesid
	stringstream idstream;
	idstream << allset.at(2);
	while (getline(idstream, branchid, ','))
	{
		output.push_back(atof(branchid.c_str()));
	}
}

// branch id&skeleton id
class PointBranchID
{
public:
	int id_in_cloud;
	int branch_id;
	PointBranchID(int idp1, int idb1) {
		this->id_in_cloud = idp1;
		this->branch_id = idb1;
	}
};

// 
class PointDirect
{
public:
	Eigen::Matrix3d Convariance;
	Eigen::Matrix3d convariance_inverse; 
	vector<double> proj_vector;
	float proj_eigen_value;
	vector<vector<double>> eigen_vectors;
	vector<double> eigen_values;
};

// 
void transTxt2Pointcloud(vector<string> &plist, pcl::PointCloud<pcl::PointXYZRGB>::Ptr& out, vector<vector<PointBranchID>> &branchlist)
{
	vector<float> data0;
	getdatafromstring(plist.at(0), data0);
	PointBranchID initp(0, data0.at(6));
	vector<PointBranchID> same_id_set;
	for (int i = 0; i < plist.size(); i++)
	{
		//0,1,2-xyz，3,4,5-rgb，6-branch_id
		vector<float> temp;
		if(plist.at(i).size()<=1)
			continue;
		getdatafromstring(plist.at(i), temp);
		pcl::PointXYZRGB ptemp;
		ptemp.x = temp.at(0), ptemp.y = temp.at(1), ptemp.z = temp.at(2);
		ptemp.r = temp.at(3)*255, ptemp.g = temp.at(4) * 255, ptemp.b = temp.at(5) * 255;
		out->push_back(ptemp);
		PointBranchID tp(i, temp.at(6));
		if (tp.branch_id == initp.branch_id)
		{
			same_id_set.push_back(tp);
		}
		else
		{
			vector<PointBranchID> same_id_set_temp;
			same_id_set_temp = same_id_set;
			branchlist.push_back(same_id_set_temp);
			same_id_set.clear();
			initp = tp;
			same_id_set.push_back(initp);
		}
	}
	if (!same_id_set.empty())
		branchlist.push_back(same_id_set);


}

// classified point cloud
void getClassedSkePointCloud(string txtpath, vector<string> &plist, pcl::PointCloud<pcl::PointXYZRGB>::Ptr& out, vector<vector<PointBranchID>> &branchlist)
{
	loadfromtxt(txtpath, plist);
	transTxt2Pointcloud(plist, out, branchlist);
}

// random color
void creatColor(float &outr, float &ob, float &og)
{
	//Sleep(300);
	//srand(time(NULL));
	float r = (rand() % (100 - 0 + 1))*256.0 / 100 +32;
	float g = (rand() % (100 - 0 + 1))*256.0 / 100 -15;
	float b = (rand() % (100 - 0 + 1))*256.0 / 100 +20;
	outr = r, ob = b, og = g;
}

// 
void ComputeGeneralConvariance(pcl::PointCloud<pcl::PointXYZRGB>::Ptr& incloud, PointDirect& tp)
{
	float ave_x = 0, ave_y = 0, ave_z = 0;
#pragma omp parallel  for
	for (int i = 0; i < incloud->size(); i++)
	{
		ave_x += incloud->at(i).x;
		ave_y += incloud->at(i).y;
		ave_z += incloud->at(i).z;
	}
	ave_x = ave_x / incloud->size();
	ave_y = ave_y / incloud->size();
	ave_z = ave_z / incloud->size();

	double xx11 = 0, xy12 = 0, xz13 = 0;
	double yy22 = 0, yz23 = 0, zz33 = 0;
	int num = incloud->size();
	for (int i = 0; i < incloud->size(); i++)
	{
		double deltax = ave_x - incloud->at(i).x;
		double deltay = ave_y - incloud->at(i).y;
		double deltaz = ave_z - incloud->at(i).z;
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
	tp.Convariance = tempMatrix;
	tp.convariance_inverse = tp.Convariance.inverse();
	return ;
}

// rerange
void getDirectionVector(pcl::PointCloud<pcl::PointXYZRGB>::Ptr& incloud, PointDirect& out)
{
	//PointDirect pdire;
	ComputeGeneralConvariance(incloud, out);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(out.Convariance);
	auto eigenvalue = es.eigenvalues();
	auto vec = es.eigenvectors();
	//cout << eigenvalue << endl;
	//cout << vec << endl;

	MatrixXd evalsReal;
	evalsReal = eigenvalue.real();
	MatrixXf::Index evalsMax;
	evalsReal.rowwise().sum().maxCoeff(&evalsMax);
	
	out.proj_eigen_value = eigenvalue(0);
	out.proj_vector = { vec.real()(0, evalsMax),vec.real()(1, evalsMax), vec.real()(2, evalsMax) };
	vector<double> vmax = { vec.real()(0, evalsMax),vec.real()(1, evalsMax), vec.real()(2, evalsMax) };
	vector<double> vmin = { vec.real()(0, 0),vec.real()(1, 0), vec.real()(2, 0) };
	vector<double>vmid = { vec.real()(0, 1),vec.real()(1, 1), vec.real()(2, 1) };
	out.eigen_vectors = {vmax,vmid,vmin};
	out.eigen_values = { eigenvalue.real()(2),eigenvalue.real()(1),eigenvalue.real()(0) };
//	cout << vec.col(evalsMax) << endl;

	return;
}

template <typename PointT> 
float computeDis(PointT p1, PointT p2)
{
	float x = p1.x - p2.x;
	float y = p1.y - p2.y;
	float z = p1.z - p2.z;
	float d2 = x * x + y * y + z * z;
	float out = sqrt(d2);
	return out;
}

float getBranchesDistance(pcl::PointCloud<pcl::PointXYZRGB>::Ptr& incloud)
{
	float dis = 0;
	for (int i = 0; i < incloud->size()-1; i++)
	{
		pcl::PointXYZRGB p1 = incloud->at(i);
		pcl::PointXYZRGB p2 = incloud->at(i+1);
		float dist = computeDis(p1, p2);
		dis += dist;
	}
	return dis;
}

float computeAngle(vector<double> v1, vector<double> v2)
{
	float v1v2 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	float v1mod = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
	float v2mod = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
	float cosAngle = v1v2 * 1.0 / (v1mod*v2mod);
	float outAngle = acos(cosAngle) * 180 / M_PI;
	return outAngle;
}

template <typename PointT>
float computeDisProj(PointT p1, Geo::ProPlane plane1)
{
	float up = p1.x*plane1.A + p1.y*plane1.B + p1.z*plane1.C + plane1.D;
	float bottom = sqrt(plane1.A*plane1.A + plane1.B*plane1.B + plane1.C*plane1.C);
	float out = abs(up)*1.0 / bottom;
	return out;
}

template <typename PointT1, typename PointT2>
float computeDisProj(PointT1 p1, PointT2 p2, Geo::ProPlane inplane)
{
	double x = p1.x - p2.x;
	double y = p1.y - p2.y;
	double z = p1.z - p2.z;
	double dis = sqrt(x*x + y * y + z * z);
	vector<double> vec_temp = { x,y,z };
	vector<double> plane_n = { inplane.A,inplane.B,inplane.C };
	float anlge = computeAngle(vec_temp, plane_n);
	if (anlge > 90)
		dis = -dis;
	return dis;
}

template <typename PointT>
PointT computePointProjPos(PointT p1, Geo::ProPlane pa)
{
	float xp, yp, zp;
	if (pa.C != 0)
	{		xp = 0, yp = 0, zp = -pa.D*1.0 / pa.C;	}
	else
	{
		if(pa.B!=0)
		{	xp = 0, yp = -pa.D*1.0 / pa.B, zp=0; }
		else
		{   xp = -pa.D / pa.A, yp = 0, zp = 0; }
	}

	//float up = p1.x*plane1.A + p1.y*plane1.B + p1.z*plane1.C + plane1.D;
	float up = pa.A*(pa.A*xp + pa.B*yp + pa.C*zp) + p1.x*pa.B2C2 - pa.A*(pa.B*p1.y + pa.C*p1.z);
	float bottom = pa.A2B2C2;
	float xout = up * 1.0 / bottom;
	float yout = pa.B / pa.A*(xout - p1.x) + p1.y;
	float zout = pa.C / pa.A*(xout - p1.x) + p1.z;
	PointT pout;
	pout.x = xout, pout.y = yout, pout.z = zout;
	return pout;
}

// length width
float computeProjLenWidth(pcl::PointCloud<pcl::PointXYZRGB>::Ptr& incloud, Geo::ProPlane &pl,
	vector<pcl::PointXYZRGB> &board)
{
	pcl::PointXYZRGB p;
	p.x = -pl.D*1.0 / pl.A, p.y = 0, p.z = 0;
	vector<float> disset;
	for (int i=0;i<incloud->size();i++)
	{
		float dis = computeDisProj(p, incloud->at(i), pl);
		disset.push_back(dis);
	}
	pcl::PointXYZRGB pmin,pmax;
	float dismin = DBL_MAX, dismax = DBL_MIN;
	for (int i = 0; i < incloud->size(); i++)
	{
		if (disset.at(i) >= dismax)
		{
			pmax = incloud->at(i);
			dismax = disset.at(i);
		}
		if (disset.at(i) <= dismin)
		{
			pmin = incloud->at(i);
			dismin = disset.at(i);
		}
			
	}
	board.push_back(pmax);
	board.push_back(pmin);
	vector<float> disset_src;
	disset_src = disset;
	sort(disset.begin(), disset.end());
	float out = abs(disset.at(0) - disset.at(disset.size() - 1));
	return out;

}

void computeVolume(pcl::PointCloud<pcl::PointXYZRGB>::Ptr& incloud,
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr& skecloud, PointDirect &dire,
	vector<float> &volume_set,bool isradiuslimit=false)
{
	float len_one = 0, width_one = 0;
	if (skecloud->size() <= 2)
		return;
	for (int i = 0; i < skecloud->size(); i++)
	{
		vector<double> d1 = dire.proj_vector;
		vector<double> d2;
		vector<double> dire;
		pcl::PointXYZRGB pt = skecloud->at(i);
		float searchdis = 0;
		float searchdis2 = 0;

		if (i == 0)
		{
			pcl::PointXYZRGB p1 = skecloud->at(0);
			pcl::PointXYZRGB p2 = skecloud->at(1);
			d2.push_back(p1.x - p2.x), d2.push_back(p1.y - p2.y), d2.push_back(p1.z - p2.z);
			searchdis = computeDis(p1, p2)*1.0 / 2;
		}
		else if (i == skecloud->size() - 1)
		{
			pcl::PointXYZRGB p1 = skecloud->at(i - 1);
			pcl::PointXYZRGB p2 = skecloud->at(i);
			d2.push_back(p1.x - p2.x), d2.push_back(p1.y - p2.y), d2.push_back(p1.z - p2.z);
			searchdis = computeDis(p1, p2)*1.0 / 2;
		}
		else
		{
			pcl::PointXYZRGB p1 = skecloud->at(i);
			pcl::PointXYZRGB p2 = skecloud->at(i - 1);
			pcl::PointXYZRGB p3 = skecloud->at(i + 1);
			d2.push_back(p2.x - p3.x), d2.push_back(p2.y - p3.y), d2.push_back(p2.z - p3.z);
			float dist1 = computeDis(p1, p2), dist2 = computeDis(p1, p3);
			searchdis = (dist1 + dist2)*1.0 / 4;
		}
		float vangle = computeAngle(d1, d2);
		if (vangle >= 15)
			dire = d2;
		else
			dire = d1;
		Geo::ProPlane plane_temp(dire[0], dire[1], dire[2], pt.x, pt.y, pt.z);
		ptt = plane_temp;
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr neipts(new pcl::PointCloud<pcl::PointXYZRGB>);
		vector<float> disset;
		for (int j = 0; j < incloud->size(); j++)
		{
			float dis = computeDisProj(incloud->at(j), plane_temp);
			disset.push_back(dis);
			if (dis <= searchdis)
			{
				neipts->push_back(incloud->at(j));
				//neipts_proj_test->push_back(incloud->at(j));
			}
		}
		float len_t, width_t;
		if (neipts->size() == 0)
		{
			len_t = len_one;
			width_t = width_one;
		}
		else
		{
			plane_temp.GetSumSquare();
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr neipts_projected(new pcl::PointCloud<pcl::PointXYZRGB>);
			vector<float> ppdist;
			for (int j = 0; j < neipts->size(); j++)
			{
				pcl::PointXYZRGB p_temp;
				p_temp = computePointProjPos(neipts->at(j), plane_temp);
				float dist = computeDis(p_temp, skecloud->at(i));
				ppdist.push_back(dist);
			}
			for (int j = 0; j < neipts->size(); j++)
			{
				pcl::PointXYZRGB p_temp;
				p_temp = computePointProjPos(neipts->at(j), plane_temp);		
				if (isradiuslimit == true)
				{
					if(ppdist.at(j)<=3)
						neipts_projected->push_back(p_temp);
				}
				else
					neipts_projected->push_back(p_temp);
				/*neipts_proj_test->push_back(p_temp);*/
			}

			/*	pcl::io::savePLYFile("C:\\Users\\zhihong\\Desktop\\programtest\\2020.12.17\\skeClassSort\\neipts2.ply", *neipts);
				pcl::io::savePLYFile("C:\\Users\\zhihong\\Desktop\\programtest\\2020.12.17\\skeClassSort\\neipts_projected2.ply", *neipts_projected);*/
			PointDirect dire_projeced;
			getDirectionVector(neipts_projected, dire_projeced);
			//ComputeGeneralConvariance(neipts_projected, dire_projeced);

			vector<double> pl1 = dire_projeced.eigen_vectors.at(0);//长
			vector<double> pl2 = dire_projeced.eigen_vectors.at(1);//宽
			vector<double> pl3 = dire_projeced.eigen_vectors.at(2);//所在平面
			Geo::ProPlane plane_len(pl1[0], pl1[1], pl1[2], pt.x, pt.y, pt.z);
			Geo::ProPlane plane_width(pl2[0], pl2[1], pl2[2], pt.x, pt.y, pt.z);
			//Geo::ProPlane plane_len3(pl3[0], pl3[1], pl3[2], pt.x, pt.y, pt.z);
			ptt1 = plane_len; ptt2 = plane_width; //ptt3 = plane_len3;
			plane_len.GetSumSquare(); plane_width.GetSumSquare();
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr neipts_len(new pcl::PointCloud<pcl::PointXYZRGB>);
			for (int j = 0; j < neipts_projected->size(); j++)
			{
				pcl::PointXYZRGB p_temp;
				p_temp = computePointProjPos(neipts_projected->at(j), plane_len);
				neipts_len->push_back(p_temp);
			}
			vector<pcl::PointXYZRGB> board_len;
			len_t = computeProjLenWidth(neipts_len, plane_len, board_len);

			pcl::PointCloud<pcl::PointXYZRGB>::Ptr neipts_width(new pcl::PointCloud<pcl::PointXYZRGB>);
			for (int j = 0; j < neipts_projected->size(); j++)
			{
				pcl::PointXYZRGB p_temp;
				p_temp = computePointProjPos(neipts_projected->at(j), plane_width);
				neipts_width->push_back(p_temp);
			}
			vector<pcl::PointXYZRGB> board_width;
			width_t = computeProjLenWidth(neipts_width, plane_width, board_width);
		}
		// correct length,width
		if (isnan(len_t))
			len_t = len_one;
		if (isnan(width_t))
			width_t = width_one;		
		if (len_t > 30 && width_t > 30)
		{
			len_t = 0, width_t = 0;
		}
		else if (len_t > 30 && width_t < 30)
			len_t = 3 * width_t;
		else if (len_t <= 30 && width_t >= 30)
			width_t = len_t / 1.5;
		else
			len_t = len_t;
		float volume_temp = len_t * width_t*searchdis / 1000;
		volume_set.push_back(volume_temp);
		if (isnan(volume_temp))
			cout << "Nan error" << endl;

	
	}
}

float getBranchesVolume(pcl::PointCloud<pcl::PointXYZRGB>::Ptr& incloud,
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr& skecloud, PointDirect &dire)
{
	vector<float> volume_set;
	computeVolume(incloud, skecloud, dire, volume_set, false);
	float outvolume = 0;
	for (int i=0;i< volume_set.size();i++)
	{
		outvolume += volume_set.at(i);
	}
	if (outvolume > 100)
	{
		volume_set.clear();
		computeVolume(incloud, skecloud, dire, volume_set, true);
		outvolume = 0;
		for (int i = 0; i < volume_set.size(); i++)
		{
			outvolume += volume_set.at(i);
		}
	}
	return outvolume;
}

vtkSmartPointer<vtkPolyData> createPlane(const pcl::ModelCoefficients& coefficients, float scale[2] = nullptr)
{
	vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();

	plane->SetNormal(coefficients.values[0], coefficients.values[1], coefficients.values[2]);
	double norm_sqr = coefficients.values[0] * coefficients.values[0]
		+ coefficients.values[1] * coefficients.values[1]
		+ coefficients.values[2] * coefficients.values[2];


	plane->Push(-coefficients.values[3] / sqrt(norm_sqr));
	plane->SetResolution(200, 200);
	plane->Update();

	double pt1[3], pt2[3], orig[3], center[3];
	plane->GetPoint1(pt1);
	plane->GetPoint2(pt2);
	plane->GetOrigin(orig);
	plane->GetCenter(center);

	double _pt1[3], _pt2[3];
	float scale1 = 3.0;
	float scale2 = 3.0;
	if (scale != nullptr)
	{
		scale1 = scale[0];
		scale2 = scale[1];
	}
	for (int i = 0; i < 3; i++) {
		_pt1[i] = scale1 * (pt1[i] - orig[i]);
		_pt2[i] = scale2 * (pt2[i] - orig[i]);
	}
	for (int i = 0; i < 3; ++i)
	{
		pt1[i] = orig[i] + _pt1[i];
		pt2[i] = orig[i] + _pt2[i];
	}
	plane->SetPoint1(pt1);
	plane->SetPoint2(pt2);

	plane->Update();
	return (plane->GetOutput());
}

void Simplify(pcl::PointCloud<pcl::PointXYZRGB>::Ptr& in,
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr& out)
{
	pcl::KdTreeFLANN<pcl::PointXYZRGB>::Ptr tree(new pcl::KdTreeFLANN<pcl::PointXYZRGB>);
	tree->setInputCloud(in);
	for (int i = 0; i < in->size(); i++)
	{
		if(in->at(i).r==255 && in->at(i).g == 255 && in->at(i).b == 255)
			continue;
		vector<int> kdid;
		vector<float> kddis;
		if (tree->radiusSearch(in->at(i), 10e-5, kdid, kddis) > 0)
		{
			for (int k = 1; k < kdid.size(); k++)
			{
				in->at(kdid.at(k)).r = 255;
				in->at(kdid.at(k)).g = 255;
				in->at(kdid.at(k)).b = 255;
			}
		}
	}
	for (int i = 0; i < in->size(); i++)
	{
		if (in->at(i).r == 255 && in->at(i).g == 255 && in->at(i).b == 255)
			continue;
		out->push_back(in->at(i));
	}
}


float AverageDis(vector<float> &in)
{
	float out=0;
	for (int i = 0; i < in.size(); i++)
		out += in[i];
	out = out / in.size();
	return out;
}

float StdDis(vector<float> &in)
{
	float avedis = AverageDis(in);
	float sum = 0;
	for (int i = 0; i < in.size(); i++)
	{
		sum += pow((in[i] - avedis), 2);
	}
	sum = sqrt(sum / in.size());
	return sum;
}

int main_0()
{
	
	/*cout << "Ready.....";
	getchar();*/
	cout << "Adjust coordinate" << endl;
	map<string, string> param;
	// Files file;
	Files::ConfigFileRead(param, "./计算体积.ini");
	istringstream incloudstr(param["singleinpath"]);
	istringstream outcloudstr(param["singleoutpath"]);
	istringstream intxtpathstr(param["singleintxtpath"]);
	string incloudpath, intxtpath, oucloudpath;
	incloudstr >> incloudpath;
	intxtpathstr >> intxtpath;
	outcloudstr >> oucloudpath;

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src(new pcl::PointCloud<pcl::PointXYZ>);
	int error = pcl::io::loadPLYFile(incloudpath, *cloud_src);
	if (error == 1)
		return -1;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_all(new pcl::PointCloud<pcl::PointXYZRGB>);
	pcl::copyPointCloud(*cloud_src, *cloud_all);
	for (int i = 0; i < cloud_src->size(); i++)
	{
		cloud_all->at(i).r = 255;
		cloud_all->at(i).g = 0;
		cloud_all->at(i).b = 0;
	}

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ske(new pcl::PointCloud<pcl::PointXYZRGB>);
	vector<string> idlist;
	vector<vector<PointBranchID>> branchlist;
	getClassedSkePointCloud(intxtpath, idlist, cloud_ske, branchlist);
	clock_t time_1 = time_stamp();

	pcl::KdTreeFLANN<pcl::PointXYZRGB>::Ptr tree(new pcl::KdTreeFLANN<pcl::PointXYZRGB>);
	tree->setInputCloud(cloud_all);
	vector< pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloudset;
	vector< pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_ske_set;
#pragma omp parallel for 
	for (int i =0;i<branchlist.size();i++)
	{
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_temp(new pcl::PointCloud<pcl::PointXYZRGB>);
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ske_temp(new pcl::PointCloud<pcl::PointXYZRGB>);
		vector<PointBranchID> temp = branchlist.at(i);
		for (int j = 0; j < temp.size(); j++)
		{
			vector<int> kdid;
			vector<float> kddis;

			pcl::PointXYZRGB tp = cloud_ske->at(temp.at(j).id_in_cloud);
			if (tree->radiusSearch(tp, 10, kdid, kddis) > 0)
			{

				for (int m = 1; m < kdid.size(); m++)
				{
					cloud_temp->push_back(cloud_all->at(kdid.at(m)));
				}
			}
			cloud_ske_temp->push_back(cloud_ske->at(temp.at(j).id_in_cloud));
		}
		cloudset.push_back(cloud_temp);
		cloud_ske_set.push_back(cloud_ske_temp);
	}


	for (int i = 0; i < cloudset.size(); i++)
	{
		float r, g, b;
		creatColor(r, g, b);
#pragma omp parallel for 
		for (int j = 0; j < cloudset.at(i)->size(); j++)
		{
			cloudset.at(i)->at(j).r = r;
			cloudset.at(i)->at(j).g = g;
			cloudset.at(i)->at(j).b = b;
		}
	}


	vector<PointDirect> dire_set;
	for (int i=0;i<cloud_ske_set.size();i++)
	{
		PointDirect out;
		getDirectionVector(cloud_ske_set.at(i), out);
		dire_set.push_back(out);
	}
	

	vector<float> dist_set;
	float dis_ave = 0;
	for (int i = 0; i < cloud_ske_set.size(); i++)
	{
		float dis_temp=getBranchesDistance(cloud_ske_set.at(i));
		dist_set.push_back(dis_temp);
		dis_ave += dis_temp;
	}

	dis_ave = (dis_ave-dist_set.at(cloud_ske_set.size()-1)) / (cloud_ske_set.size()-1);
	vector<float> dist_set_sort;
	for (int i = 0; i < dist_set.size(); i++)
	{
		if (dist_set[i] > 2 * dis_ave)
			continue;
		else
			dist_set_sort.push_back(dist_set[i]);
	}
	//dist_set_sort = dist_set;
	sort(dist_set_sort.begin(), dist_set_sort.end());
	float middis = dist_set_sort.at(int(cloudset.size() / 2));
	cout << "ave: " << dis_ave << ", mid: " << middis << endl;
	float thetadis = StdDis(dist_set_sort);

	vector<float> volume_set;
	float total_volume = 0;
	vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_out_mid;
	vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_out_mid_ske;
	int validnum1 = 0;
	double minthre = 0;
	if (middis - 20 < 0 || middis - 1.96 * thetadis < 0)
		minthre = 0.1*middis;
	else
	{
		minthre = middis - 20;
		minthre = middis - 1.96 * thetadis;
	}	
	vector<float> disout;
	float maxthreh1 = 1.5*middis + 40;
	float maxthreh2 = middis + 1.96*thetadis;
	float maxthreh = max(maxthreh1, maxthreh2);
	for (int i = 0; i < cloud_ske_set.size(); i++)
	{
		if (cloudset.at(i)->size() <= 10 && dist_set.at(i) < 15)
			continue;
		if (dist_set.at(i) <= minthre || dist_set.at(i) >= maxthreh)
			continue;
		else
		{
			float volume_temp = getBranchesVolume(cloudset.at(i), cloud_ske_set.at(i), dire_set.at(i));
			//float volume_temp = 1;
			disout.push_back(dist_set[i]);
			volume_set.push_back(volume_temp);
			total_volume += volume_temp;

			validnum1 += 1;
			//break;
		}
		cloud_out_mid.push_back(cloudset.at(i));

	}
	float ave_volume = total_volume * 1.0 / volume_set.size();
	//cout << "Total volume: " << total_volume<<endl;
	vector<float> disoutsorted;
	sort(disout.begin(), disout.end());
	disoutsorted = disout;

	vector<float> volset;
	sort(volume_set.begin(), volume_set.end());
	volset = volume_set;


	float Maxdis = disoutsorted.at(disoutsorted.size() - 1);
	
	clock_t time_2 = time_stamp();
	cout << "using time: " << (time_2 - time_1)*1.0 / 1000 << "s" << endl;
	
	vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_out_mid2;
	float tov2 = 0;
	for (int i = 0; i < cloud_out_mid.size(); i++)
	{
		float voltemp = volume_set[i];
		if (voltemp > 1.5 || voltemp <0.01*AverageDis(volume_set))
			continue;
		else
		{
			cloud_out_mid2.push_back(cloud_out_mid[i]);
			tov2 += voltemp;
		}
	}
	cout << "Max dis: " << disoutsorted.at(disoutsorted.size() - 1) << endl;
	cout << "Mid dis: " << disoutsorted.at(disoutsorted.size() / 2) << endl;
	cout << "Ave dis: " << AverageDis(disoutsorted) << endl;
	cout << "Ave volume: " << AverageDis(volset) << endl;
	cout << "Total volume: " << tov2 << endl;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_out(new pcl::PointCloud<pcl::PointXYZRGB>);
	//// 显示
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("viewer1"));
	viewer->setBackgroundColor(1, 1, 1);
	//viewer->setBackgroundColor(0, 0, 0);
	int validnum = 0;

	for (int i = 0; i < cloud_out_mid2.size(); i++)
	{
		/*if (volume_set.at(i)<=100)
			continue;*/
		string cname = "Cloud" + to_string(i);
		string cname_ske = "Cloud_ske" + to_string(i);

		for (int j = 0; j < cloud_out_mid2.at(i)->size(); j++)
		{
			cloud_out->push_back(cloud_out_mid2.at(i)->at(j));
		}
		
		validnum++;
	}

	cout << "validnum = " << validnum << endl;

	viewer->addPointCloud(cloud_out, "out");
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);

		//viewer3->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}

	getchar();
	return 0;
}

int main_1()
{
	//cout << "Ready.....";
	//getchar();
	cout << "Adjust coordinate" << endl;
	map<string, string> param;
	vector<string> files;
	vector<string> filenames;
	// Files file;
	Files::ConfigFileRead(param, "./计算体积.ini");
	istringstream incloudstr(param["inpath"]);
	istringstream outcloudstr(param["outpath"]);
	istringstream intxtpathstr(param["intxtpath"]);
	string incloudpath, intxtpath, outpath;
	incloudstr >> incloudpath;
	intxtpathstr >> intxtpath;
	outcloudstr >> outpath;
	Files::getFiles(incloudpath, files, filenames, ".ply");
	vector<string> txtfiles;
	vector<string> txtfilenames;
	Files::getFiles(intxtpath, txtfiles, txtfilenames, ".txt");
	vector<string> volumvec;
	for (int i = 0; i <files.size(); i++)// 1; i++)//
	{
		string prinout;
		string incpth = files.at(i);
		string intxtpth = txtfiles.at(i);
		
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src(new pcl::PointCloud<pcl::PointXYZ>);
		int error = pcl::io::loadPLYFile(incpth, *cloud_src);
		if (error == 1)
			return -1;
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_all(new pcl::PointCloud<pcl::PointXYZRGB>);
		pcl::copyPointCloud(*cloud_src, *cloud_all);
		for (int i = 0; i < cloud_src->size(); i++)
		{
			cloud_all->at(i).r = 255;
			cloud_all->at(i).g = 0;
			cloud_all->at(i).b = 0;
		}

		// 
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ske(new pcl::PointCloud<pcl::PointXYZRGB>);
		vector<string> idlist;
		vector<vector<PointBranchID>> branchlist;
		getClassedSkePointCloud(intxtpth, idlist, cloud_ske, branchlist);
		clock_t time_1 = time_stamp();
		
		pcl::KdTreeFLANN<pcl::PointXYZRGB>::Ptr tree(new pcl::KdTreeFLANN<pcl::PointXYZRGB>);
		tree->setInputCloud(cloud_all);
		vector< pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloudset;
		vector< pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_ske_set;
#pragma omp parallel for 
		for (int i = 0; i < branchlist.size(); i++)
		{
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_temp(new pcl::PointCloud<pcl::PointXYZRGB>);
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ske_temp(new pcl::PointCloud<pcl::PointXYZRGB>);
			vector<PointBranchID> temp = branchlist.at(i);
			for (int j = 0; j < temp.size(); j++)
			{
				vector<int> kdid;
				vector<float> kddis;

				pcl::PointXYZRGB tp = cloud_ske->at(temp.at(j).id_in_cloud);
				if (tree->radiusSearch(tp, 10, kdid, kddis) > 0)
				{

					for (int m = 1; m < kdid.size(); m++)
					{
						cloud_temp->push_back(cloud_all->at(kdid.at(m)));
					}
				}
				cloud_ske_temp->push_back(cloud_ske->at(temp.at(j).id_in_cloud));
			}
			cloudset.push_back(cloud_temp);
			cloud_ske_set.push_back(cloud_ske_temp);
		}


		for (int i = 0; i < cloudset.size(); i++)
		{
			float r, g, b;
			creatColor(r, g, b);
#pragma omp parallel for 
			for (int j = 0; j < cloudset.at(i)->size(); j++)
			{
				cloudset.at(i)->at(j).r = r;
				cloudset.at(i)->at(j).g = g;
				cloudset.at(i)->at(j).b = b;
			}
		}


		vector<PointDirect> dire_set;
		for (int i = 0; i < cloud_ske_set.size(); i++)
		{
			PointDirect out;
			getDirectionVector(cloud_ske_set.at(i), out);
			dire_set.push_back(out);
		}
		
		vector<float> dist_set;
		float dis_ave = 0;
		for (int i = 0; i < cloud_ske_set.size(); i++)
		{
			float dis_temp = getBranchesDistance(cloud_ske_set.at(i));
			dist_set.push_back(dis_temp);
			dis_ave += dis_temp;
		}

		dis_ave = (dis_ave - dist_set.at(cloud_ske_set.size() - 1)) / (cloud_ske_set.size() - 1);
		vector<float> dist_set_sort;
		for (int i = 0; i < dist_set.size(); i++)
		{
			if (dist_set[i] > 2 * dis_ave)
				continue;
			else
				dist_set_sort.push_back(dist_set[i]);
		}
		//dist_set_sort = dist_set;
		sort(dist_set_sort.begin(), dist_set_sort.end());


		int beginpos = dist_set.size()*0.15;
		int endpos = dist_set.size()*0.85-1;
		dis_ave = 0;
		vector<float> disok;
		int disoknum = 0;
		for (int i = beginpos; i < endpos; i++)
		{
			disok.push_back(dist_set_sort[i]);
			dis_ave += dist_set_sort[i];
			disoknum++;
		}
		dis_ave = dis_ave / disoknum;
		//float middis = dist_set_sort.at(int(cloudset.size() / 2));
		float middis = disok.at(int(disok.size() / 2));
		cout << "ave: " << dis_ave << ", mid: " << middis << endl;
		prinout = filenames[i] + " ave dis: " + to_string(dis_ave) + ", mid: " +to_string(middis);

	
		vector<float> volume_set;
		vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_out_mid;
		vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_out_mid_ske;
		float thetadis = StdDis(dist_set_sort);
		double minthre = 0;
		if (middis - 20 < 0 || middis - 1.96 * thetadis < 0)
			minthre = 0.1*middis;
		else
		{
			minthre = middis - 20;
			minthre = middis - 1.96 * thetadis;
		}
		vector<float> disout;
		float maxthreh1 = 1.5*middis + 40;
		float maxthreh2 = middis + 1.96*thetadis;
		float maxthreh = max(maxthreh1, maxthreh2);
		maxthreh = maxthreh2;

		for (int i = 0; i < cloud_ske_set.size(); i++)
		{
			if (cloudset.at(i)->size() <= 10 && dist_set.at(i) < 15)
				continue;
			if (dist_set.at(i) <= minthre || dist_set.at(i) <= 15 || dist_set.at(i) >= maxthreh)
				continue;
			else
			{
				float volume_temp = getBranchesVolume(cloudset.at(i), cloud_ske_set.at(i), dire_set.at(i));
				//float volume_temp = 1;
				disout.push_back(dist_set[i]);
				volume_set.push_back(volume_temp);
			}
			cloud_out_mid.push_back(cloudset.at(i));
		}


		clock_t time_2 = time_stamp();
		
		cout << "using time: " << (time_2 - time_1)*1.0 / 1000 << "s" << endl;

	
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_out(new pcl::PointCloud<pcl::PointXYZRGB>);

		vector<float> disoutsorted;
		sort(disout.begin(), disout.end());
		disoutsorted = disout;

		vector<float> volset;
		sort(volume_set.begin(), volume_set.end());
		volset = volume_set;
		
		string volume_out = outpath + "\\" + filenames[i] + "_volume.txt";
		ofstream of1(volume_out);
		vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_out_mid2;
		float tov2 = 0;
		for (int i = 0; i < cloud_out_mid.size(); i++)
		{
			float voltemp = volume_set[i];
			of1 << voltemp << " ";
			if (voltemp > 1.0 || voltemp < 0.01*AverageDis(volume_set))
			{
				of1 << endl;
				continue;
			}
			else
			{
				cloud_out_mid2.push_back(cloud_out_mid[i]);
				tov2 += voltemp;
				of1 << voltemp;
			}
			of1 << endl;
		}
		of1.close();

		int validnum = 0;

		for (int i = 0; i < cloud_out_mid2.size(); i++)
		{
			/*if (volume_set.at(i)<=100)
				continue;*/
			string cname = "Cloud" + to_string(i);
			string cname_ske = "Cloud_ske" + to_string(i);

			for (int j = 0; j < cloud_out_mid2.at(i)->size(); j++)
			{
				cloud_out->push_back(cloud_out_mid2.at(i)->at(j));
			}

			validnum++;
		}
		cout << "validnum = " << validnum << endl;
		prinout = prinout + " Volume: " + to_string(tov2);
		prinout += ", Max dis: " + to_string(disout[disout.size() - 1]);
		prinout += ", Siliques: " + to_string(validnum);
		prinout = prinout + ", Using time: " + to_string((time_2 - time_1)*1.0 / 1000) + "s";
		string outcloudpath;
		outcloudpath = outpath + "\\" + filenames[i] + "_classed.ply";
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_out_sim(new pcl::PointCloud<pcl::PointXYZRGB>);
		//cloud_out_sim = cloud_out;
		Simplify(cloud_out, cloud_out_sim);
		pcl::io::savePLYFile(outcloudpath, *cloud_out_sim);
		volumvec.push_back(prinout);


	}
	
	for (int i=0;i<volumvec.size();i++)
	{
		cout << volumvec[i] << endl;
	}

	 getchar();
	 return 0;
}

int main_2()
{
	//getchar();
	cout << "Adjust coordinate" << endl;
	map<string, string> param;
	vector<string> files;
	vector<string> filenames;
	// Files file;
	Files::ConfigFileRead(param, "./计算体积.ini");
	istringstream incloudstr(param["inpath"]);
	istringstream outcloudstr(param["outpath"]);
	istringstream intxtpathstr(param["intxtpath"]);
	string incloudpath, intxtpath, outpath;
	incloudstr >> incloudpath;
	intxtpathstr >> intxtpath;
	outcloudstr >> outpath;
	Files::getFiles(incloudpath, files, filenames, ".ply");
	vector<string> txtfiles;
	vector<string> txtfilenames;
	Files::getFiles(intxtpath, txtfiles, txtfilenames, ".txt");
	vector<string> volumvec;
	for (int i = 0; i < files.size(); i++)// 1; i++)//
	{
		string prinout;
		string incpth = files.at(i);
		string intxtpth = txtfiles.at(i);
	
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src(new pcl::PointCloud<pcl::PointXYZ>);
		int error = pcl::io::loadPLYFile(incpth, *cloud_src);
		if (error == 1)
			return -1;
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_all(new pcl::PointCloud<pcl::PointXYZRGB>);
		pcl::copyPointCloud(*cloud_src, *cloud_all);
		for (int i = 0; i < cloud_src->size(); i++)
		{
			cloud_all->at(i).r = 255;
			cloud_all->at(i).g = 0;
			cloud_all->at(i).b = 0;
		}

		// 
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ske(new pcl::PointCloud<pcl::PointXYZRGB>);
		vector<string> idlist;
		vector<vector<PointBranchID>> branchlist;
		getClassedSkePointCloud(intxtpth, idlist, cloud_ske, branchlist);
		clock_t time_1 = time_stamp();

		pcl::KdTreeFLANN<pcl::PointXYZRGB>::Ptr tree(new pcl::KdTreeFLANN<pcl::PointXYZRGB>);
		tree->setInputCloud(cloud_all);
		vector< pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloudset;
		vector< pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_ske_set;
#pragma omp parallel for 
		for (int i = 0; i < branchlist.size(); i++)
		{
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_temp(new pcl::PointCloud<pcl::PointXYZRGB>);
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ske_temp(new pcl::PointCloud<pcl::PointXYZRGB>);
			vector<PointBranchID> temp = branchlist.at(i);
			for (int j = 0; j < temp.size(); j++)
			{
				vector<int> kdid;
				vector<float> kddis;

				pcl::PointXYZRGB tp = cloud_ske->at(temp.at(j).id_in_cloud);
				if (tree->radiusSearch(tp, 10, kdid, kddis) > 0)
				{

					for (int m = 1; m < kdid.size(); m++)
					{
						cloud_temp->push_back(cloud_all->at(kdid.at(m)));
					}
				}
				cloud_ske_temp->push_back(cloud_ske->at(temp.at(j).id_in_cloud));
			}
			cloudset.push_back(cloud_temp);
			cloud_ske_set.push_back(cloud_ske_temp);
		}


		for (int i = 0; i < cloudset.size(); i++)
		{
			float r, g, b;
			creatColor(r, g, b);
#pragma omp parallel for 
			for (int j = 0; j < cloudset.at(i)->size(); j++)
			{
				cloudset.at(i)->at(j).r = r;
				cloudset.at(i)->at(j).g = g;
				cloudset.at(i)->at(j).b = b;
			}
		}


		vector<PointDirect> dire_set;
		for (int i = 0; i < cloud_ske_set.size(); i++)
		{
			PointDirect out;
			getDirectionVector(cloud_ske_set.at(i), out);
			dire_set.push_back(out);
		}

		vector<float> volume_set;
		vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_out_mid;
		vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_out_mid_ske;
	
		cout << "计算体积" << endl;
		for (int i = 0; i < cloud_ske_set.size(); i++)
		{
				float volume_temp = getBranchesVolume(cloudset.at(i), cloud_ske_set.at(i), dire_set.at(i));
				volume_set.push_back(volume_temp);
			cloud_out_mid.push_back(cloudset.at(i));
		}


		clock_t time_2 = time_stamp();

		cout << "using time: " << (time_2 - time_1)*1.0 / 1000 << "s" << endl;

		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_out(new pcl::PointCloud<pcl::PointXYZRGB>);

		vector<float> volset;
		sort(volume_set.begin(), volume_set.end());
		volset = volume_set;

		vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_out_mid2;
		float tov2 = 0;
		for (int i = 0; i < cloud_out_mid.size(); i++)
		{
			float voltemp = volume_set[i];
			if (voltemp > 1.0 || voltemp < 0.01*AverageDis(volume_set))
				continue;
			else
			{
				cloud_out_mid2.push_back(cloud_out_mid[i]);
				tov2 += voltemp;
			}
		}

		int validnum = 0;

		for (int i = 0; i < cloud_out_mid2.size(); i++)
		{
			/*if (volume_set.at(i)<=100)
				continue;*/
			string cname = "Cloud" + to_string(i);
			string cname_ske = "Cloud_ske" + to_string(i);

			for (int j = 0; j < cloud_out_mid2.at(i)->size(); j++)
			{
				cloud_out->push_back(cloud_out_mid2.at(i)->at(j));
			}

			validnum++;
		}
		cout << "validnum = " << validnum << endl;
		prinout = prinout + " Volume: " + to_string(tov2);
		prinout += ", Siliques: " + to_string(validnum);
		prinout = prinout + ", Using time: " + to_string((time_2 - time_1)*1.0 / 1000) + "s";
		string outcloudpath;
		outcloudpath = outpath + "\\" + filenames[i] + "_classed.ply";

		volumvec.push_back(prinout);

	}

	for (int i = 0; i < volumvec.size(); i++)
	{
		cout << volumvec[i] << endl;
	}
	getchar();
	return 0;
}

int main()
{
	getchar();
	cout << "Adjust coordinate" << endl;
	map<string, string> param;
	vector<string> files;
	vector<string> filenames;
	// Files file;
	Files::ConfigFileRead(param, "./computeVolume.ini");
	istringstream modelstr(param["model"]);
	string model;
	modelstr >> model;
	cout << model << endl;
	if (model == "single")
	{
		cout << "Single" << endl;
		main_0();
	}
	else if (model == "multi")
	{
		cout << "Multi" << endl;
		main_1();
	}
	else
	{
		cout << "Multi no stem" << endl;
		main_2();
	}
}