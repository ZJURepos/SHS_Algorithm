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

using namespace std;
vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_vector;
vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_delete;
pcl::PointCloud<pcl::PointNei>::Ptr c_in_ini(new pcl::PointCloud<pcl::PointNei>); 
pcl::PointCloud<pcl::PointNei>::Ptr c_down_ini(new pcl::PointCloud<pcl::PointNei>); 
pcl::PointCloud<pcl::PointXYZRGB>::Ptr c_down_rgb_ini(new pcl::PointCloud<pcl::PointXYZRGB>);
vector<ske::deltPoint> pts_set_ini; 
vector<ske::deltPoint> pts_set_down_ini; 
Files::ConfigFile parameters;

volatile bool extrac_ske_finished = 0;
float ransacdis = 0.001;

pcl::PointCloud<pcl::PointXYZRGB>::Ptr c_down_show;
boost::mutex updateModelMutex;

int *viewerRunner()
{
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("v1"));
	viewer->setBackgroundColor(0.75, 0.75, 0.75);
	while (!viewer->wasStopped())
	{
		viewer->removeAllPointClouds();

		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100));
	}
	return 0;
}

void *viewerCreat(ViewerThread *v1)
{
	boost::shared_ptr<pcl::visualization::PCLVisualizer>
		tv(new pcl::visualization::PCLVisualizer(v1->viewername));
	v1->viewer = boost::shared_ptr<pcl::visualization::PCLVisualizer>(tv);
	v1->viewer->setBackgroundColor(0.85, 0.85, 0.85);
	v1->updatePointCloud();
	while (!v1->viewer->wasStopped())
	{
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
			FOREGROUND_INTENSITY | FOREGROUND_GREEN | FOREGROUND_BLUE);
		cout << v1->cloud_show->size() << endl;
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY |
			FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
		v1->viewer->removeAllPointClouds();
		v1->viewer->addPointCloud(v1->cloud_show, "p");
		v1->viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100));
	}
	return 0;
}

void initialization(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &cloud_in,
	pcl::PointCloud<pcl::PointNei>::Ptr& c_in,
	pcl::PointCloud<pcl::PointNei>::Ptr& c_down,
	vector<ske::deltPoint>& pts_set,
	float &h0, float &rate)
{
	copySrc2PointNei(cloud_in, c_in);
	getdownsample(c_in, c_down, rate);
#pragma omp parallel  for
	for (int i = 0; i < c_in->size(); i++)
	{
		ske::deltPoint temp;
		temp.p = c_in->at(i);
		temp.id_in_down = c_in->at(i).id_in_down;
		temp.id_in_src = c_in->at(i).id_in_src;
		temp.current_h = h0;
		temp.fixed_h = h0;
		pts_set.push_back(temp);
	}
}

void RANSAC_Plane(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloudin,
	pcl::ModelCoefficients::Ptr &coefficients,
	pcl::PointIndices::Ptr &inliers,
	pcl::PointCloud<pcl::PointXYZ>::Ptr &out)
{
	pcl::SACSegmentation<pcl::PointXYZ> seg;
	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_PLANE);
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setDistanceThreshold(0.01);
	seg.setInputCloud(cloudin);
	seg.segment(*inliers, *coefficients);
	pcl::ExtractIndices<pcl::PointXYZ> extract;
	extract.setInputCloud(cloudin);
	extract.setIndices(inliers);
	extract.setNegative(false);
	extract.filter(*out);
}


void OneStep(pcl::PointCloud<pcl::PointNei>::Ptr& c_in,
	pcl::PointCloud<pcl::PointNei>::Ptr& c_down, vector<ske::deltPoint>& pts_set,
	vector<ske::deltPoint> &pts_set_down_out, float &h0, float u, int setTimes,
	ViewerThread *v1)
{
	cout << c_in->size() << endl;
	cout << c_down->size() << endl;	
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr alltemp(new pcl::PointCloud<pcl::PointXYZRGB>);
	copyDel2PCRGB(pts_set, alltemp);


	int iterationTimes = 0;
	float h = h0;
	pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr t_src(new pcl::KdTreeFLANN<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr tc_src(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr t_down(new pcl::KdTreeFLANN<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr tc_down(new pcl::PointCloud<pcl::PointXYZ>);

	float ave_move_dis = DBL_MAX;
	float ave_move_dis_before = DBL_MAX;
	float stdv = DBL_MAX;
	vector<float> tempdis;
	float maxdis = 0;
	vector<ske::deltPoint> pts_set_down;
#pragma omp parallel  for
	for (int i = 0; i < pts_set.size(); i++)
	{
		pts_set.at(i).current_h = h;
	}
	
	clock_t time_ini = time_stamp();
	ske::build_kdtree(c_in, t_src, tc_src);
	ske::build_kdtree(c_down, t_down, tc_down);

	vector<int> countsdow;
	clock_t time_findner = time_stamp();
#pragma omp parallel  for
	for (int i = 0; i < pts_set.size(); i++)
	{
		if (pts_set.at(i).has_fixed == true)
			continue;
		pts_set.at(i).FIndNeiborInSrc(t_src, tc_src, c_in);

		//Ransac
		pcl::PointCloud<pcl::PointXYZ>::Ptr neibor_point_src(new pcl::PointCloud<pcl::PointXYZ>);
		for (int j = 0; j < pts_set.at(i).neibor_point_src->size(); i++)
		{
			pcl::PointXYZ ttp;
			ttp.x = pts_set.at(i).neibor_point_src->at(j).x;
			ttp.y = pts_set.at(i).neibor_point_src->at(j).y;
			ttp.z = pts_set.at(i).neibor_point_src->at(j).z;
			neibor_point_src->push_back(ttp);
		}
		pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
		pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
		RANSAC_Plane(neibor_point_src, coefficients, inliers);
		pcl::PointCloud<pcl::PointNei>::Ptr neibor_point_srcrenew(new pcl::PointCloud<pcl::PointNei>);
		for (int j = 0; i < inliers->indices.size(); j++)
		{
			int id = inliers->indices[j];
			neibor_point_srcrenew->push_back(pts_set.at(i).neibor_point_src->at(id));
		}
		//renew src neibor
		pts_set.at(i).neibor_point_src->clear();
		for (int j = 0; i < neibor_point_srcrenew->size(); j++)
		{
			pts_set.at(i).neibor_point_src->push_back(neibor_point_srcrenew->at(j));
		}

		// down nei
		if (pts_set.at(i).id_in_down != -1)
		{
			pts_set.at(i).FIndNeiborInDown(t_down, tc_down, c_down);

		}

		// renew down nei
		float pA = coefficients->values[0];
		float pB = coefficients->values[1];
		float pC = coefficients->values[2];
		float pD = coefficients->values[3];
		float planeModel = sqrt(pA * pA + pB * pB + pC * pC);
		pcl::PointCloud<pcl::PointNei>::Ptr neibor_point_downrenew(new pcl::PointCloud<pcl::PointNei>);
		for (int j = 0; i < pts_set.at(i).neibor_point_down->size(); j++)
		{
			pcl::PointNei ttp = pts_set.at(i).neibor_point_down->at(j);
			float dis = abs(ttp.x*pA+ttp.y*pB+ttp.z*pC+pD)/planeModel;
			if (dis <= ransacdis)
				neibor_point_downrenew->push_back(ttp);
		}
		pts_set.at(i).neibor_point_down->clear();
		for (int j = 0; i < neibor_point_downrenew->size(); j++)
		{
			pts_set.at(i).neibor_point_down->push_back(neibor_point_downrenew->at(j));
		}

		pts_set.at(i).NeiborSrcDownUpate();
	}
	clock_t time_findnerend = time_stamp();
	cout << countsdow.size();
	cout << "Points find Ransac neibor using time: " << (time_findnerend - time_findner)*1.0 / 1000 << "s" << endl;
	
	copyDel2PCRGB(pts_set, alltemp);
	v1->UpdateOutRGBAPts(alltemp, "all");

	clock_t time_ini_end = time_stamp();
	cout << "Initial Kdtree time used: " << (time_ini_end - time_ini)*1.0 / 1000 << "s" << endl;

	while (iterationTimes < setTimes)
	{
		iterationTimes += 1;
		cout << "**************Iteration Times: " << iterationTimes << "**************" << endl;
		clock_t time_b = time_stamp();

#pragma omp parallel  for
		for (int i = 0; i < pts_set.size(); i++)
		{
			int id_in_src = pts_set.at(i).id_in_src;
			pts_set.at(id_in_src).valid_ornot = 0;
		}

#pragma omp parallel  for
		for (int i = 0; i < pts_set.size(); i++)
		{
			if (pts_set.at(i).neibor_point_down->size() < 2)
			{
				pts_set.at(i).sample_invalid = -1;
				pts_set.at(i).valid_ornot = -1;
			}
		}

		if (iterationTimes == 1)
		{
#pragma omp parallel  for num(8) 
			for (int i = 0; i < pts_set.size(); i++)
			{
				c_in->at(i).dj = pts_set.at(i).local_density_src(h0); //dj
				if (isnan(c_in->at(i).dj))
				{
					cout << "current " << i << " dj is: Nan" << endl;
					getchar();
				}
				if (pts_set.at(i).valid_ornot != -1)
				{
					pts_set.at(i).compute_distribution_value();  //delta
				}
				else
					pts_set.at(i).delta = -1;
			}
		}
		else
		{
#pragma omp parallel  for num(8)
			for (int i = 0; i < pts_set.size(); i++)
			{
				if (pts_set.at(i).valid_ornot != -1)
				{
					c_in->at(i).dj = pts_set.at(i).local_density(h0); //dj
					pts_set.at(i).compute_distribution_value();  //delta
					if (isnan(c_in->at(i).dj))
					{
						cout << "current " << i << " dj is: Nan" << endl;
						getchar();
					}

				}
				else
					pts_set.at(i).delta = -1;
			}
		}


		pts_set_down.clear();
#pragma omp parallel  for num(8)
		for (int i = 0; i < c_down->size(); i++)
		{
			int id_in_src = c_down->at(i).id_in_src;
			pts_set_down.push_back(pts_set.at(id_in_src));
		}
		float dmax = 0, dmin = 0;
		ske::compute_dv_range(pts_set_down, dmax, dmin);
		vector<pcl::PointNei> tempvector;
		ave_move_dis = 0;
		stdv = 0;
		tempdis.clear();
		float neibormaxdis_src = 0, neibormaxdis_down = 0;

		clock_t time_renew = time_stamp();
		cout << "Preparation time used: " << (time_renew - time_b)*1.0 / 1000 << "s" << endl;
		
#pragma omp parallel  for
		for (int i = 0; i < pts_set_down.size(); i++)
		{
			if (pts_set_down.at(i).valid_ornot != 0)
				continue;
			if (pts_set_down.at(i).sample_invalid != 0)
				continue;
			if (pts_set_down.at(i).has_fixed == true)
				continue;

			pcl::PointXYZ left, right, temp;
			float current_h = pts_set_down.at(i).current_h;

			left = pts_set_down.at(i).compute_average(pts_set_down.at(i).neibor_point_src, c_in, current_h);
			right = pts_set_down.at(i).compute_repulsion(pts_set_down.at(i).neibor_point_down, c_in, current_h);
			float delta = (pts_set_down.at(i).delta - dmin) / (dmax - dmin);
			if (delta > 0.99)
			{
				int id_in_src = pts_set_down.at(i).id_in_src;
				pts_set_down.at(i).has_fixed = true;
				pts_set.at(id_in_src).has_fixed = true;
			}
			temp.x = left.x + delta * u*right.x;
			temp.y = left.y + delta * u*right.y;
			temp.z = left.z + delta * u*right.z;
			if (isnan(temp.x) || isnan(temp.y) || isnan(temp.z))
				break;
			float tx = pts_set_down.at(i).p.x - temp.x;
			float ty = pts_set_down.at(i).p.y - temp.y;
			float tz = pts_set_down.at(i).p.z - temp.z;
			tempdis.push_back(sqrt(tx*tx + ty * ty + tz * tz));

			pcl::PointNei tempnei;
			tempnei.x = temp.x;
			tempnei.y = temp.y;
			tempnei.z = temp.z;
			tempnei.id_in_down = pts_set_down.at(i).id_in_down;
			tempnei.id_in_src = pts_set_down.at(i).id_in_src;
			tempnei.dj = c_in->at(tempnei.id_in_src).dj;
			tempvector.push_back(tempnei);
		}

#pragma omp parallel  for
		for (int i = 0; i < tempvector.size(); i++)
		{
			pcl::PointNei p = tempvector.at(i);
			int id_in_down = p.id_in_down;
			int id_in_src = p.id_in_src;
			c_in->at(id_in_src) = p;
			c_down->at(id_in_down) = p;

			pts_set.at(id_in_src).p = p;
			pts_set_down.at(id_in_down).p = p;
		}
		tc_down->clear();
		tc_src->clear();
		copyPointCloud(*c_in, *tc_src);
		copyPointCloud(*c_down, *tc_down);

		pts_set_down_out.clear();
		pts_set_down_out.assign(pts_set_down.begin(), pts_set_down.end());

		clock_t time_renewend = time_stamp();
		cout << "Renew time used: " << (time_renewend - time_renew)*1.0 / 1000 << "s" << endl;

		sort(tempdis.begin(), tempdis.end());
		for (int i = 0; i < tempdis.size(); i++)
		{
			ave_move_dis += tempdis.at(i);
		}
		ave_move_dis = ave_move_dis * 1.0 / tempdis.size();

		for (int i = 0; i < tempdis.size(); i++)
		{
			stdv += (tempdis.at(i) - ave_move_dis)*(tempdis.at(i) - ave_move_dis);
		}
		stdv = sqrt(stdv*1.0 / tempdis.size());

		clock_t time_upkdtree = time_stamp();
		if ((ave_move_dis > (h*1.0 / 100) && stdv > (h*1.0 / 100)))
		{
			tc_down->clear();
			tc_src->clear();
			ske::build_kdtree(c_in, t_src, tc_src);
			ske::build_kdtree(c_down, t_down, tc_down);
		}
		if ((ave_move_dis > (h*1.0 / 100) && stdv > (h*1.0 / 100)))
		{
#pragma omp parallel  for
			for (int i = 0; i < pts_set.size(); i++)
			{
				if (pts_set.at(i).id_in_down == -1)
					continue;
				if (pts_set.at(i).has_fixed == true)
					continue;
				if (iterationTimes <= setTimes * 0.75)
				{
					pts_set.at(i).FIndNeiborInSrc(t_src, tc_src, c_in);
				}
				else
				{
					pts_set.at(i).FIndNeiborInSrcFromPlane(t_src, tc_src, c_in);
				}
				if (iterationTimes <= setTimes * 0.75)
					pts_set.at(i).FIndNeiborInDown(t_down, tc_down, c_down);
				else
					pts_set.at(i).FIndNeiborInDownFromPlane(t_down, tc_down, c_down);
				pts_set.at(i).NeiborSrcDownUpate();
			}
		}

#pragma omp parallel  for
		for (int i = 0; i < pts_set_down.size(); i++)
		{
			ske::deltPoint tempp = pts_set_down.at(i);
			sort(tempp.disvec_down.begin(), tempp.disvec_down.end());
			reverse(tempp.disvec_down.begin(), tempp.disvec_down.end());
			int threnum = 0;
			for (int j = 0; j < tempp.disvec_down.size(); j++)
			{
				if (tempp.disvec_down.at(j) > tempp.current_h*0.5)
					threnum++;
				else
					break;
			}
			if (threnum <= 6)
			{
				int idsrc = tempp.id_in_src;
				int iddown = tempp.id_in_down;
				pts_set_down.at(iddown).has_fixed = true;
				pts_set.at(idsrc).has_fixed = true;
			}
		}

#pragma omp parallel  for num(8)
		for (int i = 0; i < pts_set_down.size(); i++)
		{
			if (pts_set_down.at(i).has_fixed == true)
				continue;
			pts_set_down.at(i).p.r = 0;
			pts_set_down.at(i).p.g = 255;
			pts_set_down.at(i).p.b = 0;
			int id = pts_set_down.at(i).id_in_src;
			pts_set.at(id).p = pts_set_down.at(i).p;
		}
		copyDel2PCRGB(pts_set, alltemp);
		v1->UpdateOutRGBAPts(alltemp, "all");

		copyDel2PCRGB(pts_set_down, v1->cloud_show);
		clock_t time_e3 = time_stamp();
		float diffe = abs(ave_move_dis - ave_move_dis_before);
		float radius_thre = h * 1.0 / 400;
		cout << " Average moving distance: " << ave_move_dis << endl;
		cout << " Average moving distance last time: " << ave_move_dis_before << endl;
		cout << " Average moving distance difference: " << abs(ave_move_dis - ave_move_dis_before) << endl;
		cout << " Current radius : " << pts_set.at(0).current_h << endl;
		cout << " Current radius threshold: " << radius_thre << endl;
		cout << "Updating Kdtree and Neiboring points time used: " << (time_e3 - time_upkdtree)*1.0 / 1000 << "s" << endl;
		cout << "Total time used: " << (time_e3 - time_b)*1.0 / 1000 << "s" << endl;
		ave_move_dis_before = ave_move_dis;


}
}

void SkeletonSteBiggerH(pcl::PointCloud<pcl::PointNei>::Ptr& c_in,
	pcl::PointCloud<pcl::PointNei>::Ptr& c_down, vector<ske::deltPoint>& pts_set,
	vector<ske::deltPoint> &pts_set_down, float &h0, float u, int setTimes,
	ViewerThread *v1)
{
	v1->UpdateMySphere(pts_set_down.at(100).p.x, pts_set_down.at(100).p.y, pts_set_down.at(100).p.z, h0, "h0");
	cout << pts_set_down.at(100).id_in_src << endl;

	int iterationTimes = 0;
	float h = h0;
	pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr t_src(new pcl::KdTreeFLANN<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr tc_src(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr t_down(new pcl::KdTreeFLANN<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr tc_down(new pcl::PointCloud<pcl::PointXYZ>);

	float ave_move_dis = DBL_MAX;
	float ave_move_dis_before = DBL_MAX;
	float stdv = DBL_MAX;
	vector<float> tempdis;
	float maxdis = 0;
	//vector<ske::deltPoint> pts_set_down;
#pragma omp parallel  for
	for (int i = 0; i < pts_set.size(); i++)
	{
		pts_set.at(i).current_h = h;
	}

	clock_t time_ini = time_stamp();
	ske::build_kdtree(c_in, t_src, tc_src);
	ske::build_kdtree(c_down, t_down, tc_down);

#pragma omp parallel  for
	for (int i = 0; i < pts_set.size(); i++)
	{
		if (pts_set.at(i).has_fixed == true)
			continue;
		if (pts_set.at(i).id_in_down == -1)
			continue;
		pts_set.at(i).FIndNeiborInSrcFromPlane(t_src, tc_src, c_in);
		pts_set.at(i).FIndNeiborInDownFromPlane(t_down, tc_down, c_down);
		pts_set.at(i).NeiborSrcDownUpate();
	}

	for (int i = 0; i < pts_set_down.size(); i++)
	{
		int id_in_src = pts_set_down.at(i).id_in_src;
		pts_set_down.at(i) = pts_set.at(id_in_src);
	}
	clock_t time_ini_end = time_stamp();
	cout << "Initial Kdtree time used: " << (time_ini_end - time_ini)*1.0 / 1000 << "s" << endl;

	while (iterationTimes < setTimes)
	{
		iterationTimes += 1;
		cout << "**************Iteration Times: " << iterationTimes << "**************" << endl;
		clock_t time_b = time_stamp();
	
#pragma omp parallel  for
		for (int i = 0; i < pts_set_down.size(); i++)
		{
			int id_in_src = pts_set_down.at(i).id_in_src;
			if (pts_set_down.at(i).neibor_point_down->size() < 2)
			{
		
				pts_set.at(id_in_src).sample_invalid = -1;
				pts_set.at(id_in_src).valid_ornot = -1;
				pts_set_down.at(i).sample_invalid = -1;
				pts_set_down.at(i).valid_ornot = -1;
			}
		}

#pragma omp parallel  for num(8)
		for (int i = 0; i < pts_set_down.size(); i++)
		{
			int id_in_src = pts_set_down.at(i).id_in_src;
				c_in->at(id_in_src).dj = pts_set.at(id_in_src).local_density(h0); //dj
			pts_set.at(id_in_src).compute_distribution_value();  //delta
			c_down->at(i).dj = c_in->at(id_in_src).dj;
			pts_set_down.at(i).delta = pts_set.at(id_in_src).delta;
		}
		clock_t time_dj_delta = time_stamp();
		cout << "dj_delta time used: " << (time_dj_delta - time_b)*1.0 / 1000 << "s" << endl;
		float dmax = 0, dmin = 0;
		ske::compute_dv_range(pts_set_down, dmax, dmin);
		vector<pcl::PointNei> tempvector;
		ave_move_dis = 0;
		stdv = 0;
		tempdis.clear();
		float neibormaxdis_src = 0, neibormaxdis_down = 0;

		clock_t time_renew = time_stamp();
		cout << "delta update time used: " << (-time_dj_delta + time_renew)*1.0 / 1000 << "s" << endl;
		cout << "Preparation time used: " << (time_renew - time_b)*1.0 / 1000 << "s" << endl;

#pragma omp parallel  for
		for (int i = 0; i < pts_set_down.size(); i++)
		{
			if (pts_set_down.at(i).valid_ornot != 0)
				continue;
			if (pts_set_down.at(i).sample_invalid != 0)
				continue;
			if (pts_set_down.at(i).has_fixed == true)
				continue;
			pcl::PointXYZ left, right, temp;
			float current_h = pts_set_down.at(i).current_h;
			left = pts_set_down.at(i).compute_average(pts_set_down.at(i).neibor_point_src, c_in, current_h);
			right = pts_set_down.at(i).compute_repulsion(pts_set_down.at(i).neibor_point_down, c_in, current_h);
			float delta = (pts_set_down.at(i).delta - dmin) / (dmax - dmin);
			if (delta > 0.98)
			{
				int id_in_src = pts_set_down.at(i).id_in_src;
				pts_set_down.at(i).has_fixed = true;
				pts_set.at(id_in_src).has_fixed = true;
			}
			temp.x = left.x + delta * u*right.x;
			temp.y = left.y + delta * u*right.y;
			temp.z = left.z + delta * u*right.z;
			if (isnan(temp.x) || isnan(temp.y) || isnan(temp.z))
				break;
			float tx = pts_set_down.at(i).p.x - temp.x;
			float ty = pts_set_down.at(i).p.y - temp.y;
			float tz = pts_set_down.at(i).p.z - temp.z;
			tempdis.push_back(sqrt(tx*tx + ty * ty + tz * tz));

			pcl::PointNei tempnei;
			tempnei.x = temp.x;
			tempnei.y = temp.y;
			tempnei.z = temp.z;
			tempnei.id_in_down = pts_set_down.at(i).id_in_down;
			tempnei.id_in_src = pts_set_down.at(i).id_in_src;
			tempnei.dj = c_in->at(tempnei.id_in_src).dj;
			tempvector.push_back(tempnei);
		}

#pragma omp parallel  for
		for (int i = 0; i < tempvector.size(); i++)
		{
			pcl::PointNei p = tempvector.at(i);
			int id_in_down = p.id_in_down;
			int id_in_src = p.id_in_src;		
			c_in->at(id_in_src) = p;		
			c_down->at(id_in_down) = p;

			pts_set.at(id_in_src).p = p;
			pts_set_down.at(id_in_down).p = p;
		}
		tc_down->clear();
		tc_src->clear();
		copyPointCloud(*c_in, *tc_src);
		copyPointCloud(*c_down, *tc_down);

		clock_t time_renewend = time_stamp();
		cout << "Renew time used: " << (time_renewend - time_renew)*1.0 / 1000 << "s" << endl;
		
		sort(tempdis.begin(), tempdis.end());
		for (int i = 0; i < tempdis.size(); i++)
		{
			ave_move_dis += tempdis.at(i);
		}
		ave_move_dis = ave_move_dis * 1.0 / tempdis.size();

		for (int i = 0; i < tempdis.size(); i++)
		{
			stdv += (tempdis.at(i) - ave_move_dis)*(tempdis.at(i) - ave_move_dis);
		}
		stdv = sqrt(stdv*1.0 / tempdis.size());

		
		clock_t time_upkdtree = time_stamp();
		if ((ave_move_dis > (h*1.0 / 100) && stdv > (h*1.0 / 100)))
		{
			tc_down->clear();
			tc_src->clear();
			ske::build_kdtree(c_in, t_src, tc_src);
			ske::build_kdtree(c_down, t_down, tc_down);
		}
		if ((ave_move_dis > (h*1.0 / 100) && stdv > (h*1.0 / 100)))
		{
#pragma omp parallel  for
			for (int i = 0; i < pts_set.size(); i++)
			{
				if (pts_set.at(i).id_in_down == -1)
					continue;
				if (pts_set.at(i).has_fixed == true)
					continue;
				if (iterationTimes <= setTimes * 0.75)
				{
					pts_set.at(i).FIndNeiborInSrc(t_src, tc_src, c_in)
				}
				else
				{
					pts_set.at(i).FIndNeiborInSrcFromPlane(t_src, tc_src, c_in);
				}
				if (iterationTimes <= setTimes * 0.75)
					pts_set.at(i).FIndNeiborInDown(t_down, tc_down, c_down);
				else
					pts_set.at(i).FIndNeiborInDownFromPlane(t_down, tc_down, c_down);
				pts_set.at(i).NeiborSrcDownUpate();
			}
		
		}
		clock_t time_finishkdtree = time_stamp();
	
#pragma omp parallel  for
		for (int i = 0; i < pts_set_down.size(); i++)
		{
			ske::deltPoint tempp = pts_set_down.at(i);
			sort(tempp.disvec_down.begin(), tempp.disvec_down.end());
			reverse(tempp.disvec_down.begin(), tempp.disvec_down.end());
			int threnum = 0;
			for (int j = 0; j < tempp.disvec_down.size(); j++)
			{
				if (tempp.disvec_down.at(j) > tempp.current_h*0.5)
					threnum++;
				else
					break;
			}
			if (threnum <= 6)
			{
				int idsrc = tempp.id_in_src;
				int iddown = tempp.id_in_down;
				pts_set_down.at(iddown).has_fixed = true;
				pts_set.at(idsrc).has_fixed = true;
			}
		}
		clock_t time_resample = time_stamp();

		copyDel2PCRGB(pts_set_down, v1->cloud_show);
		v1->UpdateMySphere(pts_set_down.at(100).p.x, pts_set_down.at(100).p.y, pts_set_down.at(100).p.z, h0, "h0");
		clock_t time_e3 = time_stamp();
		float diffe = abs(ave_move_dis - ave_move_dis_before);
		float radius_thre = h * 1.0 / 400;
		cout << " Average moving distance: " << ave_move_dis << endl;
		cout << " Average moving distance last time: " << ave_move_dis_before << endl;
		cout << " Average moving distance difference: " << abs(ave_move_dis - ave_move_dis_before) << endl;
		cout << " Current radius : " << pts_set.at(0).current_h << endl;
		cout << " Current radius threshold: " << radius_thre << endl;
		cout << " Updating Kdtree and Neiboring points time used: " << (time_e3 - time_upkdtree)*1.0 / 1000 << "s" << endl;
		cout << " Total time used: " << (time_e3 - time_b)*1.0 / 1000 << "s" << endl;
		ave_move_dis_before = ave_move_dis;
	}
}

void PointsFix(pcl::PointCloud<pcl::PointNei>::Ptr& c_in, pcl::PointCloud<pcl::PointNei>::Ptr& c_down,
	vector<ske::deltPoint>& pts_set, vector<ske::deltPoint>& pts_set_down,
	vector<ske::deltPoint>& pts_set_down_out, float &h0,
	ViewerThread *v1)
{
	cout << endl;
	cout << "**************Fixing Points" << "**************" << endl;
	clock_t time_smooth_begin = time_stamp();
	float dmax = 0, dmin = 0;
	ske::compute_dv_range(pts_set_down, dmax, dmin);
#pragma omp parallel  for
	for (int i = 0; i < pts_set_down.size(); i++)
	{
		if (pts_set_down.at(i).valid_ornot != 0)
			continue;
		float delta = (pts_set_down.at(i).delta - dmin) / (dmax - dmin);
		pts_set_down.at(i).delta = delta;
	}
	clock_t time_delete = time_stamp();
	for (int i = 0; i < pts_set_down.size(); i++)
	{
		if (pts_set_down.at(i).valid_ornot != 0)
			continue;
		if (pts_set_down.at(i).has_fixed == true)
			continue;
		vector<ske::deltPoint> temp;
		ske::deltPoint p = pts_set_down.at(i);
		if (p.has_moved == true)
			continue;
#pragma omp parallel  for
		for (int j = 0; j < p.disvec_down.size(); j++)
		{
			if (p.disvec_down.at(j) <= (h0*1.0 / 1000))
			{
				int tempid = p.neibor_point_down->at(j).id_in_down;
				pts_set_down.at(tempid).has_moved = true;
				temp.push_back(pts_set_down.at(tempid));
			}
			else
				continue;
		}
		if (temp.size() == 0)
			continue;
		float x = pts_set_down.at(i).p.x,
			y = pts_set_down.at(i).p.y,
			z = pts_set_down.at(i).p.z;
		for (int j = 0; j < temp.size(); j++)
		{
			x += temp.at(j).p.x;
			y += temp.at(j).p.y;
			z += temp.at(j).p.z;
		}
		pts_set_down.at(i).p.x = x * 1.0 / (temp.size() + 1);
		pts_set_down.at(i).p.y = y * 1.0 / (temp.size() + 1);
		pts_set_down.at(i).p.z = z * 1.0 / (temp.size() + 1);
	}
	clock_t time_delete_end = time_stamp();
	cout << "Delete Points Time used: " << (time_delete_end - time_delete)*1.0 / 1000 << "s" << endl;
	pcl::PointCloud<pcl::PointNei>::Ptr c_down_out(new pcl::PointCloud<pcl::PointNei>);
	for (int i = 0; i < pts_set_down.size(); i++)
	{
		pcl::PointNei p;
		p = pts_set_down.at(i).p;
		c_down_out->push_back(p);
	}
	pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr tree_knn(new pcl::KdTreeFLANN<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr t_knn_down(new pcl::PointCloud<pcl::PointXYZ>);
	ske::build_kdtree(c_down_out, tree_knn, t_knn_down);

	clock_t time_fix = time_stamp();
	vector<float> delta_vec;
#pragma omp parallel  for
	for (int i = 0; i < pts_set_down.size(); i++)
	{
		if (pts_set_down.at(i).has_fixed == true)
			continue;
		if (pts_set_down.at(i).has_moved == true)
			continue;
		if (pts_set_down.at(i).sample_invalid != 0)
			continue;
		if (pts_set_down.at(i).valid_ornot != 0)
			continue;
		vector<int> kdid;
		vector<float> kddis;
		int id_in_down = pts_set_down.at(i).id_in_down;
		float delta = pts_set_down.at(i).delta;
		int searchsize = 7;
		if (pts_set_down.at(i).neibor_point_down->size() < 7)
			searchsize = pts_set_down.at(i).neibor_point_down->size();
		cout << "id_in_down = " << id_in_down << endl;
		if (tree_knn->nearestKSearch(t_knn_down->at(id_in_down), searchsize, kdid, kddis) > 0)
		{
			for (int j = 1; j < kdid.size(); j++)
			{
				if (kddis.at(j) == 0)
					continue;
				delta += pts_set_down.at(kdid.at(j)).delta;
			}
			delta = delta * 1.0 / 7;
		}
		delta_vec.push_back(delta);
				if (delta > 0.9)
		{
			int id_in_src = pts_set_down.at(i).id_in_src;
			pts_set_down.at(i).has_fixed = true;
			pts_set.at(id_in_src).has_fixed = true;
		}

	}
	pts_set_down_out.clear();
	pts_set_down_out.resize(pts_set_down.size());
	pts_set_down_out.assign(pts_set_down.begin(), pts_set_down.end());
	clock_t time_smooth_end = time_stamp();
	copyDel2PCRGB(pts_set_down_out, v1->cloud_show);
	cout << "Fixint time used" << (time_smooth_end - time_fix)*1.0 / 1000 << "s" << endl;
	cout << "Smoothing and fixing time used: " << (time_smooth_end - time_smooth_begin)*1.0 / 1000 << "s" << endl;
}


void ConnectPts(pcl::PointCloud<pcl::PointNei>::Ptr& c_in, pcl::PointCloud<pcl::PointNei>::Ptr& c_down,
	vector<ske::deltPoint>& pts_set, vector<ske::deltPoint>& pts_set_down_in,
	vector<ske::deltPoint>& pts_set_down, float &resolution)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr c_octree(new pcl::PointCloud<pcl::PointXYZRGB>);
	copyPointCloud(*c_down, *c_octree);
	pcl::octree::OctreePointCloudSearch<pcl::PointXYZRGB>::Ptr
		tmp_octree(new pcl::octree::OctreePointCloudSearch<pcl::PointXYZRGB>(resolution));

	tmp_octree->setInputCloud(c_octree); 
	tmp_octree->addPointsFromInputCloud();
	cout << "Octree Depth of Original Cloud: " << tmp_octree->getTreeDepth() << endl;
	cout << "Number of Leaf Nodes in Original Cloud: " << tmp_octree->getLeafCount() << endl;

	pcl::PointXYZRGB p = c_octree->at(0);
	vector<int> pointidx_vec;

	clock_t oct1 = time_stamp();
	vector<int> pointidx_vecall;
	for (int i = 0; i < 1; i++)
	{
		pcl::PointXYZRGB p = c_octree->at(i);
		vector<int> pointidx_vec;
		if (tmp_octree->voxelSearch(p, pointidx_vec))
		{
			std::cout << "Neighbors within voxel search at (" << p.x
				<< " " << p.y
				<< " " << p.z << ")"
				<< std::endl;

			std::cout << "Searched Points Number : " << pointidx_vec.size() << endl;

			std::cout << "Leaf Count : " << tmp_octree->getLeafCount() << std::endl;             // 叶子数
			std::cout << "Tree Depth : " << tmp_octree->getTreeDepth() << std::endl;             // 八叉树深度
			std::cout << "Branch Count : " << tmp_octree->getBranchCount() << std::endl;         // 非叶子结点数
			std::cout << "Voxel Diameter : " << tmp_octree->getVoxelSquaredDiameter() << std::endl;  // ???Voxel Side Length*3
			std::cout << "Voxel Side Length : " << tmp_octree->getVoxelSquaredSideLen() << std::endl;// 分辨率的平方
			double minx, miny, minz, maxx, maxy, maxz;
			tmp_octree->getBoundingBox(minx, miny, minz, maxx, maxy, maxz);
			std::cout << "BoundingBox: " << "(" << minx << " - " << maxx << ")" << " , "
				<< "(" << miny << " - " << maxy << ")" << " , "
				<< "(" << minz << " - " << maxz << ")" << std::endl;       // 整个八叉树的范围
		}
	}
	for (int j = 0; j < pointidx_vec.size(); j++)
	{
		pointidx_vecall.push_back(pointidx_vec.at(j));
	}

	clock_t oct2 = time_stamp();
	cout << "20 points octree neiboring search time used: " << (oct2 - oct1)*1.0 / 1000 << "s" << endl;
}


void ResampleClusters(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &ske_in,
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr& skeout)
{
	pcl::KdTreeFLANN<pcl::PointXYZRGB>::Ptr tree(new pcl::KdTreeFLANN<pcl::PointXYZRGB>);

	tree->setInputCloud(ske_in);
	for (int i = 0; i < ske_in->size(); i++)
	{
		if (ske_in->at(i).x == 0 && ske_in->at(i).y == 0 &&
			ske_in->at(i).z == 0)
			continue;
		float x = ske_in->at(i).x,
			y = ske_in->at(i).y, z = ske_in->at(i).z;
		vector<int> kdid;
		vector<float> kddis;
		if (tree->radiusSearch(ske_in->at(i), 3, kdid, kddis) > 0)
		{
			int gathernum = 1;
			for (int j = 1; j < kdid.size(); j++)
			{
				if (ske_in->at(kdid.at(j)).x == 0 &&
					ske_in->at(kdid.at(j)).y == 0 &&
					ske_in->at(kdid.at(j)).z == 0)
					continue;
				if (sqrt(kddis.at(j)) < 1.5)
				{
					x += ske_in->at(kdid.at(j)).x;
					y += ske_in->at(kdid.at(j)).y;
					z += ske_in->at(kdid.at(j)).z;
					ske_in->at(kdid.at(j)).x = 0;
					ske_in->at(kdid.at(j)).y = 0;
					ske_in->at(kdid.at(j)).z = 0;
					gathernum++;
				}
			}
			ske_in->at(i).x = x * 1.0 / gathernum;
			ske_in->at(i).y = y * 1.0 / gathernum;
			ske_in->at(i).z = z * 1.0 / gathernum;
			skeout->push_back(ske_in->at(i));
		}
	}
	cout << "Finish resampling! " << endl;
}

int Processing(ViewerThread *v1)
{
	vector<string> files;
	vector<string> filenames;

	cout << "Adjust coordinate" << endl;
	map<string, string> param;
	Files::ConfigFileRead(param, "./ske_id.ini");
	parameters.configed(param);
	cout << parameters.outformat << endl;
	Files::getFiles(parameters.inpath,files,filenames,parameters.informat);
	files[0] = parameters.infile;
	pcl::PointCloud<pcl::PointXYZRGB> * ini = new pcl::PointCloud<pcl::PointXYZRGB>;
	for (int i = 0; i < 1; i++)
	{
		if (ini->empty()!=true)
		{
			ini->clear();
		}
		if (c_in_ini->empty() != true)
		{
			c_in_ini->clear();
		}
		if (c_down_ini->empty() != true)
		{
			c_down_ini->clear();
		}
		if (c_down_rgb_ini->empty() != true)
		{
			c_down_rgb_ini->clear();
		}
		if (pts_set_ini.empty()!=true)
		{
			pts_set_ini.clear();
		}
		if (pts_set_down_ini.empty() != true)
		{
			pts_set_down_ini.clear();
		}

		string outfilename;
		outfilename.append(parameters.outpaths)
			.append("\\").append(filenames[i]).
			append("-ske").append(parameters.outformat);

		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_in(new pcl::PointCloud<pcl::PointXYZRGB>());
		int error = Files::Load_Pts_format(cloud_in, files[i],parameters.informat);
		if (error == 1)
			return -1;
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_temp(new pcl::PointCloud<pcl::PointXYZRGB>());
		for (int i = 0; i < cloud_in->size(); i++)
		{
			if (cloud_in->at(i).x == 0 && cloud_in->at(i).y == 0 && cloud_in->at(i).z == 0)
				continue;
			cloud_temp->push_back(cloud_in->at(i));
		}
		cloud_in->resize(cloud_temp->size());
		copyPointCloud(*cloud_temp, *cloud_in);
		copyPointCloud(*cloud_temp, *v1->cloud_show);
		float hini = getInitialRadius(cloud_in)*1.0;
		float h0 = getInitialRadius(cloud_in)*1.0*parameters.h_times;
		cout << "Radius ini: " << getInitialRadius(cloud_in) << endl;
		cout << "Radius: " << h0 << endl;
		v1->normal_size = h0;

		if (parameters.model == 1)
			initialization(cloud_in, c_in_ini, c_down_ini, pts_set_ini, h0, parameters.rate);
		else
			initialization_2(cloud_in, c_in_ini, c_down_ini, pts_set_ini, h0, parameters.rate, parameters.leafsize);

		cout << c_down_ini->size() << endl;
		int setTimes = parameters.iterationTimes; 
		float u = parameters.u; 

		float current_h = h0;
		vector<ske::deltPoint> pts_set_down_2;
		vector<vector<ske::deltPoint>> skeall;

		OneStep(c_in_ini, c_down_ini, pts_set_ini, pts_set_down_ini, current_h, u, setTimes, v1);
		current_h += 0.25*h0;
		{
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr c_show_temp(new pcl::PointCloud<pcl::PointXYZRGB>);
			copyDel2PCRGB(pts_set_down_ini, c_show_temp);
			copyPointCloud(*c_show_temp, *v1->cloud_show);
			PointsFix(c_in_ini, c_down_ini, pts_set_ini, pts_set_down_ini, pts_set_down_2, current_h, v1);
			copyDel2PCRGB(pts_set_down_ini, c_show_temp);
			copyPointCloud(*c_show_temp, *v1->cloud_show);
			int fixednum = 0;
			int fixednum_down = 0;
#pragma omp parallel for 
			for (int i = 0; i < pts_set_ini.size(); i++)
			{
				if (pts_set_ini.at(i).has_fixed == true)
				{
					fixednum++;
				}
			}
			skeall.push_back(pts_set_down_2);
#pragma omp parallel for 
			for (int i = 0; i < pts_set_down_ini.size(); i++)
			{
				if (pts_set_down_ini.at(i).has_fixed == true)
				{
					int id = pts_set_down_2.at(i).id_in_down;
					int idsrc = pts_set_down_2.at(i).id_in_src;
					pts_set_down_ini.at(i).p.r = 255;
					pts_set_down_ini.at(i).p.g = 0;
					pts_set_down_ini.at(i).p.b = 0;
					c_down_ini->at(id).r = 255;
					c_down_ini->at(id).g = 0;
					c_down_ini->at(id).b = 0;
					c_in_ini->at(idsrc).r = 255;
					c_in_ini->at(idsrc).g = 0;
					c_in_ini->at(idsrc).b = 0;
					pts_set_ini.at(idsrc).p.r = 255;
					pts_set_ini.at(idsrc).p.g = 0;
					pts_set_ini.at(idsrc).p.b = 0;
					fixednum_down++;
				}
			}
			cout << "Fixed points number: " << fixednum << endl;
			cout << "Fixed points number down: " << fixednum_down << endl;
			cout << "------------------------------------------------" << endl;
			cout << "------------------------------------------------" << endl;
			cout << endl;
		}

		while (current_h <= (parameters.h_finals*hini))
		{
			cout << "-----------------Current H: " << current_h << "-------------------" << endl;
			SkeletonSteBiggerH(c_in_ini, c_down_ini, pts_set_ini, pts_set_down_ini, current_h, u, setTimes, v1);
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr c_show_temp(new pcl::PointCloud<pcl::PointXYZRGB>);
			copyDel2PCRGB(pts_set_down_ini, c_show_temp);
			copyPointCloud(*c_show_temp, *v1->cloud_show);
			PointsFix(c_in_ini, c_down_ini, pts_set_ini, pts_set_down_ini, pts_set_down_2, current_h, v1);
			copyDel2PCRGB(pts_set_down_ini, c_show_temp);
			copyPointCloud(*c_show_temp, *v1->cloud_show);
			int fixednum = 0;
			int fixednum_down = 0;
#pragma omp parallel for 
			for (int i = 0; i < pts_set_ini.size(); i++)
			{
				if (pts_set_ini.at(i).has_fixed == true)
				{
					fixednum++;
				}
			}
			current_h += 0.25*h0;
			skeall.push_back(pts_set_down_2);

#pragma omp parallel for 
			for (int i = 0; i < pts_set_down_ini.size(); i++)
			{
				if (pts_set_down_ini.at(i).has_fixed == true)
				{
					int id = pts_set_down_2.at(i).id_in_down;
					int idsrc = pts_set_down_2.at(i).id_in_src;
					pts_set_down_ini.at(i).p.r = 255;
					pts_set_down_ini.at(i).p.g = 0;
					pts_set_down_ini.at(i).p.b = 0;
					c_down_ini->at(id).r = 255;
					c_down_ini->at(id).g = 0;
					c_down_ini->at(id).b = 0;
					c_in_ini->at(idsrc).r = 255;
					c_in_ini->at(idsrc).g = 0;
					c_in_ini->at(idsrc).b = 0;
					pts_set_ini.at(idsrc).p.r = 255;
					pts_set_ini.at(idsrc).p.g = 0;
					pts_set_ini.at(idsrc).p.b = 0;
					fixednum_down++;
				}
			}
			cout << "Fixed points number: " << fixednum << endl;
			cout << "Fixed points number down: " << fixednum_down << endl;
			cout << "------------------------------------------------" << endl;
			cout << "------------------------------------------------" << endl;
			cout << endl;
		}
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr c_down_rgb2(new pcl::PointCloud<pcl::PointXYZRGB>);
#pragma omp parallel for 
		for (int i = 0; i < skeall.size(); i++)
		{
			vector<ske::deltPoint> temp = skeall.at(i);
			for (int j = 0; j < temp.size(); j++)
			{
				if (temp.at(j).has_fixed != true)
					continue;
				pcl::PointXYZRGB tp;
				tp.x = temp.at(j).p.x;
				tp.y = temp.at(j).p.y;
				tp.z = temp.at(j).p.z;
				tp.g = 0;
				if (temp.at(j).has_fixed == true)
				{
					tp.r = 255;
					tp.b = 0;
				}
				else
				{
					tp.r = 0; tp.b = 255;
				}
				c_down_rgb2->push_back(tp);
			}
		}
		cout << c_down_rgb2->size() << endl;

#pragma omp parallel for 
		for (int i = 0; i < c_down_rgb2->size(); i++)
		{
			pcl::PointXYZRGB p1 = c_down_rgb2->at(i);
			for (int j = 0; j < c_down_rgb2->size(); j++)
			{
				pcl::PointXYZRGB tp = c_down_rgb2->at(j);
				float x = p1.x - tp.x;
				float y = p1.y - tp.y;
				float z = p1.z - tp.z;
				float dis = sqrt(x*x + y * y + z * z);

				if (dis < h0*1.0 / 10000)
				{
					if (tp.r = 255 || p1.r == 255)
					{
						tp.r = 255;
						tp.b = 0;
						p1.r = 255;
						p1.b = 0;
					}
				}
			}
		}
		copyPointNei2PCRGB(c_down_ini, c_down_rgb_ini);

		for (int i = 0; i < pts_set_down_2.size(); i++)
		{
			if (pts_set_down_2.at(i).has_fixed == true)
			{
				int id = pts_set_down_2.at(i).id_in_down;
				int idsrc = pts_set_down_2.at(i).id_in_src;
				c_down_rgb_ini->at(i).r = 255;
				c_down_rgb_ini->at(i).b = 0;
				pts_set_down_2.at(i).p.r = 255;
				pts_set_down_2.at(i).p.g = 0;
				pts_set_down_2.at(i).p.b = 0;
				pts_set_down_ini.at(i).p.r = 255;
				pts_set_down_ini.at(i).p.g = 0;
				pts_set_down_ini.at(i).p.b = 0;
				c_down_ini->at(id).r = 255;
				c_down_ini->at(id).g = 0;
				c_down_ini->at(id).b = 0;
				c_in_ini->at(idsrc).r = 255;
				c_in_ini->at(idsrc).g = 0;
				c_in_ini->at(idsrc).b = 0;
				pts_set_ini.at(idsrc).p.r = 255;
				pts_set_ini.at(idsrc).p.g = 0;
				pts_set_ini.at(idsrc).p.b = 0;
			}
		}
	
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr ske_resample(new pcl::PointCloud<pcl::PointXYZRGB>);
		ResampleClusters(c_down_rgb_ini, ske_resample);
		pcl::copyPointCloud(*ske_resample, *c_down_rgb_ini);
		cout << ske_resample->size() << endl;

		Files::Save_Pts_format(ske_resample, outfilename, parameters.outformat);
		cout << "Finish extract ske" << endl;

	}
	return 0;
}

int GetSkeBranch(ViewerThread *v1)
{
	while (true)
	{
		if (extrac_ske_finished == 1)
			break;
	
	}
	cout << "Start Growing Branches. " << endl;
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr
		cloud_with_normals(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
	FindDireSke(c_down_rgb_ini, pts_set_down_ini, cloud_with_normals);
	copyPointCloud(*c_down_rgb_ini, *v1->cloud_show);
	copyPointCloud(*cloud_with_normals, *v1->cloud_with_normals);
	v1->addNormal();
	vector<pcl::PointCloud<pcl::PointNei>::Ptr> cs_out;
	ConnectSke(c_down_ini, pts_set_down_ini, cloud_with_normals, cs_out);
	cout << "Finish Growing Branches. " << endl;

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr branch_rgb(new pcl::PointCloud<pcl::PointXYZRGB>);
	for (int i = 0; i < cs_out.size(); i++)
	{
		float r = (rand() % 256); // *1.0 / 100;
		float g = (rand() % 256); // *1.0 / 100;
		float b = (rand() % 256); // *1.0 / 100;
		for (int j = 0; j < cs_out.at(i)->size(); j++)
		{
			pcl::PointXYZRGB temp;
			temp.x = cs_out.at(i)->at(j).x;
			temp.y = cs_out.at(i)->at(j).y;
			temp.z = cs_out.at(i)->at(j).z;
			temp.r = r;
			temp.b = b;
			temp.g = g;
			branch_rgb->push_back(temp);
		}
	}
	cout << "Enter 'Enter' and Save Branches. " << endl;
	getchar();
	copyPointCloud(*branch_rgb, *v1->cloud_show);
	pcl::io::savePLYFile(parameters.skeoutpath, *branch_rgb);
	return 0;
}


int viewerthread();

int main_ske()
{
	//Processing();
	cout << "Ready.....";
	getchar();
	ViewerThread *v1 = new ViewerThread("v1");
	boost::thread vthread(Processing, v1);
	//boost::thread connect_ske(GetSkeBranch, v1);
	/*Processing();*/
	//viewerthread();
	vtkOutputWindow::SetGlobalWarningDisplay(0);
	boost::shared_ptr<pcl::visualization::PCLVisualizer>
		tv(new pcl::visualization::PCLVisualizer(v1->viewername));
	v1->viewer = boost::shared_ptr<pcl::visualization::PCLVisualizer>(tv);
	v1->viewer->setBackgroundColor(0.85, 0.85, 0.85);
	v1->viewer->addPointCloud(v1->cloud_show, "skeleton");
	v1->updatePointCloud();
	while (!v1->viewer->wasStopped())
	{
		v1->viewer->removePointCloud("skeleton");
		//v1->viewer->removeAllPointClouds();
		v1->viewer->addPointCloud(v1->cloud_show, "skeleton");
		v1->viewer->spinOnce(100000);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
	vthread.join();
	getchar();
	return 0;
}

int main()
{
	main_ske();
}

int viewerthread()
{
	cout << "Ready.....";
	getchar();
	cout << "Adjust coordinate" << endl;
	map<string, string> param;
	Files::ConfigFileRead(param, "./ske.ini");
	Files::ConfigFile parameters;
	parameters.configed(param);
	cout << parameters.outformat << endl;
	pcl::PointCloud<pcl::PointXYZRGB> * ini = new pcl::PointCloud<pcl::PointXYZRGB>;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_in(new pcl::PointCloud<pcl::PointXYZRGB>());
	int error = Files::Load_Pts(cloud_in, parameters);
	if (error == 1)
		return -1;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_temp(new pcl::PointCloud<pcl::PointXYZRGB>());
	for (int i = 0; i < cloud_in->size(); i++)
	{
		if (cloud_in->at(i).x == 0 && cloud_in->at(i).y == 0 && cloud_in->at(i).z == 0)
			continue;
		cloud_temp->push_back(cloud_in->at(i));
	}
	cloud_in->resize(cloud_temp->size());
	copyPointCloud(*cloud_temp, *cloud_in);

	boost::thread vthread(&viewerRunner);

	for (int i = 0; i < cloud_in->size(); i++)
	{
		cloud_in->at(i).r = 0;
		cloud_in->at(i).g = 255;
		cloud_in->at(i).b = 0;
	}
	while (1)
	{
		int  value = 255 * (rand() % 100 * 1.0 / 100);
		for (int i = 0; i < cloud_in->size(); i++)
		{
			cloud_in->at(i).r = value;
			cloud_in->at(i).b = 255 * (rand() % 100 * 1.0 / 100);
			cloud_in->at(i).g = 255 * (rand() % 100 * 1.0 / 100);
		}
		copyPointCloud(*cloud_in, *c_down_show);
	}
	vthread.join();
}
