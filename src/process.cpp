#include "process.h"

void getdownsample(pcl::PointCloud<pcl::PointNei>::Ptr &input,
	pcl::PointCloud<pcl::PointNei>::Ptr &output, float rate)
{
	vector<int> downIndices;
	for (int i = 0; i < input->size(); i++)
	{
		downIndices.push_back(i);
	}
	random_shuffle(downIndices.begin(), downIndices.end());
	for (int i = 0; i < input->size()*rate; i++)
	{
		pcl::PointNei tempp;
		int point_id = downIndices.at(i);
		input->at(point_id).id_in_src = point_id;
		input->at(point_id).id_in_down = i;
		tempp = input->at(point_id);
		output->push_back(tempp);
	}
}

void pointCloudAdd(pcl::PointCloud<pcl::PointNei>::Ptr& c_in,
	pcl::PointCloud<pcl::PointNei>::Ptr& c_down)
{
	int in_size = c_in->size();
	for (int i = 0; i < c_down->size(); i++)
	{
		c_down->at(i).id_in_src = in_size + i;
		c_in->push_back(c_down->at(i));
	}

}


void getdownsample_tra(pcl::PointCloud<pcl::PointNei>::Ptr &input,
	pcl::PointCloud<pcl::PointNei>::Ptr &output, float rate, float leavesize)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr temp(new pcl::PointCloud<pcl::PointXYZRGB>);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr tempout(new pcl::PointCloud<pcl::PointXYZRGB>);
	copyPointNei2PCRGB(input, temp);

	pcl::VoxelGrid<pcl::PointXYZRGB> vg;
	vg.setInputCloud(temp);
	vg.setLeafSize(leavesize, leavesize, leavesize);//change leaf size into 0.5cm
	vg.filter(*tempout);
	
	copySrc2PointNei(tempout, output);
	for (int i = 0; i < output->size(); i++)
	{
		output->at(i).id_in_down = i;
	}
	pointCloudAdd(input, output);

}

float getInitialRadius(const pcl::PointCloud<pcl::PointXYZRGB>::Ptr &cloudin)
{
	float temxmax = -INT_MAX, temymax = -INT_MAX, temzmax = -INT_MAX,
		temxmin = INT_MAX, temymin = INT_MAX, temzmin = INT_MAX;
	pcl::PointXYZRGB zmaxpoint, zminpoint, xminpoint, xmaxpoint, ymaxpoint, yminpoint;
	omp_set_num_threads(8);
#pragma omp parallel for
	for (int i = 0; i < cloudin->size(); i++) {
		auto point = cloudin->at(i);
		if (point.x == 0 || point.y == 0 || point.z == 0)
			continue;
		else
		{
			if (point.x <= temxmin)
			{
				temxmin = point.x;
				xminpoint = point;
			}
			if (point.x >= temxmax)
			{
				temxmax = point.x;
				xmaxpoint = point;
			}
			if (point.y <= temymin)
			{
				temymin = point.y;
				yminpoint = point;
			}
			if (point.y >= temymax)
			{
				temymax = point.y;
				ymaxpoint = point;
			}
			if (point.z <= temzmin)
			{
				temzmin = point.z;
				zminpoint = point;
			}
			if (point.z >= temzmax)
			{
				temzmax = point.z;
				zmaxpoint = point;
			}
		}
	}
	float xdis = temxmax - temxmin;
	float ydis = temymax - temymin;
	float zdis = temzmax - temzmin;
	float dis = xdis * xdis + ydis * ydis
		+ zdis * zdis;
	float pointSize = cloudin->size();
	return sqrt(dis) / pow(pointSize, 1.0 / 3);
}

bool isbigger(ske::deltPoint &p1, ske::deltPoint&p2)
{
	if (p1.delta >= p2.delta)
		return true;
	else
		return false;
}

void findMaxDeltaPoint(const vector<ske::deltPoint> &Vec, vector<int> &idvec)
{
	vector<ske::deltPoint> temp;
	temp = Vec;
	sort(temp.begin(), temp.end(), isbigger);
	for (int i=0;i<temp.size();i++)
	{
		idvec.push_back(temp.at(i).id_in_down);
	}
}

void copySrc2PointNei(const pcl::PointCloud<pcl::PointXYZRGB>::Ptr &in,
	pcl::PointCloud<pcl::PointNei>::Ptr &out)
{
	pcl::copyPointCloud(*in, *out);
	for (int i = 0; i < in->size(); i++)
	{
		out->at(i).id_in_src = i;
		out->at(i).id_in_down = -1;
	}
}

void copyPointNei2PCRGB(const pcl::PointCloud<pcl::PointNei>::Ptr &in,
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr &out)
{
	pcl::copyPointCloud(*in, *out);
	for (int i = 0; i < out->size(); i++)
	{
		out->at(i).r = in->at(i).r;
		out->at(i).g = in->at(i).g;
		out->at(i).b = in->at(i).b;
	}
}

void copyDelt2PCNei(const vector<ske::deltPoint> &Vec,
	pcl::PointCloud<pcl::PointNei>::Ptr &out)
{
	out->resize(Vec.size());
	for (int i = 0; i < Vec.size(); i++)
	{
		pcl::PointNei p;
		p = Vec.at(i).p;
		out->at(i) = p;
	}
}

void copyDel2PCRGB(const vector<ske::deltPoint> &Vec,
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr &out)
{
	
	out->resize(Vec.size());
	for (int i = 0; i < Vec.size(); i++)
	{
		pcl::PointXYZRGB p;
		p.x = Vec.at(i).p.x;
		p.y = Vec.at(i).p.y;
		p.z = Vec.at(i).p.z;
		if (Vec.at(i).has_fixed == true)
		{
			p.r = 255;
			p.g = 0;
			p.b = 0;
		}	
		else
		{
			p.r = Vec.at(i).p.r;
			p.g = Vec.at(i).p.g;
			p.b = Vec.at(i).p.b;
		}		
		out->at(i) = p;
	}
}