#pragma once
#include<iostream>
#include <time.h>
#include "PCLlibrary.h"
#include <string>
#include <map>
#include "geometry.h"
clock_t time_stamp();

class ViewerThread
{
public:
	std::string viewername = "default";
	float normal_size=0.01;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_show;
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud_with_normals;
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;
	std::map<std::string, pcl::PointXYZ> textvec;

	ViewerThread()
	{
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr tc(new pcl::PointCloud<pcl::PointXYZRGB>);
		pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr nt(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
		this->cloud_show = boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>>(tc);
		this->cloud_with_normals= boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGBNormal>>(nt);
	}
	ViewerThread(std::string viewername)
	{
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr tc(new pcl::PointCloud<pcl::PointXYZRGB>);
		pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr nt(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
		this->cloud_show = boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>>(tc);
		this->cloud_with_normals = boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGBNormal>>(nt);
		this->viewername = viewername;
		//boost::shared_ptr<pcl::visualization::PCLVisualizer> 
		//	tv(new pcl::visualization::PCLVisualizer(viewername));
		//this->viewer = boost::shared_ptr<pcl::visualization::PCLVisualizer>(tv);
	}

	void setBackGroud(float r, float g, float b);

	void addPointCLoud();

	void updatePointCloud();

	void addNormal();

	void UpdateNormal();

	void addLine(float p1x, float p1y, float p1z, float p2x, float p2y, float p2z, std::string name);

	void AddMySphere(float x, float y, float z, float radius, std::string name);

	void UpdateMySphere(float x, float y, float z, float radius, std::string name);

	void AddPlane(Geo::ProPlane plane_in, float x, float y, float z, std::string name);

	void UpdatePlane(Geo::ProPlane plane_in, float x, float y, float z, std::string name);

	void AddOutRGBAPts(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &pts, std::string name);

	void UpdateOutRGBAPts(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &pts, std::string name);

	template <class T>
	void AddText3D(int id, T &pin)
	{
		std::string idstr = std::to_string(id);
		std::string idstr2 = std::string(1, '*') + std::to_string(id);
		pcl::PointXYZ p;
		p.x = pin.x;
		p.y = pin.y;
		p.z = pin.z;	
		textvec[idstr]= p;
		
		//this->viewer->addText3D(idstr, p, 1.0, 0, 255, 255, idstr2);
	}

	void showtext();


	template<class T>
	void addLine(T p1, T p2, std::string name)
	{
		this->viewer->addLine(p1, p2, name);
	}
};