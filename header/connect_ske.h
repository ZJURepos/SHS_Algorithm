#pragma once
#include "PCLlibrary.h"
#include "skeleton.h"
#include "geometry.h"
#include <vector>
//#include <vtkAutoInit.h>
//VTK_MODULE_INIT(vtkRenderingOpenGL2_AutoInit);
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingFreeType);
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);

struct PCurvature
{
	int index;
	float curvature;
};

void FindDireSke(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &c_in);
void FindDireSke(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &c_in,
	vector<ske::deltPoint>& pts_set_down,
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud_with_normals);

void AddDirection(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &c_in,
	vector<ske::deltPoint>& pts_set,
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud_with_normals);

void ConnectSke(pcl::PointCloud<pcl::PointNei>::Ptr &c_in,
	vector<ske::deltPoint>& pts_set_down,
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud_with_normals,
	vector<pcl::PointCloud<pcl::PointNei>::Ptr> &cs_out);
