#pragma once
#include "PCLlibrary.h"
#include "skeleton.h"

void getdownsample(pcl::PointCloud<pcl::PointNei>::Ptr &input,
	pcl::PointCloud<pcl::PointNei>::Ptr &output, float rate);

void getdownsample_tra(pcl::PointCloud<pcl::PointNei>::Ptr &input,
	pcl::PointCloud<pcl::PointNei>::Ptr &output, float rate, float leavsize);

float getInitialRadius(const pcl::PointCloud<pcl::PointXYZRGB>::Ptr &cloudin);

void findMaxDeltaPoint(const vector<ske::deltPoint> &Vec, vector<int> &idvec);


void copySrc2PointNei(const pcl::PointCloud<pcl::PointXYZRGB>::Ptr &in,
	pcl::PointCloud<pcl::PointNei>::Ptr &out);

void copyPointNei2PCRGB(const pcl::PointCloud<pcl::PointNei>::Ptr &in,
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr &out);

void copyDelt2PCNei(const vector<ske::deltPoint> &Vec,
	pcl::PointCloud<pcl::PointNei>::Ptr &out);

void copyDel2PCRGB(const vector<ske::deltPoint> &Vec,
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr &out);
