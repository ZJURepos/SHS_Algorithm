#include "tools.h"


clock_t time_stamp()
{
	clock_t  clocknow;
	clocknow = clock();
	return clocknow;
}

void ViewerThread::setBackGroud(float r, float g, float b)
{
	this->viewer->setBackgroundColor(r, g, b);
}

void ViewerThread::addPointCLoud()
{
	this->viewer->addPointCloud(this->cloud_show,"p");
}

void ViewerThread::addNormal()
{
	this->viewer->addPointCloudNormals<pcl::PointXYZRGBNormal>(this->cloud_with_normals, 1, this->normal_size*0.5, "normal");
}

void ViewerThread::UpdateNormal()
{
	this->viewer->removeShape("normal");
	this->viewer->addPointCloudNormals<pcl::PointXYZRGBNormal>(this->cloud_with_normals, 1, this->normal_size*0.5, "normal");
}

void ViewerThread::updatePointCloud()
{
	this->viewer->updatePointCloud(this->cloud_show,"p");
}

void ViewerThread::addLine(float p1x, float p1y, float p1z, 
	float p2x, float p2y, float p2z, std::string name)
{
	pcl::PointXYZ p1(p1x,p1y,p1z);
	pcl::PointXYZ p2(p2x,p2y,p2z);
	this->viewer->addLine(p1, p2, name);
}

void ViewerThread::AddMySphere(float x, float y, float z, float radius,std::string name)
{
	pcl::PointXYZ t;
	t.x = x;
	t.y = y;
	t.z = z;
	this->viewer->addSphere(t, radius, name);
}

void ViewerThread::UpdateMySphere(float x, float y, float z, float radius, std::string name)
{
	pcl::PointXYZ t;
	t.x = x;
	t.y = y;
	t.z = z;
	this->viewer->updateSphere(t, radius,0,255,0,name);
}

void ViewerThread::AddPlane(Geo::ProPlane plane_in,float x, float y, float z,std::string name)
{
	pcl::ModelCoefficients coff;
	coff.values.push_back(plane_in.A);
	coff.values.push_back(plane_in.B);
	coff.values.push_back(plane_in.C);
	coff.values.push_back(plane_in.D);
	this->viewer->addPlane(coff,x,y,z,name);
}

void ViewerThread::UpdatePlane(Geo::ProPlane plane_in, float x, float y, float z, std::string name)
{
	this->viewer->removeShape(name);
	pcl::ModelCoefficients coff;
	coff.values.push_back(plane_in.A);
	coff.values.push_back(plane_in.B);
	coff.values.push_back(plane_in.C);
	coff.values.push_back(plane_in.D);
	this->viewer->addPlane(coff, x, y, z, name);
	
}

void ViewerThread::AddOutRGBAPts(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &pts, std::string name)
{
	this->viewer->addPointCloud(pts, name);
}

void ViewerThread::UpdateOutRGBAPts(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &pts, std::string name)
{
	this->viewer->updatePointCloud(pts, name);
}

void ViewerThread::showtext()
{
	for (auto i=this->textvec.begin();i!=this->textvec.end();i++)
	{
		//this->viewer->addText3D(i->first, i->second);
		this->viewer->addText(i->first, i->second.x, i->second.y,0,255,0);
	}
}