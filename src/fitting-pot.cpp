#include "PCLlibrary.h"
#include <iostream>
#include <string>
#include <map>
#include <omp.h>
#include <pcl/ModelCoefficients.h>

using namespace std;
float xrange = 0;
float yrange = 0;
float zrange = 0;
float delete_threshold;
//r-x,g-y轴，b-z；
float xmin = 0, ymin = 0, zmin = 0;
float xmax = 0, ymax = 0, zmax = 0;

//#define Pi 3.141592657;
typedef pcl::PointXYZRGB PointT;
typedef pcl::PointCloud<PointT> PointCloud;
typedef pcl::PointNormal PointNormalT;
typedef pcl::PointCloud<PointNormalT> PointCloudWithNormals;
typedef pcl::PointCloud<pcl::PointXYZRGB>::Ptr PointXYZRGBPtr;
typedef pcl::PointCloud<pcl::PointXYZ>::Ptr PointXYZPtr;


void ConfigFileRead(map<string, string>& m_mapConfigInfo)
{
	ifstream configFile;
	string path = "./toZero.ini";
	configFile.open(path.c_str());
	string str_line;
	if (configFile.is_open())
	{
		while (!configFile.eof())
		{
			getline(configFile, str_line);
			if (str_line.find('#') == 0) 
			{
				continue;
			}
			size_t pos = str_line.find('=');
			string str_key = str_line.substr(0, pos);
			string str_value = str_line.substr(pos + 1);
			m_mapConfigInfo.insert(pair<string, string>(str_key, str_value));
		}
	}
	else
	{
		cout << "Cannot open config file setting.ini, path: ";
		exit(-1);
	}
}

class ConfigFile
{
public:
	string inpath;
	string outpath;
	string inname;
	string informat;
	string pre_post; // if yes,name+i,if not i+name
	int num;
	float delete_threshold;
	void configed(map<string, string> &param);
};

void ConfigFile::configed(map<string, string> &param)
{
	this->inpath = param["inpath"];
	this->outpath = param["outpath"];
	this->inname = param["name"];
	this->informat = param["format"];
	this->pre_post = param["pre"]; // if yes,name+i,if not i+name
	istringstream numstr(param["number"]);
	istringstream delete_thresholdstr(param["delete_threshold"]);
	numstr >> this->num;
	delete_thresholdstr >> this->delete_threshold;
	delete_threshold = this->delete_threshold;
}

template <class T>
int Load_Pts(T &c_in, ConfigFile &para, int& i)
{
	stringstream infile;
	if (para.pre_post == "yes")
		infile << para.inpath << "\\" << para.inname << i << para.informat;
	else
		infile << para.inpath << "\\" << i << para.inname << para.informat;
	if (para.informat == ".ply")
	{
		pcl::io::loadPLYFile(infile.str(), *c_in);
	}
	else if (para.informat == ".pcd")
	{
		pcl::io::loadPCDFile(infile.str(), *c_in);
	}
	else
	{
		cout << "Unsupported Format!" << endl;
		getchar();
		return 1;
	}
	return 0;
}

template <class T>
int Save_Pts(T &c_in, ConfigFile &para, int& i)
{
	stringstream outfile;
	//outfile << parameters.outpath << "\\" << i << parameters.informat;
	if (para.pre_post == "yes")
		outfile << para.outpath << "\\" << para.inname << i << para.informat;
	else
		outfile << para.outpath << "\\" << i << para.inname << para.informat;
	if (para.informat == ".ply")
	{
		pcl::io::savePLYFile(outfile.str(), *c_in);
	}
	else if (para.informat == ".pcd")
	{
		pcl::io::savePCDFile(outfile.str(), *c_in);
	}
	else
	{
		cout << "Unsupported Format!" << endl;
		getchar();
		return 1;
	}
	return 0;
}

class CircleFit
{
public:
	float X, Y;
	float height;
	float radius;
	template <class T>
	inline void fittingCircle(T &cloud)
	{
		float x1 = 0, y1 = 0;
		float x3 = 0, y3 = 0, x2 = 0, y2 = 0, xy = 0, x2y1 = 0, x1y2 = 0;
		int N = cloud->size();
		for (int i = 0; i < cloud->size(); i++)
		{
			float x = cloud->at(i).x; float y = cloud->at(i).z;
			x1 += x;
			y1 += y;
			x3 += x * x*x;
			y3 += y * y*y;
			x2 += x * x;
			y2 += y * y;
			xy += x * y;
			x2y1 += x * x*y;
			x1y2 += x * y*y;
		}


		float C, D, E, G, H;
		float a, b, c;
		C = N * x2 - x1*x1;
		D = N * xy - x1 * y1;
		E = N * x3 + N * x1y2 - (x2 + y2)*x1;
		G = N * y2 - y1 * y1;
		H = N * x2y1 + N * y3 - (x2 + y2)*y1;
		a = (H*D - E * G) / (C*G - D * D);
		b = (H*C - E * D) / (D*D - G * C);
		c = -(a*x1 + b * y1 + x2 + y2) / N;
		//cout << "a: " << a << " b: " << b << " c: " << c << endl;
		this->X = a * 1.0 / (-2);
		this->Y = b * 1.0 / (-2);
		this->radius = sqrt(a*a + b * b - 4 * c)/2;
	}
};

template <class T>
void compute_range(T &cloud)
{
	for (int i = 0; i < cloud->size(); i++)
	{
		if (cloud->at(i).x >= xmax)
			xmax = cloud->at(i).x;
		if (cloud->at(i).y >= ymax)
			ymax = cloud->at(i).y;
		if (cloud->at(i).z >= zmax)
			zmax = cloud->at(i).z;
		if (cloud->at(i).x <= xmin)
			xmin = cloud->at(i).x;
		if (cloud->at(i).y <= ymin)
			ymin = cloud->at(i).y;
		if (cloud->at(i).z <= zmin)
			zmin = cloud->at(i).z;
	}
	xrange = xmax - xmin;
	yrange = ymax - ymin;
	zrange = zmax - zmin;
}

template <class T1, class T2, class T3>
void get_pts(T1 &in, T2 &out, T3 &p_planec, float threshold)
{
#pragma omp parallel  for
	for (int i = 0; i < in->size(); i++)
	{
		float y_value = in->at(i).y;
		if (y_value <= (p_planec + threshold))
			out->push_back(in->at(i));
	}
}

template <class T>//, class T2>
void DrawCirclePoints(T &out, float x, float z, float y, float radius)
{
	float t = 0;
	float angle = (t / 180.0)*M_PI;
	while (t < 360.0)
	{
		pcl::PointXYZRGB p;
		p.x = x + radius * cos(angle);
		p.z = z + radius * sin(angle);
		p.y = y;
		p.r = 255;
		p.g = 0;
		p.b = 0;
		out.push_back(p);
		t = t + 1;
		angle = (t / 180.0)*M_PI;
	}
}

void TransToZeros(pcl::PointCloud<PointXYZRGB>::Ptr &c_in,
	pcl::PointCloud<PointXYZRGB>::Ptr &c_out)
{
	compute_range(c_in);
	vector<CircleFit> circles;
	vector<pcl::PointCloud<PointXYZRGB>::Ptr> every_xyzrgb;
	for (int i = 0; i < 500; i++)
	{
		pcl::PointCloud<PointXYZRGB>::Ptr circle_pts(new pcl::PointCloud<PointXYZRGB>);
		get_pts(c_in, circle_pts, ymin, i*yrange*1.0 / 1000);
		if (circle_pts->size() < 20)
			continue;
		CircleFit centerP;
		centerP.height = ymin + (i * yrange*1.0 / 1000) *1.0 / 2;
		centerP.fittingCircle(circle_pts);
		every_xyzrgb.push_back(circle_pts);
		circles.push_back(centerP);
	}
	pcl::PointCloud<PointXYZRGB>::Ptr centerpts(new pcl::PointCloud<PointXYZRGB>);
	for (int i = 0; i < circles.size(); i++)
	{
		pcl::PointXYZRGB p;
		p.x = circles.at(i).X;
		p.y = circles.at(i).height;
		p.z = circles.at(i).Y;
		p.r = 255;
		centerpts->push_back(p);
	}

	vector<pcl::PointCloud<PointXYZ>::Ptr> every_xyz;
	cout << every_xyzrgb.at(0)->size() << endl;
	for (int i = 0; i < every_xyzrgb.size(); i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr temp(new pcl::PointCloud<pcl::PointXYZ>);

		pcl::copyPointCloud(*every_xyzrgb.at(i), *temp);
		every_xyz.push_back(temp);
	}
	Eigen::Matrix4f trans;
	trans << 1, 0, 0, -circles.at(0).X,
		0, 1, 0, -circles.at(0).heght,
		0, 0, 1, -circles.at(0).Y,
		0, 0, 0, 1;
	Eigen::Transform<float, 3, Eigen::Affine> a3f_transform(trans);
	pcl::transformPointCloud(*c_in, *c_out, trans);
}

int main()
{
	cout << "Ready.....";
	getchar();
	cout << "Adjust coordinate" << endl;
	map<string, string> param;
	ConfigFileRead(param);
	ConfigFile parameters;
	parameters.configed(param);
	
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr c_in(new pcl::PointCloud<pcl::PointXYZRGB>);
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("result"));
	viewer->setBackgroundColor(1, 1, 1);
	viewer->addCoordinateSystem(0.2);// , p.x, p.y, p.z);
	viewer->addPointCloud(c_in, "1");
	for (int i=1;i<parameters.num;i++)
	{
		int error = Load_Pts(c_in, parameters, i);
		if (error == 1)
			break;
		pcl::PointCloud<PointXYZRGB>::Ptr c_out(new pcl::PointCloud<PointXYZRGB>);
		TransToZeros(c_in, c_out);
		error = Save_Pts(c_out, parameters, i);
		if (error == 1)
			break;
		viewer->updatePointCloud(c_out, "1");
	}

	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
	

	getchar();
	return 0;

}

