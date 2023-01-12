#pragma once
#include "PCLlibrary.h"
#include <vector>
#include <math.h>
#define DOT(a,b,c, m,n,q) (a*m+b*n+c*q) 
#define CROSS(dest,a,b,c, m,n,q)\
dest[0]=b*q-c*n;\
dest[1]=c*m-a*q;\
dest[2]=a*n-b*m; 
#define norm_of_vector(a,b,c) (sqrt(a*a+b*b+c*c));

namespace Geo
{
	class Quaternion_s
	{
	public:
		float a;
		float b;
		float c;
		float d;
		Eigen::Affine3f trans_matrix;
		Eigen::Affine3f trans_matrix_reverse;
		void cal_quaternion(float inagle, float axle[3]);
		void toTransMatrix();
	};


	class RotationMatrixs
	{
	public:
		Quaternion_s r1_around_x, r2_around_z, r3_around_y;
	};

	class CircleFit
	{
	public:
		float X, Y;
		float height;
		float radius;
		template <class T>
		inline void fittingCircle(T &cloud);
	};

	template<class T>
	float YoZAngle_2Z(T v_in, float axle[3]);// the angle between input axle with projected vector in plane yoz

	template<class T>
	float YoXAngle_2Z(T v_in, float axle[3]);// the angle between z axls with projected vector in plane yox
	
	template<class T>
	float ZoXAngle_2Z(T v_in, float axle[3]);// the angle between z axls with projected vector in plane yox

	void Compute_Rotation(RotationMatrixs &r_matrixs, pcl::PointXYZ &norp);

	template <class T, class T2>
	void Compute_Range(T &cloud, T2 & para);

	pcl::PointXYZ Compute_center(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &cloud);

	template <class T>//, class T2>
	void DrawCirclePoints(T &out, float x, float z, float y, float radius);

	float ComputeL1Dis(std::vector<double> v);

	float ComputeAngleOfVectors(std::vector<double> v1, std::vector<double> v2);

	class ProPlane
	{
	public:
		float A;
		float B;
		float C;
		float D;
		float A2B2;
		float B2C2;
		float A2C2;
		float A2B2C2;
		ProPlane() {};
		ProPlane(float a, float b, float c, float d);
		ProPlane(float a, float b, float c, float x, float y, float z);
		void GetSumSquare()
		{
			this->A2B2 = this->A*this->A + this->B*this->B;
			this->A2C2 = this->A*this->A + this->C*this->C;
			this->B2C2 = this->C*this->C + this->B*this->B;
			this->A2B2C2 = this->A*this->A + this->B*this->B + this->C*this->C;
		}
	};

	template <class T>
	T ComputeProjectPoint(const ProPlane &plane, T p_in)
	{
		float a2_b2_c2 = plane.A*plane.A + plane.B*plane.B + plane.C*plane.C;
		float a2_b2 = plane.A*plane.A + plane.B*plane.B;
		float a2_c2 = plane.A*plane.A + plane.C*plane.C;
		float b2_c2 = plane.B*plane.B + plane.C*plane.C;
		float x = (b2_c2*p_in.x - plane.A*(plane.B*p_in.y + plane.C*p_in.z + plane.D))
			*1.0 / a2_b2_c2;
		float y = (a2_c2*p_in.y - plane.B*(plane.A*p_in.x + plane.C*p_in.z + plane.D))
			*1.0 / a2_b2_c2;
		float z = (a2_b2*p_in.z - plane.C*(plane.A*p_in.x + plane.B*p_in.y + plane.D))
			*1.0 / a2_b2_c2;
		T temp;
		temp.x = x;
		temp.y = y;
		temp.z = z;
		temp.id_in_src = p_in.id_in_src;
		temp.id_in_down = p_in.id_in_down;
		return temp;
	}

	template <class T1, class T2>
	float ComputeDistance(T1 &p1, T2 &p2)
	{
		float x = p1.x - p2.x;
		float y = p1.y - p2.y;
		float z = p1.z - p2.z;
		float dis = sqrt(x*x + y * y + z * z);
		return dis;
	}

	template <class T>
	bool ComputeInPlaneRange(const ProPlane &plane1, T p1, float d)
	{
		T project_p = ComputeProjectPoint(plane1, p1);
		
	}

	class MyVector
	{
	public:
		pcl::PointXYZ pStart;
		pcl::PointXYZ pEnd;
		float x;
		float y;
		float z;
		float dis;
		MyVector() {};
		MyVector(pcl::PointXYZ p1, pcl::PointXYZ p2):pStart(p1),pEnd(p2)
		{
			this->x = p1.x - p2.x;
			this->y = p1.y - p2.y;
			this->z = p1.z - p2.z;
		}
		template <class T1, class T2>
		MyVector(T1 _p1, T2 _p2)
		{
			pcl::PointXYZ p1,p2;
			p1.x = _p1.x;
			p1.y = _p1.y;
			p1.z = _p1.z;
			p2.x = _p2.x;
			p2.y = _p2.y;
			p2.z = _p2.z;
			this->pStart = p1;
			this->pEnd = p2;
			this->x = p1.x - p2.x;
			this->y = p1.y - p2.y;
			this->z = p1.z - p2.z;
		}
		MyVector(float _x, float _y, float _z) :x(_x),y(_y),z(_z){};
		
		void myVectorNormalize();

		void averageVec(int avesize);

		void operator +=(MyVector p2);

		MyVector operator *(float _times);
		MyVector operator -();
	};

	Eigen::Matrix3d ComputeGeneralConvariance(pcl::PointCloud<pcl::PointXYZ>::Ptr &Pts);

	float ComputeMahalanobisDistance(pcl::PointXYZ &p1, pcl::PointXYZ &p2, Eigen::Matrix3d &Convariance);
}


