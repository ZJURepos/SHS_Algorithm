#include "geometry.h"

namespace Geo
{
	void Quaternion_s::cal_quaternion(float inagle, float axle[3])
	{
		this->a = cos(inagle*1.0 / 2);
		this->b = sin(inagle*1.0 / 2)*axle[0];
		this->c = sin(inagle*1.0 / 2)*axle[1];
		this->d = sin(inagle*1.0 / 2)*axle[2];
		this->toTransMatrix();
	}

	void Quaternion_s::toTransMatrix()
	{
		float a = this->a, b = this->b, c = this->c, d = this->d;
		float a11 = 1 - 2 * c*c - 2 * d*d;
		float a12 = 2 * b*c - 2 * a*d;
		float a13 = 2 * a*c + 2 * b*d;
		float a21 = 2 * b*c + 2 * a*d;
		float a22 = 1 - 2 * b*b - 2 * d*d;
		float a23 = 2 * c*d - 2 * a*b;
		float a31 = 2 * b*d - 2 * a*c;
		float a32 = 2 * a*b + 2 * c*d;
		float a33 = 1 - 2 * b*b - 2 * c*c;
		Eigen::Matrix3f R_trans;
		R_trans << a11, a12, a13,
			a21, a22, a23,
			a31, a32, a33;
		//R_trans.reverse();
		Eigen::Matrix4f trans;
		trans << R_trans(0), R_trans(1), R_trans(2), 0,
			R_trans(3), R_trans(4), R_trans(5), 0,
			R_trans(6), R_trans(7), R_trans(8), 0,
			0, 0, 0, 1;
		Eigen::Transform<float, 3, Eigen::Affine> a3f_transform(trans);
		this->trans_matrix = a3f_transform;

		Eigen::Matrix3f R_trans_reverse;
		R_trans_reverse = R_trans.reverse();
		Eigen::Matrix4f trans_reverse;
		trans_reverse << R_trans_reverse(0), R_trans_reverse(1), R_trans_reverse(2), 0,
			R_trans_reverse(3), R_trans_reverse(4), R_trans_reverse(5), 0,
			R_trans_reverse(6), R_trans_reverse(7), R_trans_reverse(8), 0,
			0, 0, 0, 1;
		Eigen::Transform<float, 3, Eigen::Affine> a3f_transform_reverse(trans_reverse);
		this->trans_matrix_reverse = a3f_transform_reverse;
	}

	template <class T>
	inline void CircleFit::fittingCircle(T &cloud)
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
		C = N * x2 - x1 * x1;
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
		this->radius = sqrt(a*a + b * b - 4 * c) / 2;
	}

	template<class T>
	float YoZAngle_2Z(T v_in, float axle[3])// the angle between input axle with projected vector in plane yoz
	{
		float norm_in = norm_of_vector(0, v_in.y, v_in.z);
		float dot_vec = DOT(0, v_in.y, v_in.z, axle[0], axle[1], axle[2]);
		float angle = acos(dot_vec / norm_in);
		float out;
		if (angle > (M_PI / 2))
		{
			out = M_PI - angle;
		}
		else
			out = angle;
		return out;
	}

	template<class T>
	float YoXAngle_2Z(T v_in, float axle[3])// the angle between z axls with projected vector in plane yox
	{
		float norm_in = norm_of_vector(v_in.x, v_in.y, 0);
		float dot_vec = DOT(v_in.x, v_in.y, 0, axle[0], axle[1], axle[2]);
		float angle = acos(dot_vec / norm_in);
		float out;
		if (angle > (M_PI / 2))
			out = M_PI - angle;
		else
			out = angle;
		return out;
	}

	template<class T>
	float ZoXAngle_2Z(T v_in, float axle[3])// the angle between z axls with projected vector in plane yox
	{
		float norm_in = norm_of_vector(v_in.x, 0, v_in.z);
		float dot_vec = DOT(v_in.x, 0, v_in.z, axle[0], axle[1], axle[2]);
		float angle = acos(dot_vec / norm_in);
		float out;
		if (angle > (M_PI / 2))
			out = M_PI - angle;
		else
			out = angle;
		return out;
	}

	void Compute_Rotation(RotationMatrixs &r_matrixs, pcl::PointXYZ &norp)
	{
		// define coordinate;
		float x_axis[3] = { 1,0,0 };
		float y_axis[3] = { 0,1,0 };
		float z_axis[3] = { 0,0,-1 };
		// compute rotation angle around the X-axis
		float angle_X = YoZAngle_2Z(norp, y_axis);
		// compute rotation angle around the Z-axis
		float angle_Z = YoXAngle_2Z(norp, y_axis);
		// compute rotation angle around the y-axis
		float angle_Y = ZoXAngle_2Z(norp, x_axis);
		r_matrixs.r1_around_x.cal_quaternion(angle_X, x_axis);
		r_matrixs.r2_around_z.cal_quaternion(angle_Z, z_axis);
		r_matrixs.r3_around_y.cal_quaternion(angle_Y, y_axis);
	}

	template <class T, class T2>
	void Compute_Range(T &cloud, T2 & para)
	{
		float xmin = 0, ymin = 0, zmin = 0;
		float xmax = 0, ymax = 0, zmax = 0;
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
		para.xrange = xmax - xmin;
		para.yrange = ymax - ymin;
		para.zrange = zmax - zmin;
	}

	pcl::PointXYZ Compute_center(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &cloud)
	{
		float x = 0, y = 0, z = 0;
		for (int i = 0; i < cloud->size(); i++)
		{
			x += cloud->at(i).x;
			y += cloud->at(i).y;
			z += cloud->at(i).z;
		}
		pcl::PointXYZ p;
		p.x = x * 1.0 / cloud->size();
		p.y = y * 1.0 / cloud->size();
		p.z = z * 1.0 / cloud->size();
		return p;
	}

	template <class T>//, class T2>
	void DrawCirclePoints(T &out, float x, float z, float y, float radius)
		//xoz平面
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

	float ComputeL1Dis(std::vector<double> v)
	{
		float sum= v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
		return sqrt(sum);
	}

	float ComputeAngleOfVectors(std::vector<double> v1, std::vector<double> v2)
	{
		float ab = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
		float aL1 = ComputeL1Dis(v1);
		float bL1 = ComputeL1Dis(v2);
		float cosAngle = ab * 1.0 / (aL1*bL1);
		float angle = acos(cosAngle);
		angle = angle * 180.0 / M_PI;
		if (angle > 90)
			angle = 180 - angle;
		return angle;
	}

	ProPlane::ProPlane(float a, float b, float c, float d)
	{
		this->A = a;
		this->B = b;
		this->C = c;
		this->D = d;
	}
	ProPlane::ProPlane(float a, float b, float c, float x, float y, float z)
	{
		this->A = a;
		this->B = b;
		this->C = c;
		this->D = -(a*x + b * y + c * z);
	}
	// 
	void MyVector::myVectorNormalize()
	{
		float tx = this->x;
		float ty = this->y;
		float tz = this->z;
		this->dis = sqrt(tx*tx + ty * ty + tz * tz);
		this->x = tx / this->dis;
		this->y = ty / this->dis;
		this->z = tz / this->dis;
	}

	//
	void MyVector::averageVec(int avesize)
	{
		this->x = this->x*1.0 / avesize;
		this->y = this->y*1.0 / avesize;
		this->z = this->z*1.0 / avesize;
	}

	void MyVector::operator +=(MyVector p2)
	{
		this->x += p2.x;
		this->y += p2.y;
		this->z += p2.z;
	}

	MyVector MyVector::operator *(float _times)
	{
		this->x = this->x*_times;
		this->y = this->y*_times;
		this->z = this->z*_times;
		return *this;
	}

	MyVector MyVector::operator -()
	{
		this->x = -this->x;
		this->y = -this->y;
		this->z = -this->z;
		return *this;
	}

	Eigen::Matrix3d ComputeGeneralConvariance(pcl::PointCloud<pcl::PointXYZ>::Ptr &Pts)
	{
		float ave_x = 0, ave_y = 0, ave_z = 0;
#pragma omp parallel  for
		for (int i = 0; i < Pts->size(); i++)
		{
			ave_x += Pts->at(i).x;
			ave_y += Pts->at(i).y;
			ave_z += Pts->at(i).z;
		}
		ave_x = ave_x / Pts->size();
		ave_y = ave_y / Pts->size();
		ave_z = ave_z / Pts->size();

		double xx11 = 0, xy12 = 0, xz13 = 0;
		double yy22 = 0, yz23 = 0, zz33 = 0;
		int num = Pts->size();
		for (int i = 0; i < Pts->size(); i++)
		{
			double deltax = ave_x - Pts->at(i).x;
			double deltay = ave_y - Pts->at(i).y;
			double deltaz = ave_z - Pts->at(i).z;
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
		cout << "---------" << endl;

		double yx21 = xy12, zx31 = xz13, zy32 = yz23;
		Eigen::Matrix3d tempMatrix;
		tempMatrix << xx11, xy12, xz13,
			yx21, yy22, yz23,
			zx31, zy32, zz33;
		return tempMatrix;
	}

	float ComputeMahalanobisDistance(pcl::PointXYZ &p1, pcl::PointXYZ &p2, Eigen::Matrix3d &Convariance)
	{
		Eigen::Matrix3d convariance_inverse = Convariance.inverse();
		float x = p1.x - p2.x;
		float y = p1.y - p2.y;
		float z = p1.z - p2.z;
		Eigen::Vector3d diff_temp, out;
		diff_temp << x, y, z;
		out = diff_temp.transpose() * convariance_inverse;
		float outnum = out.dot(diff_temp);
		return sqrt(outnum);
	}
}

