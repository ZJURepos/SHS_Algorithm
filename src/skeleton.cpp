#include "skeleton.h"
#include "ClusterAnalysis.h"
#include "tools.h"


void ske::build_kdtree(const pcl::PointCloud<pcl::PointNei>::Ptr & c_in,
	pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
	pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud)
{
	copyPointCloud(*c_in, *treeCloud);
	tree->setInputCloud(treeCloud);
}

inline pcl::PointNei & pcl::PointNei::operator=(const pcl::PointXYZRGB &p2)
{
	this->x = p2.x;
	this->y = p2.y;
	this->z = p2.z;
	this->r = p2.r;
	this->g = p2.g;
	this->b = p2.g;
	return *this;
}

inline pcl::PointNei & pcl::PointNei::operator=(const pcl::PointNei &p2)
{
	this->x = p2.x;
	this->y = p2.y;
	this->z = p2.z;
	this->r = p2.r;
	this->g = p2.g;
	this->b = p2.g;
	this->id_in_down = p2.id_in_down;
	this->id_in_src=p2.id_in_src;
	this->dj = p2.dj;
	return *this;
}


void ske::deltPoint::FIndNeiborInSrc(const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
	const pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud,
	const pcl::PointCloud<pcl::PointNei>::Ptr &c_in)
{
	vector<int> kdid;
	vector<float> kddis;
	if (tree->radiusSearch(treeCloud->at(this->id_in_src), this->current_h, kdid, kddis) > 0)
	{
		this->neibor_point_src->clear();
		this->neibor_pts_id_src.clear();
		this->disvec_src.clear();
		for (int j = 1; j < kdid.size(); j++)
		{
			pcl::PointNei p2;
			p2 = c_in->at(kdid.at(j));
			p2.id_in_src = c_in->at(kdid.at(j)).id_in_src;
			if (p2.id_in_src == this->id_in_src)
				continue;
			if(this->compute_distance(p2)<=0.0001)
				continue;
			//this->neibor_point_src->push_back(p2);
			ske::PointID tempID(p2.id_in_src,p2.id_in_down);
			this->neibor_pts_id_src.push_back(tempID);
			this->neibor_point_src->push_back(p2);
			this->disvec_src.push_back(this->compute_distance(p2));
		}
	}

}


void ske::deltPoint::FIndNeiborInDown(const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
	const pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud,
	const pcl::PointCloud<pcl::PointNei>::Ptr &c_down_in)
{
	vector<int> kdid;
	vector<float> kddis;
	if (tree->radiusSearch(treeCloud->at(this->id_in_down), this->current_h, kdid, kddis) > 0)
	{
		//float tempmaxdis = 0;
		this->neibor_point_down->clear();
		this->disvec_down.clear();
		this->neibor_pts_id_down.clear();
		for (int j = 1; j < kdid.size(); j++)
		{
			pcl::PointNei p2;
			p2 = c_down_in->at(kdid.at(j));
			//p2.id_in_down = c_down->at(kdid.at(j)).id_in_down;
			p2.id_in_down = kdid.at(j);
			p2.id_in_src = c_down_in->at(kdid.at(j)).id_in_src;
			if (p2.id_in_src == this->id_in_src || p2.id_in_down==this->id_in_down)
				continue;
			if (this->compute_distance(p2) <= 0.0001)
				continue;
			ske::PointID tempID(p2.id_in_src, p2.id_in_down);
			this->neibor_pts_id_down.push_back(tempID);
			this->neibor_point_down->push_back(p2);
			this->disvec_down.push_back(this->compute_distance(p2));
		}
	}	
}

void ske::deltPoint::UpdateDbsacn()
{
	ClusterAnalysis test;
	test.transPCL2DBpoint(*this, 0.37, 3, false);
	test.preProcess(false);
	test.DoDBscan();
	test.DataAddClusterID();
	vector<int> DBscanNeibor;
	DBscanNeibor = test.GetKernelCluster();
	pcl::PointCloud<pcl::PointNei>::Ptr neibor_point_src
	(new pcl::PointCloud<pcl::PointNei>);
	vector<PointID> tempIDvec;
	for (int i = 0; i < DBscanNeibor.size(); i++)
	{
		neibor_point_src->push_back(this->neibor_point_src->at(DBscanNeibor.at(i)));
		ske::PointID tempID
		(this->neibor_point_src->at(DBscanNeibor.at(i)).id_in_src, 
			this->neibor_point_src->at(DBscanNeibor.at(i)).id_in_down);
		tempIDvec.push_back(tempID);
	}
	this->neibor_pts_id_src.clear();
	this->neibor_point_src->clear();
	this->neibor_point_src = neibor_point_src;
	this->neibor_pts_id_src = tempIDvec;
	
}

void ske::deltPoint::FIndNeiborInSrcFromPlane(const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
	const pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud,
	const pcl::PointCloud<pcl::PointNei>::Ptr &c_in)
{
	vector<int> kdid;
	vector<float> kddis;
	if (tree->radiusSearch(treeCloud->at(this->id_in_src), this->current_h, kdid, kddis) > 0)
	{
		this->neibor_point_src->clear();
		this->neibor_pts_id_src.clear();
		this->disvec_src.clear();
		for (int j = 1; j < kdid.size(); j++)
		{
			pcl::PointNei p2;
			p2 = c_in->at(kdid.at(j));
			p2.id_in_src = c_in->at(kdid.at(j)).id_in_src;
			if (p2.id_in_src == this->id_in_src)
				continue;
		
			pcl::PointNei p2_proj = Geo::ComputeProjectPoint(this->direction_plane, p2);
			float dis = Geo::ComputeDistance(p2, p2_proj);
			if(dis>=this->fixed_h*0.5)
				continue;
			if (this->compute_distance(p2) <= 0.0001)
				continue;
			//this->neibor_point_src->push_back(p2);
			ske::PointID tempID(p2.id_in_src, p2.id_in_down);
			this->neibor_pts_id_src.push_back(tempID);
			this->neibor_point_src->push_back(p2);
			this->disvec_src.push_back(this->compute_distance(p2));
		}
	}
	//this->disvec_src = kddis;
	//maxdis = sqrt(kddis.at(kddis.size() - 1));
}

void ske::deltPoint::FIndNeiborInDownFromPlane(const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
	const pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud,
	const pcl::PointCloud<pcl::PointNei>::Ptr &c_down_in)
{
	vector<int> kdid;
	vector<float> kddis;
	if (tree->radiusSearch(treeCloud->at(this->id_in_down), this->current_h, kdid, kddis) > 0)
	{
		//float tempmaxdis = 0;
		this->neibor_point_down->clear();
		this->disvec_down.clear();
		this->neibor_pts_id_down.clear();
		for (int j = 1; j < kdid.size(); j++)
		{
			pcl::PointNei p2;
			p2 = c_down_in->at(kdid.at(j));
			//p2.id_in_down = c_down->at(kdid.at(j)).id_in_down;
			p2.id_in_down = kdid.at(j);
			p2.id_in_src = c_down_in->at(kdid.at(j)).id_in_src;
			if (p2.id_in_src == this->id_in_src || p2.id_in_down == this->id_in_down)
				continue;
			pcl::PointNei p2_proj = Geo::ComputeProjectPoint(this->direction_plane, p2);
			float dis = Geo::ComputeDistance(p2, p2_proj);
			if (dis >= this->fixed_h*0.5)
				continue;
			if (this->compute_distance(p2) <= 0.0001)
				continue;
			ske::PointID tempID(p2.id_in_src, p2.id_in_down);
			this->neibor_pts_id_down.push_back(tempID);
			this->neibor_point_down->push_back(p2);
			this->disvec_down.push_back(this->compute_distance(p2));

		}
	}
}

int ske::deltPoint::GetIdInSrc()
{
	return this->id_in_src;
}

int ske::deltPoint::GetIdInDown()
{
	return this->id_in_down;
}

void ske::deltPoint::NeiborSrcDownUpate()
{
#pragma omp parallel  for num(8)
	for (int i = 0; i < this->neibor_pts_id_down.size(); i++)
	{
		int id_in_src_temp = this->neibor_pts_id_down.at(i).id_in_src;
		int id_in_down_temp = this->neibor_pts_id_down.at(i).id_in_down;
		for (int j = 0; j < this->neibor_pts_id_src.size(); j++)
		{
			if (id_in_src_temp == this->neibor_pts_id_src.at(j).id_in_src)
			{
				this->neibor_pts_id_src.at(j).id_in_down = id_in_down_temp;
				break;
			}
		}
	}
}


float ske::deltPoint::compute_distance(const pcl::PointNei &p2)
{
	float x = this->p.x - p2.x;
	float y = this->p.y - p2.y;
	float z = this->p.z - p2.z;
	return sqrt(x*x + y * y + z * z);
}


float ske::deltPoint::comput_theta(const pcl::PointNei &p2)
{
	float d = compute_distance(p2);
	float h = this->current_h;

	float out = exp(-(d*d) * 4 * 1.0 / (h*h));
	if (isnan(out))
	{
		cout << p.x << ", " << p.y << ", " << p.z << endl;
		cout << p2.x << ", " << p2.y << ", " << p2.z << endl;
		cout << p.id_in_down << ", " << p2.id_in_down << endl;
		cout << "theta= " << out << endl;
		cout << "d= " << d << ", h= " << h << endl;
		//getchar();
	}
	return exp(-(d*d) * 4 * 1.0 / (h*h));
}


float ske::deltPoint::local_density(const float h_radius)
{
	double dj = 0;
	vector<double>djall;
	for (int i = 0; i < this->neibor_point_src->size(); i++)
	{
		dj += comput_theta(this->neibor_point_src->at(i));
		djall.push_back(dj);
	}
	if (isnan(dj))
	{
		cout << "point down"<<this->id_in_down<<" down_dj " << dj << endl;
		cout << "point src" << this->id_in_src << " down_dj " << dj << endl;
		cout << djall.size() << endl;
		for (int i = 0; i < djall.size(); i++)
			cout << djall.at(i) << endl;
		//getchar();
	}
	this->local_dj = dj + 1;
	return (dj + 1);
}


float ske::deltPoint::local_density_src(const float h_radius)
{
	double dj = 0;
	for (int i = 0; i < this->neibor_point_src->size(); i++)
	{
		dj += comput_theta(this->neibor_point_src->at(i));
	}
	if (isnan(dj))
	{
		cout << "sir_dj " << dj << endl;
		//getchar();
	}
	this->local_dj_src = dj + 1;
	return (dj + 1);
}

Eigen::Matrix3d compute_convariance(ske::deltPoint &p)
{
	float ave_x = 0, ave_y = 0, ave_z = 0;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	copyPointCloud(*p.neibor_point_down, *cloud);
#pragma omp parallel  for
	for (int i = 0; i < cloud->size(); i++)
	{
		ave_x += cloud->at(i).x;
		ave_y += cloud->at(i).y;
		ave_z += cloud->at(i).z;
	}
	ave_x = ave_x / cloud->size();
	ave_y = ave_y / cloud->size();
	ave_z = ave_z / cloud->size();

	double xx11 = 0, xy12 = 0, xz13 = 0;
	double yy22 = 0, yz23 = 0, zz33 = 0;
	int num = cloud->size();
	for (int i = 0; i < cloud->size(); i++)
	{
		double theta = p.comput_theta(p.neibor_point_down->at(i));
		double deltax = ave_x - cloud->at(i).x;
		double deltay = ave_y - cloud->at(i).y;
		double deltaz = ave_z - cloud->at(i).z;
		xx11 += deltax * deltax *theta;
		xy12 += deltax * deltay *theta;
		xz13 += deltax * deltaz *theta;
		yy22 += deltay * deltay*theta;
		yz23 += deltay * deltaz*theta;
		zz33 += deltaz * deltaz*theta;
	}
	/*xx11 = xx11 * 1.0 / num;
	xy12 = xy12 * 1.0 / num;
	xz13 = xz13 * 1.0 / num;
	yy22 = yy22 * 1.0 / num;
	yz23 = yz23 * 1.0 / num;
	zz33 = zz33 * 1.0 / num;*/
	double yx21 = xy12, zx31 = xz13, zy32 = yz23;
	Eigen::Matrix3d tempMatrix;
	tempMatrix << xx11, xy12, xz13,
		yx21, yy22, yz23,
		zx31, zy32, zz33;
	return tempMatrix;
}



void compute_eigens(const Eigen::Matrix3d &matrix_in,
	vector<vector<double>> & eigen_vector, vector<double> & eigen_value)
{
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(matrix_in);
	double lambda1 = es.eigenvalues().data()[0];
	double lambda2 = es.eigenvalues().data()[1];
	double lambda3 = es.eigenvalues().data()[2];
	auto vec1 = es.eigenvectors().col(0);
	auto vec2 = es.eigenvectors().col(1);
	auto vec3 = es.eigenvectors().col(2);
	eigen_value.push_back(lambda1);
	eigen_value.push_back(lambda2);
	eigen_value.push_back(lambda3);
	
	eigen_vector.push_back(vector<double> {vec1(0), vec1(1), vec1(2), lambda1});
	eigen_vector.push_back(vector<double> {vec2(0), vec2(1), vec2(2), lambda2});
	eigen_vector.push_back(vector<double> {vec3(0), vec3(1), vec3(2), lambda3});
}

bool isbigger_vec(vector<double> &vec1, vector<double> &vec2)
{
	if (vec1.at(3) >= vec2.at(3))
		return true;
	else
		return false;
}


float ske::deltPoint::compute_distribution_value()
{
	Eigen::Matrix3d convariance = compute_convariance(*this);
	vector<vector<double>> eigen_vector;
	vector<double> eigen_value;
	compute_eigens(convariance, eigen_vector, eigen_value);
	this->delta_order.clear();
	this->delta_order.push_back(eigen_value[0]);
	this->delta_order.push_back(eigen_value[1]);
	this->delta_order.push_back(eigen_value[2]);

	sort(eigen_vector.begin(), eigen_vector.end(), isbigger_vec);
	this->delta_vec_value = eigen_vector;
	sort(eigen_value.begin(), eigen_value.end());
	this->delta_max = eigen_value[2];
	this->delta_mid = eigen_value[1];
	this->delta_min = eigen_value[0];
	this->delta = eigen_value[2] * 1.0 /
		(eigen_value[2] + eigen_value[1] + eigen_value[0]);
	this->direction_plane = Geo::ProPlane(
		delta_vec_value.at(0).at(0), delta_vec_value.at(0).at(1), delta_vec_value.at(0).at(2),
		this->p.x, this->p.y, this->p.z);
	
	return this->delta;
}

inline float ske::deltPoint::compute_aij(const pcl::PointNei &p2)
{
	float aij_theta = comput_theta(p2);
	float dis = compute_distance(p2);
	return (aij_theta*1.0 / dis);
}

pcl::PointXYZ ske::deltPoint::compute_average(const pcl::PointCloud<pcl::PointNei>::Ptr &c_nei,
	const pcl::PointCloud<pcl::PointNei>::Ptr &c_in, float &h0)
{
	
	double x = 0, y = 0, z = 0;
	pcl::PointXYZ temp2;
	double aij_all = 0;
	for (int i = 0; i < c_nei->size()-1; i++)
	{
		int id = c_nei->at(i).id_in_src;
		if (c_nei->at(i).id_in_src == this->id_in_src)
			continue;
		double dis= this->compute_distance(c_in->at(c_nei->at(i).id_in_src));
		if (dis==0 || dis > h0 || dis<=0.0001)
			continue;
		double aij = this->compute_aij(c_in->at(c_nei->at(i).id_in_src));
		double xt = c_in->at(c_nei->at(i).id_in_src).x*aij / c_in->at(c_nei->at(i).id_in_src).dj;
		double yt = c_in->at(c_nei->at(i).id_in_src).y*aij / c_in->at(c_nei->at(i).id_in_src).dj;
		double zt = c_in->at(c_nei->at(i).id_in_src).z*aij / c_in->at(c_nei->at(i).id_in_src).dj;
		x += c_in->at(c_nei->at(i).id_in_src).x*aij / c_in->at(c_nei->at(i).id_in_src).dj;
		y += c_in->at(c_nei->at(i).id_in_src).y*aij / c_in->at(c_nei->at(i).id_in_src).dj;
		z += c_in->at(c_nei->at(i).id_in_src).z*aij / c_in->at(c_nei->at(i).id_in_src).dj;
		aij_all += aij / c_in->at(c_nei->at(i).id_in_src).dj;

	}
	if (aij_all == 0)
	{
		temp2.x = this->p.x;
		temp2.y = this->p.y;
		temp2.z = this->p.z;
	}
	else
	{
		temp2.x = x * 1.0 / aij_all;
		temp2.y = y * 1.0 / aij_all;
		temp2.z = z * 1.0 / aij_all;
	}
	
	return temp2;

}

inline float ske::deltPoint::compute_bii(const pcl::PointNei &p2)
{
	float aij_theta = comput_theta(p2);
	float dis = compute_distance(p2);
	dis = dis * dis;
	return (aij_theta*1.0 / dis);
}

pcl::PointXYZ ske::deltPoint::compute_repulsion(const pcl::PointCloud<pcl::PointNei>::Ptr &c_nei,
	const pcl::PointCloud<pcl::PointNei>::Ptr &c_in, float &h0)
{
	float x = 0, y = 0, z = 0;
	pcl::PointXYZ temp;
	float bii_all = 0;
	for (int i = 0; i < c_nei->size(); i++)
	{
		if (c_nei->at(i).id_in_src == this->id_in_src)
			continue;
		float dis = this->compute_distance(c_in->at(c_nei->at(i).id_in_src));
		if (dis == 0 || dis >h0)
			continue;
		float bii = this->compute_bii(c_in->at(c_nei->at(i).id_in_src));
		float xtemp= (this->p.x - c_in->at(c_nei->at(i).id_in_src).x)*bii;
		float ytemp = (this->p.y - c_in->at(c_nei->at(i).id_in_src).y)*bii;
		float ztemp = (this->p.z - c_in->at(c_nei->at(i).id_in_src).z)*bii;
		x += xtemp;
		y += ytemp;
		z += ztemp;
		bii_all += bii;
	}
	if (bii_all==0)
	{
		temp.x = 0;
		temp.y = 0;
		temp.z = 0;
	}
	else
	{
		temp.x = x * 1.0 / bii_all;
		temp.y = y * 1.0 / bii_all;
		temp.z = z * 1.0 / bii_all;
	}
	
	return temp;
}


void ske::compute_dv_range(const vector<ske::deltPoint> &pts, float &omax, float &omin)
{
	vector<float> temp;
	for (int i = 0; i < pts.size(); i++)
	{
		if(pts.at(i).valid_ornot==0)
			temp.push_back(pts.at(i).delta);
	}
	sort(temp.begin(), temp.end());
	double out1 = temp.at(0);
	double out2 = temp.at(temp.size() - 1);
	if (out1 >= out2)
	{
		omax = out1;
		omin = out2;
	}
	else
	{
		omax = out2;
		omin = out1;
	}
}


void ske::compute_dv_range_region(const vector<ske::deltPoint> &pts, ske::deltPoint &p,
	float &omax, float &omin)
{
	vector<float> temp;
	if (p.valid_ornot != 0)
		return;
	if (p.has_fixed == true)
		return;

	for (int i = 0; i < p.neibor_point_down->size(); i++)
	{
		int id_in_down_temp = p.neibor_point_down->at(i).id_in_down;
		if (pts.at(id_in_down_temp).valid_ornot == 0)
			temp.push_back(pts.at(id_in_down_temp).delta);
	}

	sort(temp.begin(), temp.end());
	double out1 = temp.at(0);
	double out2 = temp.at(temp.size() - 1);
	if (out1 >= out2)
	{
		omax = out1;
		omin = out2;
	}
	else
	{
		omax = out2;
		omin = out1;
	}
}


void ske::findEdgePoint(const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
	const pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud,
	pcl::PointCloud<pcl::PointNei>::Ptr &c_in,
	vector<ske::deltPoint> &pts_set)
{
	float search_raduis = 0.15;
	vector<Geo::MyVector> vectorout;
#pragma omp parallel  for
	for (int i = 0; i < pts_set.size(); i++)
	{
		vector<int> kdid;
		vector<float> kddis;
		vector<Geo::MyVector> vectorTemp;
		int id_in_src = pts_set.at(i).id_in_src;
		Geo::MyVector tempvec(0, 0, 0);
		pts_set.at(id_in_src).vector_direction = tempvec;
		if (tree->radiusSearch(treeCloud->at(id_in_src), pts_set.at(i).fixed_h*search_raduis, kdid, kddis) > 0)
		{
			if (kdid.size()<3)
				continue;
			for (int j = 1; j < kdid.size(); j++)
			{
				if (c_in->at(kdid.at(j)).id_in_src == id_in_src)
					continue;
				pcl::PointNei p2;
				p2 = c_in->at(kdid.at(j));
				p2.id_in_src = c_in->at(kdid.at(j)).id_in_src;
				Geo::MyVector temp(treeCloud->at(id_in_src), p2);
				temp.myVectorNormalize();
				vectorTemp.push_back(temp);
				tempvec += temp;// *(temp.dis);
			}
			tempvec.averageVec(kdid.size());
			pts_set.at(id_in_src).vector_direction = -tempvec;
			vectorout.push_back(tempvec);
			if (abs(tempvec.x) >= 0.2 || abs(tempvec.y) >= 0.2 || abs(tempvec.z >= 0.2))
			{
				pts_set.at(id_in_src).is_edge = true;

			}
		}
	}
	std::cout << "Edge points: " << vectorout.size() << endl;
}
