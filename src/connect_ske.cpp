#include "connect_ske.h"
#include "geometry.h"

void FindDireSke(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &c_in)
{
	pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZRGB>);
	tree->setInputCloud(c_in);

	pcl::NormalEstimation<pcl::PointXYZRGB, pcl::Normal> n;
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);

	n.setInputCloud(c_in);
	n.setSearchMethod(tree);
	n.setKSearch(20);
	n.compute(*normals);

	std::vector<PCurvature> tempCV;
	pcl::PointXYZ tempPoint;
	float curvature = 0.0;

	tempPoint.x = tempPoint.y = tempPoint.z = 0.0;
	for (int i = 0; i < normals->size(); i++) {
		struct PCurvature pv;
		pv.index = i;
		pv.curvature = normals->at(i).curvature;
		tempCV.insert(tempCV.end(), pv);
	}
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud_with_normals(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
	pcl::concatenateFields(*c_in, *normals, *cloud_with_normals);


	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D viewer"));
	viewer->setBackgroundColor(0, 0, 0);
	viewer->addCoordinateSystem(0.5, "v1");


	viewer->addPointCloud<pcl::PointXYZRGB>(c_in, "sample cloud");

	viewer->addPointCloudNormals<pcl::PointXYZRGB, pcl::Normal>(c_in, normals, 1, 10, "normals");

	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "sample cloud");

	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100));
	}
	getchar();

}


void FindDireSke(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &c_in, 
	vector<ske::deltPoint>& pts_set_down,
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud_with_normals)
{
	pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZRGB>);
	tree->setInputCloud(c_in);

	pcl::NormalEstimation<pcl::PointXYZRGB, pcl::Normal> n;
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);

	n.setInputCloud(c_in);
	n.setSearchMethod(tree);

	n.setKSearch(20);


	n.compute(*normals);
	
	for (int i = 0; i < pts_set_down.size(); i++)
	{
		if (pts_set_down.at(i).valid_ornot != 0)
			continue;
		int id = pts_set_down.at(i).id_in_down;
		//cout << id << endl;
		vector<double > delta_vec;
		delta_vec = pts_set_down.at(id).delta_vec_value.at(0);
		float tx = delta_vec.at(0);
		float ty = delta_vec.at(1);
		float tz = delta_vec.at(2);
		//float tx = pts_set_down.at(i).delta_order.at(0);
		//float ty = pts_set_down.at(i).delta_order.at(1);
		//float tz = pts_set_down.at(i).delta_order.at(2);
		float all = sqrt(tx*tx + ty * ty + tz * tz);
		tx = tx / all;
		ty = ty / all;
		tz = tz / all;
		normals->at(id).normal_x = tx;
		normals->at(id).normal_y = ty;
		normals->at(id).normal_z = tz;
	}
	pcl::concatenateFields(*c_in, *normals, *cloud_with_normals);
}

void AddDirection(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &c_in,
	vector<ske::deltPoint>& pts_set,
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud_with_normals)
{
	pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZRGB>);
	tree->setInputCloud(c_in);

	pcl::NormalEstimation<pcl::PointXYZRGB, pcl::Normal> n;
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);


	n.setInputCloud(c_in);
	n.setSearchMethod(tree);

	n.setKSearch(20);


	n.compute(*normals);

	for (int i = 0; i < pts_set.size(); i++)
	{
		if (pts_set.at(i).is_edge==false)
		{
			normals->at(i).normal_x = 0;
			normals->at(i).normal_y = 0;
			normals->at(i).normal_z = 0;
			continue;
		}
		//cout << id << endl;
		float tx = pts_set.at(i).vector_direction.x;
		float ty = pts_set.at(i).vector_direction.y;
		float tz = pts_set.at(i).vector_direction.z;
		//float tx = pts_set_down.at(i).delta_order.at(0);
		//float ty = pts_set_down.at(i).delta_order.at(1);
		//float tz = pts_set_down.at(i).delta_order.at(2);
		float all = sqrt(tx*tx + ty * ty + tz * tz);
		tx = tx / all;
		ty = ty / all;
		tz = tz / all;
		normals->at(i).normal_x = tx;
		normals->at(i).normal_y = ty;
		normals->at(i).normal_z = tz;
	}
	pcl::concatenateFields(*c_in, *normals, *cloud_with_normals);
}


void VecNormalize(vector<double> &delta_vec)
{
	float tx = delta_vec.at(0);
	float ty = delta_vec.at(1);
	float tz = delta_vec.at(2);
	float all = sqrt(tx*tx + ty * ty + tz * tz);
	tx = tx / all;
	ty = ty / all;
	tz = tz / all;
	delta_vec.at(0) = tx;
	delta_vec.at(1) = ty;
	delta_vec.at(2) = tz;
}


void growBranch(pcl::PointNei p, int id_c,
	vector<ske::deltPoint>& pts_set_down,
	pcl::PointCloud<pcl::PointNei>::Ptr &neipts,
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud_with_normals,
	pcl::PointCloud<pcl::PointNei>::Ptr &out)
{
	
	vector<double> delta_vec; 
	vector<int> neipts_temp;
	for (int i = 0; i < neipts->size(); i++)
	{
		neipts_temp.push_back(neipts->at(i).id_in_down);
	}
	int ini = 0;
	
	vector<int> cenid(pts_set_down.size() * pts_set_down.size()*0.5, -1);
	for (int i = 0; i < neipts_temp.size(); i++)
		cenid[i] = id_c;
	int currentIDc = p.id_in_down;
	while (ini < neipts_temp.size())
	{
		int id_tempc = cenid[ini];
		int id = neipts_temp.at(ini);
		if (pts_set_down.at(id).is_branch == true || id_tempc == -1 
			|| pts_set_down.at(id).valid_ornot != 0 || pts_set_down.at(id).sample_invalid != 0)
		{
			ini++;
			continue;
		}
		if (pts_set_down.at(id_tempc).valid_ornot != 0 || pts_set_down.at(id_tempc).sample_invalid != 0)
		{
			ini++;
			continue;
		}
		
		vector<double> delta_vec_nei;
		delta_vec_nei = pts_set_down.at(id).delta_vec_value.at(0);
		VecNormalize(delta_vec_nei);
	
		vector<double> delta_vec_cen;
		delta_vec_cen = pts_set_down.at(id_tempc).delta_vec_value.at(0);
		float angle = Geo::ComputeAngleOfVectors(delta_vec_cen, delta_vec_nei);
		ini++;
	
		int currentsize = neipts_temp.size();
		if (angle <= 10)
		{
			pts_set_down.at(id).is_branch = true;
			pcl::PointCloud<pcl::PointNei>::Ptr neipts_nei = pts_set_down.at(id).neibor_point_down;
			out->push_back(pts_set_down.at(id).p);;
#pragma omp parallel for 
			for (int j = 0; j < neipts_nei->size(); j++)
			{
				int id_nei_temp = neipts_nei->at(j).id_in_down;
				if (pts_set_down.at(id_nei_temp).is_branch == true)
					continue;
				else
				{
					cenid[currentsize + j] = id;
					
					neipts_temp.push_back(id_nei_temp);
				}

			}
		}
		else
		{
			continue;
		}
	}
}

void growBranch_id(pcl::PointNei p, int id_c,
	vector<ske::deltPoint>& pts_set_down,
	vector<ske::PointID> &neipts,
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud_with_normals,
	vector<int> &out)
{

	vector<double> delta_vec; /
	vector<int> neipts_temp;
	for (int i = 0; i < neipts.size(); i++)
	{
		neipts_temp.push_back(neipts.at(i).id_in_down);
	}
	int ini = 0;
	
	vector<int> cenid(pts_set_down.size() * 200, -1);
	for (int i = 0; i < neipts_temp.size(); i++)
		cenid[i] = id_c;
	int currentIDc = p.id_in_down;
	while (ini < neipts_temp.size())
	{
		int id_tempc = cenid[ini];
		int id = neipts_temp.at(ini);
		if (pts_set_down.at(id).is_branch == true || id_tempc == -1)
		{
			ini++;
			continue;
		}
		vector<double> delta_vec_nei;
		delta_vec_nei = pts_set_down.at(id).delta_vec_value.at(0);
		VecNormalize(delta_vec_nei);
		
		vector<double> delta_vec_cen;
		delta_vec_cen = pts_set_down.at(id_tempc).delta_vec_value.at(0);
		float angle = Geo::ComputeAngleOfVectors(delta_vec_cen, delta_vec_nei);
		ini++;
		
		int currentsize = neipts_temp.size();
		if (angle <= 15)
		{
			pts_set_down.at(id).is_branch = true;
			
			vector<ske::PointID> neipts_nei_id;
			neipts_nei_id = pts_set_down.at(id).neibor_pts_id_down;
			out.push_back(id);
#pragma omp parallel for 
			for (int j = 0; j < neipts_nei_id.size(); j++)
			{
				int id_nei_temp = neipts_nei_id.at(j).id_in_down;
				if (pts_set_down.at(id_nei_temp).is_branch == true)
					continue;
				else
				{
					cenid[currentsize + j] = id;
					neipts_temp.push_back(id_nei_temp);
				}

			}
			/*		cout << id << endl;*/
		}
		else
		{
			continue;
		}
	}
}


void ConnectSke(pcl::PointCloud<pcl::PointNei>::Ptr &c_in,
	vector<ske::deltPoint>& pts_set_down,
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr &cloud_with_normals,
	vector<pcl::PointCloud<pcl::PointNei>::Ptr> &cs_out)
{
	for (int i = 0; i < pts_set_down.size(); i++)
	{
		if (pts_set_down.at(i).is_branch==true)
			continue;
		pcl::PointCloud<pcl::PointNei>::Ptr temp(new pcl::PointCloud<pcl::PointNei>);
		vector<int> temp_id;
		int id = pts_set_down.at(i).id_in_down;
		temp->push_back(c_in->at(id));
		temp_id.push_back(id);
		auto neipts_id = pts_set_down.at(i).neibor_pts_id_down;
		auto neipts = pts_set_down.at(i).neibor_point_down;
		growBranch(pts_set_down.at(id).p, id, pts_set_down, neipts, cloud_with_normals, temp);
	
		cs_out.push_back(temp);
	}
	cout << "Branches number is: " << cs_out.size() << endl;
}