#pragma once
#include "PCLlibrary.h"
#include "geometry.h"
#include <cmath>
#include "Warning.h"
#include <vector>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>

using namespace std;

namespace pcl
{
	struct EIGEN_ALIGN16 PointNei
	{
		PCL_ADD_POINT4D;                  /
		PCL_ADD_RGB;
		int id_in_src;
		int id_in_down=-1;
		float dj = 0;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
			inline PointNei(const _PointXYZRGB &p)
		{
			x = p.x; y = p.y; z = p.z; data[3] = 1.0f;
			rgb = p.rgb;
		}
		inline PointNei()
		{
			x = y = z = 0.0f;
			data[3] = 1.0f;
			r = g = b = 0;
			a = 255;
		}
		inline PointNei(uint8_t _r, uint8_t _g, uint8_t _b)
		{
			x = y = z = 0.0f;
			data[3] = 1.0f;
			r = _r;
			g = _g;
			b = _b;
			a = 255;
		}
		PointNei & operator=(const PointXYZRGB &p2);
		PointNei & operator=(const PointNei &p2);
	}EIGEN_ALIGN16;
	//POINT_CLOUD_REGISTER_POINT_WRAPPER(PointNei, _PointNei)
}
POINT_CLOUD_REGISTER_POINT_STRUCT(pcl::PointNei,
(float, x, x)
(float, y, y)
(float, z, z)
(float, r, r)
(float, g, g)
(float, b, b)
(int, id_in_src, id_in_src)
(int, id_in_down, id_in_down)
(float, dj, dj)
)
typedef pcl::PointCloud<pcl::PointNei> PointCloudNei;
typedef boost::shared_ptr<pcl::PointCloud<pcl::PointNei>> Ptr;
typedef boost::shared_ptr<const pcl::PointCloud<pcl::PointNei>> ConstPtr;
PCL_INSTANTIATE_PointCloud(pcl::PointNei)

//}

namespace ske
{

	void build_kdtree(const pcl::PointCloud<pcl::PointNei>::Ptr & c_in,
		pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
		pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud);

	class PointID
	{
	public:
		int id_in_src;
		int id_in_down;
		PointID(int id1, int id2)
		{
			this->id_in_src = id1;
			this->id_in_down = id2;
		}
	};

	class deltPoint
	{
	public:
	
		pcl::PointNei p; // current，color pos
		int id_in_src=-1; //nei id;
		int id_in_down=-1; //nei down sample id

		int valid_ornot = 0; 
		int sample_invalid = 0; 
		float fixed_h = 0.25;
		float current_h=0.25; 
		float next_h=0.5; 
		float local_dj;
		float local_dj_src;  
		float delta; // distribution value
		float delta_max;  
		float delta_mid;  
		float delta_min; 
		float delta_max_region; 
		float delta_min_region;  
		Geo::ProPlane direction_plane; 
		Geo::MyVector vector_direction; 
		vector<int> delta_order;
		vector<float> disvec_down; 
		vector<float> disvec_src; 
		
		vector<vector<double>> delta_vec_value;
		
		vector<double> delta_vec_final;
		float ave_dis_src; 
		float ave_dis_down; 
		pcl::PointCloud<pcl::PointNei>::Ptr neibor_point_src; 
		pcl::PointCloud<pcl::PointNei>::Ptr neibor_point_down; 
		vector<PointID> neibor_pts_id_src; 
		vector<PointID> neibor_pts_id_down; 
		
		bool has_moved = false; 
		bool has_smooth = false; 
		bool has_fixed = false; 
		bool is_connection = false; 
		bool is_edge = false; 

	
		bool is_branch = false; 
		bool is_terminal = false; 

		deltPoint() {
			this->id_in_src = -1;
			this->id_in_down = -1;
			this->current_h = 0.25;
			this->next_h = 0.5;
			pcl::PointCloud<pcl::PointNei> * down = new pcl::PointCloud<pcl::PointNei>;
			pcl::PointCloud<pcl::PointNei> * src = new pcl::PointCloud<pcl::PointNei>;
			this->neibor_point_down = boost::shared_ptr<pcl::PointCloud<pcl::PointNei>>(down);
			this->neibor_point_src = boost::shared_ptr<pcl::PointCloud<pcl::PointNei>>(src);
		}
		deltPoint(pcl::PointXYZRGB &p1) {
			this->p = p1;
			this->id_in_src = -1; 
			this->id_in_down = -1; 
			this->current_h = 0.25; 
			this->next_h = 0.5; 
			pcl::PointCloud<pcl::PointNei> * down = new pcl::PointCloud<pcl::PointNei>;
			pcl::PointCloud<pcl::PointNei> * src = new pcl::PointCloud<pcl::PointNei>;
			this->neibor_point_down = boost::shared_ptr<pcl::PointCloud<pcl::PointNei>>(down);
			this->neibor_point_src = boost::shared_ptr<pcl::PointCloud<pcl::PointNei>>(src);
		}

		void FIndNeiborInSrc(const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
			const pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud,
			const pcl::PointCloud<pcl::PointNei>::Ptr &c_in);
		void FIndNeiborInDown(const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
				const pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud,
				const pcl::PointCloud<pcl::PointNei>::Ptr &c_in);
	
		void NeiborSrcDownUpate();

		void UpdateDbsacn();

		void FIndNeiborInSrcFromPlane(const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
			const pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud,
			const pcl::PointCloud<pcl::PointNei>::Ptr &c_in);
		void FIndNeiborInDownFromPlane(const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
			const pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud,
			const pcl::PointCloud<pcl::PointNei>::Ptr &c_in);
	
		int GetIdInSrc();
		int GetIdInDown();

		float compute_distance(const pcl::PointNei &p2);
	
		//inline float computTheta(pcl::PointXYZRGB &p2);
		float comput_theta(const pcl::PointNei &p2);
	
		float local_density(const float h_radius);

		float local_density_src(const float h_radius);
		
		float compute_distribution_value();

		// aij
		inline float compute_aij(const pcl::PointNei &p);
		inline float compute_bii(const pcl::PointNei &p);
		//pcl::PointXYZ &compute_average(const pcl::PointCloud<pcl::PointNei>::Ptr &c_nei,
		//	const pcl::PointCloud<pcl::PointNei>::Ptr &c_in, float &h0);
		pcl::PointXYZ compute_average(const pcl::PointCloud<pcl::PointNei>::Ptr &c_nei,
			const pcl::PointCloud<pcl::PointNei>::Ptr &c_in, float &h0);
		//pcl::PointXYZ &compute_repulsion(const pcl::PointCloud<pcl::PointNei>::Ptr &c_nei,
		//	const pcl::PointCloud<pcl::PointNei>::Ptr &c_in, float &h0);
		pcl::PointXYZ compute_repulsion(const pcl::PointCloud<pcl::PointNei>::Ptr &c_nei,
			const pcl::PointCloud<pcl::PointNei>::Ptr &c_in, float &h0);


	};

	void compute_dv_range(const vector<ske::deltPoint> &pts, float &omax, float &omin);

	void compute_dv_range_region(const vector<ske::deltPoint> &pts, ske::deltPoint &p,
		float &omax, float &omin);

	void findEdgePoint(const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr &tree,
		const pcl::PointCloud<pcl::PointXYZ>::Ptr &treeCloud,
		pcl::PointCloud<pcl::PointNei>::Ptr &c_in,
		vector<ske::deltPoint> &pts_set);
}
