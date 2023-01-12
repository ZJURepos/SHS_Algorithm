#pragma once
#include "skeleton.h"

class BranchPoint
{
public:
	float x, y, z;
	float r, g, b;
	int id_in_src;
	int id_in_down;
	float delta;
	bool has_moved = false; 
	vector<BranchPoint> nei_in_src; 
	vector<BranchPoint> nei_in_down; 

	inline float compute_dis(BranchPoint &p2);

	BranchPoint& transPoint(const ske::deltPoint &p);

	BranchPoint& operator=(const BranchPoint &p);

};

void KnnSearch(const BranchPoint& p, vector<BranchPoint> &vec,)