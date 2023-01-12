#include "Branches.h"

inline float BranchPoint::compute_dis(BranchPoint &p2)
{
	float x = this->x - p2.x;
	float y = this->y - p2.y;
	float z = this->z - p2.z;
	return (sqrt(x*x + y * y + z * z));
}


BranchPoint& BranchPoint::transPoint(const ske::deltPoint &p)
{
	this->x = p.p.x;
	this->y = p.p.y;
	this->z = p.p.z;
	this->id_in_src = p.id_in_src;
	this->id_in_down = p.id_in_down;
	this->delta = p.delta;
}

BranchPoint& BranchPoint::operator=(const BranchPoint &p)
{
	this->x = p.x;
	this->y = p.y;
	this->z = p.z;
	this->id_in_src = p.id_in_src;
	this->id_in_down = p.id_in_down;
	this->delta = p.delta;	
	this->has_moved = p.has_moved;
}

