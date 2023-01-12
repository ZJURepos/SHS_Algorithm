#include "DBscan.h"

DataPoint::DataPoint(){}

DataPoint::DataPoint(int dpID, float* dimension, bool isKey)
{
	this->pId = dpID;
	for (int i=0;i<Dim_num;i++)
	{
		this->Dimension[i] = dimension[i];
	}
	this->isKey = isKey;
}

void DataPoint::SetDimension(float* dimension)
{
	for (int i = 0; i < Dim_num; i++)
	{
		this->Dimension[i] = dimension[i];
	}
}

float* DataPoint::GetDimension()
{
	return this->Dimension;
}

float* DataPoint::GetDimensionWithClusterId()
{
	return this->DimensionWithClusterId;
}

bool DataPoint::IsKey()
{
	return this->isKey;
}

void DataPoint::SetKey(bool isKey)
{
	this->isKey = isKey;
}

bool DataPoint::GetIsKernel()                           
{
	return this->iskernel;
}


void DataPoint::SetKernel(bool iskernel)
{
	this->iskernel = iskernel;
}

int DataPoint::GetDpId()
{
	return this->pId;
}

void DataPoint::SetDpId(int dpID)
{
	this->pId = dpID;
}

int DataPoint::GetPidinSrc()
{
	return this->IdinSrc;
}

void DataPoint::SetPidinSrc(int pid)
{
	this->IdinSrc = pid;
}

bool DataPoint::isVisited()
{
	return  this->visited;
}

void DataPoint::SetVisited(bool visited)
{
	this->visited = visited;
}

int DataPoint::GetClusterId()
{
	return this->clusterId;
}

void DataPoint::SetClusterId(int pid)
{
	this->clusterId = pid;
}

vector<int>& DataPoint::GetArrivalPoints()
{
	return this->arrivalPoints;
}

void DataPoint::DemisionAddClusterId()
{
	for (int i = 0; i < Dim_num; i++)
	{
		this->DimensionWithClusterId[i] = this->Dimension[i];
	}
	this->DimensionWithClusterId[Dim_num] = this->clusterId;
}