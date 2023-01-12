#pragma once
#include <vector>
using namespace std;
const int Dim_num = 3;

class DataPoint
{
private:
	int pId; 
	int IdinSrc; 
	float Dimension[Dim_num]; 
	float DimensionWithClusterId[Dim_num + 1]; 
	int clusterId;                   
	bool isKey;                      
	bool visited;                    
	bool iskernel;				
	vector<int> arrivalPoints;    
public:
	DataPoint();                                                   
	DataPoint(int dpID, float* dimension, bool isKey);    
	int GetDpId();                
	void SetDpId(int dpID);        
	int GetPidinSrc();
	void SetPidinSrc(int dpID);
	float* GetDimension(); 
	float* GetDimensionWithClusterId();
	void SetDimension(float* dimension);    
	bool IsKey();                           
	void SetKey(bool isKey);                
	bool isVisited();                      
	void SetVisited(bool visited);           
	bool GetIsKernel();                            
	void SetKernel(bool iskernel);               
	int GetClusterId();                    
	void SetClusterId(int classId);        
	vector<int>& GetArrivalPoints();    
	void DemisionAddClusterId();

};