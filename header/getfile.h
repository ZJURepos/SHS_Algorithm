#pragma once
#include "PCLlibrary.h"

using namespace std;

namespace Files
{
	class ConfigFile
	{
	public:
		string inpath;
		string infile
		string outpaths;
		string outpath;
		string inname;
		string informat;
		string pre_post; 
		string inipath; 
		string outformat;
		string skeoutpath;
		int num;
		float delete_threshold = 0.1;
		float xrange = 0;
		float yrange = 0;
		float zrange = 0;
		float u = 0.35; 
		float h_times = 0.25; 
		float h_finals = 2.0; 
		float iterationTimes = 30;
		float rate = 0.2; 
		float leafsize = 1.0; 
		int model = 2; 
		void configed(map<string, string> &param);
	};

	int getNumberInString(string &istring, bool &hasnumbr);

	template <class T>
	int reRangeFileName(vector<T>& files, vector<int> &SerialNumber);

	int GetAllFiles_CertainFormat(string path, vector<string>& files, string format);

	string getpostfixname(string file_name);

	void ConfigFileRead(map<string, string>& m_mapConfigInfo, string inipath);

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
		else if (para.informat == ".obj")
		{
			pcl::io::loadOBJFile(infile.str(), *c_in);
		}
		else
		{
			cout << "Unsupported Format!" << endl;
			getchar();
			return -1;
		}
		return 0;
	}

	template <class T>
	int Save_Pts(T &c_in, ConfigFile &para, int& i)
	{
		stringstream outfile;
		if (para.pre_post == "yes")
			outfile << para.outpath << "\\" << para.inname << i << para.outformat;
		else
			outfile << para.outpath << "\\" << i << para.inname << para.outformat;
		if (para.outformat == ".ply")
		{
			pcl::io::savePLYFile(outfile.str(), *c_in);
		}
		else if (para.outformat == ".pcd")
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

	template <class T>
	int Load_Pts(T &c_in, ConfigFile &para)
	{
		stringstream infile;
		infile << para.inpath << "\\" << para.inname << para.informat;
		if (para.informat == ".ply")
		{
			pcl::io::loadPLYFile(infile.str(), *c_in);
		}
		else if (para.informat == ".pcd")
		{
			pcl::io::loadPCDFile(infile.str(), *c_in);
		}
		else if (para.informat == ".obj")
		{
			pcl::io::loadOBJFile(infile.str(), *c_in);
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
	int Load_Pts_format(T &c_in, string &infile, string &informat)
	{

		if (informat == ".ply")
		{
			pcl::io::loadPLYFile(infile, *c_in);
		}
		else if (informat == ".pcd")
		{
			pcl::io::loadPCDFile(infile, *c_in);
		}
		else if (informat == ".obj")
		{
			pcl::io::loadOBJFile(infile, *c_in);
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
	int Save_Pts(T &c_in, ConfigFile &para)
	{
		stringstream outfile;

		outfile << para.outpath;//<< "\\" << para.inname << para.informat;

		if (para.outformat == ".ply")
		{
			pcl::io::savePLYFile(outfile.str(), *c_in);
		}
		else if (para.outformat == ".pcd")
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


	template <class T>
	int Save_Pts_format(T &c_in, string &outfile, string &outformat)
	{
		if (outformat == ".ply")
		{
			pcl::io::savePLYFile(outfile, *c_in);
		}
		else if (outformat == ".pcd")
		{
			pcl::io::savePCDFile(outfile, *c_in);
		}

		else
		{
			cout << "Unsupported Format!" << endl;
			getchar();
			return 1;
		}
		return 0;
	}

	void getFiles(string path, vector<string>& files, 
		vector<string>& filenames, string informat);

};
