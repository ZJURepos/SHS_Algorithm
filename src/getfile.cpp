#include"getfile.h"

void Files::ConfigFile::configed(map<string, string> &param)
{
	this->inpath = param["inpath"];
	this->infile = param["infile"];
	this->outpaths = param["outpaths"];
	this->outpath = param["outpath"];
	this->inname = param["name"];
	this->informat = param["format"];
	this->pre_post = param["pre"]; // if yes,name+i,if not i+name
	this->outformat = param["outformat"];
	this->skeoutpath = param["skeoutpath"];
	istringstream numstr(param["number"]);
	numstr >> this->num;
	istringstream ustr(param["u"]);
	istringstream h_timesstr(param["h_times"]);
	istringstream h_finalsstr(param["h_finals"]);
	istringstream iterationTimesstr(param["iterationTimes"]);
	istringstream ratestr(param["rate"]);
	istringstream leafsizestr(param["leafsize"]);
	istringstream modelstr(param["model"]);
	ustr >> this->u;
	h_timesstr >> this->h_times;
	h_finalsstr >> this->h_finals;
	iterationTimesstr >> this->iterationTimes;
	ratestr >> this->rate;
	leafsizestr >> this->leafsize;
	modelstr >> this->model;
}

int Files::getNumberInString(string &istring, bool &hasnumbr)
{
	int number = 0;
	string filterstring;

	for (int i = istring.size(); i > 0 && istring[i - 1] != '_'&& istring[i - 1] != '/'&&istring[i - 1] != '\\'; i--)
	{
		if (istring[i - 1] >= '0'&&istring[i - 1] <= '9')
			filterstring.insert(filterstring.begin(), istring[i - 1]);
	}
	number = atoi(filterstring.c_str());
	if (number == 0)
		hasnumbr = false;
	return number;
}

template <class T>
int Files::reRangeFileName(vector<T>& files, vector<int> &SerialNumber)
{
	if (files.size() != SerialNumber.size())
	{
		//cout << "The number of Files and Serial Number is wrong" << endl;
		return -1;
	}
	for (int i = 0; i < files.size(); i++)
	{
		for (int j = i; j < files.size(); j++)
		{
			if (SerialNumber[i] >= SerialNumber[j])
			{
				//swap(SerialNumber[i], SerialNumber[j]);

				int tmp;
				tmp = SerialNumber[j];
				SerialNumber[j] = SerialNumber[i];
				SerialNumber[i] = tmp;
				//swap(files[i], files[j]);
				T tmpname;
				tmpname = files[j];
				files[j] = files[i];
				files[i] = tmpname;
			}
		}

	}
	return 0;
}

int Files::GetAllFiles_CertainFormat(string path, vector<string>& files, string format)
{
	intptr_t hFile = 0;
	struct  _finddata_t fileinfo;
	vector<int> SerialNumber;
	bool hasnumber = true;
	int nonum = 0;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("/*" + format).c_str(), &fileinfo)) != -1)
	{
		do
		{
			if ((fileinfo.attrib & _A_SUBDIR))
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, ".") != 0)
				{
					GetAllFiles_CertainFormat((p.assign(path).append("/")).append(fileinfo.name), files, format);

				}
			}
			else
			{
				string serialnumber;
				serialnumber.append(fileinfo.name);
				int num = getNumberInString(serialnumber, hasnumber);
				{
					SerialNumber.push_back(num);

				}
				else
				{
					SerialNumber.push_back(nonum);
				}
				p.assign(path).append("/").append(fileinfo.name);
				files.push_back(p);
			}
			nonum++;

		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
	reRangeFileName(files, SerialNumber);
	//cout << "-------------------------------------------------" << endl;
	return 0;
}

string Files::getpostfixname(string file_name)
{
	string postname;
	for (auto i = file_name.end() - 1; *i != '/'; i--)
	{
		postname.insert(postname.begin(), *i);
	}
	
	return postname;
}

void Files::ConfigFileRead(map<string, string>& m_mapConfigInfo, string inipath)
{
	ifstream configFile;
	string path = inipath;

	configFile.open(path.c_str());
	string str_line;
	if (configFile.is_open())
	{
		while (!configFile.eof())
		{
			getline(configFile, str_line);
			if (str_line.find('#') == 0) {
				continue;
			}
			size_t pos = str_line.find('=');
			string str_key = str_line.substr(0, pos);
			string str_value = str_line.substr(pos + 1);
			m_mapConfigInfo.insert(pair<string, string>(str_key, str_value));
		}
	}
	else
	{
		cout << "Cannot open config file setting.ini, path: ";
		exit(-1);
	}

}

void Files::getFiles(string path, vector<string>& files,
	vector<string>& filenames,string informat)
{
	intptr_t  hFile = 0;
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").append(informat).c_str(), &fileinfo)) != -1)
	{
		do
		{
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				cout << "No Files" << endl;
				continue;
			}
			else
			{
				string filename = p.assign(path).append("\\").append(fileinfo.name);
				cout << filename << endl;
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
				string purename; 
				purename.append(fileinfo.name);
				int dotpos = purename.find_last_of(".");
				filenames.push_back(purename.substr(0, dotpos));
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}

