/********************************************************
文件名：		InitSearch.cpp
创建人：		Mako Wan
加注日期：	2015-7-23
描述：		读取search.param文件
********************************************************/

#include "Head.h"

//================================================
//对应search.param参数

string input_path_spectra;
string input_path_index;			//上个作业建立的索引
string input_path_pFind;			//pFind搜索结果
string output_path_result;			//输出文件，有详细的格式要求
int mod_fix_total;
vector<string> mod_fix_1;
int mod_var_total;					//本作业必须考虑可变修饰，但本作业不要求多于3个
vector<string> mod_var_1;			//示例中只放了这一个
int mod_var_max_number;				//一个肽段上几个修饰位点
string tol_precursor_type;          //ppm和Da两种
double tol_precursor_value;	        //正负20ppm，严格大于或小于
string tol_fragment_type;
double tol_fragment_value;			//正负0.5Da，严格大于或小于

//对应search.param参数
//================================================


/*
函数名：	InitSearchParam
功能：	初始化search.param参数
输入：	无
输出：	无
*/
void InitSearchParam(){

	string tempStr;

	ifstream searchFile("search.param", ios::in);
	if (!searchFile.is_open()){

		cerr << "不能打开 search.param." << endl;
		exit(1);

	}

	searchFile >> input_path_spectra;
	auto pos = input_path_spectra.find('=');
	input_path_spectra = input_path_spectra.substr(pos + 1);

	for (int i = 0; i < input_path_spectra.size(); i++){

		if (input_path_spectra[i] == '\\')
		{
			input_path_spectra[i] = '/';
		}

	}

	searchFile >> input_path_index;
	pos = input_path_index.find('=');
	input_path_index = input_path_index.substr(pos+1);

	for (int i = 0; i < input_path_index.size(); i++){
	
		if (input_path_index[i] == '\\')
			input_path_index[i] = '/';
	
	}

	searchFile >> input_path_pFind;
	pos = input_path_pFind.find('=');
	input_path_pFind = input_path_pFind.substr(pos+1);

	for (int i = 0; i < input_path_pFind.size(); i++){
	
		if (input_path_pFind[i] == '\\')
			input_path_pFind[i] = '/';
	
	}

	searchFile >> output_path_result;
	pos = output_path_result.find('=');
	output_path_result = output_path_result.substr(pos+1);

	for (int i = 0; i < output_path_result.size(); i++){
	
		if (output_path_result[i] == '\\')
			output_path_result[i] = '/';
	
	}

	searchFile >> tempStr;
	pos = tempStr.find('=');
	tempStr = tempStr.substr(pos+1);
	mod_fix_total = atoi(tempStr.c_str());

	searchFile >> tempStr;
	pos = tempStr.find('=');
	tempStr = tempStr.substr(pos + 1);
	mod_fix_1.push_back(tempStr);

	searchFile >> tempStr;
	pos = tempStr.find('=');
	tempStr = tempStr.substr(pos+1);
	mod_var_total = atoi(tempStr.c_str());

	searchFile >> tempStr;
	pos = tempStr.find('=');
	tempStr = tempStr.substr(pos + 1);
	mod_var_1.push_back(tempStr);			//TODO：暂时只push_back一个修饰，记得改成for读取

	searchFile >> tempStr;
	pos = tempStr.find('=');
	tempStr = tempStr.substr(pos+1);
	mod_var_max_number = atoi(tempStr.c_str());

	searchFile >> tol_precursor_type;
	pos = tol_precursor_type.find('=');
	tol_precursor_type = tol_precursor_type.substr(pos+1);

	searchFile >> tempStr;
	pos = tempStr.find('=');
	tempStr = tempStr.substr(pos+1);
	tol_precursor_value = atof(tempStr.c_str());

	searchFile >> tol_fragment_type;
	pos = tol_fragment_type.find('=');
	tol_fragment_type = tol_fragment_type.substr(pos+1);

	searchFile >> tempStr;
	pos = tempStr.find('=');
	tempStr = tempStr.substr(pos+1);
	tol_fragment_value = atof(tempStr.c_str());


	searchFile.close();
	//cout << input_path_spectra << endl;
	//cout << input_path_index << endl;
	//cout << input_path_pFind << endl;
	//cout << output_path_result << endl;
	//cout << mod_fix_total << endl;
	//cout << mod_fix_1[0] << endl;
	//cout << mod_var_total << endl;
	//cout << mod_var_1[0] << endl;
	//cout << mod_var_max_number << endl;
	//cout << tol_precursor_type << endl;
	//cout << tol_precursor_value << endl;
	//cout << tol_fragment_type << endl;
	//cout << tol_fragment_value << endl;

}