/********************************************************
�ļ�����		InitSearch.cpp
�����ˣ�		Mako Wan
��ע���ڣ�	2015-7-23
������		��ȡsearch.param�ļ�
********************************************************/

#include "Head.h"

//================================================
//��Ӧsearch.param����

string input_path_spectra;
string input_path_index;			//�ϸ���ҵ����������
string input_path_pFind;			//pFind�������
string output_path_result;			//����ļ�������ϸ�ĸ�ʽҪ��
int mod_fix_total;
vector<string> mod_fix_1;
int mod_var_total;					//����ҵ���뿼�ǿɱ����Σ�������ҵ��Ҫ�����3��
vector<string> mod_var_1;			//ʾ����ֻ������һ��
int mod_var_max_number;				//һ���Ķ��ϼ�������λ��
string tol_precursor_type;          //ppm��Da����
double tol_precursor_value;	        //����20ppm���ϸ���ڻ�С��
string tol_fragment_type;
double tol_fragment_value;			//����0.5Da���ϸ���ڻ�С��

//��Ӧsearch.param����
//================================================


/*
��������	InitSearchParam
���ܣ�	��ʼ��search.param����
���룺	��
�����	��
*/
void InitSearchParam(){

	string tempStr;

	ifstream searchFile("search.param", ios::in);
	if (!searchFile.is_open()){

		cerr << "���ܴ� search.param." << endl;
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
	mod_var_1.push_back(tempStr);			//TODO����ʱֻpush_backһ�����Σ��ǵøĳ�for��ȡ

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