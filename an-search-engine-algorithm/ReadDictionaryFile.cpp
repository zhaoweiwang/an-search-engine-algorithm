/********************************************************
文件名：		ReadDictionaryFile.cpp
创建人：		Mako Wan
加注日期：	2015-7-26
描述：		读取DictionaryFile.txt文件，以肽段为中心进行打分
********************************************************/

#include "Head.h"

extern string input_path_index;			//上个作业建立的索引

vector<peptideInfo> peptideBuffer;

/*
函数名：	ReadPeptide
功能：	读取酶切后的所有肽段进内存
输入：	无
输出：	肽段个数
*/
int ReadDictionaryFile(){

	char * DictionaryBuffer;

	int peptideCount = 0;

	FILE * PeptideFile = fopen(input_path_index.c_str(), "r");
	if (!PeptideFile){

		cerr << "打开不了 " << input_path_index << " 文件." << endl;
		exit(1);

	}

	fseek(PeptideFile, 0, SEEK_END);
	int sizePeptideFile = ftell(PeptideFile);
	fseek(PeptideFile, 0, SEEK_SET);
	DictionaryBuffer = new char[sizePeptideFile + 1];
	sizePeptideFile = fread(DictionaryBuffer, sizeof(char), sizePeptideFile, PeptideFile);
	//cout << sizePeptideFile << endl;
	fclose(PeptideFile);

	peptideInfo tempPeptide;
	string s1;
	int i = 0;
	while (i < sizePeptideFile){

		tempPeptide.squence = "";
		s1 = "";

		while (DictionaryBuffer[i] != ' ') tempPeptide.squence += DictionaryBuffer[i++]; i++;
		while (DictionaryBuffer[i] != '\n') s1 += DictionaryBuffer[i++]; i++;
		tempPeptide.mass = atof(s1.c_str());

		peptideBuffer.push_back(tempPeptide);

	}

	//for (auto i = 0; i < peptideBuffer.size(); i++)
	//	cout << peptideBuffer[i].squence << " " << peptideBuffer[i].mass << endl;

	delete[] DictionaryBuffer;

	return peptideBuffer.size();
}
