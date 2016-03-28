/********************************************************
�ļ�����		ReadDictionaryFile.cpp
�����ˣ�		Mako Wan
��ע���ڣ�	2015-7-26
������		��ȡDictionaryFile.txt�ļ������Ķ�Ϊ���Ľ��д��
********************************************************/

#include "Head.h"

extern string input_path_index;			//�ϸ���ҵ����������

vector<peptideInfo> peptideBuffer;

/*
��������	ReadPeptide
���ܣ�	��ȡø�к�������Ķν��ڴ�
���룺	��
�����	�Ķθ���
*/
int ReadDictionaryFile(){

	char * DictionaryBuffer;

	int peptideCount = 0;

	FILE * PeptideFile = fopen(input_path_index.c_str(), "r");
	if (!PeptideFile){

		cerr << "�򿪲��� " << input_path_index << " �ļ�." << endl;
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
