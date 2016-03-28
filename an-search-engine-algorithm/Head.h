/******************************************
�ļ�����		Head.h
�����ˣ�		Mako Wan
��ע���ڣ�	2015-7-23
������		�������趫��
******************************************/

#ifndef IHEAD_H
#define IHEAD_H

#include <iostream>
#include <string>
#include <vector>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <algorithm>

using namespace std;

typedef struct peak{

	double mass2chargeRatio;	//�ʺɱ�
	double abundance;			//���

}peaksInfo;

typedef struct record{

	vector<int> a;
	vector<int> b;

}recordInfo;

typedef struct spectra{

	string title = "";				//��ͼȫ��
	int charge;						//ĸ���ӵ��
	double pepMass;					//ĸ�����ʺɱ�
	//peaksInfo peaks[1000];		//��¼�׷���Ϣ
	vector<peaksInfo> peaks;		//�ĳ�������Ҫ�Ƿ�������ȡǰ200ǿ�ȵķ�
	int peaksCount = 0;				//�׷�����
	string peptide = "";			//�÷���߶�Ӧ���Ķ�
	double scoreMax = 0.0;			//��߷���
	int top = 0;					//top��

	recordInfo modifyFix;			//��¼�̶�����λ��ͺ�������
	recordInfo modifyVar;			//��¼�ɱ�����λ��ͺ�������

	vector<int> isMatch;

}spectraInfo;

typedef struct peptide{

	string squence = "";
	double mass;

}peptideInfo;

typedef struct modify{

	string modifyName = "";
	char modifyAminoacid;			//�������εİ�����
	double modifyMass = 0.0;		//��������

}modifyInfo;


typedef struct pFindResult{

	string squence;
	double score;

}pFindResultInfo;

/*
��������	InitSearchParam
���ܣ�	��ʼ��search.param����
���룺	��
�����	��
*/
void InitSearchParam();


/*
��������	ReadMgfFile
���ܣ�	��ȡ��ͼ����
���룺	��
�����	��ͼ��
*/
int ReadMgfFile();

/*
��������	ReadDictionaryFile
���ܣ�	��ȡø�к�������Ķν��ڴ�
���룺	��
�����	�Ķθ���
*/
int ReadDictionaryFile();

/*
��������	ScoreRank
���ܣ�	���Ķ�Ϊ���Ľ��д��
���룺	��
�����	��ѯ���ؽ������pFind�ԱȽ��
*/
void ScoreRank();


/*
��������	ReadModify
���ܣ�	��ȡmodify.ini�ļ�����ȡ���ζ�Ӧ����
���룺	��
�����	��ȡ��Ӧ�������飬�̶��Ϳɱ���������vector
*/
void ReadModify();

/*
��������	dpsPep
���ܣ�	�ݹ�������пɱ������Ķ�
���룺	
�����	
*/
void dpsPep(int index, recordInfo record, string pep);

/*
��������	checkpFind
���ܣ�	�Ա�pFind���
���룺
�����
*/
void checkpFind();

#endif