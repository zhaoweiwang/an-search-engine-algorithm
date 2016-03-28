/******************************************
文件名：		Head.h
创建人：		Mako Wan
加注日期：	2015-7-23
描述：		放置所需东西
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

	double mass2chargeRatio;	//质荷比
	double abundance;			//丰度

}peaksInfo;

typedef struct record{

	vector<int> a;
	vector<int> b;

}recordInfo;

typedef struct spectra{

	string title = "";				//谱图全称
	int charge;						//母离子电荷
	double pepMass;					//母离子质荷比
	//peaksInfo peaks[1000];		//记录谱峰信息
	vector<peaksInfo> peaks;		//改成容器主要是方便排序取前200强度的峰
	int peaksCount = 0;				//谱峰数量
	string peptide = "";			//得分最高对应的肽段
	double scoreMax = 0.0;			//最高分数
	int top = 0;					//top几

	recordInfo modifyFix;			//记录固定修饰位点和何种修饰
	recordInfo modifyVar;			//记录可变修饰位点和何种修饰

	vector<int> isMatch;

}spectraInfo;

typedef struct peptide{

	string squence = "";
	double mass;

}peptideInfo;

typedef struct modify{

	string modifyName = "";
	char modifyAminoacid;			//发生修饰的氨基酸
	double modifyMass = 0.0;		//修饰质量

}modifyInfo;


typedef struct pFindResult{

	string squence;
	double score;

}pFindResultInfo;

/*
函数名：	InitSearchParam
功能：	初始化search.param参数
输入：	无
输出：	无
*/
void InitSearchParam();


/*
函数名：	ReadMgfFile
功能：	读取谱图数据
输入：	无
输出：	谱图数
*/
int ReadMgfFile();

/*
函数名：	ReadDictionaryFile
功能：	读取酶切后的所有肽段进内存
输入：	无
输出：	肽段个数
*/
int ReadDictionaryFile();

/*
函数名：	ScoreRank
功能：	以肽段为中心进行打分
输入：	无
输出：	查询返回结果及与pFind对比结果
*/
void ScoreRank();


/*
函数名：	ReadModify
功能：	读取modify.ini文件，获取修饰对应质量
输入：	无
输出：	获取对应修饰数组，固定和可变修饰两个vector
*/
void ReadModify();

/*
函数名：	dpsPep
功能：	递归求出所有可变修饰肽段
输入：	
输出：	
*/
void dpsPep(int index, recordInfo record, string pep);

/*
函数名：	checkpFind
功能：	对比pFind结果
输入：
输出：
*/
void checkpFind();

#endif