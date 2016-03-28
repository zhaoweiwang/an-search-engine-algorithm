/********************************************************
文件名：		Search.cpp
创建人：		Mako Wan
加注日期：	2015-7-23
描述：		搜索引擎的主函数
********************************************************/

#include "Head.h"

int numSpectra = 0;			//谱图数
int numPeptide = 0;			//肽段数

extern vector<spectraInfo> spectraBuffer;

int main(){

	//================================================
	//================================================
	//读取search.param文件 && modify.ini修饰文件（参数详见于InitSearchParam.cpp && ReadModify.cpp 全局变量）
	clock_t readParamS, readParamE;
	double readParamT;
	readParamS = clock();

	InitSearchParam();							//读取.para参数文件
	ReadModify();								//读取.ini修饰文件

	readParamE = clock();
	readParamT = (double)(readParamE - readParamS) / CLOCKS_PER_SEC;
	cout << "Step1. 读取search.param配置文件 && modify.ini修饰文件的时间为:" << setprecision(6) << setiosflags(ios::fixed | ios::showpoint) << readParamT << "秒." << endl;
	//读取search.param文件 && modify.ini修饰文件（参数详见于InitSearchParam.cpp && ReadModify.cpp 全局变量）
	//================================================
	//================================================

	

	//================================================
	//================================================
	//读取spectra.mgf谱图文件（返回谱图数量）
	clock_t readMgfS, readMgfE;
	double readMgfT;
	readMgfS = clock();

	numSpectra = ReadMgfFile();					//读取.mgf文件，每张二级谱取Top-200的谱峰

	readMgfE = clock();
	readMgfT = (double)(readMgfE - readMgfS) / CLOCKS_PER_SEC;
	cout << endl << "Step2. 读取spectra.mgf谱图文件的时间为:" << setprecision(6) << setiosflags(ios::fixed | ios::showpoint) << readMgfT << "秒. && 总共读取 " << numSpectra << " 个谱图." << endl;
	//读取spectra.mgf谱图文件（返回谱图数量）
	//================================================
	//================================================



	//================================================
	//================================================
	//读取DictionaryFile.txt肽段文件
	clock_t readDicS, readDicE;
	double readDicT;
	readDicS = clock();

	numPeptide = ReadDictionaryFile();
	
	readDicE = clock();
	readDicT = (double)(readDicE - readDicS) / CLOCKS_PER_SEC;
	cout << endl << "Step3. 读取DictionaryFile.txt肽段文件的时间为:" << setprecision(6) << setiosflags(ios::fixed | ios::showpoint) << readDicT << "秒. && 总共读取 " << numPeptide << " 个肽段." << endl;
	//读取DictionaryFile.txt肽段文件
	//================================================
	//================================================



	//================================================
	//================================================
	//以肽段为中心进行打分
	clock_t scoreS, scoreE;
	double scoreT;
	scoreS = clock();

	ScoreRank();

	scoreE = clock();
	scoreT = (double)(scoreE - scoreS) / CLOCKS_PER_SEC;
	cout << endl << "Step4. 打分的时间为:" << setprecision(6) << setiosflags(ios::fixed | ios::showpoint) << scoreT << "秒." << endl;
	//以肽段为中心进行打分
	//================================================
	//================================================


	//================================================
	//================================================
	//和pFind的结果进行对比
	clock_t checkS, checkE;
	double checkT;
	checkS = clock();

	checkpFind();

	checkE = clock();
	checkT = (double)(checkE - checkS) / CLOCKS_PER_SEC;
	cout << endl << "Step5. 和pFind对比的时间为:" << setprecision(6) << setiosflags(ios::fixed | ios::showpoint) << checkT << "秒." << endl;
	//和pFind的结果进行对比
	//================================================
	//================================================

	getchar();
	return 0;
}