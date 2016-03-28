/********************************************************
�ļ�����		Search.cpp
�����ˣ�		Mako Wan
��ע���ڣ�	2015-7-23
������		���������������
********************************************************/

#include "Head.h"

int numSpectra = 0;			//��ͼ��
int numPeptide = 0;			//�Ķ���

extern vector<spectraInfo> spectraBuffer;

int main(){

	//================================================
	//================================================
	//��ȡsearch.param�ļ� && modify.ini�����ļ������������InitSearchParam.cpp && ReadModify.cpp ȫ�ֱ�����
	clock_t readParamS, readParamE;
	double readParamT;
	readParamS = clock();

	InitSearchParam();							//��ȡ.para�����ļ�
	ReadModify();								//��ȡ.ini�����ļ�

	readParamE = clock();
	readParamT = (double)(readParamE - readParamS) / CLOCKS_PER_SEC;
	cout << "Step1. ��ȡsearch.param�����ļ� && modify.ini�����ļ���ʱ��Ϊ:" << setprecision(6) << setiosflags(ios::fixed | ios::showpoint) << readParamT << "��." << endl;
	//��ȡsearch.param�ļ� && modify.ini�����ļ������������InitSearchParam.cpp && ReadModify.cpp ȫ�ֱ�����
	//================================================
	//================================================

	

	//================================================
	//================================================
	//��ȡspectra.mgf��ͼ�ļ���������ͼ������
	clock_t readMgfS, readMgfE;
	double readMgfT;
	readMgfS = clock();

	numSpectra = ReadMgfFile();					//��ȡ.mgf�ļ���ÿ�Ŷ�����ȡTop-200���׷�

	readMgfE = clock();
	readMgfT = (double)(readMgfE - readMgfS) / CLOCKS_PER_SEC;
	cout << endl << "Step2. ��ȡspectra.mgf��ͼ�ļ���ʱ��Ϊ:" << setprecision(6) << setiosflags(ios::fixed | ios::showpoint) << readMgfT << "��. && �ܹ���ȡ " << numSpectra << " ����ͼ." << endl;
	//��ȡspectra.mgf��ͼ�ļ���������ͼ������
	//================================================
	//================================================



	//================================================
	//================================================
	//��ȡDictionaryFile.txt�Ķ��ļ�
	clock_t readDicS, readDicE;
	double readDicT;
	readDicS = clock();

	numPeptide = ReadDictionaryFile();
	
	readDicE = clock();
	readDicT = (double)(readDicE - readDicS) / CLOCKS_PER_SEC;
	cout << endl << "Step3. ��ȡDictionaryFile.txt�Ķ��ļ���ʱ��Ϊ:" << setprecision(6) << setiosflags(ios::fixed | ios::showpoint) << readDicT << "��. && �ܹ���ȡ " << numPeptide << " ���Ķ�." << endl;
	//��ȡDictionaryFile.txt�Ķ��ļ�
	//================================================
	//================================================



	//================================================
	//================================================
	//���Ķ�Ϊ���Ľ��д��
	clock_t scoreS, scoreE;
	double scoreT;
	scoreS = clock();

	ScoreRank();

	scoreE = clock();
	scoreT = (double)(scoreE - scoreS) / CLOCKS_PER_SEC;
	cout << endl << "Step4. ��ֵ�ʱ��Ϊ:" << setprecision(6) << setiosflags(ios::fixed | ios::showpoint) << scoreT << "��." << endl;
	//���Ķ�Ϊ���Ľ��д��
	//================================================
	//================================================


	//================================================
	//================================================
	//��pFind�Ľ�����жԱ�
	clock_t checkS, checkE;
	double checkT;
	checkS = clock();

	checkpFind();

	checkE = clock();
	checkT = (double)(checkE - checkS) / CLOCKS_PER_SEC;
	cout << endl << "Step5. ��pFind�Աȵ�ʱ��Ϊ:" << setprecision(6) << setiosflags(ios::fixed | ios::showpoint) << checkT << "��." << endl;
	//��pFind�Ľ�����жԱ�
	//================================================
	//================================================

	getchar();
	return 0;
}