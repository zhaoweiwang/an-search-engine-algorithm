/********************************************************
�ļ�����		ReadModify.cpp
�����ˣ�		ZhaoweiWang
��ע���ڣ�	2015-7-27
������		��ȡmodify.ini�ļ�����ȡ���ζ�Ӧ����
********************************************************/

#include "Head.h"


extern int mod_fix_total;
extern vector<string> mod_fix_1;
extern int mod_var_total;					//����ҵ���뿼�ǿɱ����Σ�������ҵ��Ҫ�����3��
extern vector<string> mod_var_1;			//ʾ����ֻ������һ��

vector<modifyInfo> modifyFix;
vector<modifyInfo> modifyVar;

/*
��������	ReadModify
���ܣ�	��ȡmodify.ini�ļ�����ȡ���ζ�Ӧ����
���룺	��
�����	��ȡ��Ӧ�������飬�̶��Ϳɱ���������vector
*/
void ReadModify(){

	char * ModifyBuffer;

	FILE * ModifyFile = fopen("modify.ini", "r");
	if (!ModifyFile){
	
		cerr << "�򿪲��� modify.ini �ļ�." << endl;
		exit(1);
	
	}

	fseek(ModifyFile, 0, SEEK_END);
	int sizeModifyFile = ftell(ModifyFile);
	fseek(ModifyFile, 0, SEEK_SET);
	ModifyBuffer = new char[sizeModifyFile+1];
	sizeModifyFile = fread(ModifyBuffer, sizeof(char), sizeModifyFile, ModifyFile);
	fclose(ModifyFile);

	int i = 0;
	string tempName, tempMass;
	modifyInfo tempModify;
	while (i < sizeModifyFile){
	
		tempName = "";
		tempMass = "";

		while (ModifyBuffer[i] != '=') i++;
		while (ModifyBuffer[++i] != '\n') tempName += ModifyBuffer[i];
		if (i >= sizeModifyFile) break;
		while (ModifyBuffer[i] != '=') i++; 
		i++; i++;
		while (ModifyBuffer[i] != 'L') i++;
		i++; i++;
		while (ModifyBuffer[i] != ' ') tempMass += ModifyBuffer[i++];
			
		for (auto i = 0; i < mod_fix_1.size(); i++){
		
			if (tempName == mod_fix_1[i]){
			
				tempModify.modifyName = mod_fix_1[i];
				tempModify.modifyMass = atof(tempMass.c_str());
				tempModify.modifyAminoacid = mod_fix_1[i][mod_fix_1[i].size()-2];
				modifyFix.push_back(tempModify);
				break;
			
			}
		
		}

		for (auto i = 0; i < mod_var_1.size(); i++){

			if (tempName == mod_var_1[i]){

				tempModify.modifyName = mod_var_1[i];
				tempModify.modifyMass = atof(tempMass.c_str());
				tempModify.modifyAminoacid = mod_var_1[i][mod_var_1[i].size() - 2];
				modifyVar.push_back(tempModify);
				break;

			}

		}	
	
	}
	//tempModify.modifyName = "Oxidation[M]";
	//tempModify.modifyMass = -63.998285;
	//tempModify.modifyAminoacid = 'M';
	//modifyVar.push_back(tempModify);

	delete[] ModifyBuffer;

	//for (int i = 0; i < modifyFix.size(); i++){
	//
	//	cout << modifyFix[i].modifyName << " " << modifyFix[i].modifyMass << " " << modifyFix[i].modifyAminoacid << endl;

	//}

	//for (int i = 0; i < modifyVar.size(); i++){

	//	cout << modifyVar[i].modifyName << " " << modifyVar[i].modifyMass << " " << modifyVar[i].modifyAminoacid << endl;

	//}

}
