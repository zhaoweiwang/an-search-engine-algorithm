/********************************************************
文件名：		ReadModify.cpp
创建人：		ZhaoweiWang
加注日期：	2015-7-27
描述：		读取modify.ini文件，获取修饰对应质量
********************************************************/

#include "Head.h"


extern int mod_fix_total;
extern vector<string> mod_fix_1;
extern int mod_var_total;					//本作业必须考虑可变修饰，但本作业不要求多于3个
extern vector<string> mod_var_1;			//示例中只放了这一个

vector<modifyInfo> modifyFix;
vector<modifyInfo> modifyVar;

/*
函数名：	ReadModify
功能：	读取modify.ini文件，获取修饰对应质量
输入：	无
输出：	获取对应修饰数组，固定和可变修饰两个vector
*/
void ReadModify(){

	char * ModifyBuffer;

	FILE * ModifyFile = fopen("modify.ini", "r");
	if (!ModifyFile){
	
		cerr << "打开不了 modify.ini 文件." << endl;
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
