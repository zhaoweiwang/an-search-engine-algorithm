/********************************************************
文件名：		dpsPep.cpp
创建人：		ZhaoweiWang
加注日期：	2015-7-27
描述：		递归求出所有可变修饰肽段
********************************************************/
#include "Head.h"

extern vector<modifyInfo> modifyVar;

vector<recordInfo> allRecord;

extern int mod_var_max_number;				//一个肽段上允许发生的最大可变修饰位点数

/*
函数名：	dpsPep
功能：	递归求出所有可变修饰肽段
输入：
输出：
*/

void dpsPep(int index, recordInfo record, string pep){

	if (index == pep.size()){

		allRecord.push_back(record);
		return;

	}
		
	for (int i = 0; i < modifyVar.size(); i++){
	
		if (pep[index] == modifyVar[i].modifyAminoacid){
			
			//cur += pep[index];

			dpsPep(index+1, record, pep);
			
			if (record.a.size() == mod_var_max_number)
				return;

			record.a.push_back(index);
			record.b.push_back(i);
			dpsPep(index+1, record, pep);
			return;
		}
	
	}

	//cur += pep[index];

	dpsPep(index+1, record, pep);


}