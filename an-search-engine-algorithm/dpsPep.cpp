/********************************************************
�ļ�����		dpsPep.cpp
�����ˣ�		ZhaoweiWang
��ע���ڣ�	2015-7-27
������		�ݹ�������пɱ������Ķ�
********************************************************/
#include "Head.h"

extern vector<modifyInfo> modifyVar;

vector<recordInfo> allRecord;

extern int mod_var_max_number;				//һ���Ķ��������������ɱ�����λ����

/*
��������	dpsPep
���ܣ�	�ݹ�������пɱ������Ķ�
���룺
�����
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