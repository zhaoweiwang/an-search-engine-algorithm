/********************************************************
文件名：		ReadMgfFile.cpp
创建人：		Mako Wan
加注日期：	2015-7-23
描述：		读取spectra.mgf文件
********************************************************/

#include "Head.h"

extern string input_path_spectra;

vector<spectraInfo> spectraBuffer;

bool cmpAbun(const peaksInfo &a, const peaksInfo &b){

	return a.abundance > b.abundance;

}

bool cmpM2C(const peaksInfo &a, const peaksInfo &b){

	return a.mass2chargeRatio < b.mass2chargeRatio;

}

/*
函数名：	ReadMgfFile
功能：	读取谱图数据
输入：	无
输出：	谱图数
*/
int ReadMgfFile(){

	char * MgfBuffer;

	int spectraCount = 0;

	FILE * MgfFile = fopen(input_path_spectra.c_str(), "r");
	if (!MgfFile){
	
		cerr << "打开不了 " << input_path_spectra << " 文件." << endl;
		exit(1);
	
	}

	fseek(MgfFile, 0, SEEK_END);
	int sizeMgfFile = ftell(MgfFile);
	fseek(MgfFile, 0, SEEK_SET);
	MgfBuffer = new char[sizeMgfFile+1];
	sizeMgfFile = fread(MgfBuffer, sizeof(char), sizeMgfFile, MgfFile);
	fclose(MgfFile);
	//cout << sizeMgfFile << endl;

	spectraInfo tempSpectra;
	string tempPepMass, s1, s2;
	int i = 0;
	peaksInfo tempPeak;
	while (i < sizeMgfFile){
	
		while (MgfBuffer[i] != 'T'){

			i++;

		}
		if (i >= sizeMgfFile) break;

		spectraCount++;

		tempSpectra.title = "";
		while (MgfBuffer[i] != '\n') tempSpectra.title += MgfBuffer[i++];
		while (MgfBuffer[i] != '+') i++;
		tempSpectra.charge = MgfBuffer[i-1]-'0';
		while (MgfBuffer[i] != '=') i++;
		tempPepMass = ""; i++;
		while (MgfBuffer[i] != '\n') tempPepMass += MgfBuffer[i++];
		tempSpectra.pepMass = atof(tempPepMass.c_str());

		//tempSpectra.peaksCount = 0;
		tempSpectra.peaks.clear();
		tempSpectra.isMatch.clear();
		while (MgfBuffer[++i] != 'E'){

			if (MgfBuffer[i] == '\n') continue; //有些mgf中的END IONS与谱峰信息中间有回车行，需要过掉

			s1 = ""; s2 = "";
			while (MgfBuffer[i] != ' ' && MgfBuffer[i] != '\t'){
			
				s1 += MgfBuffer[i++];
			
			}
			i++;
			while (MgfBuffer[i] != '\n'){

				s2 += MgfBuffer[i++];

			}

			//cout << s1 << " " << s2 << endl;
			tempPeak.mass2chargeRatio = atof(s1.c_str());
			tempPeak.abundance = atof(s2.c_str());
			//tempSpectra.peaks[tempSpectra.peaksCount].mass2chargeRatio = atof(s1.c_str());
			//tempSpectra.peaks[tempSpectra.peaksCount++].abundance = atof(s2.c_str());
			tempSpectra.peaks.push_back(tempPeak);
		
		}
		//cout << "haha" << endl;

		//保留强度前200的峰
		sort(tempSpectra.peaks.begin(), tempSpectra.peaks.end(), cmpAbun);
		
		if (tempSpectra.peaks.size() > 200){
		
			int count = 0;
			for (int t = 199; t < tempSpectra.peaks.size(); t++){
			
				if (tempSpectra.peaks[t + 1].abundance == tempSpectra.peaks[t].abundance)	count++;
			
			}
			tempSpectra.peaks.erase(tempSpectra.peaks.begin()+199+count, tempSpectra.peaks.end());
		
		}
		//if (tempSpectra.title == "TITLE=cncp2012.4174.4174.4.dta")
		//	for (int z = 0; z < tempSpectra.peaks.size(); z++)
		//		cout << tempSpectra.peaks[z].mass2chargeRatio << " " << tempSpectra.peaks[z].abundance << endl;

		//tempSpectra.isMatch.resize(tempSpectra.peaks.size(), 0);

		//保留强度前200的峰
		//sort(tempSpectra.peaks.begin(), tempSpectra.peaks.end(), cmpM2C);

		spectraBuffer.push_back(tempSpectra);
		
		
	
	}

	//for (int i = 0; i < spectraBuffer.size(); i++){
	//
	//	cout << spectraBuffer[i].title << " " << spectraBuffer[i].charge << " " << spectraBuffer[i].pepMass << endl;
	//	for (auto j = 0; j < spectraBuffer[i].peaksCount; j++)
	//		cout << spectraBuffer[i].peaks[j].mass2chargeRatio << " " << spectraBuffer[i].peaks[j].abundance << endl;
	//
	//}

	//cout << "总共读取 " << spectraCount << " 个谱图." << endl;

	delete[] MgfBuffer;

	return spectraCount;
}