/********************************************************
文件名：		ScoreRank.cpp
创建人：		Mako Wan
加注日期：	2015-7-26
描述：		读取DictionaryFile.txt文件，以肽段为中心进行打分
********************************************************/

#include "Head.h"

extern vector<spectraInfo> spectraBuffer;
extern vector<peptideInfo> peptideBuffer;

extern vector<recordInfo> allRecord;

extern vector<modifyInfo> modifyFix;
extern vector<modifyInfo> modifyVar;


extern string tol_precursor_type;           //ppm和Da两种
extern double tol_precursor_value;	        //正负20ppm，严格大于或小于
extern string tol_fragment_type;
extern double tol_fragment_value;			//正负0.5Da，严格大于或小于

extern string output_path_result;

const double pMass = 1.00727647012;			//质子质量

vector<double> MassTable = { 71.03711, 166.99836, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891, 113.08406, 181.01401, 128.09496, 113.08406, 131.04048, 114.04293, 243.02965, 97.05276, 128.05858, 156.10111, 87.03203, 101.04768, 0.00000, 99.06841, 186.07931, 113.08406, 163.06332, 128.55059 };
vector<double> MassTableTemp;


double detaPrecursor = 0.0;
double detaFragment = 0.0;
double gama = -0.9;
double theta = 0.5;


//ofstream FragmentMassFile("FragmentMassFile.txt");

//针对可变修饰进行打分
void ScoreKSDP_Var(peptideInfo pep, int k, int z, double addMass){

	double mzIon[4][200];				//200对应肽段的最长长度
	double tempMass = 0.0;
	
	//for (int k = 0; k < allRecord[z].a.size(); k++){

	//	pep.mass += modifyVar[allRecord[z].b[k]].modifyMass;

	//}
	pep.mass += addMass;
	
	spectraBuffer[k].isMatch.assign(spectraBuffer[k].peaks.size(), 0);

	//如果母离子的电荷<=2,只考虑计算b+、y+离子，不考虑b++、y++离子
	if (spectraBuffer[k].charge <= 2){

		//求出肽段的b、y离子
		int len = pep.squence.size() - 1;
		for (int i = 0; i < len; i++){

			//for (int j = 0; j < allRecord[z].a.size(); j++){
			//
			//	if (i == allRecord[z].a[j])
			//		tempMass += (MassTableTemp[pep.squence[i] - 'A'] + modifyVar[allRecord[z].b[j]].modifyMass);
			//	else
			//		tempMass += MassTableTemp[pep.squence[i] - 'A'];
			//
			//}
			tempMass += MassTableTemp[pep.squence[i] - 'A'];
			for (int j = 0; j < allRecord[z].a.size(); j++){

				if (i == allRecord[z].a[j])
					tempMass += modifyVar[allRecord[z].b[j]].modifyMass;

			}

			mzIon[0][i] = tempMass + pMass;
			//mzIon[1][i] = (mzIon[0][i] + pMass) / 2;
			mzIon[1][len - i - 1] = pep.mass - tempMass + pMass;
			//mzIon[3][len - i - 1] = (mzIon[2][len - i - 1] + pMass) / 2;

		}



		//FragmentMassFile << pep.squence << "\t" << pep.mass << endl;
		//for (int i = 0; i < 4; i++){
		//	for (int j = 0; j < len; j++)
		//		FragmentMassFile << mzIon[i][j] << "\t";
		//	FragmentMassFile << endl;
		//}


		//初始化和计算得分
		int speIon[2][200], pepIon[2][200];		//实验谱图和理论谱图的b、y离子
		double sumAbundance = 0.0;
		double maxAbundance = 0.0;
		int l1 = 2, l2 = 2;							//对应窗口5时的l1和l2参数

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < len + l2 + 1; j++){

				speIon[i][j] = 0;
				pepIon[i][j] = 0;

			}

		if (tol_fragment_type == "ppm"){

			for (int i = 0; i < 2; i++){

				for (int j = 0; j < len; j++){

					pepIon[i][j] = 1;

					for (int l = 0; l < spectraBuffer[k].peaks.size(); l++){

						if (maxAbundance < spectraBuffer[k].peaks[l].abundance)
							maxAbundance = spectraBuffer[k].peaks[l].abundance;

						//if (abs( (int)(spectraBuffer[k].peaks[l].mass2chargeRatio * 10000) - (int)(mzIon[i][j] * 10000)) < (detaFragment*(int)(mzIon[i][j] * 10000)) ){
						//if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < (detaFragment*spectraBuffer[k].peaks[l].mass2chargeRatio)){
						//if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < (detaFragment*1000)){
						if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < detaFragment*mzIon[i][j]){

							if (spectraBuffer[k].isMatch[l] == 0){
							
								sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
								spectraBuffer[k].isMatch[l] == 1;
							
							}
							//cout << spectraBuffer[k].peaks[l].mass2chargeRatio << endl;
							speIon[i][j] = 1;
							break;

						}

					}

				}

			}

		}
		else if (tol_fragment_type == "Da"){

			for (int i = 0; i < 4; i++){

				for (int j = 0; j < len; j++){

					pepIon[i][j] = 1;

					for (int l = 0; l < spectraBuffer[k].peaks.size(); l++){

						if (maxAbundance < spectraBuffer[k].peaks[l].abundance)
							maxAbundance = spectraBuffer[k].peaks[l].abundance;

						if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < detaFragment){

							//sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);

							if (spectraBuffer[k].isMatch[l] == 0){

								sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
								spectraBuffer[k].isMatch[l] == 1;

							}
							speIon[i][j] = 1;
							break;

						}

					}

				}

			}

		}

		double tempScore0 = 0.0;

		double temp1 = 0.0;

		for (int i = 0; i < 2; i++){

			int temp2 = 0;

			for (int j = 0; j <= l2; j++)
				temp2 += ((pepIon[i][j] - speIon[i][j]) * (pepIon[i][j] - speIon[i][j]));
			temp1 += exp(gama*temp2);

			for (int j = l2 + 1; j <= l1 + l2; j++){

				temp2 += (speIon[i][j] - pepIon[i][j])*(speIon[i][j] - pepIon[i][j]);
				temp1 += exp(gama*temp2);

			}

			for (int j = l1 + 1; j < len; j++){

				temp2 = temp2 + (speIon[i][j + l2] - pepIon[i][j + l2])*(speIon[i][j + l2] - pepIon[i][j + l2]) - (speIon[i][j - l1 - 1] - pepIon[i][j - l1 - 1])*(speIon[i][j - l1 - 1] - pepIon[i][j - l1 - 1]);
				temp1 += exp(gama*temp2);

			}

		}

		//tempScore0 = sumAbundance * pow(temp1, theta) / (len + 1);
		tempScore0 = (sumAbundance * 100 / sqrt(maxAbundance)) * pow(temp1, theta) / (len << 1);
		//cout << "分数:" << tempScore0 << endl;

		if (tempScore0 > spectraBuffer[k].scoreMax){

			spectraBuffer[k].scoreMax = tempScore0;
			spectraBuffer[k].peptide = pep.squence;
			//待记录修饰信息
			//加固定修饰信息
			spectraBuffer[k].modifyFix.a.clear();
			spectraBuffer[k].modifyFix.b.clear();

			for (int i = 0; i < modifyFix.size(); i++)
				for (int j = 0; j < pep.squence.size(); j++)
					if (pep.squence[j] == modifyFix[i].modifyAminoacid){

						spectraBuffer[k].modifyFix.a.push_back(j + 1);
						spectraBuffer[k].modifyFix.b.push_back(i);

					}
			//加可变修饰信息
			spectraBuffer[k].modifyVar.a.clear();
			spectraBuffer[k].modifyVar.b.clear();
			for (int i = 0; i < allRecord[z].a.size(); i++){
			
				spectraBuffer[k].modifyVar.a.push_back(allRecord[z].a[i] + 1);
				spectraBuffer[k].modifyVar.b.push_back(allRecord[z].b[i]);
			
			}

			//cout << spectraBuffer[k].modifyFix.a.size() << endl;
			//cout << spectraBuffer[k].modifyFix.b.size() << endl;
			//cout << spectraBuffer[k].modifyVar.a.size() << endl;
			//cout << spectraBuffer[k].modifyVar.b.size() << endl;
			//getchar();

		}

	}
	else{

		//求出肽段的b、y离子
		int len = pep.squence.size() - 1;
		for (int i = 0; i < len; i++){

			tempMass += MassTableTemp[pep.squence[i] - 'A'];
			for (int j = 0; j < allRecord[z].a.size(); j++){

				if (i == allRecord[z].a[j])
					tempMass += modifyVar[allRecord[z].b[j]].modifyMass;

			}

			mzIon[0][i] = tempMass + pMass;
			mzIon[1][i] = (mzIon[0][i] + pMass) / 2;
			mzIon[2][len - i - 1] = pep.mass - tempMass + pMass;
			mzIon[3][len - i - 1] = (mzIon[2][len - i - 1] + pMass) / 2;

		}



		//FragmentMassFile << pep.squence << "\t" << pep.mass << endl;
		//for (int i = 0; i < 4; i++){
		//	for (int j = 0; j < len; j++)
		//		FragmentMassFile << mzIon[i][j] << "\t";
		//	FragmentMassFile << endl;
		//}


		//初始化和计算得分
		int speIon[4][200], pepIon[4][200];		//实验谱图和理论谱图的b、y离子
		double sumAbundance = 0.0;
		double maxAbundance = 0.0;
		int l1 = 2, l2 = 2;							//对应窗口5时的l1和l2参数

		for (int i = 0; i < 4; i++)
			for (int j = 0; j < len + l2 + 1; j++){

				speIon[i][j] = 0;
				pepIon[i][j] = 0;

			}

		if (tol_fragment_type == "ppm"){

			for (int i = 0; i < 4; i++){

				for (int j = 0; j < len; j++){

					pepIon[i][j] = 1;

					for (int l = 0; l < spectraBuffer[k].peaks.size(); l++){

						if (maxAbundance < spectraBuffer[k].peaks[l].abundance)
							maxAbundance = spectraBuffer[k].peaks[l].abundance;

						if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < (detaFragment*mzIon[i][j])){
						//if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < (detaFragment*spectraBuffer[k].peaks[l].mass2chargeRatio)){
						//if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < (detaFragment*1000)){

							//sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);

							if (spectraBuffer[k].isMatch[l] == 0){

								sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
								spectraBuffer[k].isMatch[l] == 1;

							}
							speIon[i][j] = 1;
							break;

						}

					}

				}

			}

		}
		else if (tol_fragment_type == "Da"){

			for (int i = 0; i < 4; i++){

				for (int j = 0; j < len; j++){

					pepIon[i][j] = 1;

					for (int l = 0; l < spectraBuffer[k].peaks.size(); l++){

						if (maxAbundance < spectraBuffer[k].peaks[l].abundance)
							maxAbundance = spectraBuffer[k].peaks[l].abundance;

						if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < detaFragment){

							//sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);

							if (spectraBuffer[k].isMatch[l] == 0){

								sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
								spectraBuffer[k].isMatch[l] == 1;

							}
							speIon[i][j] = 1;
							break;

						}

					}

				}

			}

		}

		double tempScore = 0.0;

		double temp1 = 0.0;

		for (int i = 0; i < 4; i++){

			int temp2 = 0;

			for (int j = 0; j <= l2; j++)
				temp2 += ((pepIon[i][j] - speIon[i][j]) * (pepIon[i][j] - speIon[i][j]));
			temp1 += exp(gama*temp2);

			for (int j = l2 + 1; j <= l1 + l2; j++){

				temp2 += (speIon[i][j] - pepIon[i][j])*(speIon[i][j] - pepIon[i][j]);
				temp1 += exp(gama*temp2);

			}

			for (int j = l1 + 1; j < len; j++){

				temp2 = temp2 + (speIon[i][j + l2] - pepIon[i][j + l2])*(speIon[i][j + l2] - pepIon[i][j + l2]) - (speIon[i][j - l1 - 1] - pepIon[i][j - l1 - 1])*(speIon[i][j - l1 - 1] - pepIon[i][j - l1 - 1]);
				temp1 += exp(gama*temp2);

			}

		}

		tempScore = (sumAbundance*100 / sqrt(maxAbundance)) * pow(temp1, theta) / (len << 2);

		if (tempScore > spectraBuffer[k].scoreMax){

			spectraBuffer[k].scoreMax = tempScore;
			spectraBuffer[k].peptide = pep.squence;
			//待记录修饰信息
			//加固定修饰信息
			spectraBuffer[k].modifyFix.a.clear();
			spectraBuffer[k].modifyFix.b.clear();

			for (int i = 0; i < modifyFix.size(); i++)
				for (int j = 0; j < pep.squence.size(); j++)
					if (pep.squence[j] == modifyFix[i].modifyAminoacid){

						spectraBuffer[k].modifyFix.a.push_back(j + 1);
						spectraBuffer[k].modifyFix.b.push_back(i);

					}
			//加可变修饰信息
			spectraBuffer[k].modifyVar.a.clear();
			spectraBuffer[k].modifyVar.b.clear();
			for (int i = 0; i < allRecord[z].a.size(); i++){

				spectraBuffer[k].modifyVar.a.push_back(allRecord[z].a[i] + 1);
				spectraBuffer[k].modifyVar.b.push_back(allRecord[z].b[i]);

			}

			//cout << spectraBuffer[k].modifyFix.a.size() << endl;
			//cout << spectraBuffer[k].modifyFix.b.size() << endl;
			//cout << spectraBuffer[k].modifyVar.a.size() << endl;
			//cout << spectraBuffer[k].modifyVar.b.size() << endl;
			//getchar();
		}
	}
}


//针对无修饰和只有固定修饰下的打分
void ScoreKSDP(peptideInfo pep, int k, double addMass){

	double mzIon[4][200];				//200对应肽段的最长长度
	double tempMass = 0.0;
	pep.mass += addMass;

	spectraBuffer[k].isMatch.assign(spectraBuffer[k].peaks.size(), 0);

	//如果母离子的电荷<=2,只考虑计算b+、y+离子，不考虑b++、y++离子
	if (spectraBuffer[k].charge <= 2){
	
		//求出肽段的b、y离子
		int len = pep.squence.size() - 1;
		for (int i = 0; i < len; i++){

			tempMass += MassTable[pep.squence[i] - 'A'];

			mzIon[0][i] = tempMass + pMass;
			//mzIon[1][i] = (mzIon[0][i] + pMass) / 2;
			mzIon[1][len - i - 1] = pep.mass - tempMass + pMass;
			//mzIon[3][len - i - 1] = (mzIon[2][len - i - 1] + pMass) / 2;

		}



		//FragmentMassFile << pep.squence << "\t" << pep.mass << endl;
		//for (int i = 0; i < 4; i++){
		//	for (int j = 0; j < len; j++)
		//		FragmentMassFile << mzIon[i][j] << "\t";
		//	FragmentMassFile << endl;
		//}


		//初始化和计算得分
		int speIon[2][200], pepIon[2][200];		//实验谱图和理论谱图的b、y离子
		double sumAbundance = 0.0;
		double maxAbundance = 0.0;
		int l1 = 2, l2 = 2;							//对应窗口5时的l1和l2参数

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < len + l2 + 1; j++){

				speIon[i][j] = 0;
				pepIon[i][j] = 0;

			}

		if (tol_fragment_type == "ppm"){

			for (int i = 0; i < 2; i++){

				for (int j = 0; j < len; j++){

					pepIon[i][j] = 1;

					for (int l = 0; l < spectraBuffer[k].peaks.size(); l++){

						if (maxAbundance < spectraBuffer[k].peaks[l].abundance)
							maxAbundance = spectraBuffer[k].peaks[l].abundance;

						if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < (detaFragment*mzIon[i][j])){
							//cout << spectraBuffer[k].peaks[l].mass2chargeRatio << endl;
							//sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
							if (spectraBuffer[k].isMatch[l] == 0){
								sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
								spectraBuffer[k].isMatch[l] = 1;
							}

							speIon[i][j] = 1;
							break;

						}

					}

				}

			}

		}
		else if (tol_fragment_type == "Da"){

			for (int i = 0; i < 4; i++){

				for (int j = 0; j < len; j++){

					pepIon[i][j] = 1;

					for (int l = 0; l < spectraBuffer[k].peaks.size(); l++){

						if (maxAbundance < spectraBuffer[k].peaks[l].abundance)
							maxAbundance = spectraBuffer[k].peaks[l].abundance;

						if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < detaFragment){

							//sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
							if (spectraBuffer[k].isMatch[l] == 0){
								sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
								spectraBuffer[k].isMatch[l] = 1;
							}
							speIon[i][j] = 1;
							break;

						}

					}

				}

			}

		}

		double tempScore0 = 0.0;

		double temp1 = 0.0;

		for (int i = 0; i < 2; i++){

			int temp2 = 0;

			for (int j = 0; j <= l2; j++)
				temp2 += ((pepIon[i][j] - speIon[i][j]) * (pepIon[i][j] - speIon[i][j]));
			temp1 += exp(gama*temp2);

			for (int j = l2 + 1; j <= l1 + l2; j++){

				temp2 += (speIon[i][j] - pepIon[i][j])*(speIon[i][j] - pepIon[i][j]);
				temp1 += exp(gama*temp2);

			}

			for (int j = l1 + 1; j < len; j++){

				temp2 = temp2 + (speIon[i][j + l2] - pepIon[i][j + l2])*(speIon[i][j + l2] - pepIon[i][j + l2]) - (speIon[i][j - l1 - 1] - pepIon[i][j - l1 - 1])*(speIon[i][j - l1 - 1] - pepIon[i][j - l1 - 1]);
				temp1 += exp(gama*temp2);

			}

		}

		tempScore0 = ((sumAbundance * 100) / sqrt(maxAbundance)) * pow(temp1, theta) / (len << 1);
		//cout << tempScore0 << endl;

		if (tempScore0 > spectraBuffer[k].scoreMax){

			spectraBuffer[k].scoreMax = tempScore0;
			spectraBuffer[k].peptide = pep.squence;
			//待记录修饰信息
			spectraBuffer[k].modifyFix.a.clear();
			spectraBuffer[k].modifyFix.b.clear();

			for (int i = 0; i < modifyFix.size(); i++)
				for (int j = 0; j < pep.squence.size(); j++)
					if (pep.squence[j] == modifyFix[i].modifyAminoacid){

						spectraBuffer[k].modifyFix.a.push_back(j + 1);
						spectraBuffer[k].modifyFix.b.push_back(i);

					}


		}

	}else{
	
		//求出肽段的b、y离子
		int len = pep.squence.size() - 1;
		for (int i = 0; i < len; i++){

			tempMass += MassTable[pep.squence[i] - 'A'];

			mzIon[0][i] = tempMass + pMass;
			mzIon[1][i] = (mzIon[0][i] + pMass) / 2;
			mzIon[2][len - i - 1] = pep.mass - tempMass + pMass;
			mzIon[3][len - i - 1] = (mzIon[2][len - i - 1] + pMass) / 2;

		}



		//FragmentMassFile << pep.squence << "\t" << pep.mass << endl;
		//for (int i = 0; i < 4; i++){
		//	for (int j = 0; j < len; j++)
		//		FragmentMassFile << mzIon[i][j] << "\t";
		//	FragmentMassFile << endl;
		//}


		//初始化和计算得分
		int speIon[4][200], pepIon[4][200];		//实验谱图和理论谱图的b、y离子
		double sumAbundance = 0.0;
		double maxAbundance = 0.0;
		int l1 = 2, l2 = 2;							//对应窗口5时的l1和l2参数

		for (int i = 0; i < 4; i++)
			for (int j = 0; j < len + l2 + 1; j++){

				speIon[i][j] = 0;
				pepIon[i][j] = 0;

			}

		if (tol_fragment_type == "ppm"){

			for (int i = 0; i < 4; i++){

				for (int j = 0; j < len; j++){

					pepIon[i][j] = 1;

					for (int l = 0; l < spectraBuffer[k].peaks.size(); l++){

						if (maxAbundance < spectraBuffer[k].peaks[l].abundance)
							maxAbundance = spectraBuffer[k].peaks[l].abundance;

						if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < (detaFragment*mzIon[i][j])){

							//sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
							if (spectraBuffer[k].isMatch[l] == 0){
								sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
								spectraBuffer[k].isMatch[l] = 1;
							}
							speIon[i][j] = 1;
							break;

						}

					}

				}

			}

		}
		else if (tol_fragment_type == "Da"){

			for (int i = 0; i < 4; i++){

				for (int j = 0; j < len; j++){

					pepIon[i][j] = 1;

					for (int l = 0; l < spectraBuffer[k].peaks.size(); l++){

						if (maxAbundance < spectraBuffer[k].peaks[l].abundance)
							maxAbundance = spectraBuffer[k].peaks[l].abundance;

						if (fabs(spectraBuffer[k].peaks[l].mass2chargeRatio - mzIon[i][j]) < detaFragment){

							//sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
							if (spectraBuffer[k].isMatch[l] == 0){
								sumAbundance += sqrt(spectraBuffer[k].peaks[l].abundance);
								spectraBuffer[k].isMatch[l] = 1;
							}
							speIon[i][j] = 1;
							break;

						}

					}

				}

			}

		}

		double tempScore = 0.0;

		double temp1 = 0.0;

		for (int i = 0; i < 4; i++){

			int temp2 = 0;

			for (int j = 0; j <= l2; j++)
				temp2 += ((pepIon[i][j] - speIon[i][j]) * (pepIon[i][j] - speIon[i][j]));
			temp1 += exp(gama*temp2);

			for (int j = l2 + 1; j <= l1 + l2; j++){

				temp2 += (speIon[i][j] - pepIon[i][j])*(speIon[i][j] - pepIon[i][j]);
				temp1 += exp(gama*temp2);

			}

			for (int j = l1 + 1; j < len; j++){

				temp2 = temp2 + (speIon[i][j + l2] - pepIon[i][j + l2])*(speIon[i][j + l2] - pepIon[i][j + l2]) - (speIon[i][j - l1 - 1] - pepIon[i][j - l1 - 1])*(speIon[i][j - l1 - 1] - pepIon[i][j - l1 - 1]);
				temp1 += exp(gama*temp2);

			}

		}

		tempScore = ((sumAbundance * 100) / sqrt(maxAbundance)) * pow(temp1, theta) / (len << 2);

		if (tempScore > spectraBuffer[k].scoreMax){

			spectraBuffer[k].scoreMax = tempScore;
			spectraBuffer[k].peptide = pep.squence;
			//待记录修饰信息
			//加固定修饰信息
			spectraBuffer[k].modifyFix.a.clear();
			spectraBuffer[k].modifyFix.b.clear();

			for (int i = 0; i < modifyFix.size(); i++)
				for (int j = 0; j < pep.squence.size(); j++)
					if (pep.squence[j] == modifyFix[i].modifyAminoacid){

						spectraBuffer[k].modifyFix.a.push_back(j + 1);
						spectraBuffer[k].modifyFix.b.push_back(i);

					}
		}
	}
}


/*
函数名：	ScoreRank
功能：	以肽段为中心进行打分
输入：	无
输出：	查询返回结果及与pFind对比结果
*/
void ScoreRank(){

	double tempMass = 0.0;

	//求出母离子Frag和碎片离子的Frag
	if (tol_precursor_type == "ppm")
		detaPrecursor = tol_precursor_value * 0.000001;
	else
		detaPrecursor = tol_precursor_value;

	if (tol_fragment_type == "ppm")
		detaFragment = tol_fragment_value * 0.000001;
	else
		detaFragment = tol_fragment_value;
	//求出母离子deta和碎片离子的deta


	//进行遍历查询

	if (tol_precursor_type == "ppm"){

		if (modifyFix.size() != 0){		//固定修饰不为零
		
			//修改MassTable，将发生固定修饰的氨基酸变为另一种氨基酸，加质量
			for (int l = 0; l < modifyFix.size(); l++){
				MassTable[modifyFix[l].modifyAminoacid - 'A'] += modifyFix[l].modifyMass;
			}
			
			if (modifyVar.size() != 0){
			
				//可变修饰情况下搜索
				for (int i = 0; i < peptideBuffer.size(); i++){

					recordInfo r;
					allRecord.clear();			//清零
					dpsPep(0, r, peptideBuffer[i].squence);

					//加固定修饰总质量
					double addMass = 0.0;

					//加一层循环实现遍历动态生成的加了可变修饰肽段集合
					for (int l = 0; l < allRecord.size(); l++){
					
						//MassTableTemp.clear();			
						MassTableTemp.assign(MassTable.begin(), MassTable.end());		//assign函数自身带有删除功能

						addMass = 0.0;

						//加固定修饰总质量
						for (int k = 0; k < modifyFix.size(); k++)
							for (int t = 0; t < peptideBuffer[i].squence.size(); t++)
								if (peptideBuffer[i].squence[t] == modifyFix[k].modifyAminoacid)
									addMass += modifyFix[k].modifyMass;

						//加可变修饰总质量
						double VarMass = 0.0;
						for (int k = 0; k < allRecord[l].a.size(); k++){
						
							//这边不能修改质量表，因为存在没有发生可变修饰的氨基酸，待会碎片离子会计算进去
							//MassTableTemp[peptideBuffer[i].squence[allRecord[l].a[k]] - 'A'] += modifyVar[allRecord[l].b[k]].modifyMass;
							
							//cout << MassTableTemp[peptideBuffer[i].squence[allRecord[l].a[k]] - 'A'] << endl;
							addMass += modifyVar[allRecord[l].b[k]].modifyMass;
							VarMass += modifyVar[allRecord[l].b[k]].modifyMass;
						
						}
						//cout << addMass << endl;
						//getchar();

						//cout << addMass << endl;
						for (int j = 0; j < spectraBuffer.size(); j++){			//TODO: 这里可以用二分优化

							tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;//TODO: 这边可以加速，存在结构体里

							if (fabs(tempMass - (peptideBuffer[i].mass + addMass)) > (detaPrecursor * (peptideBuffer[i].mass + addMass)))
							//if (fabs(tempMass - (peptideBuffer[i].mass + addMass)) > (detaPrecursor * (tempMass)))
								continue;

							//cout << peptideBuffer[i].squence << endl;
							//getchar();
							ScoreKSDP_Var(peptideBuffer[i], j, l, addMass);

						}

					
					}

				}
				cout << "固定修饰和可变修饰都有的情况结束." << endl;


			}else if (modifyVar.size() == 0){
				//固定修饰不为零，可变修饰为零 的情况下搜索
				cout << "固定修饰不为零，可变修饰为零 的情况下搜索." << endl;

				for (int i = 0; i < peptideBuffer.size(); i++){

					double addMass = 0.0;			//所有修饰的质量
					for (int k = 0; k < modifyFix.size(); k++)
						for (int t = 0; t < peptideBuffer[i].squence.size(); t++)
							if (peptideBuffer[i].squence[t] == modifyFix[k].modifyAminoacid)
								addMass += modifyFix[k].modifyMass;

					for (int j = 0; j < spectraBuffer.size(); j++){

						tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;//这边可以加速，存在结构体里

						if (fabs(tempMass - (peptideBuffer[i].mass + addMass)) > (detaPrecursor * (peptideBuffer[i].mass + addMass)))
							continue;
						//cout << peptideBuffer[i].squence << endl;
						//getchar();
						ScoreKSDP(peptideBuffer[i], j, addMass);

					}

				}

			}
		
		}

		
		if (modifyFix.size() == 0){

			if (modifyVar.size() == 0){			//固定修饰为零，可变修饰为零的情况

				for (int i = 0; i < peptideBuffer.size(); i++){

					for (int j = 0; j < spectraBuffer.size(); j++){

						tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;

						if (fabs(tempMass - peptideBuffer[i].mass) > (detaPrecursor * peptideBuffer[i].mass))
							continue;

						ScoreKSDP(peptideBuffer[i], j, 0.0);

					}

				}
				cout << "即没有固定修饰，也没有可变修饰." << endl;

			}
			else if (modifyVar.size() != 0){
			
				for (int i = 0; i < peptideBuffer.size(); i++){

					recordInfo r;
					allRecord.clear();			//清零
					dpsPep(0, r, peptideBuffer[i].squence);

					////加固定修饰总质量
					double addMass = 0.0;

					//加一层循环实现遍历动态生成的加了可变修饰肽段集合
					for (int l = 0; l < allRecord.size(); l++){

						MassTableTemp.clear();			//assign函数自身其实有删除功能，为保险
						MassTableTemp.assign(MassTable.begin(), MassTable.end());

						addMass = 0.0;

						//加可变修饰总质量
						for (int k = 0; k < allRecord[l].a.size(); k++){

							MassTableTemp[peptideBuffer[i].squence[allRecord[l].a[k]] - 'A'] += modifyVar[allRecord[l].b[k]].modifyMass;

							addMass += modifyVar[allRecord[l].b[k]].modifyMass;

						}

						for (int j = 0; j < spectraBuffer.size(); j++){			//TODO: 这里可以用二分优化

							tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;//TODO: 这边可以加速，存在结构体里

							if (fabs(tempMass - peptideBuffer[i].mass - addMass) >(detaPrecursor * (peptideBuffer[i].mass + addMass)))
							//if (fabs(tempMass - (peptideBuffer[i].mass + addMass)) > (detaPrecursor * (tempMass)))
							//if (fabs(tempMass - (peptideBuffer[i].mass + addMass)) > (detaPrecursor * (1000)))
								continue;

							ScoreKSDP_Var(peptideBuffer[i], j, l, addMass);

						}
					}
				}
				cout << "固定修饰没有，可变修饰有的情况结束." << endl;
			}
		}


	}else if (tol_precursor_type == "Da"){
	
		if (modifyFix.size() != 0){

			//修改MassTable，将发生固定修饰的氨基酸变为另一种氨基酸，加质量
			for (int l = 0; l < modifyFix.size(); l++){

				MassTable[modifyFix[l].modifyAminoacid - 'A'] += modifyFix[l].modifyMass;

			}

			if (modifyVar.size() != 0){

				//可变修饰情况下搜索
				for (int i = 0; i < peptideBuffer.size(); i++){

					recordInfo r;
					allRecord.clear();
					dpsPep(0, r, peptideBuffer[i].squence);

					double addMass = 0.0;
					double VarMass = 0.0;
					for (int l = 0; l < allRecord.size(); l++){

						MassTableTemp.clear();			//assign函数自身其实有删除功能，为保险
						MassTableTemp.assign(MassTable.begin(), MassTable.end());

						addMass = 0.0;
						VarMass = 0.0;
						//加可变修饰总质量
						for (int k = 0; k < allRecord[l].a.size(); k++){

							//MassTableTemp[peptideBuffer[i].squence[allRecord[l].a[k]] - 'A'] += modifyVar[allRecord[l].b[k]].modifyMass;

							VarMass += modifyVar[allRecord[l].b[k]].modifyMass;
							addMass += modifyVar[allRecord[l].b[k]].modifyMass;

						}

						//加固定修饰总质量
						for (int k = 0; k < modifyFix.size(); k++)
							for (int t = 0; t < peptideBuffer[i].squence.size(); t++)
								if (peptideBuffer[i].squence[t] == modifyFix[k].modifyAminoacid)
									addMass += modifyFix[k].modifyMass;
						//cout << addMass << endl;
						for (int j = 0; j < spectraBuffer.size(); j++){			//TODO: 这里可以用二分优化

							tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;//TODO: 这边可以加速，存在结构体里

							if (fabs(tempMass - peptideBuffer[i].mass - addMass) > detaPrecursor)
								continue;

							//cout << peptideBuffer[i].squence << endl;
							//getchar();
							ScoreKSDP_Var(peptideBuffer[i], j, l, VarMass);

						}

					}

				}
				cout << "固定修饰和可变修饰都有的情况结束." << endl;

			}else if (modifyVar.size() == 0){
				//固定修饰情况下搜索

				for (int i = 0; i < peptideBuffer.size(); i++){

					double addMass = 0.0;			//所有修饰的质量
					for (int k = 0; k < modifyFix.size(); k++)
						for (int t = 0; t < peptideBuffer[i].squence.size(); t++)
							if (peptideBuffer[i].squence[t] == modifyFix[k].modifyAminoacid)
								addMass += modifyFix[k].modifyMass;

					for (int j = 0; j < spectraBuffer.size(); j++){

						tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;

						if (fabs(tempMass - peptideBuffer[i].mass - addMass) > detaPrecursor)
							continue;

						ScoreKSDP(peptideBuffer[i], j, addMass);

					}

				}

			}

		}

		//固定修饰和可变修饰都没有
		if (modifyFix.size() == 0){

			if (modifyVar.size() == 0){

				for (int i = 0; i < peptideBuffer.size(); i++){

					for (int j = 0; j < spectraBuffer.size(); j++){

						tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;

						if (fabs(tempMass - peptideBuffer[i].mass) > detaPrecursor)
							continue;

						ScoreKSDP(peptideBuffer[i], j, 0.0);

					}

				}
				cout << "既没有固定修饰，也没有可变修饰." << endl;

			}
			else if (modifyVar.size() != 0){
			
				for (int i = 0; i < peptideBuffer.size(); i++){

					recordInfo r;
					allRecord.clear();			//清零
					dpsPep(0, r, peptideBuffer[i].squence);

					////加固定修饰总质量
					double addMass = 0.0;

					//加一层循环实现遍历动态生成的加了可变修饰肽段集合
					for (int l = 0; l < allRecord.size(); l++){

						MassTableTemp.clear();			//assign函数自身其实有删除功能，为保险
						MassTableTemp.assign(MassTable.begin(), MassTable.end());

						addMass = 0.0;

						//加可变修饰总质量
						for (int k = 0; k < allRecord[l].a.size(); k++){

							//MassTableTemp[peptideBuffer[i].squence[allRecord[l].a[k]] - 'A'] += modifyVar[allRecord[l].b[k]].modifyMass;

							addMass += modifyVar[allRecord[l].b[k]].modifyMass;

						}

						for (int j = 0; j < spectraBuffer.size(); j++){			//TODO: 这里可以用二分优化

							tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;//TODO: 这边可以加速，存在结构体里

							if (fabs(tempMass - peptideBuffer[i].mass - addMass) > detaPrecursor)
								continue;

							ScoreKSDP_Var(peptideBuffer[i], j, l, addMass);

						}


					}

				}
				cout << "固定修饰没有，可变修饰有的情况结束." << endl;
			
			}

		}

	}


	//输出结果
	//ofstream resultOut(output_path_result);
	ofstream checkOut("myResult.txt");
	for (int i = 0; i < spectraBuffer.size(); i++){

		//resultOut << spectraBuffer[i].title << " ";
		if (spectraBuffer[i].peptide.size() == 0){
			//resultOut << "NULL" << endl;
			checkOut << "NULL" << endl;
		}
		else{
			//resultOut << spectraBuffer[i].peptide << " " << spectraBuffer[i].scoreMax << endl;
			checkOut << spectraBuffer[i].peptide << endl;
		}


	}

	//resultOut.close();
	checkOut.close();

}