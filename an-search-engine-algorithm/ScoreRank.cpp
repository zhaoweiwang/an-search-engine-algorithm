/********************************************************
�ļ�����		ScoreRank.cpp
�����ˣ�		Mako Wan
��ע���ڣ�	2015-7-26
������		��ȡDictionaryFile.txt�ļ������Ķ�Ϊ���Ľ��д��
********************************************************/

#include "Head.h"

extern vector<spectraInfo> spectraBuffer;
extern vector<peptideInfo> peptideBuffer;

extern vector<recordInfo> allRecord;

extern vector<modifyInfo> modifyFix;
extern vector<modifyInfo> modifyVar;


extern string tol_precursor_type;           //ppm��Da����
extern double tol_precursor_value;	        //����20ppm���ϸ���ڻ�С��
extern string tol_fragment_type;
extern double tol_fragment_value;			//����0.5Da���ϸ���ڻ�С��

extern string output_path_result;

const double pMass = 1.00727647012;			//��������

vector<double> MassTable = { 71.03711, 166.99836, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891, 113.08406, 181.01401, 128.09496, 113.08406, 131.04048, 114.04293, 243.02965, 97.05276, 128.05858, 156.10111, 87.03203, 101.04768, 0.00000, 99.06841, 186.07931, 113.08406, 163.06332, 128.55059 };
vector<double> MassTableTemp;


double detaPrecursor = 0.0;
double detaFragment = 0.0;
double gama = -0.9;
double theta = 0.5;


//ofstream FragmentMassFile("FragmentMassFile.txt");

//��Կɱ����ν��д��
void ScoreKSDP_Var(peptideInfo pep, int k, int z, double addMass){

	double mzIon[4][200];				//200��Ӧ�Ķε������
	double tempMass = 0.0;
	
	//for (int k = 0; k < allRecord[z].a.size(); k++){

	//	pep.mass += modifyVar[allRecord[z].b[k]].modifyMass;

	//}
	pep.mass += addMass;
	
	spectraBuffer[k].isMatch.assign(spectraBuffer[k].peaks.size(), 0);

	//���ĸ���ӵĵ��<=2,ֻ���Ǽ���b+��y+���ӣ�������b++��y++����
	if (spectraBuffer[k].charge <= 2){

		//����Ķε�b��y����
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


		//��ʼ���ͼ���÷�
		int speIon[2][200], pepIon[2][200];		//ʵ����ͼ��������ͼ��b��y����
		double sumAbundance = 0.0;
		double maxAbundance = 0.0;
		int l1 = 2, l2 = 2;							//��Ӧ����5ʱ��l1��l2����

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
		//cout << "����:" << tempScore0 << endl;

		if (tempScore0 > spectraBuffer[k].scoreMax){

			spectraBuffer[k].scoreMax = tempScore0;
			spectraBuffer[k].peptide = pep.squence;
			//����¼������Ϣ
			//�ӹ̶�������Ϣ
			spectraBuffer[k].modifyFix.a.clear();
			spectraBuffer[k].modifyFix.b.clear();

			for (int i = 0; i < modifyFix.size(); i++)
				for (int j = 0; j < pep.squence.size(); j++)
					if (pep.squence[j] == modifyFix[i].modifyAminoacid){

						spectraBuffer[k].modifyFix.a.push_back(j + 1);
						spectraBuffer[k].modifyFix.b.push_back(i);

					}
			//�ӿɱ�������Ϣ
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

		//����Ķε�b��y����
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


		//��ʼ���ͼ���÷�
		int speIon[4][200], pepIon[4][200];		//ʵ����ͼ��������ͼ��b��y����
		double sumAbundance = 0.0;
		double maxAbundance = 0.0;
		int l1 = 2, l2 = 2;							//��Ӧ����5ʱ��l1��l2����

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
			//����¼������Ϣ
			//�ӹ̶�������Ϣ
			spectraBuffer[k].modifyFix.a.clear();
			spectraBuffer[k].modifyFix.b.clear();

			for (int i = 0; i < modifyFix.size(); i++)
				for (int j = 0; j < pep.squence.size(); j++)
					if (pep.squence[j] == modifyFix[i].modifyAminoacid){

						spectraBuffer[k].modifyFix.a.push_back(j + 1);
						spectraBuffer[k].modifyFix.b.push_back(i);

					}
			//�ӿɱ�������Ϣ
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


//��������κ�ֻ�й̶������µĴ��
void ScoreKSDP(peptideInfo pep, int k, double addMass){

	double mzIon[4][200];				//200��Ӧ�Ķε������
	double tempMass = 0.0;
	pep.mass += addMass;

	spectraBuffer[k].isMatch.assign(spectraBuffer[k].peaks.size(), 0);

	//���ĸ���ӵĵ��<=2,ֻ���Ǽ���b+��y+���ӣ�������b++��y++����
	if (spectraBuffer[k].charge <= 2){
	
		//����Ķε�b��y����
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


		//��ʼ���ͼ���÷�
		int speIon[2][200], pepIon[2][200];		//ʵ����ͼ��������ͼ��b��y����
		double sumAbundance = 0.0;
		double maxAbundance = 0.0;
		int l1 = 2, l2 = 2;							//��Ӧ����5ʱ��l1��l2����

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
			//����¼������Ϣ
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
	
		//����Ķε�b��y����
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


		//��ʼ���ͼ���÷�
		int speIon[4][200], pepIon[4][200];		//ʵ����ͼ��������ͼ��b��y����
		double sumAbundance = 0.0;
		double maxAbundance = 0.0;
		int l1 = 2, l2 = 2;							//��Ӧ����5ʱ��l1��l2����

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
			//����¼������Ϣ
			//�ӹ̶�������Ϣ
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
��������	ScoreRank
���ܣ�	���Ķ�Ϊ���Ľ��д��
���룺	��
�����	��ѯ���ؽ������pFind�ԱȽ��
*/
void ScoreRank(){

	double tempMass = 0.0;

	//���ĸ����Frag����Ƭ���ӵ�Frag
	if (tol_precursor_type == "ppm")
		detaPrecursor = tol_precursor_value * 0.000001;
	else
		detaPrecursor = tol_precursor_value;

	if (tol_fragment_type == "ppm")
		detaFragment = tol_fragment_value * 0.000001;
	else
		detaFragment = tol_fragment_value;
	//���ĸ����deta����Ƭ���ӵ�deta


	//���б�����ѯ

	if (tol_precursor_type == "ppm"){

		if (modifyFix.size() != 0){		//�̶����β�Ϊ��
		
			//�޸�MassTable���������̶����εİ������Ϊ��һ�ְ����ᣬ������
			for (int l = 0; l < modifyFix.size(); l++){
				MassTable[modifyFix[l].modifyAminoacid - 'A'] += modifyFix[l].modifyMass;
			}
			
			if (modifyVar.size() != 0){
			
				//�ɱ��������������
				for (int i = 0; i < peptideBuffer.size(); i++){

					recordInfo r;
					allRecord.clear();			//����
					dpsPep(0, r, peptideBuffer[i].squence);

					//�ӹ̶�����������
					double addMass = 0.0;

					//��һ��ѭ��ʵ�ֱ�����̬���ɵļ��˿ɱ������Ķμ���
					for (int l = 0; l < allRecord.size(); l++){
					
						//MassTableTemp.clear();			
						MassTableTemp.assign(MassTable.begin(), MassTable.end());		//assign�����������ɾ������

						addMass = 0.0;

						//�ӹ̶�����������
						for (int k = 0; k < modifyFix.size(); k++)
							for (int t = 0; t < peptideBuffer[i].squence.size(); t++)
								if (peptideBuffer[i].squence[t] == modifyFix[k].modifyAminoacid)
									addMass += modifyFix[k].modifyMass;

						//�ӿɱ�����������
						double VarMass = 0.0;
						for (int k = 0; k < allRecord[l].a.size(); k++){
						
							//��߲����޸���������Ϊ����û�з����ɱ����εİ����ᣬ������Ƭ���ӻ�����ȥ
							//MassTableTemp[peptideBuffer[i].squence[allRecord[l].a[k]] - 'A'] += modifyVar[allRecord[l].b[k]].modifyMass;
							
							//cout << MassTableTemp[peptideBuffer[i].squence[allRecord[l].a[k]] - 'A'] << endl;
							addMass += modifyVar[allRecord[l].b[k]].modifyMass;
							VarMass += modifyVar[allRecord[l].b[k]].modifyMass;
						
						}
						//cout << addMass << endl;
						//getchar();

						//cout << addMass << endl;
						for (int j = 0; j < spectraBuffer.size(); j++){			//TODO: ��������ö����Ż�

							tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;//TODO: ��߿��Լ��٣����ڽṹ����

							if (fabs(tempMass - (peptideBuffer[i].mass + addMass)) > (detaPrecursor * (peptideBuffer[i].mass + addMass)))
							//if (fabs(tempMass - (peptideBuffer[i].mass + addMass)) > (detaPrecursor * (tempMass)))
								continue;

							//cout << peptideBuffer[i].squence << endl;
							//getchar();
							ScoreKSDP_Var(peptideBuffer[i], j, l, addMass);

						}

					
					}

				}
				cout << "�̶����κͿɱ����ζ��е��������." << endl;


			}else if (modifyVar.size() == 0){
				//�̶����β�Ϊ�㣬�ɱ�����Ϊ�� �����������
				cout << "�̶����β�Ϊ�㣬�ɱ�����Ϊ�� �����������." << endl;

				for (int i = 0; i < peptideBuffer.size(); i++){

					double addMass = 0.0;			//�������ε�����
					for (int k = 0; k < modifyFix.size(); k++)
						for (int t = 0; t < peptideBuffer[i].squence.size(); t++)
							if (peptideBuffer[i].squence[t] == modifyFix[k].modifyAminoacid)
								addMass += modifyFix[k].modifyMass;

					for (int j = 0; j < spectraBuffer.size(); j++){

						tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;//��߿��Լ��٣����ڽṹ����

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

			if (modifyVar.size() == 0){			//�̶�����Ϊ�㣬�ɱ�����Ϊ������

				for (int i = 0; i < peptideBuffer.size(); i++){

					for (int j = 0; j < spectraBuffer.size(); j++){

						tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;

						if (fabs(tempMass - peptideBuffer[i].mass) > (detaPrecursor * peptideBuffer[i].mass))
							continue;

						ScoreKSDP(peptideBuffer[i], j, 0.0);

					}

				}
				cout << "��û�й̶����Σ�Ҳû�пɱ�����." << endl;

			}
			else if (modifyVar.size() != 0){
			
				for (int i = 0; i < peptideBuffer.size(); i++){

					recordInfo r;
					allRecord.clear();			//����
					dpsPep(0, r, peptideBuffer[i].squence);

					////�ӹ̶�����������
					double addMass = 0.0;

					//��һ��ѭ��ʵ�ֱ�����̬���ɵļ��˿ɱ������Ķμ���
					for (int l = 0; l < allRecord.size(); l++){

						MassTableTemp.clear();			//assign����������ʵ��ɾ�����ܣ�Ϊ����
						MassTableTemp.assign(MassTable.begin(), MassTable.end());

						addMass = 0.0;

						//�ӿɱ�����������
						for (int k = 0; k < allRecord[l].a.size(); k++){

							MassTableTemp[peptideBuffer[i].squence[allRecord[l].a[k]] - 'A'] += modifyVar[allRecord[l].b[k]].modifyMass;

							addMass += modifyVar[allRecord[l].b[k]].modifyMass;

						}

						for (int j = 0; j < spectraBuffer.size(); j++){			//TODO: ��������ö����Ż�

							tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;//TODO: ��߿��Լ��٣����ڽṹ����

							if (fabs(tempMass - peptideBuffer[i].mass - addMass) >(detaPrecursor * (peptideBuffer[i].mass + addMass)))
							//if (fabs(tempMass - (peptideBuffer[i].mass + addMass)) > (detaPrecursor * (tempMass)))
							//if (fabs(tempMass - (peptideBuffer[i].mass + addMass)) > (detaPrecursor * (1000)))
								continue;

							ScoreKSDP_Var(peptideBuffer[i], j, l, addMass);

						}
					}
				}
				cout << "�̶�����û�У��ɱ������е��������." << endl;
			}
		}


	}else if (tol_precursor_type == "Da"){
	
		if (modifyFix.size() != 0){

			//�޸�MassTable���������̶����εİ������Ϊ��һ�ְ����ᣬ������
			for (int l = 0; l < modifyFix.size(); l++){

				MassTable[modifyFix[l].modifyAminoacid - 'A'] += modifyFix[l].modifyMass;

			}

			if (modifyVar.size() != 0){

				//�ɱ��������������
				for (int i = 0; i < peptideBuffer.size(); i++){

					recordInfo r;
					allRecord.clear();
					dpsPep(0, r, peptideBuffer[i].squence);

					double addMass = 0.0;
					double VarMass = 0.0;
					for (int l = 0; l < allRecord.size(); l++){

						MassTableTemp.clear();			//assign����������ʵ��ɾ�����ܣ�Ϊ����
						MassTableTemp.assign(MassTable.begin(), MassTable.end());

						addMass = 0.0;
						VarMass = 0.0;
						//�ӿɱ�����������
						for (int k = 0; k < allRecord[l].a.size(); k++){

							//MassTableTemp[peptideBuffer[i].squence[allRecord[l].a[k]] - 'A'] += modifyVar[allRecord[l].b[k]].modifyMass;

							VarMass += modifyVar[allRecord[l].b[k]].modifyMass;
							addMass += modifyVar[allRecord[l].b[k]].modifyMass;

						}

						//�ӹ̶�����������
						for (int k = 0; k < modifyFix.size(); k++)
							for (int t = 0; t < peptideBuffer[i].squence.size(); t++)
								if (peptideBuffer[i].squence[t] == modifyFix[k].modifyAminoacid)
									addMass += modifyFix[k].modifyMass;
						//cout << addMass << endl;
						for (int j = 0; j < spectraBuffer.size(); j++){			//TODO: ��������ö����Ż�

							tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;//TODO: ��߿��Լ��٣����ڽṹ����

							if (fabs(tempMass - peptideBuffer[i].mass - addMass) > detaPrecursor)
								continue;

							//cout << peptideBuffer[i].squence << endl;
							//getchar();
							ScoreKSDP_Var(peptideBuffer[i], j, l, VarMass);

						}

					}

				}
				cout << "�̶����κͿɱ����ζ��е��������." << endl;

			}else if (modifyVar.size() == 0){
				//�̶��������������

				for (int i = 0; i < peptideBuffer.size(); i++){

					double addMass = 0.0;			//�������ε�����
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

		//�̶����κͿɱ����ζ�û��
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
				cout << "��û�й̶����Σ�Ҳû�пɱ�����." << endl;

			}
			else if (modifyVar.size() != 0){
			
				for (int i = 0; i < peptideBuffer.size(); i++){

					recordInfo r;
					allRecord.clear();			//����
					dpsPep(0, r, peptideBuffer[i].squence);

					////�ӹ̶�����������
					double addMass = 0.0;

					//��һ��ѭ��ʵ�ֱ�����̬���ɵļ��˿ɱ������Ķμ���
					for (int l = 0; l < allRecord.size(); l++){

						MassTableTemp.clear();			//assign����������ʵ��ɾ�����ܣ�Ϊ����
						MassTableTemp.assign(MassTable.begin(), MassTable.end());

						addMass = 0.0;

						//�ӿɱ�����������
						for (int k = 0; k < allRecord[l].a.size(); k++){

							//MassTableTemp[peptideBuffer[i].squence[allRecord[l].a[k]] - 'A'] += modifyVar[allRecord[l].b[k]].modifyMass;

							addMass += modifyVar[allRecord[l].b[k]].modifyMass;

						}

						for (int j = 0; j < spectraBuffer.size(); j++){			//TODO: ��������ö����Ż�

							tempMass = (spectraBuffer[j].pepMass - pMass) * spectraBuffer[j].charge;//TODO: ��߿��Լ��٣����ڽṹ����

							if (fabs(tempMass - peptideBuffer[i].mass - addMass) > detaPrecursor)
								continue;

							ScoreKSDP_Var(peptideBuffer[i], j, l, addMass);

						}


					}

				}
				cout << "�̶�����û�У��ɱ������е��������." << endl;
			
			}

		}

	}


	//������
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