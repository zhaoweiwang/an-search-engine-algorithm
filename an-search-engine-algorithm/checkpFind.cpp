/********************************************************
文件名：		checkpFind.cpp
创建人：		ZhaoweiWang
加注日期：	2015-7-29
描述：		和pFind的结果进行对比
********************************************************/

#include "Head.h"

extern vector<spectraInfo> spectraBuffer;
extern string output_path_result;

extern string input_path_pFind;			//pFind搜索结果

extern vector<modifyInfo> modifyFix;
extern vector<modifyInfo> modifyVar;

bool comp(const spectraInfo &a, const spectraInfo &b){

	return a.scoreMax > b.scoreMax;

}

/*
函数名：	checkpFind
功能：	对比pFind结果
输入：
输出：
*/
void checkpFind(){

	ifstream myResult("myResult.txt", ios::in);
	if (!myResult.is_open()){

		cerr << "不能打开 myResult.txt 文件." << endl;
		exit(1);

	}

	vector<string> my_result;

	string temp = "";
	//while (!myResult.eof()){
	for (int i = 0; i < spectraBuffer.size(); i++){
		myResult >> temp;
		//if (temp != "NULL")
		//	my_result.push_back(temp);
		//else{
		//
		//	temp = "";
		//	my_result.push_back(temp);
		//
		//}
		my_result.push_back(temp);
	}
	myResult.close();

	//ofstream myResultout("Myout.txt");
	//myResultout << my_result.size() << endl;
	//for (int i = 0; i < my_result.size(); i++){

	//		myResultout << my_result[i] << endl;

	//}
	//myResultout.close();


	//my_result.pop_back();
	//for (int i = 0; i < my_result.size(); i++)
	//	cout << my_result[i] << endl;

	//cout << my_result[my_result.size() - 1] << endl;
	//cout << my_result[my_result.size() - 2] << endl;
	//cout << my_result.size() << endl;
	//getchar();

	//cout << my_result.size() << endl;
	//getchar();

	//vector<vector<string>> pFind_result;
	vector<vector<pFindResultInfo>> pFind_result;
	pFind_result.resize(my_result.size());
	//cout << "测试 " << pFind_result[0].size() << endl;

	ifstream pFindFile(input_path_pFind, ios::in);
	if (!pFindFile.is_open()){

		cerr << "不能打开 " << input_path_pFind << " 文件." << endl;
		exit(1);

	}

	int index = 1;
	char t[10];
	sprintf(t, "%d", index);
	string s1 = t;
	string str = "[Spectrum";
	str = str + s1 + "]";
	
	int count = 0;

	while (!pFindFile.eof()){
	
		pFindFile >> temp;
		//cout << temp << endl;
		//getchar();

		if (temp == "[Protein1]") {
			//cout << temp << endl;
			break;
		}

		//cout << temp << endl;
		//cout << str << endl;
		//getchar();

		if (temp == str){

			//cout << "haha " << temp << endl;

			while (count != 7){

				count++;
				pFindFile >> temp;

			}
			count = 0;

			//cout << temp << endl; //Valicand

			int numCand = 0;
			auto pos = temp.find('=');
			temp = temp.substr(pos + 1);
			numCand = atoi(temp.c_str());

			if (numCand == 0){
			
				//string s3 = "NULL";
				pFindResultInfo temppFind;
				temppFind.squence = "NULL";
				temppFind.score = 0;
				pFind_result[index-1].push_back(temppFind);

				index++;
			
			}else{
			
				int count1 = 15;
				pFindResultInfo temppFind;

				while (numCand != 0){

					while (count1--){

						pFindFile >> temp;

						if (count1 == 14){
						
							string s2 = temp.substr(temp.find('=') + 1);
							temppFind.score = atof(s2.c_str());
						
						}
						//cout << "读No " << temp << endl;
						if (count1 == 9){

							//cout << temp << endl;
							string s2 = temp.substr(temp.find('=') + 1);
							//cout << s2 << endl;
							//pFind_result[index-1].push_back(s2);
							//cout << pFind_result[index][0] << endl;
							temppFind.squence = s2;
							pFind_result[index - 1].push_back(temppFind);
						}

					}
					//cout << temp << endl;
					count1 = 15;
					numCand--;

				}

				//cout << pFind_result[index].size() << endl;
				//for (auto i = 0; i < pFind_result[index].size(); i++)
				//	cout << pFind_result[index][i] << endl;

				index++;

			}

		}

		//stringstream ss1;
		//ss1 << index;
		//s1 = ss1.str();

		sprintf(t, "%d", index);
		s1 = t;
		str = "[Spectrum";
		str = str + s1 + "]";
		//cout << str << endl;
		//getchar();
		//pFindFile >> temp;
	
	}
	pFindFile.close();

	//ofstream pFindout("pFindout.txt");
	//pFindout << pFind_result.size() << endl;
	//for (int i = 0; i < pFind_result.size(); i++){
	//	pFindout << "第 " << i+1 << " 个谱图：" << endl;
	//	for (int j = 0; j < pFind_result[i].size(); j++)
	//		pFindout << pFind_result[i][j] << endl;

	//}
	//pFindout.close();

	//cout << pFind_result.size() << endl;
	//cout << pFind_result[pFind_result.size() - 2].size() << endl;
	//cout << pFind_result[pFind_result.size() - 2][pFind_result[pFind_result.size() - 2].size() - 1] << endl;
	//getchar();

	int sameCount = 0;
	int No1Count = 0;
	vector<int> Cant;
	vector<int> topX;
	for (int i = 0; i < my_result.size(); i++){
		
		//if (my_result[i] == "NULL") my_result[i] = "";
		int j;
		for (j = 0; j < pFind_result[i].size(); j++){
		
			if (my_result[i] == pFind_result[i][j].squence){

				topX.push_back(j+1);
				spectraBuffer[i].top = j + 1;
				sameCount++;
				if (j == 0)
					No1Count++;
				else{
				
					if (pFind_result[i][0].score == pFind_result[i][j].score)
						No1Count++;
				
				}
				break;

			}
				
		}

		if (j == pFind_result[i].size()){

			Cant.push_back(i);
			topX.push_back(0);
			spectraBuffer[i].top = 0;

		}
			
	
	}
	//cout << topX.size() << endl;

	cout << "对比结果: " << sameCount << " " << (double)sameCount / (double)my_result.size() << "%" << endl;
	cout << "Top1对比结果: " << No1Count << " " << (double)No1Count / (double)my_result.size() << "%" << endl;

	cout << "没有匹配上的有 " << Cant.size() << " 个." << endl;
	//getchar();
	for (int i = 0; i < Cant.size(); i++)
		cout << "No: " << Cant[i]+1 << " " << my_result[Cant[i]] << endl;


	//输出结果
	sort(spectraBuffer.begin(), spectraBuffer.end(), comp);
	ofstream resultOut(output_path_result);

	for (int i = 0; i < spectraBuffer.size(); i++){

		resultOut << spectraBuffer[i].title << "\t";
		if (spectraBuffer[i].peptide.size() == 0){
			resultOut << "NULL\tnull\t" << spectraBuffer[i].top << "\t" << endl;
		}
		else{
			resultOut << spectraBuffer[i].peptide << "\t";
			if (spectraBuffer[i].modifyFix.a.size() == 0 && spectraBuffer[i].modifyVar.a.size() == 0)	 resultOut << "null\t";
			else if (spectraBuffer[i].modifyFix.a.size() != 0 && spectraBuffer[i].modifyVar.a.size() == 0){
			
				resultOut << spectraBuffer[i].modifyFix.a.size();
				for (int j = 0; j < spectraBuffer[i].modifyFix.a.size(); j++){
				
					resultOut << "|" << spectraBuffer[i].modifyFix.a[j] << "," << modifyFix[spectraBuffer[i].modifyFix.b[j]].modifyName;
				
				}
				resultOut << "\t";
			
			}
			else if (spectraBuffer[i].modifyFix.a.size() == 0 && spectraBuffer[i].modifyVar.a.size() != 0){
			
				resultOut << spectraBuffer[i].modifyVar.a.size();
				for (int j = 0; j < spectraBuffer[i].modifyVar.a.size(); j++){
				
					resultOut << "|" << spectraBuffer[i].modifyVar.a[j] << "," << modifyVar[spectraBuffer[i].modifyVar.b[j]].modifyName;
				
				}
				resultOut << "\t";
			
			}
			else if (spectraBuffer[i].modifyFix.a.size() != 0 && spectraBuffer[i].modifyVar.a.size() != 0){
			
				resultOut << (spectraBuffer[i].modifyFix.a.size() + spectraBuffer[i].modifyVar.a.size());
				for (int j = 0; j < spectraBuffer[i].modifyFix.a.size(); j++){

					resultOut << "|" << spectraBuffer[i].modifyFix.a[j] << "," << modifyFix[spectraBuffer[i].modifyFix.b[j]].modifyName;

				}
				for (int j = 0; j < spectraBuffer[i].modifyVar.a.size(); j++){

					resultOut << "|" << spectraBuffer[i].modifyVar.a[j] << "," << modifyVar[spectraBuffer[i].modifyVar.b[j]].modifyName;

				}
				resultOut << "\t";
			
			}

			resultOut << spectraBuffer[i].top << "\t" << spectraBuffer[i].scoreMax <<endl;
		}


	}

	resultOut.close();

}