#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<vector>
#include<string>
#include<math.h>
#include<vector>
#include<map>
#include<set>
#pragma once
#include "parsers.h"
using namespace std;
struct edge
{
	int varOnEdge, nextVert;
	edge(int _varOnEdge, int _nextVert)
	{
		varOnEdge = _varOnEdge;
		nextVert = _nextVert;
	}
};

class ClassificationTree
{
	int freeVert;
	int targetFunctionIndex;  //индекс признака, в котором записана целевая функция
	vector<int> vars; //признаки в вершинах
	vector<vector<edge> > g;
	map<int, int> leafVerts;

	double CalculateEntropy(Storage & s, vector<int> & availableObjects)
	{
		vector<int> proportions(s.parInfo[s.parNum - 1], 0); //размер вектора равен числу различных возможных значений целевой функции

		for (int j = 0; j < availableObjects. size(); j++)
			proportions[s.storage[availableObjects[j]].parameters[s.parNum - 1]]++;       //прибавляем 1 в ячейку массива, номером которой является значение целевой функции для данного объекта

		double entropy = 0;
		double all = availableObjects.size();

		for (int j = 0; j < proportions.size(); j++)
			if (proportions[j] != 0)
			{
				double p = proportions[j];
				entropy += (p / all) * (log(p / all) / log(2.0));
			}
			entropy *= -1;

			return entropy;
	}
	vector<vector<int> > SplitByParameter(Storage & s, vector<int> & availableObjects, int targVar)
	{
		vector<vector<int> > splittedObjects(s.parInfo[targVar]);

		for (int j = 0; j < availableObjects.size(); j++)
			splittedObjects[s.storage[availableObjects[j]].parameters[targVar]].push_back(availableObjects[j]);

		return splittedObjects;
	}
	double CalculateIG(Storage & s, vector<int> & availableObjects, int targVar, double H)
	{
		vector<vector<int> > splittedObjects = SplitByParameter(s, availableObjects, targVar);
		vector<int> proportions(s.parInfo[targVar], 0);

		for (int j = 0; j < availableObjects.size(); j++)
			proportions[s.storage[availableObjects[j]].parameters[targVar]]++;

		double IG = H;
		double all = availableObjects.size();
		for (int j = 0; j < s.parInfo[targVar]; j++)
		{
			double p = proportions[j];
			IG -= (p / all) * CalculateEntropy(s, splittedObjects[j]);
		}
		return IG;
	}
	int CreateNewVertex(int var) //возвращает номер новой вершины
	{
		vars.push_back(var);
		g.resize(g.size() + 1);

		return vars.size() - 1;
	}

	int AddVertex(int from, int toVar, int edgeVar) //from - индекс вершины в массиве g, к которой надо приклеивать, toVar - номер признака, который приклеиваем 
	{
		int newv = CreateNewVertex(toVar);

		if(from != -1)
			g[from].push_back(edge(edgeVar, newv));

		return newv; //возвращаем индекс новой вершины
	}
	void AddLeafVertex(int from, int targetFunctionValue, int edgeVar)
	{
		int newv = CreateNewVertex(targetFunctionIndex);

		if(from != -1)
			g[from].push_back(edge(edgeVar, newv));

		//cout << targetFunctionValue << endl;
		leafVerts[newv] = targetFunctionValue;
	}
	bool isBelongToOneClass(Storage & s, vector<int> & availableObjects)
	{
		if (availableObjects.size() == 0)
			return true;

		int sample = s.storage[availableObjects[0]].parameters[targetFunctionIndex];
		for (int j = 1; j < availableObjects.size(); j++)
		{
			if (s.storage[availableObjects[j]].parameters[targetFunctionIndex] != sample)
				return false;
		}
		//cout << "good on " << sample << endl;
		return true;
	}
	void GenerateDOTFile(Storage & s, string file)
	{
		ofstream dot (file);

		dot << "digraph Decisiontree {\n";
		for (int j = 0; j < vars.size(); j++)
		{
			if (g[j].size() != 0)
			{
				dot << j << " [shape=box][label=\"" << s.varNames[vars[j]] << "\"]" << ";\n";           // обработка вершины
			}
			else
			{
				dot << j << " [shape=none][label=\"" << s.assignNames[targetFunctionIndex][leafVerts[j]] << "\"]" << ";\n";      //обработка листа
			}
		}
		for (int j = 0; j < g.size(); j++)
		{
			for (int i = 0; i < g[j].size(); i++)
			{
				dot << j << " -> " << g[j][i].nextVert << " [label=\"" << s.assignNames[vars[j]][g[j][i].varOnEdge] << "\"]" << ";\n"; 
			}
		}
		dot << "}\n";
	}
public:
	ClassificationTree()
	{
		freeVert = 0;
	}

	void PrintCrossValidationTestResults(Storage & testSt, string & testingResultsFile, vector<int> & TP, vector<int> & FP, vector<int> & FN, vector<int> & TN)
	{
		ofstream fout (testingResultsFile);
		for (int j = 0; j < testSt.parInfo[targetFunctionIndex]; j++)
		{
			fout << "Testing results for class " << testSt.assignNames[targetFunctionIndex][j] << ":\n";
			fout << "TP = " << TP[j] << " | FP = " << FP[j] << "\n";
			fout << "FN = " << FN[j] << " | TN = " << TN[j] << "\n";
			fout << "Precision = " << 100 * TP[j] / (TP[j] + FP[j]) << "%\n";
			fout << "Recall = " << 100 * TP[j] / (TP[j] + FN[j]) << "%\n";
			fout << "\n";
		}
		fout.close();
	}
	void CrossValidationTest(Storage & testSt, string & testingResultsFile)
	{
		vector<int> TP(testSt.parInfo[targetFunctionIndex], 0);
		vector<int> FP(testSt.parInfo[targetFunctionIndex], 0);
		vector<int> FN(testSt.parInfo[targetFunctionIndex], 0);
		vector<int> TN(testSt.parInfo[targetFunctionIndex], 0);
		for  (int typeOfClass = 0; typeOfClass < testSt.parInfo[targetFunctionIndex]; typeOfClass++)
		{
			for (int j = 0; j < testSt.objNum; j++)
			{
				int testObjClass = testSt.storage[j].parameters[targetFunctionIndex];  //класс, к которому принадлежит текущий объект
				int predictedClass = CalculateTargetFunction(testSt.storage[j].parameters);

				if (testObjClass == typeOfClass)
				{
					if (predictedClass == typeOfClass)
					{
						TP[typeOfClass]++;
					}
					else
					{
						FN[typeOfClass]++;
					}
				}
				else
				{
					if(predictedClass != typeOfClass)
					{
						TN[typeOfClass]++;
					}
					else
					{
						FP[typeOfClass]++;
					}
				}
			}
		}
		PrintCrossValidationTestResults(testSt, testingResultsFile, TP, FP, FN, TN);
	}

	int CalculateTargetFunction(vector<int> & objParameters)
	{
		int curv = 0;
		while(true)
		{
			if (g[curv].size() == 0) //попали в лист, возвращаем значение в нем
			{
				return leafVerts[curv];
			}
			int currentParVal = objParameters[vars[curv]];
			bool isFound = 0;

			for (int j = 0; j < g[curv].size(); j++)
			{
				if (g[curv][j].varOnEdge == currentParVal)
				{
					curv = g[curv][j].nextVert;
					isFound = 1;
					break;
				}
			}
			if (!isFound) //не нашли подходящее ребро
			{
				return 0; //возвращаем любое значение целевой функции, например 0
			}
		}
	}
	void LearnID3(Storage & s, vector<int> availableObjects, int parentVert, int edgeVar, set<int> freeVars)
	{
		vector<double> IG(freeVars.size()); //information gain
		double H = CalculateEntropy(s, availableObjects); //проверить, вдруг H == 0
		int p = 0, bestVar = -1;
		double bestIG = 0;

		
		if (isBelongToOneClass(s, availableObjects))
		{
			AddLeafVertex(parentVert, s.storage[availableObjects[0]].parameters[targetFunctionIndex], edgeVar);
			return;	
		}
		else
		{
			if (freeVars.size() == 0)
				cout << "Alert\n";
		}
		for (auto it:freeVars)
		{

			double curIG = CalculateIG(s, availableObjects, it, H);
			//cout << H << ' ' << curIG << endl;
			if (curIG > bestIG)
			{
				bestIG = curIG;
				bestVar = it;
			}
		}
		vector<vector<int> > splittedObjects = SplitByParameter(s, availableObjects, bestVar);
	
		freeVars.erase(bestVar);

		int curVertIndex = AddVertex(parentVert, bestVar, edgeVar);

		for (int j = 0; j < splittedObjects.size(); j++)
		{
			if (splittedObjects[j].size() != 0)
			{
				LearnID3(s, splittedObjects[j], curVertIndex, j, freeVars);
			}
		}
		//return bestVar;
	}

	void RunID3Algorithm(string learningFile, string testingFile, string testingResultsFile, string dotFile)
	{
		Storage mainSt;

		ScanLearningSet(learningFile, mainSt); //learningFile = "c:\\coursework\\cars\\test.txt"
		targetFunctionIndex = mainSt.parNum - 1;

		vector<int> availableObjects(mainSt.objNum);
		for (int j = 0; j < mainSt.objNum; j++)
			availableObjects[j] = j;

		set<int> freeVars;
		for (int j = 0; j < mainSt.parNum; j++)
			if (j != targetFunctionIndex)
				freeVars.insert(j);

		cout << freeVars.size() << endl;
		LearnID3(mainSt, availableObjects, -1, -1, freeVars);
		GenerateDOTFile(mainSt, dotFile);

		/*** Тестирование  ***/

		Storage testSt;

		ScanLearningSet(testingFile, testSt);

		CrossValidationTest(testSt, testingResultsFile);
	}
	
};
