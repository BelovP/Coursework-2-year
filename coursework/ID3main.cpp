#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<vector>
#include<string>
#include<math.h>
#include<vector>
#include<set>
#include "parsers.h"
using namespace std;
struct edge
{
	int varOnEdge, nextVert;
};

class ClassificationTree
{
	int freeVert;
	vector<vector<edge> > g;
	
	double CalculateEntropy(Storage & s, vector<int> & availableObjects)
	{
		vector<int> proportions(s.parInfo[s.parNum - 1], 0); //размер вектора равен числу различных возможных значений целевой функции

		for (int j = 0; j < availableObjects. size(); j++)
			proportions[s.storage[availableObjects[j]].parameters[s.parNum - 1]];       //прибавл€ем 1 в €чейку массива, номером которой €вл€етс€ значение целевой функции дл€ данного объекта
		
		double entropy = 0;
		
		for (int j = 0; j < proportions.size(); j++)
			if (proportions[j] != 0)
			{
				double p = proportions[j], all = availableObjects.size();
				entropy += (p / all) * (log(p / all) / log(2.0));
			}
		entropy *= -1;

		return entropy;
	}
public:
	ClassificationTree()
	{
		freeVert = 0;
	}

	void LearnID3(Storage & s, vector<int> availableObjects, int targVar, set<int> freeVars)
	{
		vector<double> IG(freeVars.size()); //entropy and information gain
		double H = CalculateEntropy(s, availableObjects); //проверить, вдруг H == 0

		for (auto it:freeVars)
		{
			
		}
	}

	void RunID3Algorithm()
	{
		Storage mainSt;

		ScanLearningSet("c:\\coursework\\cars\\test.txt", mainSt);

	}
};
