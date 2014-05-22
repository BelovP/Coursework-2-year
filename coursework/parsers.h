#include <iostream>
#include <fstream>
#include <string>
#pragma once
#include "algdata.h"
using namespace std;
void ScanLearningSet(string file, Storage & s) //последний параметр - значение целевой функции для данного объекта
{
	int i, j;
	ifstream fin (file);
	
	fin >> s.objNum >> s.parNum;
	s.storage.resize(s.objNum);
	for (j = 0; j < s.objNum; j++)
	{
		s.storage[j].parameters.resize(s.parNum);
		for (i = 0; i < s.parNum; i++)
		{
			fin >> s.storage[j].parameters[i];
		}
	}
	s.parInfo.resize(s.parNum);
	for (j = 0; j < s.parNum; j++)
		fin >> s.parInfo[j];
	s.varNames.resize(s.parNum);
	s.assignNames.resize(s.parNum);
	for (j = 0; j < s.parNum; j++)
	{
		fin >> s.varNames[j];
		s.assignNames[j].resize(s.parInfo[j]);
		for (i = 0; i < s.parInfo[j]; i++)
			fin >> s.assignNames[j][i];

	}
	fin.close();
}