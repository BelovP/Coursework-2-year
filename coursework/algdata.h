#include<vector>
#include<string>
#pragma once
using namespace std;
struct Object
{
	vector<int> parameters;
};
struct Storage
{
	vector<Object> storage;
	vector<int> parInfo;
	vector<string> varNames;
	vector<vector<string> > assignNames; 
	int objNum, parNum;
};
struct SimpleStorage
{
	vector<vector<string> > storage;
	int objNum, parNum;
};