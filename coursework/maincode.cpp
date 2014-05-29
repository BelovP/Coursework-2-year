#include<vector>
#include<stdio.h>
#include<iostream>
#include<string>
#include<time.h>
#include<fstream>
#include<set>
#include<queue>
#include<map>
#include "ID3main.h"
#define INF 100000
using namespace std;

#define NUMBER_OF_TESTS 1000
ifstream fin("c:\\coursework\\cars\\test.txt"); //формат входных данных: объекты, далее число параметров каждого признака, далее имена признаков и параметров, далее предпочтения
ifstream testfin("c:\\coursework\\cars\\objects.txt");
ofstream flog("c:\\coursework\\logs\\log.txt");
ofstream fout("c:\\coursework\\output\\output.txt");
ofstream testout("c:\\coursework\\testing\\results.txt");
ofstream objtest("c:\\coursework\\testing\\bigtestresults.txt");
ofstream watchsets ("c:\\coursework\\testing\\sets.txt");

int globalIsInfoLoaded = 0; //переделать


struct Preferences
{
	vector<pair<int, int> > g;
	int vertNum;
};
class Node
{
	vector<int> parameters;
	vector<vector<pair< pair<int, int>, int> > > CPT; //1 == x > !x, firstarg, secondarg, result
	map<vector<int>, int> CPTlines;
	int curln;

	void ConstructCPTlines(int curpos,const vector<int> & parLegend, vector<int> curvec)
	{
		
		for (int j = 0; j < parLegend[curpos]; j++)
		{
			curvec.push_back(j);
			if (curpos < parLegend.size() - 1)
				ConstructCPTlines(curpos + 1, parLegend, curvec);
			else
				CPTlines[curvec] = curln++;
		}
	}
public:
	Node()
	{
		parameters.clear();
		CPT.clear();
	}
	Node(const vector<int> & par, const vector<int> & parLegend)
	{
		parameters.clear();
		CPT.clear();
		parameters = par;
		curln = 0;
		if(parLegend.size() == 0)
		{
			vector<int> tmp;
			CPTlines[tmp] = 0;
			CPT.resize(1);
			return;
		}
		else
		{
			ConstructCPTlines(0, parLegend, vector<int>(0));
		}
		CPT.resize(curln);
		
	}
	/*Node(const vector<int> & par, const vector<int> solution)
	{
		parameters = par;
		CPT.resize(1 << par.size());
		for (int j = 0; j < solution.size(); j++)
			CPT[solution[j] / 2] = 1 - solution[j] % 2;
	}*/
	int FindCycle(int curv, vector<vector<int> > & graph, vector<int> & col)
	{
		col[curv] = 1;
		for(int j = 0; j < graph[curv].size(); j++)
		{
			int to = graph[curv][j];

			if (col[to] == 0)
			{
				if (FindCycle(to, graph, col))
					return true;	
			}
			else if (col[to] == 1)
			{
				return true;
			}
		}
		col[curv] = 2;
		return false;
	}
	bool isCycle()         //checking if there is cycle in CPT lines
	{
		for (map<vector<int>, int>::iterator it = CPTlines.begin(); it != CPTlines.end(); it++)
		{
			int line = it->second;

			vector<vector<int> > graph;                             //building oriented graph of pairwise preferences in CPT line
			for (int i = 0; i < CPT[line].size(); i++)
			{
				int firstv = CPT[line][i].first.second, secondv = CPT[line][i].first.first;		//arguments of pairwise comparison
				int compresult = CPT[line][i].second;												// 1 - >, 0 - <

				if (firstv >= graph.size() || secondv >= graph.size())										 //resizing graph 
					graph.resize(max(firstv + 1, secondv + 1));

				if (compresult)
					graph[firstv].push_back(secondv);
				else
					graph[secondv].push_back(firstv);
			}
			vector<int> col(graph.size(), 0); //vector of vertices' colors for searching cycle

			for (int j = 0; j < graph.size(); j++)                //dfs for searching cycle
				if (col[j] == 0)
					if (FindCycle(j, graph, col))
						return true;
		}
		return false;
	}
	void AddSolution(const vector<int> & par, const vector<int> solution, int firstarg, int secondarg)
	{
		for(int j = 0; j < solution.size(); j++)
		{
			cout << CPT.size() << ' ' << solution[j] / 2 << endl;	
			CPT[solution[j] / 2].push_back(make_pair(make_pair(firstarg, secondarg), 1 - solution[j] % 2));
		}
	}
	int getLns()
	{
		return curln;
	}
	void Print(vector<string> & varNames, vector<vector<string> > & names, int num)
	{
		vector<bool> used(parameters.size(), 0);
			//fout << CPT[j] << endl;
			/*for (int i = 0; i < parameters.size(); i++)
			{
				if (!used[i])
					fout << "not ";
				fout << names[parameters[i]] << " | ";
			}
			if (CPT[j])
			{
				fout << ": " << names[num] << " > not " << names[num] << "\n";
			}
			else
			{
				fout << ": not " << names[num] << " > " << names[num] << "\n";
			}
			for (int i = 0; i < used.size(); i++)
			{
				if (!used[i])
				{
					used[i] = 1;
					break;
				}
				else
					used[i] = 0;
			}
			*/
			for (map<vector<int>, int>::iterator it = CPTlines.begin(); it != CPTlines.end(); it++)
			{
				vector<int> curvec = it->first;
				int line = it->second;
				for (int i = 0; i < curvec.size(); i++)
				{
					fout << varNames[parameters[i]] << " = " << names[parameters[i]][curvec[i]] << " | ";
				}
				fout << ": ";
				for (int i = 0; i < CPT[line].size(); i++)
				{
					if (CPT[line][i].second == 1)
					{
						fout << names[num][CPT[line][i].first.second] << " > " << names[num][CPT[line][i].first.first] << " & ";
					}
					else
					{
						fout << names[num][CPT[line][i].first.second] << " < " << names[num][CPT[line][i].first.first] << " & ";
					}
				}
				fout << endl;
			}
	}
	int CPTlinenum(Object & curobj)
	{
		int j, result = 0;

		/*for (j = 0; j < parameters.size(); j++)
		{
			result *= 2;
			result += curobj.parameters[parameters[j]];
		}
		return result;
		*/
		vector<int> toSearch;
		for (j = 0; j < parameters.size(); j++)
		{
			toSearch.push_back(curobj.parameters[parameters[j]]);
		}
		return CPTlines[toSearch];
	}
	bool CouldBeCompared(Object & first, Object & second)
	{
		for (int j = 0; j < parameters.size(); j++)
		{
			if (first.parameters[parameters[j]] != second.parameters[parameters[j]])
				return false;
		}

		return true;
	}
	int CompareObjects(Object & first, Object & second, int checkedPar) //0 - same, -1  - <, 1 - >
	{
		/*if (first.parameters[checkedPar] == second.parameters[checkedPar])
			return 0;

		int targetLine = CPTlinenum(first);
		if (CPT[targetLine] == 1)
		{
			if (first.parameters[checkedPar])
				return 1;
			else
				return -1;
		}
		else
		{
			if (first.parameters[checkedPar])
				return -1;
			else
				return 1;
		}
		*/
		return 0;//заглушка, потом убрать
	}
};
class CPnet
{
	int VertNum;
	vector<string> parnames;
	vector<Node> CPTs;
	vector<pair<int, int> > g;
    

	void formulatographs(vector<pair<int, int> > & f, vector<vector<int> > & g, vector<vector<int> > & revg)
	{
		for (int j = 0; j < f.size(); j++)
		{
			int from = f[j].first, to = f[j].second;

			g[(from^1)].push_back(to);
			g[(to^1)].push_back(from);
			revg[to].push_back((from^1));
			revg[from].push_back((to^1));
		}
	}
	void dfs(int v, vector<vector<int> > & g, vector<int> & used, vector<int> & order) 
	{
		used[v] = true;
		//cerr << v << endl;
		for (size_t i=0; i<g[v].size(); ++i) 
		{
			int to = g[v][i];
			if (!used[to])
				dfs(to, g, used, order);
		}
		order.push_back (v);
	}
	void dfsrev(int v, int cl, vector<vector<int> > & revg, vector<int> & comp)
	{
		comp[v] = cl;
		for (size_t i=0; i<revg[v].size(); ++i) 
		{
			int to = revg[v][i];
			if (comp[to] == -1)
				dfsrev(to, cl, revg, comp);
		}
	}
	vector<int> SATsolve(int n, vector<pair<int, int> > & f)
	{
		vector<vector<int> > g;
		vector<vector<int> > revg;
		vector<int> used;
		vector<int> order;
		vector<int> comp;

		//flog << "Solving SAT..." << endl;
		used.assign (n, false);
		g.resize(n);
		revg.resize(n);
		formulatographs(f, g, revg);

		for (int i=0; i<n; ++i)
			if (!used[i])
				dfs(i, g, used, order);

		//flog << "Direct dfs completed" << endl;
		comp.assign (n, -1);
		for (int i=0, j=0; i<n; ++i) 
		{
			int v = order[n-i-1];
			if (comp[v] == -1)
				dfsrev(v, j++, revg, comp);
		}
		//flog << "Reverse dfs completed" << endl;
		vector<int> solution;
		for (int i=0; i<n; ++i)
			if (comp[i] == comp[i^1])
			{
				//flog << "2-SAT instance found: false" << endl;
				return solution;
			}

			//flog << "2-SAT instance found: true" << endl;
			for (int i=0; i<n; i += 2)
			{
				int ans = comp[i] > comp[i^1] ? i : i^1;
				solution.push_back(ans);
			}
			//flog << "Solution: ";
			//for (int j = 0; j < solution.size(); j++)
			//	flog << solution[j] << ' ';
			//flog << endl;
			return solution;
	}
	/*int CPTlinenum(Object curobj, vector<int> & u)
	{
		int j, result = 0;

		for (j = 0; j < u.size(); j++)
		{
			result *= 2;
			result += curobj.parameters[u[j]];
		}
		return result;
	}
	*/
	
	bool CreateCPT(int curx, vector<int> & u, Preferences & p, Storage & s, vector<int> & parLegend)
	{
		int j, selectedNum = 0;
		Object objfirst, objsecond;
		map<pair<int, int>, int> testedpairs; //стоит ли переделать на двумерный массив?

		//flog << "Creating CPT..." << endl;
		vector<pair<int, int> > satformula; //if 2x, then !curx > curx
		vector<vector<pair<int, int> > > satFormulas; //разбиваем на  пары и строим формулы для них
		vector<pair<int, int> > pairArgs;
		Node tmp(u, parLegend);
		for (j = 0; j < p.g.size(); j++)
		{
			int line1, line2, formulaID;
			objfirst = s.storage[p.g[j].first];
			objsecond = s.storage[p.g[j].second];
			
			
			if (objfirst.parameters[curx] == objsecond.parameters[curx])
			{
				continue;
			}

			pair<int, int> targpair = make_pair(min(objfirst.parameters[curx], objsecond.parameters[curx]),max(objfirst.parameters[curx], objsecond.parameters[curx]) );
			

			if (testedpairs.find(targpair)==testedpairs.end())
			{
				testedpairs[targpair] = selectedNum;
				formulaID = selectedNum;
				pairArgs.push_back(targpair);
				selectedNum++;
				satFormulas.resize(selectedNum);
			}
			else
				formulaID = testedpairs[targpair];
			
			line1 = tmp.CPTlinenum(objfirst);
			line2 = tmp.CPTlinenum(objsecond);
			
			line1 *= 2; //line *= s.parInfo[j];
			line2 *= 2;
			int flag = 1;
		
			/*if (!objfirst.parameters[curx]) 
				satformula.push_back(make_pair(line1 ^ 1, line2 ^ 1));
			else
				satformula.push_back(make_pair(line1, line2));*/
			if (objfirst.parameters[curx] < objsecond.parameters[curx]) 
				satFormulas[formulaID].push_back(make_pair(line1 ^ 1, line2 ^ 1));
			else
				satFormulas[formulaID].push_back(make_pair(line1, line2));
		}
		//flog << "Formula of " << (1 << (u.size() + 1)) << " variables\n";
		/*for (j = 0; j < satformula.size(); j++)
		{
			if (j != 0)
				flog << " and ";
			flog << "(" << satformula[j].first << " or " << satformula[j].second << ")";
		}
		flog << endl;*/
		int targverts = max(1, tmp.getLns());
	//	cout << targverts << endl;
		for (j = 0; j < satFormulas.size(); j++)
		{
			vector<int> solution = SATsolve(2 * targverts, satFormulas[j]);
			if (solution.size() == 0)
				return false;

			tmp.AddSolution(u, solution, pairArgs[j].first, pairArgs[j].second);
			//Node tmp(u, solution);
		}

		if (tmp.isCycle())             //check if there is a cycle in resulting pairwise order
			return false;

		CPTs[curx] = tmp;
		for (j = 0; j < u.size(); j++)
			g.push_back(make_pair(u[j], curx));
		return true;
	}
	bool GetNextSubset(vector<int> & v1, int k, vector<int> & u, vector<int> & selected)
	{
		bool f = 0;
		int j;

		for (j = 0; j < selected.size(); j++)
		{
			if (selected[j] == 1)
				selected[j] = 0;
			else if (selected[j] == 0)
			{
				selected[j] = 1;
				f = 1;
				break;
			}
		}
		if (!f)
		{
			if (selected.size() == k || selected.size() == v1.size())
				return false;
			selected.push_back(1);
		}
		u.clear();
		for (j = 0; j < selected.size(); j++)
			if (selected[j])
				u.push_back(v1[j]);
		return true;
	}
	int ExtendNetwork(set<int> x, vector<int> & v1, Preferences & p , Storage & s, int k)
	{
		while(x.size() != 0)
		{
			int curx = *(x.begin());
			vector<int> u;
			vector<int> selected;	
			bool CPTcreated = false, subsetExists;
			do 
			{
			/*	flog << "Trying to build new CPT for " << curx << ":\n";
				flog << "Parents: ";
				for (int i = 0; i < u.size(); i++)
					flog << u[i] << ' ';
				flog << endl; */

				if (selected.size() <= curx || (selected.size() > curx && !selected[curx]))
				{
				  vector<int> parLegend;
				  for (int t = 0; t < u.size(); t++)
					  parLegend.push_back(s.parInfo[u[t]]);

				  CPTcreated = CreateCPT(curx, u, p, s, parLegend);
				}
				subsetExists = GetNextSubset(v1, k, u, selected);
			//	flog << "subset sizes: " << u.size() << ' ' << selected.size() << endl;
			} while (!CPTcreated && subsetExists);

			if (CPTcreated)
				return curx;

			x.erase(x.begin());
		}
		return -1;
	}
public:
	int GetVertsNum()
	{
		return CPTs.size();
	}
	vector<vector<int> > GetAdjacencyList()
	{
		vector<vector<int> > tmpg(VertNum);

		for (int j = 0; j < g.size(); j++)
			tmpg[g[j].first].push_back(g[j].second);

		return tmpg;
	}
	Node GetCPT(int num)
	{
		return CPTs[num];
	}
	bool Learn(Storage & mainSt, Preferences & mainP)
	{
		int allowedParents = 0, i, j;

		set<int> v, x;

		VertNum = mainSt.parNum;
		CPTs.resize(mainSt.parNum);

		for (j = 0; j < mainSt.parNum; j++)
			v.insert(j);
		x = v;
		//cfflog << "Building CP-net for " << mainSt.parNum << " parameters..." << endl;
		for (allowedParents = 0; allowedParents != v.size(); allowedParents++)
		{
			int added;
			//flog << "Trying " << allowedParents << " parents...\n";
			do 
			{
				vector<int> dif;
				for (auto it = v.begin(); it != v.end(); it++)
					if (x.find(*it) == x.end())
						dif.push_back(*it);
				added = ExtendNetwork(x, dif, mainP, mainSt, allowedParents);
				if (added != -1)
					x.erase(added);
			} while (added != -1);
			
			if (x.size() == 0)
			{
				//flog << "CP-net learned: true" << endl;
				return true;	
			}

		}
	//	flog << "CP-net learned: false" << endl;
		return false;
	}
	void ScanParametersLegend()
	{
		int cnt;

		ifstream parfin("c:\\coursework\\cars\\parnames.txt");
		parfin >> cnt; //+++++++++++++++++++++++++++++++++++++++++++++++++++++ временно сломано !!!!
		parnames.resize(cnt);
		getline(fin, parnames[0]);
		for (int j = 0; j < cnt; j++)
			getline(fin, parnames[j]);
	}
	void print(Storage & s)
	{
		//flog << "Printing CPnet...\n";
		fout << "CPT vertices legend:\n";
		for (int j = 0; j < parnames.size(); j++)
			fout << j << ": " << parnames[j] << endl;
		fout << "CPTs:\n";
		for (int j = 0; j < CPTs.size(); j++)
		{
			fout << "CPT #" << j << " for " << s.varNames[j] << ":\n";
			CPTs[j].Print(s.varNames, s.assignNames, j);
		}
		fout << "Edges:\n";
		for (int j = 0; j < g.size(); j++)
		{
			fout << parnames[g[j].first] << " -> " << parnames[g[j].second] << endl;
		}
	}
	void GenerateDOTFile()
	{
		ofstream dot ("c:\\coursework\\system\\graph.dot");

		dot << "digraph CP-net {\n";
		for (int j = 0; j < parnames.size(); j++)
			dot << j << " [shape=circle][label=\"" << parnames[j] << "\"]" << ";\n";
		for (int j = 0; j < g.size(); j++)
		{
			dot << g[j].first << " -> " << g[j].second << ";\n"; 
		}
		dot << "}\n";
	}
};
Storage ScanStorage()
{
  Storage In;

  fin >> In.objNum >> In.parNum;
  In.storage.resize(In.objNum);
  for (int j = 0; j < In.objNum; j++)
  {
	  In.storage[j].parameters.resize(In.parNum);
	  for (int i = 0; i < In.parNum; i++)
		  fin >> In.storage[j].parameters[i];
  }
  In.parInfo.resize(In.parNum);
  for (int j = 0; j < In.parNum; j++)
	  fin >> In.parInfo[j];
  In.varNames.resize(In.parNum);
  In.assignNames.resize(In.parNum);
  for(int j = 0; j < In.parNum; j++)
  {
	  string s;
	  fin >> s;
	  //cout << s << ' ';
	  In.varNames[j] = s;
	  In.assignNames[j].resize(In.parInfo[j]);
	  for (int i = 0; i < In.parInfo[j]; i++)
	  {
		  fin >> In.assignNames[j][i];
		 // cout << In.assignNames[j][i] << ' ';
	  }
	//  cout << endl;
  }

  //flog << "Objects loaded";
  return In;
}
Storage CreateStorageFromVec(vector<vector<int> > & vec)
{
	Storage In;

	In.objNum = vec.size();
	In.parNum = vec[0].size();
	In.storage.resize(In.objNum);
	for (int j = 0; j < In.objNum; j++)
	{
		In.storage[j].parameters.resize(In.parNum);
		for (int i = 0; i < In.parNum; i++)
			In.storage[j].parameters[i] = vec[j][i];
	}
	return In;
}
Preferences ScanPreferences(int numObj)
{
	Preferences In;
	int prefNum;

	fin >> prefNum;
	for (int j = 0; j < prefNum; j++)
	{
		int st, fn;

		fin >> st >> fn;
		In.g.push_back(make_pair(st, fn)); //st > fn
	}
	//flog << "Preferences loaded";
	return In;
}
Preferences CreatePreferencesFromVec(vector<pair<int, int> > & vec)
{
	Preferences In;
	
	In.g = vec;
	/*for (int j = 0; j < vec.size(); j++)
	{
		int st, fn;

		fin >> st >> fn;
		In.g.push_back(make_pair(st, fn)); //st > fn
	}*/
	
	return In;
}
class TestManager
{
	Storage objects;
	Preferences pref;
	vector<vector<int> > ParametersGrades;
	int total, equal, agreed, notagreed, testnum;

	void SetVariables()
	{
		total = 0;
		equal = 0;
		agreed = 0;
		notagreed = 0;
	}
	void ProcessSingleTestResult(int result)
	{
		total++;
		if (result == 0)
		{
			equal++;
			notagreed++;
		}
		if (result == 1)
		{
			agreed++;
		}
		if (result == -1)
		{
			notagreed++;
		}
	}
public:
	void PrintTestingResults()
	{
		testout << "Test #" << testnum << ":\n";
		testout << "Tested on " << NUMBER_OF_TESTS << " tests.\n";
		testout << "Agreed preferences: " << agreed << " (" << 100 * agreed / total << "%)\n";
		testout << "Not agreed preferences: " << notagreed << " (" << 100 * notagreed / total << "%)\n";
		testout << "Equal preferences: " << equal << " (" << 100 * equal / total << "%)\n";
		testout << "----------------------------\n";
	}
	void BuildGrades(CPnet & curNet)   //bfs для разделения признаков на слои "важности"
	{
		int i, j;

		vector<vector<int> > CPgraph = curNet.GetAdjacencyList();
		vector<int> Roots(CPgraph.size(), 0);
		queue<int> q;
		vector<int> dist(CPgraph.size(), INF);

		for(j = 0; j < CPgraph.size(); j++)
			for (i = 0; i < CPgraph[j].size(); i++)
				Roots[CPgraph[j][i]]++;

		for (j = 0; j < Roots.size(); j++)
			if (Roots[j] == 0)
			{
				q.push(j);
				dist[j] = 0;
			}

		int maxdist = 0;

		while(!q.empty())
		{
			int curv = q.front();

			q.pop();
			for (j = 0; j < CPgraph[curv].size(); j++)
			{
				if (dist[CPgraph[curv][j]] > dist[curv] + 1)
				{
					dist[CPgraph[curv][j]] = dist[curv] + 1;
					maxdist = max(maxdist, dist[CPgraph[curv][j]]);
					q.push(CPgraph[curv][j]);
				}
			}
		}

		ParametersGrades.resize(maxdist + 1);
		for (j = 0; j < dist.size(); j++)
		{
			ParametersGrades[dist[j]].push_back(j);
		}
	}
	void LoadInfo()
	{
		globalIsInfoLoaded = 1;
		fin = ifstream("c:\\coursework\\cars\\datavar.txt");

		objects = ScanStorage();
		pref = ScanPreferences(objects.objNum);

		fin = ifstream("c:\\coursework\\cars\\databin.txt");
	}
	int TestObjectPair(int firstobj, int secondobj, CPnet & curNet) //0 - =, 1 - >, -1 - <
	{
		int i, j;

		for (j = 0; j < ParametersGrades.size(); j++)
		{
			for (i = 0; i < ParametersGrades[j].size(); i++)
			{

				int curPar = ParametersGrades[j][i];

				/*if (objects.storage[firstobj].parameters[curPar] == objects.storage[secondobj].parameters[curPar])
					continue;*/

				Node activeCPT = curNet.GetCPT(curPar);
				if (activeCPT.CouldBeCompared(objects.storage[firstobj], objects.storage[secondobj]))
				{
					int tmpres = activeCPT.CompareObjects(objects.storage[firstobj], objects.storage[secondobj], curPar);

					if (tmpres != 0)
						return tmpres;
				}

			}
		}
		return 0;
	}
	void TestOnLoadedInfo(CPnet & curNet)
	{

		for (int test = 0; test < NUMBER_OF_TESTS; test++)
		{
			int firstobj, secondobj;
			int testedPref = rand() % pref.g.size();
			
			firstobj = pref.g[testedPref].first;
			secondobj = pref.g[testedPref].second;
			ProcessSingleTestResult(TestObjectPair(firstobj, secondobj, curNet));
		}
		testnum++;
	}
	/*void FullTest(CPnet & curNet)
	{
		for (j = 0; j < 50)
	}*/
	void AddDataNodeToFile(int objnum, int percentage)
	{
		objtest << objnum << ' ' << percentage << endl; 
	}
	void TestCPnet(CPnet & curNet)
	{
		testnum = 0;
		//SetVariables();
		if (!globalIsInfoLoaded) //сделать нормально!
		  LoadInfo();

		BuildGrades(curNet);
		//cerr << "Info loaded...\n";
		int aggr = 0;
		for (int j = 0; j < 10; j++)
		{
		//	cerr << "Test #" << j << endl;

			SetVariables();
			TestOnLoadedInfo(curNet);
			aggr += agreed * 100 / total;
		//	PrintTestingResults();
		}
		AddDataNodeToFile(curNet.GetVertsNum(), aggr / 10);
	}
	void FormDataVectors(vector<vector<int> > & setofobj, vector<pair<int, int> > & pr, int objnum)
	{
		set<int> selected;
		map<int, int> numtie;
		setofobj.resize(objnum);
		for(int j = 0; j < objnum; j++)
		{
			int targ = rand() % objects.objNum;

			selected.insert(targ);
			numtie[targ] = j;
			setofobj[j].resize(objects.parNum);
			for (int i = 0; i < objects.parNum; i++)
			{
				setofobj[j][i] = objects.storage[targ].parameters[i];
			}
		}
		for (int j = 0; j < pref.g.size(); j++)
		{
			if(selected.find(pref.g[j].first) != selected.end() && selected.find(pref.g[j].second) != selected.end())
				pr.push_back(make_pair(numtie[pref.g[j].first], numtie[pref.g[j].second]));
		}
		for (auto it = selected.begin(); it != selected.end(); it++)

		{
			watchsets << *it << ' ';
		}
		watchsets << endl;
	}
	bool CreateDataNode(int objnum)
	{
		Storage st;
		Preferences p;
		vector<vector<int> > setofobj;
		vector<pair<int, int> > pr;

		FormDataVectors(setofobj, pr, objnum);
		st = CreateStorageFromVec(setofobj);
		p = CreatePreferencesFromVec(pr);

		CPnet nodeNet;
		nodeNet.ScanParametersLegend();
		if(nodeNet.Learn(st, p))
		{
			TestCPnet(nodeNet);

//			AddDataNodeToFile(objnum);
			return true;
		}
		return false;
	}
	void ObjectsSetsTest(int minobj, int maxobj)
	{
		int i, j, step = 3;
		LoadInfo();
		for(int objnum = minobj; objnum <= maxobj; objnum+=step)
		{
			while(!CreateDataNode(objnum))
			{

			}

			cerr << objnum << endl;
		}
	}
};
int main()
{
	srand(time(0));

	//TestManager BigTest;
//	BigTest.ObjectsSetsTest(10, 60);
	//exit(0);



	/*** код для CP-net ***/

	/*Storage mainSt = ScanStorage();
	Preferences mainP = ScanPreferences(mainSt.objNum);
	
	CPnet mainNet;
//	mainNet.ScanParametersLegend();

	if(mainNet.Learn(mainSt, mainP))
	{
		mainNet.print(mainSt);
		mainNet.GenerateDOTFile();

		TestManager MainTest;
		MainTest.TestCPnet(mainNet);
	}
	else
		fout << "No CPnet";	
	*/

	/*** завершение кода для CP-net ***/

	ClassificationTree DecisionTree;
	DecisionTree.RunID3Algorithm(string("c:/coursework/ID3/learningdata.txt"), string("c:/coursework/ID3/testingdata.txt"), string("c:/coursework/ID3/testingresults.txt"), string("c:/coursework/ID3/dot/decisiontree.dot"));
	//system("PAUSE");
}