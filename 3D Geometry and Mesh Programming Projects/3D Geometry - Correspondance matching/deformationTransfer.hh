#include "DeformationViewer.hh"
#include "correspondance.hh"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <list>
using namespace std;

extern vector<CorrespondanceEntry> corresList;
extern unordered_map<int, list<int>> corresMap;
void LoadCorrespondanceList(char* filename);
void PerformDeformationTransfer();