#include "deformationTransfer.hh"

vector<CorrespondanceEntry> corresList;
unordered_map<int, list<int>> corresMap;

void LoadCorrespondanceList(char* filename)
{
	string line;
	ifstream file;
	file.open(filename);
	while(!file.eof())
	{
		CorrespondanceEntry corres;
		corres.distance = 0;
		getline(file, line);
		sscanf(line.c_str(),"%d %d",&corres.src_triangle_id, &corres.target_triangle_id);
		corresList.push_back(corres);
		unordered_map<int,list<int>>::iterator it = corresMap.find(corres.src_triangle_id);
		if(it == corresMap.end())
		{
			list<int> l;
			l.push_back(corres.target_triangle_id);
			corresMap.insert(make_pair<int, list<int>>(corres.src_triangle_id, l));
		}
		else
		{
			list<int> l = it->second;
			l.push_back(corres.target_triangle_id);
			corresMap.insert(make_pair<int, list<int>>(corres.src_triangle_id, l));
		}
	}
	file.close();
}

void PerformDeformationTransfer()
{
	Initialize();
	cholmod_dense * deformation;
	targetWindow->DTransfer();
	Deinitialize();
}