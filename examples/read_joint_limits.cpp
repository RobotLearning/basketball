#include "stdlib.h"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "constants.h"

bool read_joint_limits(double *lb, double *ub) {

	using namespace std;
	int idx;
	string line;
	vector<string> lines;
	string foldername = "/usr/home/cbiwork/workspace/SL_prog/cbUser/config/";
	string name = "SensorOffset.cf";
	string filename = foldername + name;
	ifstream myfile(filename);
	if (myfile.is_open()) {
		while (myfile.good()) {
			getline(myfile,line);
			lines.push_back(line);
		}
	}
	else {
		cout << "Error: cannot open file: " << filename << " !\n";
		return false;
	}
	for (unsigned i = 0; i < lines.size(); i++) {
		istringstream iss(lines[i]);
		for (unsigned j = 0; j < joint_names.size(); j++) {
			idx = lines[i].find(joint_names[j]);
			if (idx != lines[i].npos) { // get the next two doubles
				cout << "Reading joint limits for " << j << endl;
				iss.seekg(idx + 5);
				iss >> lb[j];
				iss >> ub[j];
				break;
			}
		}
	}
	return true;
}

int main() {
	
	double lb[NDOF];
	double ub[NDOF];
	if (read_joint_limits(lb,ub)) {
		for (int i = 0; i < NDOF; i++) {
			std::cout << "Lower limit for joint " << i << " : " << lb[i] << std::endl;
			std::cout << "Upper limit for joint " << i << " : " << ub[i] << std::endl;
			std::cout << std::endl;
		}
	}
	return 1;
}
