#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "Grid.h"

using namespace std;

int main() 
{
	double H, L, k, c, ro, alfa, ambient_temp, initial_temp, simulation_time, time_step;
	int nH, nL;
	ifstream file("plik.txt");
	string value;

	getline(file, value);
	H = atof(value.c_str());
	getline(file, value);
	L = atof(value.c_str());
	getline(file, value);
	nH = atoi(value.c_str());
	getline(file, value);
	nL = atoi(value.c_str());
	getline(file, value);
	k = atof(value.c_str());
	getline(file, value);
	c = atof(value.c_str());
	getline(file, value);
	ro = atof(value.c_str());
	getline(file, value);
	alfa = atof(value.c_str());
	getline(file, value);
	ambient_temp = atof(value.c_str());
	getline(file, value);
	initial_temp = atof(value.c_str());
	getline(file, value);
	simulation_time = atof(value.c_str());
	getline(file, value);
	time_step = atof(value.c_str());

	cout << "Time" << "\t" << "MinTemp" << "\t" << "MaxTemp" << endl;

	Grid grid(nH, nL, initial_temp); 
	grid.fulfill_grid(H / (nH - 1), L / (nL - 1));

	for (double time = time_step; time <= simulation_time; time += time_step)
	{
		for (int i = 0; i < nL - 1; i++) {
			for (int j = 0; j < nH - 1; j++) {
				grid.elements[i][j].calculate_matrix_H(k);
				grid.elements[i][j].calculate_matrix_C(c, ro);
				grid.elements[i][j].calculate_BC(alfa, ambient_temp);
			}
		}

		grid.agregate();
		double* t = grid.calculate_temp(time_step);
		grid.clear_data();

		cout << time << "\t" << *min_element(t, t + grid.nH*grid.nL) << "\t" << *max_element(t, t + grid.nH*grid.nL) << endl;
	}

	//grid.elements[0][0].calculate_matrix_H(k);

	//grid.elements[0][0].calculate_matrix_C(c, ro);
	////cout << c << "\t" << ro << endl;
	//grid.elements[0][0].calculate_BC(alfa, ambient_temp);

	//for (int i = 0; i < grid.nL * grid.nH; i++)
	//{
	//	for (int j = 0; j < grid.nL * grid.nH; j++)
	//	{
	//		cout << setprecision(4) << grid.GH[i][j] << "\t";
	//	}
	//	cout << endl;
	//}

	system("pause");
	return 0;
}
