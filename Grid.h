#pragma once
#include "Node.h"
#include "Element.h"
#include <iomanip>

class Grid
{
public:
	int nH, nL;
	Node** nodes;
	Element** elements;
	double** GH;
	double** GC;
	double* GP;
	double* t;
	
	Grid(int, int, double);
	~Grid();
	void fulfill_grid(double, double);
	void agregate();
	double* calculate_temp(double);
	void clear_data();

private:
	double* gauss_elimination(double**, double*);
};

