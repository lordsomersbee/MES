#pragma once
#include "Node.h"
#include <cmath>
#include <iostream>

using namespace std;

class Element
{
public:
	int** id;
	Node* nodes;
	double** H;
	double** C;
	double* ksi;
	double* eta;
	int* BC;
	double* P;
	double* l;

	Element();
	~Element();
	void setValues(Node, Node, Node, Node);
	void calculate_matrix_H(double);
	void calculate_matrix_C(double, double);
	void calculate_BC(double, double);

	//------
	void print_matrix(double**);
};

