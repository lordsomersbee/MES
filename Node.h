#pragma once
class Node
{
public:
	int id;
	double x, y, det_J;
	double* jacobian;
	double* jacobian_invert;
	Node(int, double, double);
	Node();
	~Node();
};

