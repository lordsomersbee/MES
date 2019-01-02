#include "Node.h"



Node::Node(int id, double x, double y)
{
	this->id = id;
	this->x = x;
	this->y = y;
	this->jacobian = new double[4];
	this->jacobian_invert = new double[4];
	this->det_J = 0;
}

Node::Node()
{
}

Node::~Node()
{
}
