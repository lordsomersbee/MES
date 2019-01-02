#include "Grid.h"



Grid::Grid(int nH, int nL, double initial_temp)
{
	this->nH = nH;
	this->nL = nL;
	this->nodes = new Node*[nL];
	for (int i = 0; i < nL; i++)
	{
		this->nodes[i] = new Node[nH];
	}

	this->elements = new Element*[nL-1];
	for (int i = 0; i < nL-1; i++)
	{
		this->elements[i] = new Element[nH-1];
	}

	this->t = new double[this->nL * this->nH];
	for (int i = 0; i < this->nL * this->nH; i++) t[i] = initial_temp;
}

Grid::~Grid()
{
}

void Grid::fulfill_grid(double H, double L) 
{
	for (int i = 0; i < this->nL; i++)
	{
		for (int j = 0; j < this->nH; j++)
		{
			this->nodes[i][j].id = i * this->nH + j;
			this->nodes[i][j].x = i * L;
			this->nodes[i][j].y = j * H;

			//cout << this->nodes[i][j].id << " (" << this->nodes[i][j].x << ", " << this->nodes[i][j].y << ") ";
		}
		//cout << endl;
	}

	for (int i = 0; i < this->nL - 1; i++)
	{
		for (int j = 0; j < this->nH - 1; j++)
		{
			this->elements[i][j].setValues(this->nodes[i][j], this->nodes[i + 1][j], this->nodes[i + 1][j + 1], this->nodes[i][j + 1]);

			// Set edges
			if (i == 0) this->elements[i][j].BC[3] = 1;
			if (i == this->nL - 2) this->elements[i][j].BC[1] = 1;
			if (j == 0) this->elements[i][j].BC[0] = 1;
			if (j == this->nH - 2) this->elements[i][j].BC[2] = 1;

			//cout << this->elements[i][j].id[0][0] << " " << this->elements[i][j].id[0][1] << " "
			//<< this->elements[i][j].id[1][1] << " " << this->elements[i][j].id[1][0] << " "
			//<< this->elements[i][j].BC[0] << this->elements[i][j].BC[1] << this->elements[i][j].BC[2] << this->elements[i][j].BC[3] << endl;
		}
	}
}

void Grid::agregate()
{
	const int len = 4;

	double** GH = new double*[this->nL * this->nH];
	double** GC = new double*[this->nL * this->nH];
	double* GP = new double[this->nL * this->nH];

	for (int i = 0; i < this->nL * this->nH; i++)
	{
		GP[i] = 0;

		GH[i] = new double[this->nL * this->nH];
		GC[i] = new double[this->nL * this->nH];
		for (int j = 0; j < this->nL * this->nH; j++)
		{
			GH[i][j] = 0;
			GC[i][j] = 0;
		}
	}

	for (int i = 0; i < this->nL-1; i++)
	{
		for (int j = 0; j < this->nH - 1; j++)
		{
			for (int k = 0; k < len; k++)
			{
				int index = this->elements[i][j].nodes[k].id;
				GP[index] += this->elements[i][j].P[k];
				//cout << this->elements[i][j].P[k] << endl;

				for (int l = 0; l < len; l++)
				{
					int index1 = this->elements[i][j].nodes[l].id;
					int index2 = this->elements[i][j].nodes[k].id;

					GH[index1][index2] += this->elements[i][j].H[k][l];
					GC[index1][index2] += this->elements[i][j].C[k][l];
					//cout << this->elements[i][j].H[k][l] << endl;
				}
			}
			//cout << endl;
		}
	}

	//for (int i = 0; i < this->nL * this->nH; i++)
	//{
	//	cout << std::setprecision(7) << GP[i] << endl;
	//}

	/*for (int i = 0; i < this->nL * this->nH; i++)
	{
		for (int j = 0; j < this->nL * this->nH; j++)
		{
			cout << std::setprecision(5) << GC[i][j] << "\t";
		}
		cout << endl;
	}*/

	this->GH = GH;
	this->GC = GC;
	this->GP = GP;
}

double* Grid::calculate_temp(double time_step)
{
	double** factor_a = new double*[this->nL * this->nH];
	double* factor_b = new double[this->nL * this->nH];
	for (int i = 0; i < this->nL * this->nH; i++)
	{
		factor_b[i] = this->GP[i];
		
		factor_a[i] = new double[this->nL * this->nH];
		for (int j = 0; j < this->nL * this->nH; j++)
		{
			factor_b[i] += (this->GC[i][j] / time_step) * this->t[j];
			factor_a[i][j] = this->GH[i][j] + (this->GC[i][j] / time_step);
		}
	}

	this->t = gauss_elimination(factor_a, factor_b);

	//for (int i = 0; i < this->nL * this->nH; i++)
	//{
	//	cout << std::setprecision(8) << factor_b[i] << "\t";
	//	for (int j = 0; j < this->nL * this->nH; j++)
	//	{
	//		cout << std::setprecision(4) << factor_a[i][j] << "\t";
	//	}
	//	cout << endl;
	//}
	//cout << endl;

	return this->t;
}

void Grid::clear_data()
{
	int len = 4;
	for (int i = 0; i < this->nL - 1; i++)
	{
		for (int j = 0; j < this->nH - 1; j++)
		{
			for (int k = 0; k < len; k++)
			{
				this->elements[i][j].P[k] = 0;
				for (int l = 0; l < len; l++)
				{
					this->elements[i][j].H[k][l] = 0;
					this->elements[i][j].C[k][l] = 0;
				}
			}
		}
	}

	for (int i = 0; i < this->nL * this->nH; i++)
	{
		this->GP[i] = 0;

		for (int j = 0; j < this->nL * this->nH; j++)
		{
			this->GH[i][j] = 0;
			this->GC[i][j] = 0;
		}
	}
}

//--------------------------------------------------------------------------
double* Grid::gauss_elimination(double** L, double* P)
{
	double **AB = new double*[this->nL * this->nH];
	for (int i = 0; i < this->nL * this->nH; i++) AB[i] = new double[this->nL * this->nH + 1];
	
	for (int i = 0; i < this->nL * this->nH; i++)
	{
		for (int j = 0; j <= this->nL * this->nH; j++)
		{
			AB[i][j] = L[i][j];
		}

		AB[i][this->nL * this->nH] = P[i];
	}

	double m, s;
	double* X = new double[this->nL * this->nH];

	for (int i = 0; i < this->nL * this->nH - 1; i++)
	{
		for (int j = i + 1; j < this->nL * this->nH; j++)
		{
			//if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (int k = i + 1; k <= this->nL * this->nH; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych

	for (int i = this->nL * this->nH - 1; i >= 0; i--)
	{
		s = AB[i][this->nL * this->nH];
		for (int j = this->nL * this->nH - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		//if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}

	/*for (int i = 0; i < this->nL * this->nH; i++)
	{
		for (int j = 0; j <= this->nL * this->nH; j++)
		{
			cout << std::setprecision(4) << AB[i][j] << "\t";
		}
		cout << endl;
	}*/

	return X;
}