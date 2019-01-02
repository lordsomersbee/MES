#include "Element.h"

void Element::setValues(Node a, Node b, Node c, Node d)
{
	this->id = new int*[2];
	this->id[0] = new int[2];
	this->id[1] = new int[2];

	this->id[0][0] = a.id;
	this->id[0][1] = b.id;
	this->id[1][1] = c.id;
	this->id[1][0] = d.id;

	this->nodes = new Node[4];

	this->nodes[0] = a;
	this->nodes[1] = b;
	this->nodes[2] = c;
	this->nodes[3] = d;

	this->H = new double*[4];
	this->H[0] = new double[4];
	this->H[1] = new double[4];
	this->H[2] = new double[4];
	this->H[3] = new double[4];

	this->C = new double*[4];
	this->C[0] = new double[4];
	this->C[1] = new double[4];
	this->C[2] = new double[4];
	this->C[3] = new double[4];

	this->BC = new int[4]{0,0,0,0};

	this->P = new double[4]{ 0,0,0,0 };

	this->l = new double[4];
	for (int i = 0; i < 4; i++)
	{
		l[i] = pow(pow(this->nodes[(i + 1)%4].x - this->nodes[i%4].x, 2) + pow(this->nodes[(i + 1) % 4].y - this->nodes[i % 4].y, 2), (1.0 / 2.0));
		//cout << l[i] << endl;
	}
}

Element::Element()
{
}

Element::~Element()
{
}

void Element::calculate_matrix_H(double k)
{
	double sqrt3 = 1.0 / pow(3.0, (1.0 / 2.0));
	double eta[] = { -sqrt3, -sqrt3, sqrt3, sqrt3 };
	double ksi[] = { -sqrt3, sqrt3, sqrt3, -sqrt3 };
	
	this->eta = new double[4]{ -sqrt3, -sqrt3, sqrt3, sqrt3 };
	this->ksi = new double[4]{ -sqrt3, sqrt3, sqrt3, -sqrt3 };

	//---------------------------------------------------
	// d_N/d_ksi and d_N/d_eta

	const int len = 4;
	double** dN_dksi = new double*[len];
	double** dN_deta = new double*[len];
	for (int i = 0; i < len; i++)
	{
		dN_dksi[i] = new double[len];
		dN_deta[i] = new double[len];

		int factor1, factor2;
		(i == 0 || i == len - 1) ? factor1 = -1 : factor1 = 1;
		(i < 2) ? factor2 = -1 : factor2 = 1;

		for (int j = 0; j < len; j++)
		{
			dN_dksi[i][j] = factor1 * (1.0 / 4.0) * (1 + factor2 * ksi[j]);
			dN_deta[i][j] = factor2 * (1.0 / 4.0) * (1 + factor1 * eta[j]);

			//cout << dN_deta[i][j] << " ";
		}
		//cout << endl;

	}

	//---------------------------------------------------
	// Jacobian

	double** jacobian = new double*[len];
	jacobian[0] = new double[len];
	jacobian[1] = new double[len];
	jacobian[2] = new double[len];
	jacobian[3] = new double[len];

	for (int j = 0; j < len; j++)
	{
		double sums[] = { 0, 0, 0, 0 };
		for (int k = 0; k < len; k++)
		{
			sums[0] += this->nodes[k].x * dN_dksi[k][j];
			sums[1] += this->nodes[k].y * dN_dksi[k][j];
			sums[2] += this->nodes[k].x * dN_deta[k][j];
			sums[3] += this->nodes[k].y * dN_deta[k][j];
		}
		jacobian[0][j] = sums[0];
		jacobian[1][j] = sums[1];
		jacobian[2][j] = sums[2];
		jacobian[3][j] = sums[3];

		//Todooooooo
		this->nodes[0].jacobian = jacobian[0];
		this->nodes[1].jacobian = jacobian[1];
		this->nodes[2].jacobian = jacobian[2];
		this->nodes[3].jacobian = jacobian[3];



	}

	//---------------------------------------------------
	// Det J

	double* det_J = new double[len];
	for (int i = 0; i < len; i++)
	{
		det_J[i] = jacobian[0][i] * jacobian[3][i] - jacobian[1][i] * jacobian[2][i];
		this->nodes[i].det_J = det_J[i];
	}

	//---------------------------------------------------
	// Jacobian^-1
	double** jacobian_invert = new double*[len];
	jacobian_invert[0] = new double[len];
	jacobian_invert[1] = new double[len];
	jacobian_invert[2] = new double[len];
	jacobian_invert[3] = new double[len];

	for (int i = 0; i < len; i++)
	{
		jacobian_invert[0][i] = jacobian[3][i] / det_J[i];
		jacobian_invert[1][i] = -1.0 * jacobian[1][i] / det_J[i];
		jacobian_invert[2][i] = -1.0 * jacobian[2][i] / det_J[i];
		jacobian_invert[3][i] = jacobian[0][i] / det_J[i];
	}

	this->nodes[0].jacobian_invert = jacobian_invert[0];
	this->nodes[1].jacobian_invert = jacobian_invert[1];
	this->nodes[2].jacobian_invert = jacobian_invert[2];
	this->nodes[3].jacobian_invert = jacobian_invert[3];

	//---------------------------------------------------
	// d_N/d_x and d_N/d_y

	double** dN_dx = new double*[len];
	double** dN_dy = new double*[len];
	for (int i = 0; i < len; i++)
	{
		dN_dx[i] = new double[len];
		dN_dy[i] = new double[len];
		for (int j = 0; j < len; j++)
		{
			dN_dx[i][j] = jacobian_invert[0][i] * dN_dksi[j][i] + jacobian_invert[1][i] * dN_deta[j][i];
			dN_dy[i][j] = jacobian_invert[2][i] * dN_dksi[j][i] + jacobian_invert[3][i] * dN_deta[j][i];
		}
	}

	//---------------------------------------------------
	// factors_x - array of matrixes containing "x" factors ( (d_N/d_x)*(d_N/d_x)T * det_J )
	// factors_y - array of matrixes containing "y" factors ( (d_N/d_y)*(d_N/d_y)T * det_J )

	double*** factors_x = new double**[len];
	double*** factors_y = new double**[len];
	for (int factor_iterator = 0; factor_iterator < len; factor_iterator++)
	{
		factors_x[factor_iterator] = new double*[len];
		factors_y[factor_iterator] = new double*[len];

		// actual matrix containing "x" factor and matrix containing "y" factor
		for (int i = 0; i < len; i++)
		{
			factors_x[factor_iterator][i] = new double[len];
			factors_y[factor_iterator][i] = new double[len];
			for (int j = 0; j < len; j++)
			{
				factors_x[factor_iterator][i][j] = dN_dx[factor_iterator][j] * dN_dx[factor_iterator][i] * det_J[factor_iterator];
				factors_y[factor_iterator][i][j] = dN_dy[factor_iterator][j] * dN_dy[factor_iterator][i] * det_J[factor_iterator];
			}
		}
	}

	//---------------------------------------------------
	// factors_sum - array of matrixes containing sum of "x" + "y" factors ( (d_N/d_x)*(d_N/d_x)T * det_J + (d_N/d_y)*(d_N/d_y)T * det_J )

	double*** factors_sum = new double**[len];
	for (int factor_iterator = 0; factor_iterator < len; factor_iterator++)
	{
		factors_sum[factor_iterator] = new double*[len];

		// actual matrix containing sum of "x" + "y" factors
		for (int i = 0; i < len; i++)
		{
			factors_sum[factor_iterator][i] = new double[len];
			for (int j = 0; j < len; j++)
			{
				factors_sum[factor_iterator][i][j] = k * (factors_x[factor_iterator][i][j] + factors_y[factor_iterator][i][j]);
			}
		}
	}

	//---------------------------------------------------
	// Calculate matrix H

	double** h = new double*[len];
	double sum;

	for (int i = 0; i < len; i++)
	{
		h[i] = new double[len];

		for (int j = 0; j < len; j++)
		{
			sum = 0;
			for (int factor_iterator = 0; factor_iterator < len; factor_iterator++)
			{
				sum += factors_sum[factor_iterator][i][j];
			}
			h[i][j] = sum;
		}
	}

	//this->print_matrix(h);
	//cout << endl;

	this->H = h;
}

void Element::calculate_matrix_C(double c, double ro)
{
	int len = 4;

	// N(ksi, eta)

	double** N = new double*[len];
	for (int i = 0; i < len; i++) N[i] = new double[len];

	for (int i = 0; i < len; i++)
	{
		N[0][i] = (1.0 / 4.0) * (1 - ksi[i]) * (1 - eta[i]);
		N[1][i] = (1.0 / 4.0) * (1 + ksi[i]) * (1 - eta[i]);
		N[2][i] = (1.0 / 4.0) * (1 + ksi[i]) * (1 + eta[i]);
		N[3][i] = (1.0 / 4.0) * (1 - ksi[i]) * (1 + eta[i]);
	}

	//---------------------------------------------------
	//factors {Ni}{Ni}*ro*c*det_J

	double*** factors = new double**[len];
	for (int factor_iterator = 0; factor_iterator < len; factor_iterator++)
	{
		factors[factor_iterator] = new double*[len];

		// actual matrix containing "intergral points" factor
		for (int i = 0; i < len; i++)
		{
			factors[factor_iterator][i] = new double[len];
			for (int j = 0; j < len; j++)
			{
				factors[factor_iterator][i][j] = N[factor_iterator][j] * N[factor_iterator][i] * this->nodes[factor_iterator].det_J * ro * c;
			}
		}
	}

	//---------------------------------------------------
	// calculate mtrix C

	double sum = 0;
	for (int i = 0; i < len; i++)
	{
		this->C[i] = new double[len];

		for (int j = 0; j < len; j++)
		{
			sum = 0;
			for (int factor_iterator = 0; factor_iterator < len; factor_iterator++)
			{
				sum += factors[factor_iterator][i][j];
			}
			this->C[i][j] = sum;
		}
	}

	//this->print_matrix(this->C);

	//for (int i = 0; i < len; i++)
	//{
	//	cout << nodes[i].x << " " << nodes[i].y << "\t" << endl;
	//}
}

void Element::calculate_BC(double alfa, double ambient_temp)
{
	if (this->BC[0] == 0 && this->BC[1] == 0 && this->BC[2] == 0 && this->BC[3] == 0) return;

	const int len = 4;

	int edge_ksi[] = { -1, 1, 1, -1 };
	int edge_eta[] = { -1, -1, 1, 1 };

	double* N = new double[len];
	double* N2 = new double[len];

	// Calcualte matrix for every surface

	double*** surfaces = new double**[len];
	double** surfaces_P = new double*[len];

	for (int factor_iterator = 0; factor_iterator < len; factor_iterator++)
	{
		// calculate Ni for 1 and 2 intergral point

		if (factor_iterator % 2)
		{
			N[0] = (1.0 / 4.0) * (1 - edge_ksi[factor_iterator]) * (1 - eta[factor_iterator]);
			N[1] = (1.0 / 4.0) * (1 + edge_ksi[factor_iterator]) * (1 - eta[factor_iterator]);
			N[2] = (1.0 / 4.0) * (1 + edge_ksi[factor_iterator]) * (1 + eta[factor_iterator]);
			N[3] = (1.0 / 4.0) * (1 - edge_ksi[factor_iterator]) * (1 + eta[factor_iterator]);

			N2[0] = (1.0 / 4.0) * (1 - edge_ksi[(factor_iterator + 1) % len]) * (1 - eta[(factor_iterator + 1) % len]);
			N2[1] = (1.0 / 4.0) * (1 + edge_ksi[(factor_iterator + 1) % len]) * (1 - eta[(factor_iterator + 1) % len]);
			N2[2] = (1.0 / 4.0) * (1 + edge_ksi[(factor_iterator + 1) % len]) * (1 + eta[(factor_iterator + 1) % len]);
			N2[3] = (1.0 / 4.0) * (1 - edge_ksi[(factor_iterator + 1) % len]) * (1 + eta[(factor_iterator + 1) % len]);
		}
		else
		{
			N[0] = (1.0 / 4.0) * (1 - ksi[factor_iterator]) * (1 - edge_eta[factor_iterator]);
			N[1] = (1.0 / 4.0) * (1 + ksi[factor_iterator]) * (1 - edge_eta[factor_iterator]);
			N[2] = (1.0 / 4.0) * (1 + ksi[factor_iterator]) * (1 + edge_eta[factor_iterator]);
			N[3] = (1.0 / 4.0) * (1 - ksi[factor_iterator]) * (1 + edge_eta[factor_iterator]);

			N2[0] = (1.0 / 4.0) * (1 - ksi[(factor_iterator + 1) % len]) * (1 - edge_eta[(factor_iterator + 1) % len]);
			N2[1] = (1.0 / 4.0) * (1 + ksi[(factor_iterator + 1) % len]) * (1 - edge_eta[(factor_iterator + 1) % len]);
			N2[2] = (1.0 / 4.0) * (1 + ksi[(factor_iterator + 1) % len]) * (1 + edge_eta[(factor_iterator + 1) % len]);
			N2[3] = (1.0 / 4.0) * (1 - ksi[(factor_iterator + 1) % len]) * (1 + edge_eta[(factor_iterator + 1) % len]);
		}

		/*for (int i = 0; i < len; i++) cout << N[i] << "\t";
		cout << endl;
		for (int i = 0; i < len; i++) cout << N2[i] << "\t";
		cout << endl;*/

		surfaces[factor_iterator] = new double*[len];
		for (int i = 0; i < len; i++)
		{
			surfaces[factor_iterator][i] = new double[len];
			for (int j = 0; j < len; j++)
			{
				surfaces[factor_iterator][i][j] = (alfa * (N[i] * N[j] + N2[i] * N2[j])) * (l[factor_iterator] / 2.0);
			}

		}

		surfaces_P[factor_iterator] = new double[len];
		for (int i = 0; i < len; i++)
		{
			surfaces_P[factor_iterator][i] = (alfa * ((N[i] + N2[i]) * ambient_temp)) * l[factor_iterator] / 2.0;
		}
	}

	double** H = new double*[len];
	double* P = new double[len] {0.0, 0.0, 0.0, 0.0};
	for (int i = 0; i < len; i++)
	{
		H[i] = new double[len] {0.0, 0.0, 0.0, 0.0};
		for (int j = 0; j < len; j++)
		{
			P[i] += this->BC[j] * surfaces_P[j][i];

			//if(this->BC[j] == 1) cout << surfaces_P[j][i] << "\t";

			for (int factor_iterator = 0; factor_iterator < len; factor_iterator++)
			{
				H[i][j] += this->BC[factor_iterator] * surfaces[factor_iterator][i][j];
			}
		}
	}

	this->P = P;
	//this->print_matrix(H);
	//this->print_matrix(this->H);
	//cout << endl;


	for (int i = 0; i < len; i++)
	{
		for (int j = 0; j < len; j++)
		{
			this->H[i][j] += H[i][j];
		}
	}

	//for (int i = 0; i < len; i++)
	//{
	//	this->print_matrix(surfaces[i]);
	//	cout << endl;
	//}

	//for (int i = 0; i < len; i++)
	//{
	//	cout << P[i] << "\t";
	//}
	//cout << endl;

	//for (int i = 0; i < len; i++)
	//{
		//cout << this->nodes[i].x << " " << this->nodes[i].y << "\t";
		//cout << this->BC[i] << "\t";

		/*for (int j = 0; j < len; j++)
		{
			cout << this->nodes[i].jacobian[j] << "\t";
		}
		cout << endl;*/
	//}
	//cout << endl;
}

//------------------
void Element::print_matrix(double** m) 
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << m[i][j] << '\t';
		}
		cout << endl;
	}
}