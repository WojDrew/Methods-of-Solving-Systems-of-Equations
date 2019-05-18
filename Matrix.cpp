#include "Matrix.h"
#include <iostream>
#include <math.h>


Matrix::Matrix(int n, int m)
{
	this->n = n;
	this->m = m;
	this->tab = new double*[n];
	for (int i = 0; i < n; i++)
		this->tab[i] = new double[m];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			this->tab[i][j] = 0.0;
}

void Matrix::print()
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			std::cout << this->tab[i][j] << " ";
		}
		std::cout << std::endl;
	}

}

Matrix* Matrix::getL()
{
	Matrix* wsk = new Matrix(this->n, this->m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if(i > j)
				(*wsk)[i][j] = this->tab[i][j];
		}
	}
	return wsk;
}

Matrix* Matrix::getU()
{
	Matrix* wsk = new Matrix(this->n, this->m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i < j)
				(*wsk)[i][j] = this->tab[i][j];
		}
	}
	return wsk;
}

Matrix* Matrix::getD()
{
	Matrix* wsk = new Matrix(this->n, this->m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j)
				(*wsk)[i][j] = this->tab[i][j];
		}
	}
	return wsk;
}

void Matrix::setOnes()
{
	for (int i = 0; i < this->n; i++)
		for (int j = 0; j < this->m; j++)
			this->tab[i][j] = 1;
}

double Matrix::norm() {
	if (this->m == 1)
	{
		double sum = 0;
		for (int i = 0; i < this->n; i++)
			sum += pow(this->tab[i][0],2);
		return sqrt(sum);
	}
	else
		return NULL;
	
}

int Matrix::getN()
{
	return this->n;
}
int Matrix::getM()
{
	return this->m;
}

void Matrix::setValueAt(int i, int j, double value)
{
	this->tab[i][j] = value;
}

double Matrix::getValueAt(int i, int j)
{
	return this->tab[i][j];
}

Matrix* Matrix::operator/(Matrix* B)
{
	if (this->n == this->m)
	{
		Matrix* wsk = new Matrix(this->n, 1);
		double* x = new double[this->n];
		if (this->tab[1][1] != 0 && this->tab[2][1] == 0)
		{
			for (int i = this->n - 1; i >= 0; i--)
			{
				double a = 0;
				for (int j = this->n; j > i; j--)
				{
					a += this->tab[i][j] * x[j];
				}
				x[i] = ((*B)[i][0] - a) / this->tab[i][i];
			}

			for (int i = this->n - 1; i >= 0; i--)
				(*wsk)[i][0] = x[i];
		}
		else
		{
			for (int i = 0; i < this->n; i++)
			{
				double a = 0;
				for (int j = 0; j < i; j++)
				{
					a += this->tab[i][j] * x[j];
				}
				x[i] = ((*B)[i][0] - a) / this->tab[i][i];
			}
			for (int i = 0; i < this->n; i++)
				(*wsk)[i][0] = x[i];
		}

		//delete x;
		return wsk;
	}
}

Matrix* Matrix::operator*(Matrix* B) 
{
	if (this->m == B->getN())
	{
		Matrix* wsk = new Matrix(this->n, B->getM());
		for (int i = 0; i < this->n; i++)
		{
			for (int j = 0; j < B->getM(); j++)
			{
				double sum = 0;
				for (int k = 0; k < this->m; k++)
					sum += this->tab[i][k] * (*B)[k][j];
				(*wsk)[i][j] =sum;
			}
		}
		return wsk;
	}
	else
		return nullptr;
}

Matrix* Matrix::operator*(double v)
{

	Matrix* wsk = new Matrix(n, m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
				(*wsk)[i][j] = this->tab[i][j]*v;
		}
	}
	return wsk;
}

Matrix* Matrix::operator-(Matrix* B)
{
	if (this->n == B->getN() && this->m == B->getM())
	{
		Matrix* wsk = new Matrix(n, m);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				(*wsk)[i][j] = this->tab[i][j] - B->getValueAt(i, j);
			}
		}
		return wsk;
	}
	else
		return nullptr;
}
Matrix* Matrix::operator+(Matrix* B)
{
	if (this->n == B->getN() && this->m == B->getM())
	{
		Matrix* wsk = new Matrix(n, m);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				(*wsk)[i][j] = this->tab[i][j] + B->getValueAt(i, j);
			}
		}
		return wsk;
	}
	else
		return nullptr;
}

void Matrix::operator=(Matrix* B)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			this->tab[i][j] = (*B)[i][j];
		}
	}
}

Matrix::~Matrix()
{
	for (int i = 0; i < n; i++)
		delete this->tab[i];
	delete this->tab;
}

double* Matrix::operator[](size_t i)
{
	return this->tab[i];
}

