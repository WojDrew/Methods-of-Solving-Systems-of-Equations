#include "Matrix.h"
#include <math.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <stdlib.h>
using namespace std;
Matrix* jacobi(Matrix* M, Matrix* b);
Matrix* gauss_seidl(Matrix* M, Matrix* b);
Matrix* LU_decomposition(Matrix* M, Matrix* b);
void makeEquation1(Matrix* A, Matrix* b);
void makeEquation2(Matrix* A, Matrix* b);

int c;
int d;
int e;
int f;
int N;

int main() {


	char index[6];
	cin >> index;
	d = index[5] - '0';
	c = index[4] - '0';
	e = index[3] - '0';
	f = index[2] - '0';
	N = 900 + c * 10 + d;
	cout << "d: " << d << endl;
	cout << "c: " << c << endl;
	cout << "e: " << e << endl;
	cout << "f: " << f << endl;
	cout << "N: " << N << endl;
	Matrix* A1 = new Matrix(N,N);
	Matrix* b1 = new Matrix(N,1);

	//ZAD A
	makeEquation1(A1,b1);
	//ZAD B
	std::cout << "Jacobi 1:" << std::endl << "----------------------------------------------" << std::endl;
	Matrix* rJ1 = jacobi(A1, b1);
	rJ1->print();
	std::cout << "Gauss-Seidl 1:" << std::endl << "----------------------------------------------" << std::endl;
	Matrix* rGS1 = gauss_seidl(A1, b1);
	rGS1->print();
	//ZAD C
	/*makeEquation2(A1, b1);
	std::cout << "Jacobi 2:" << std::endl << "----------------------------------------------" << std::endl;
	Matrix* rJ2 = jacobi(A1, b1);
	rJ2->print();
	std::cout << "Gauss-Seidl 2:" << std::endl << "----------------------------------------------" << std::endl;
	Matrix* rGS2 = gauss_seidl(A1, b1);
	rGS2->print();*/
	//ZAD D
	std::cout << "LU 2:" << std::endl << "----------------------------------------------" << std::endl;
	Matrix* rLU2 = LU_decomposition(A1, b1);
	rLU2->print(); 
	Matrix* res = (*(*A1*rLU2) - b1);
	cout << "Norma z res LU: " << res->norm() << endl;
	delete rJ1;
	delete rGS1;
	//delete rJ2;
	//delete rGS2;
	delete rLU2;

	//ZAD E
	int n[] = {100,500,1000,2000,3000,4000,5000,6000};
	ofstream jacobiFile("jacobi.txt");
	ofstream gSFile("gauss_seidl.txt");
	ofstream lUFile("LU.txt");
	for (int i = 0; i < 8; i++)
	{
		N = n[i];
		Matrix* A = new Matrix(N, N);
		Matrix* b = new Matrix(N, 1);
		makeEquation1(A, b);

		clock_t begin_time = clock();
		Matrix* rJ1 = jacobi(A, b);
		jacobiFile << float(clock() - begin_time) / CLOCKS_PER_SEC << endl;

		begin_time = clock();
		Matrix* rGS1 = gauss_seidl(A, b);
		gSFile << float(clock() - begin_time) / CLOCKS_PER_SEC << endl;

		if (i <= 4)
		{
			begin_time = clock();
			Matrix* rLU1 = LU_decomposition(A, b);
			lUFile << float(clock() - begin_time) / CLOCKS_PER_SEC << endl;
			delete rLU1;
		}

		delete rJ1;
		delete rGS1;

		delete A;
		delete b;
	}
	jacobiFile.close();
	gSFile.close();
	lUFile.close();

	return 0;
}

Matrix* jacobi(Matrix* M, Matrix* b)
{
	Matrix *L = M->getL();
	Matrix *U = M->getU();
	Matrix *D = M->getD();

	Matrix* r = new Matrix(M->getN(), 1);
	r->setOnes();
	Matrix* res = new Matrix(M->getN(), 1);
	res = *(*M * r) - b;
	int i = 0;
	while (res->norm() > 1e-9)
	{
		Matrix* Right = (*D) / b;
		Matrix* LU = *L + U;
		Matrix* LUr = *LU * (r);
		Matrix* minusD = (*D)*(-1.0);
		Matrix* Left = *minusD / LUr;
		delete r;
		r = *Left + Right;
		delete Right;
		delete LU;
		delete LUr;
		delete minusD;
		delete Left;
		delete res;
		res = *(*M * r) - b;

		i++;
	}
	delete res;
	std::cout << "Liczba teracji Jacobi: " << i << std::endl;
	return r;
}

Matrix* gauss_seidl(Matrix* M, Matrix* b)
{
	Matrix *L = M->getL();
	Matrix *U = M->getU();
	Matrix *D = M->getD();


	Matrix* r = new Matrix(M->getN(), 1);
	r->setOnes();
	Matrix* res = new Matrix(M->getN(), 1);
	res = *(*M * r) - b;
	int i = 0;
	while (res->norm() > 1e-9)
	{
		Matrix* DplusL = (*D) + L;
		Matrix* minusDplusL = *DplusL * -1;
		Matrix* Right = (*DplusL) / b;
		Matrix* Ur = (*U)*r;
		Matrix* Left = (*minusDplusL) / Ur;
		delete r;
		r = *Left + Right;

		delete DplusL;
		delete minusDplusL;
		delete Right;
		delete Ur;
		delete Left;
		delete res;
		res = *(*M * r) - b;


		i++;
	}
	delete res;
	std::cout << "Liczba teracji Gauss-Seidl: " << i << std::endl;
	return r;
}

Matrix* LU_decomposition(Matrix* M, Matrix* b)
{
	Matrix* U = new Matrix(M->getN(), M->getM());
	Matrix* L = new Matrix(M->getN(), M->getM());
	for (int i = 0; i < M->getN(); i++)
		(*L)[i][i] = 1.0;
	for (int i = 0; i < M->getN(); i++) {
		for (int j = i; j < M->getN(); j++) {
			(*U)[i][j] = (*M)[i][j];
			for (int k = 0; k < i; k++)
			{
				(*U)[i][j] -= (*L)[i][k] * (*U)[k][j];
			}
		}

		for (int j = i + 1; j < M->getN(); j++) {
			(*L)[j][i] = (*M)[j][i];
			for (int k = 0; k < i; k++)
				(*L)[j][i] -= (*L)[j][k] * (*U)[k][i];
			(*L)[j][i] /= (*U)[i][i];
		}
	}

	Matrix* z = *L / b;
	Matrix* x = *U / z;
	delete U;
	delete L;
	delete z;
	return x;
}

void makeEquation1(Matrix* A, Matrix* b)
{
	for (int i = 0; i < A->getN(); i++)
	{
		for (int j = 0; j < A->getM(); j++)
		{
			if (i == j)
			{
				(*A)[i][j] = 5 + e;
			}
			else if (i == j - 1 || j == i - 1 || i == j - 2 || j == i - 2)
			{
				(*A)[i][j] = -1.0;
			}
		}
	}
	for (int i = 0; i < b->getN(); i++)
	{
		for (int j = 0; j < b->getM(); j++)
		{
			(*b)[i][j] = sin(i*(f + 1));
		}
	}
}

void makeEquation2(Matrix* A, Matrix* b)
{
	for (int i = 0; i < A->getN(); i++)
	{
		for (int j = 0; j < A->getM(); j++)
		{
			if (i == j)
			{
				(*A)[i][j] = 3.0;
			}
			if ( i == j - 1 || j == i - 1 || i == j - 2 || j == i - 2)
			{
				(*A)[i][j] = -1.0;
			}
		}
	}
	for (int i = 0; i < b->getN(); i++)
	{
		for (int j = 0; j < b->getM(); j++)
		{
			(*b)[i][j] = sin(i*(f + 1));
		}
	}
}