#pragma once
class Matrix
{
	double** tab;
	int n;
	int m;
public:
	Matrix(int n, int m);

	void print();

	Matrix* getL();
	Matrix* getU();
	Matrix* getD();

	void setOnes();
	double norm();

	int getN();
	int getM();
	void setValueAt(int i, int j, double value);
	double getValueAt(int i, int j);

	Matrix* operator/(Matrix* B);
	Matrix* operator*(Matrix* B);
	Matrix* operator*(double v);
	Matrix* operator-(Matrix* B);
	Matrix* operator+(Matrix* B);
	void operator=(Matrix* B);

	double* operator[](size_t i);
	~Matrix();
};

