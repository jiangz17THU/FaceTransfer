#pragma once
class Matrix
{
public:
	Matrix();//Ĭ�Ϲ��캯��
	//Matrix();
	//�������Ĺ��캯��
	Matrix(const int&, const int&);
	Matrix(const int&, const int&, double[71][2]);
	Matrix(const int&, const int&, double[71][71]);
	//Matrix(const int&, const int&, double[1][4]);
	Matrix(const int&, const int&, double[4][4]);
	//���ƹ��캯��
	Matrix(const Matrix&);
	//��ֵ������
	Matrix& operator=(const Matrix&);
	Matrix& operator+=(const Matrix&);
	Matrix& operator-=(const Matrix&);
	//��������
	~Matrix();

	//��Ȳ�����
	friend bool operator==(const Matrix&, const Matrix&);
	//���Ȳ�����
	friend bool operator!=(const Matrix& lhs, const Matrix& rhs)
	{
		return !(lhs == rhs);
	}
	//�±������
	double* operator[] (int index);
	//�±������
	double* operator[] (int index)const;
	//��ӡ���
	void print()const;
	//�������������

	// *(�������)
	friend Matrix operator*(const Matrix&, const Matrix&);
	// *(ʵ��*����)
	friend Matrix operator*(const double&, const Matrix&);
	// *(����*ʵ��)
	friend Matrix operator*(const Matrix& M, const double& c)
	{
		return c * M;
	}

	//ת��
	friend Matrix trv(const Matrix&);
	//������
	friend Matrix inv(Matrix);
	//�Ծ�������
	void modify();

	double* p;//�������
	int row;//����
	int col;//����
	int n;//Ԫ�ظ���
};
// +
Matrix operator+(const Matrix&, const Matrix&);
// -
Matrix operator-(const Matrix&, const Matrix&);

#include "cmath"
#include <iostream>
using namespace std;

//start Ĭ�Ϲ��캯��
Matrix::Matrix()
{
	this->row = 0; this->col = 0;
	this->n = 0;   this->p = 0;
}//end Ĭ�Ϲ��캯��

//start �������Ĺ��캯��
Matrix::Matrix(const int& r, const int& c)
{
	this->row = r; this->col = c; this->n = r * c;
	this->p = new double[n];
	for (int i = 0; i < n; ++i)
	{
		this->p[i] = 0.0;//Ĭ��ֵ��Ϊ0
	}
}

Matrix::Matrix(const int& r, const int& c, double p[71][2])
{
	this->row = r; this->col = c; this->n = r * c;
	this->p = new double[n];
	for (int i = 0; i < r; ++i)
	{
		for(int j=0;j<c;j++)
		this->p[i*c+j] = p[i][j];//��p��ֵ
	}
	this->modify();
}//end �������Ĺ��캯��

Matrix::Matrix(const int& r, const int& c, double p[71][71])
{
	this->row = r; this->col = c; this->n = r * c;
	this->p = new double[n];
	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; j++)
			this->p[i * c + j] = p[i][j];//��p��ֵ
	}
	this->modify();
}//end �������Ĺ��캯��

Matrix::Matrix(const int& r, const int& c, double p[4][4])
{
	this->row = r; this->col = c; this->n = r * c;
	this->p = new double[n];
	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; j++)
			this->p[i * c + j] = p[i][j];//��p��ֵ
	}
	this->modify();
}//end �������Ĺ��캯��

//start ���ƹ��캯��
Matrix::Matrix(const Matrix& other)
{
	this->row = other.row; this->col = other.col;
	this->n = row * col;
	this->p = new double[n];
	for (int i = 0; i < n; ++i)
	{
		this->p[i] = other.p[i];
	}
	this->modify();
}//end ���ƹ��캯��

//start ��ֵ������
Matrix& Matrix::operator=(const Matrix& other)
{
	if (*this == other) return *this;
	this->row = other.row;
	this->col = other.col;
	this->n = row * col;
	delete[] this->p;
	this->p = new double[n];
	for (int i = 0; i < n; ++i)
	{
		this->p[i] = other.p[i];
	}
	this->modify();
	return *this;
}//end ��ֵ������

//start +=����
Matrix& Matrix::operator+=(const Matrix& other)
{
	if (this->row != other.row
		|| this->col != other.col)
	{
		cout << "������ȵľ�����ܽ���\"+=\"������" << endl;
		return *this;
	}
	for (int i = 0; i < this->n; ++i)
	{
		this->p[i] += other.p[i];
	}
	this->modify();
	return *this;
}//end +=����
//start -=����
Matrix& Matrix::operator-=(const Matrix& other)
{
	if (this->row != other.row
		|| this->col != other.col)
	{
		cout << "������ȵľ�����ܽ���\"-=\"������" << endl;
		return *this;
	}
	for (int i = 0; i < this->n; ++i)
	{
		this->p[i] -= other.p[i];
	}
	this->modify();
	return *this;

}//end -=����

//start +����
Matrix operator+(const Matrix& lhs, const Matrix& rhs)
{
	Matrix lhsTemp = lhs, rhsTemp = rhs;
	return lhsTemp += rhsTemp;
}//end +����
//start -���� 
Matrix operator-(const Matrix& lhs, const Matrix& rhs)
{
	Matrix lhsTemp = lhs, rhsTemp = rhs;
	return lhsTemp -= rhsTemp;
}//end -����

//start *����(�������)
Matrix operator*(const Matrix& lhs, const Matrix& rhs)
{
	if (rhs.row != lhs.col)
	{
		cout << "��˾��󲻷���������" << endl;
		return Matrix();
	}
	Matrix temp(lhs.row, rhs.col);
	for (int r = 0; r < temp.row; ++r)
		for (int c = 0; c < temp.col; ++c)
		{
			for (int k = 0; k < lhs.col; ++k)
				temp[r][c] += lhs[r][k] * rhs[k][c];
		}
	return temp;
}//end *����(�������)
// start *����(ʵ����������)
Matrix operator*(const double& c, const Matrix& M)
{
	Matrix temp(M.row, M.col);
	for (int i = 0; i < temp.row * temp.col; ++i)
	{
		temp.p[i] = M.p[i] * c;
	}
	return temp;
}// end *����(ʵ����������)
//start ת��
Matrix trv(const Matrix& M)
{
	Matrix temp(M.col, M.row);
	for (int r = 0; r < M.row; ++r)
		for (int c = 0; c < M.col; ++c)
			temp[c][r] = M[r][c];
	return temp;
}//end ת��
//start ��ӡ���
void Matrix::print()const
{
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
			cout << (*this)[i][j] << "  ";
		cout << endl;
	}
}//end ��ӡ���

//start ��Ȳ�����(��Ԫ����)
bool operator==(const Matrix& lhs, const Matrix& rhs)
{
	if (lhs.row != rhs.row)
		return false;
	if (lhs.row != rhs.col)
		return false;
	for (int i = 0; i < lhs.row * lhs.col; ++i)
	{
		if (lhs.p[i] != rhs.p[i])
			return false;
	}
	return true;
}//end ��Ȳ�����(��Ԫ����)

//start �±������
double* Matrix::operator[](int indexRow)const
{
	return this->p + indexRow * col;
}//end �±������
//start �±������
double* Matrix::operator[](int indexRow)
{
	return this->p + indexRow * col;
}//end �±������



//start ������
Matrix inv(Matrix M)
{
	//start �жϾ����Ƿ����
	if (M.row != M.col)
	{
		cout << "ֻ�з���������棡" << endl;
		return Matrix();
	}
	//end �жϾ����Ƿ����
	//start ���쵥λ��
	Matrix E(M.row, M.col);
	for (int i = 0; i < M.row; ++i)
		E[i][i] = 1.0;
	//end ���쵥λ��
	//start ������Ϊ������
	for (int c = 0, r = 0; c < M.col && r < M.row; ++c, ++r)
	{
		int ind;//��¼��Ԫ�к�
		//start Ѱ����Ԫ
		for (int i = r; i < M.row; ++i)
		{
			if (M[i][c] != 0)
			{
				ind = i;
				break;
			}
			if (i == M.row - 1)
			{
				cout << "���󲻿��棡" << endl;
				return Matrix();
			}
		}
		//end Ѱ����Ԫ
		//start ��������
		for (int j = 0; j < M.col; ++j)
		{
			double temp;
			temp = M[ind][j];
			M[ind][j] = M[r][j];
			M[r][j] = temp;

			temp = E[ind][j];
			E[ind][j] = E[r][j];
			E[r][j] = temp;
		}
		//end ��������
		//start ��Ԫ��Ԫ�ػ�0
		for (int i = r + 1; i < M.row; ++i)
		{
			double temp;
			temp = M[i][c] / M[r][c];
			for (int j = 0; j < M.col; ++j)
			{
				M[i][j] -= M[r][j] * temp;
				E[i][j] -= E[r][j] * temp;
			}
		}
		//end ��Ԫ��Ԫ�ػ�0   
	}//end ������Ϊ������
	//start ��Ԫ��1
	for (int r = 0; r < M.row; ++r)
	{
		double temp = M[r][r];
		for (int c = r; c < M.col; ++c)
		{
			M[r][c] /= temp;
		}
		for (int c = 0; c < M.col; ++c)
		{
			E[r][c] /= temp;
		}
	}
	//end ��Ԫ��1
	//start ���Խ���
	for (int r = M.row - 2; r >= 0; --r)
	{
		for (int i = M.row - 1; i > r; --i)
		{
			for (int j = 0; j < M.col; ++j)
			{
				E[r][j] -= E[i][j] * M[r][i];
			}
			for (int j = 0; j < M.col; ++j)
			{
				M[r][j] -= M[i][j] * M[r][i];
			}
		}
	}//end ���Խ���
	return E;

}//end ������

//start �Ծ�������
void Matrix::modify()
{
	for (int k = 0; k < this->n; ++k)
		if (fabs(this->p[k]) < 1.0e-10)
			this->p[k] = 0;
}
//start ��������
Matrix::~Matrix()
{
	delete[] this->p;
}//end ��������