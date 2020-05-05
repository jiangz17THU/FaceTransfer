#pragma once
class Matrix
{
public:
	Matrix();//默认构造函数
	//Matrix();
	//带参数的构造函数
	Matrix(const int&, const int&);
	Matrix(const int&, const int&, double[71][2]);
	Matrix(const int&, const int&, double[71][71]);
	//Matrix(const int&, const int&, double[1][4]);
	Matrix(const int&, const int&, double[4][4]);
	//复制构造函数
	Matrix(const Matrix&);
	//赋值操作符
	Matrix& operator=(const Matrix&);
	Matrix& operator+=(const Matrix&);
	Matrix& operator-=(const Matrix&);
	//析构函数
	~Matrix();

	//相等操作符
	friend bool operator==(const Matrix&, const Matrix&);
	//不等操作符
	friend bool operator!=(const Matrix& lhs, const Matrix& rhs)
	{
		return !(lhs == rhs);
	}
	//下标操作符
	double* operator[] (int index);
	//下标操作符
	double* operator[] (int index)const;
	//打印输出
	void print()const;
	//算数运算符重载

	// *(矩阵相乘)
	friend Matrix operator*(const Matrix&, const Matrix&);
	// *(实数*矩阵)
	friend Matrix operator*(const double&, const Matrix&);
	// *(矩阵*实数)
	friend Matrix operator*(const Matrix& M, const double& c)
	{
		return c * M;
	}

	//转置
	friend Matrix trv(const Matrix&);
	//求逆阵
	friend Matrix inv(Matrix);
	//对矩阵修整
	void modify();

	double* p;//存放数组
	int row;//行数
	int col;//列数
	int n;//元素个数
};
// +
Matrix operator+(const Matrix&, const Matrix&);
// -
Matrix operator-(const Matrix&, const Matrix&);

#include "cmath"
#include <iostream>
using namespace std;

//start 默认构造函数
Matrix::Matrix()
{
	this->row = 0; this->col = 0;
	this->n = 0;   this->p = 0;
}//end 默认构造函数

//start 带参数的构造函数
Matrix::Matrix(const int& r, const int& c)
{
	this->row = r; this->col = c; this->n = r * c;
	this->p = new double[n];
	for (int i = 0; i < n; ++i)
	{
		this->p[i] = 0.0;//默认值设为0
	}
}

Matrix::Matrix(const int& r, const int& c, double p[71][2])
{
	this->row = r; this->col = c; this->n = r * c;
	this->p = new double[n];
	for (int i = 0; i < r; ++i)
	{
		for(int j=0;j<c;j++)
		this->p[i*c+j] = p[i][j];//给p赋值
	}
	this->modify();
}//end 带参数的构造函数

Matrix::Matrix(const int& r, const int& c, double p[71][71])
{
	this->row = r; this->col = c; this->n = r * c;
	this->p = new double[n];
	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; j++)
			this->p[i * c + j] = p[i][j];//给p赋值
	}
	this->modify();
}//end 带参数的构造函数

Matrix::Matrix(const int& r, const int& c, double p[4][4])
{
	this->row = r; this->col = c; this->n = r * c;
	this->p = new double[n];
	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; j++)
			this->p[i * c + j] = p[i][j];//给p赋值
	}
	this->modify();
}//end 带参数的构造函数

//start 复制构造函数
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
}//end 复制构造函数

//start 赋值操作符
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
}//end 赋值操作符

//start +=重载
Matrix& Matrix::operator+=(const Matrix& other)
{
	if (this->row != other.row
		|| this->col != other.col)
	{
		cout << "行列相等的矩阵才能进行\"+=\"操作！" << endl;
		return *this;
	}
	for (int i = 0; i < this->n; ++i)
	{
		this->p[i] += other.p[i];
	}
	this->modify();
	return *this;
}//end +=重载
//start -=重载
Matrix& Matrix::operator-=(const Matrix& other)
{
	if (this->row != other.row
		|| this->col != other.col)
	{
		cout << "行列相等的矩阵才能进行\"-=\"操作！" << endl;
		return *this;
	}
	for (int i = 0; i < this->n; ++i)
	{
		this->p[i] -= other.p[i];
	}
	this->modify();
	return *this;

}//end -=重载

//start +重载
Matrix operator+(const Matrix& lhs, const Matrix& rhs)
{
	Matrix lhsTemp = lhs, rhsTemp = rhs;
	return lhsTemp += rhsTemp;
}//end +重载
//start -重载 
Matrix operator-(const Matrix& lhs, const Matrix& rhs)
{
	Matrix lhsTemp = lhs, rhsTemp = rhs;
	return lhsTemp -= rhsTemp;
}//end -重载

//start *重载(矩阵相乘)
Matrix operator*(const Matrix& lhs, const Matrix& rhs)
{
	if (rhs.row != lhs.col)
	{
		cout << "相乘矩阵不符合条件！" << endl;
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
}//end *重载(矩阵相乘)
// start *重载(实数与矩阵相乘)
Matrix operator*(const double& c, const Matrix& M)
{
	Matrix temp(M.row, M.col);
	for (int i = 0; i < temp.row * temp.col; ++i)
	{
		temp.p[i] = M.p[i] * c;
	}
	return temp;
}// end *重载(实数与矩阵相乘)
//start 转置
Matrix trv(const Matrix& M)
{
	Matrix temp(M.col, M.row);
	for (int r = 0; r < M.row; ++r)
		for (int c = 0; c < M.col; ++c)
			temp[c][r] = M[r][c];
	return temp;
}//end 转置
//start 打印输出
void Matrix::print()const
{
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
			cout << (*this)[i][j] << "  ";
		cout << endl;
	}
}//end 打印输出

//start 相等操作符(友元函数)
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
}//end 相等操作符(友元函数)

//start 下标操作符
double* Matrix::operator[](int indexRow)const
{
	return this->p + indexRow * col;
}//end 下标操作符
//start 下标操作符
double* Matrix::operator[](int indexRow)
{
	return this->p + indexRow * col;
}//end 下标操作符



//start 求逆阵
Matrix inv(Matrix M)
{
	//start 判断矩阵是否可逆
	if (M.row != M.col)
	{
		cout << "只有方阵才能求逆！" << endl;
		return Matrix();
	}
	//end 判断矩阵是否可逆
	//start 构造单位阵
	Matrix E(M.row, M.col);
	for (int i = 0; i < M.row; ++i)
		E[i][i] = 1.0;
	//end 构造单位阵
	//start 将矩阵化为阶梯型
	for (int c = 0, r = 0; c < M.col && r < M.row; ++c, ++r)
	{
		int ind;//记录主元行号
		//start 寻找主元
		for (int i = r; i < M.row; ++i)
		{
			if (M[i][c] != 0)
			{
				ind = i;
				break;
			}
			if (i == M.row - 1)
			{
				cout << "矩阵不可逆！" << endl;
				return Matrix();
			}
		}
		//end 寻找主元
		//start 交换两行
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
		//end 交换两行
		//start 主元下元素化0
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
		//end 主元下元素化0   
	}//end 将矩阵化为阶梯型
	//start 主元化1
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
	//end 主元化1
	//start 化对角型
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
	}//end 化对角型
	return E;

}//end 求逆阵

//start 对矩阵修整
void Matrix::modify()
{
	for (int k = 0; k < this->n; ++k)
		if (fabs(this->p[k]) < 1.0e-10)
			this->p[k] = 0;
}
//start 析构函数
Matrix::~Matrix()
{
	delete[] this->p;
}//end 析构函数