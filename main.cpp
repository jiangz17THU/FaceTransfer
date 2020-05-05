#include<opencv2/opencv.hpp>
#include<iostream>
//#include <time.h>
#include <cmath>
#include<sstream>
#include"matrix.h"
#include <opencv2\highgui\highgui.hpp>  //opencv����
#include <opencv2\imgproc\imgproc.hpp>
//#include <opencv2\core\core.hpp>

using namespace cv;
using namespace std;

# define pi 3.145926

class point
{
public:
	double x;
	double y;
};
point p[68];  //ȫ�ֶ���68�����Ƶ��Ŀ���
point p_[68];
double dis(point p1, point p2) //��������
{
	double dis;
	double dis_square = (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
	dis = sqrt(long double (dis_square));
	return dis;
}
float max(float a, float b)
{
	if (a > b)
		return a;
	else
		return b;
}

float* roll(int x, int y, float anglemax, int R)  //Ť���任����
{
	float radianmax = anglemax * pi / (float)180;  //�����ת�Ƕ�
	float dis = sqrt(x * x + y * y); 
	float posi[2]; 
	float a = radianmax * max((1 - dis / (float)R), 0);
	posi[0] = x * cos(a) - y * sin(a);
	posi[1] = x * sin(a) + y * cos(a);
	return posi;    //����ԭͼ��������
}
float* distort1(int x, int y,int R)  //Ͱ�ͻ��亯����͹�澵��
{
	float dis = sqrt(x * x + y * y);
	float posi[2];
	if (dis <= R&& dis>0)
	{
		posi[0] = R / (float)dis * asin(dis / (float)R) * x;  
		posi[1] = R / (float)dis * asin(dis / (float)R) * y;
	}
	else if (dis > R)    //�������Radiusʱ��ͼ����Ϊ��ɫ��posi�������ر߽磩
	{
		posi[0] = 1000;
		posi[1] = 1000;
	}
	else          //dis==0���������ֹ����0�ж�
	{
		posi[0] = x;
		posi[1] = y;
	}
	
	return posi;
}

float* distort2(int x, int y, int R)  //���ͻ��亯�������澵��
{
	float dis = sqrt(x * x + y * y);
	float posi[2];
	if  (dis > 0)
	{
		posi[0] = R / (float)dis * sin(dis / (float)R) * x;
		posi[1] = R / (float)dis * sin(dis / (float)R) * y;
	}
	
	else  //dis==0���������ֹ����0�ж�
	{
		posi[0] = x;
		posi[1] = y;
	}

	return posi;
}
int* near(float xposi, float yposi) //����ڲ�ֵ���������������������
{
	int a[2];
	a[0] = int(xposi + 0.5) ;
	a[1] = int(yposi + 0.5) ;
	return a;
}

int* bilinear(double xposi, double yposi,Mat mat)  //˫���Բ�ֵ�����ݸ������������֪λͼ�����ظõ�����RGBֵ
{
	double u = xposi - (int)xposi ;  
	double v = yposi - (int)yposi ;
	int a0=mat.at<Vec3b>((int)xposi, (int)yposi)[0];
	int b0 = mat.at<Vec3b>((int)xposi,1+((int)yposi ))[0];
	int c0= mat.at<Vec3b>(1+ (int)xposi , (int)yposi )[0];
	int d0= mat.at<Vec3b>(1+ (int)xposi , 1+ (int)yposi )[0];
	int a1 = mat.at<Vec3b>((int)xposi , (int)yposi )[1];
	int b1 = mat.at<Vec3b>((int)xposi, 1 + (int)yposi)[1];
	int c1 = mat.at<Vec3b>(1 + (int)xposi , (int)yposi )[1];
	int d1 = mat.at<Vec3b>(1 + (int)xposi , 1 + (int)yposi )[1];
	int a2 = mat.at<Vec3b>((int)xposi , (int)yposi )[2];
	int b2 = mat.at<Vec3b>((int)xposi , 1 + (int)yposi )[2];
	int c2 = mat.at<Vec3b>(1 + (int)xposi  , (int)yposi )[2];
	int d2 = mat.at<Vec3b>(1 + (int)xposi  , 1 + (int)yposi )[2];
	int t[3];
	t[0] = (1 - v) * ((1 - u) * a0 + u * c0) + v * ((1 - u) * b0 + u * d0); 
	t[1] = (1 - v) * ((1 - u) * a1 + u * c1) + v * ((1 - u) * b1 + u * d1);
	t[2] = (1 - v) * ((1 - u) * a2 + u * c2) + v * ((1 - u) * b2 + u * d2);
	return t;
}

double S(long double x)   //˫����S��x������
{
	if (fabs(x) <= 1)
		return 1 - 2 * x * x + pow(fabs(x), 3);
	else if (fabs(x) > 1 && fabs(x) < 2)
		return 4 - 8 * fabs(x) + 5 * x * x - pow(fabs(x), 3);
	else
		return 0;
}

double* bicubic(double x, double y,Mat mat)    //˫���β�ֵ�����ݸ������������֪λͼ�����ظõ�����RGBֵ
{
	double rgb[3];
	double u = x - (int)x;
	double v = y - (int)y;
	double a[1][4];
	double c[1][4];
	double b0[4][4];
	double b1[4][4];
	double b2[4][4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i == 0)
			{
				a[i][j] = S(u + 1 - j);
				c[i][j] = S(v + 1 - j);
			}
			b0[i][j]= mat.at<Vec3b>((int)x+i-1, (int)y+j-1)[0];
			b1[i][j] = mat.at<Vec3b>((int)x + i - 1, (int)y + j - 1)[1];
			b2[i][j] = mat.at<Vec3b>((int)x + i - 1, (int)y + j - 1)[2];
		}
	}
	Matrix A(1, 4);
	for (int i = 0; i < 4; i++)  A[0][i] = a[0][i];
	Matrix C(1, 4);
	for (int i = 0; i < 4; i++)  C[0][i] = c[0][i];
	Matrix B0(4, 4, b0);
	Matrix B1(4, 4, b1);
	Matrix B2(4, 4, b2);
	Matrix cT = trv(C);
	Matrix temp0 = A * B0 * cT;
	Matrix temp1 = A * B1 * cT;
	Matrix temp2 = A * B2 * cT;
	rgb[0] = temp0[0][0];
	if (rgb[0] < 0) rgb[0] = 0;   //��RGBֵ���ƣ���ֹ����0-255
	if (rgb[0] > 255) rgb[0] = 255;

	rgb[1] = temp1[0][0];
	if (rgb[1] < 0) rgb[1] = 0;
	if (rgb[1] > 255) rgb[1] = 255;

	rgb[2] = temp2[0][0];
	if (rgb[2] < 0) rgb[2] = 0;
	if (rgb[2] > 255) rgb[2] = 255;
	return rgb;
}

double U(point p1, point p2)   //tps�е�U����
{
	long double y = (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
	long double t = sqrt(y);
	double s= log(t * t);
	if (t == 0)
		return 0;
	else
		return t * t * s;
}

double* f(int x, int y,Matrix ans)  //tps�еı��κ���f��x,y��
{
	    double posi[2];
		point point;
		point.x = x; point.y = y;
		double m, n;
		m = n = 0;
		double result[68];
		for (int i = 0; i < 68; i++)
		{
			result[i] = U(p[i], point);
		}
		for (int i = 0; i < 68; i++)
		{
			m += ans[0][i] * result[i];
			n += ans[1][i] * result[i];
		}
		
		posi[0] = ans[0][68] + ans[0][69] * x + ans[0][70] * y + m;
		posi[1] = ans[1][68] + ans[1][69] * x + ans[1][70] * y + n;

	   return posi;
}

int main(int argc, char** argv)
{
	while (1)
	{
		cout << "����������򽫼���ѭ��������д��ĸ��Q���˳�" << endl;
		char c;
		cin >> c;
		if (c == 'Q')
			break;
		cout << "�������������롰1����ѡ�����������롰2����" << endl;
		int instruction;
		cin >> instruction;
		if (instruction == 1)
		{
			int mode;
			cout << "��ѡ��ͼƬ����ģʽ�����롰1������ͼ��Ť�������롰2������ͼ�����" << endl;
			cin >> mode;
			int way;
			int amax, Radius, radius;
			if (mode == 1)
			{
				cout << "��ѡ�������ת�Ƕ�amax����ת�뾶Radius:" << endl;
				cout << "��ת�Ƕ�amax:";
				cin >> amax;
				cout << "��ת�뾶Radius:";
				cin >> Radius;
			}
			if (mode == 2)
			{
				cout << "��ѡ����䷽ʽ�����롰1������Ͱ�ͻ��䣨͹�������롰2���������ͻ��䣨����" << endl;
				cin >> way;
				cout << "���������뾶radius��";
				cin >> radius;
			}
			int method;
			cout << "��ѡ���ֵ���������롰1����������ڲ�ֵ�����롰2������˫���Բ�ֵ�����롰3������˫���β�ֵ" << endl;
			cin >> method;
			if (method == 3)
				cout << "ʱ���Գ������ĵȴ�Ŷ~" << endl;
			Mat image = imread("./data/THU.jpg"); //ԭͼ��
			Mat copy_im;  //Ť��ͼ��
			image.copyTo(copy_im);
			Mat copy_im1; //����ͼ��
			image.copyTo(copy_im1);

			if (mode == 1)  //Ť��
			{
				for (int i = 1; i < copy_im.rows; i++)
				{

					for (int j = 1; j < copy_im.cols; j++)
					{

						float m = roll((i - 256), (j - 256), -amax, Radius)[0];
						float n = roll((i - 256), (j - 256), -amax, Radius)[1];
						float x = m + 256;
						float y = n + 256;
						int x0, y0;
						if (method == 1)
						{
							x0 = near(x, y)[0];
							y0 = near(x, y)[1];
						}

						if (x >= 510 || y >= 510 || x <= 1 || y <= 1)
						{
							copy_im.at<Vec3b>(i, j)[0] = image.at<Vec3b>(i, j)[0];
							copy_im.at<Vec3b>(i, j)[1] = image.at<Vec3b>(i, j)[1];
							copy_im.at<Vec3b>(i, j)[2] = image.at<Vec3b>(i, j)[2];
						}
						else
						{
							int t;
							int s;
							int q;
							if (method == 1)  //�����
							{
								copy_im.at<Vec3b>(i, j)[0] = image.at<Vec3b>(x0, y0)[0];
								copy_im.at<Vec3b>(i, j)[1] = image.at<Vec3b>(x0, y0)[1];
								copy_im.at<Vec3b>(i, j)[2] = image.at<Vec3b>(x0, y0)[2];
							}
							else if (method == 2)  //˫����
							{
								q = bilinear(x, y, image)[2];
								s = bilinear(x, y, image)[1];
								t = bilinear(x, y, image)[0];
								copy_im.at<Vec3b>(i, j)[0] = t;
								copy_im.at<Vec3b>(i, j)[1] = s;
								copy_im.at<Vec3b>(i, j)[2] = q;
							}
							else if (method == 3) //˫����
							{
								q = bicubic(x, y, image)[2];
								s = bicubic(x, y, image)[1];
								t = bicubic(x, y, image)[0];
								copy_im.at<Vec3b>(i, j)[0] = t;
								copy_im.at<Vec3b>(i, j)[1] = s;
								copy_im.at<Vec3b>(i, j)[2] = q;
							}
						}

					}
				}

			}


			else if (mode == 2)  //����
			{
				for (int i = 0; i < copy_im1.rows; i++)
				{
					for (int j = 0; j < copy_im1.cols; j++)
					{
						float m, n;
						if (way == 1)
						{
							m = distort1((i - 256), (j - 256), radius)[0];
							n = distort1((i - 256), (j - 256), radius)[1];
						}
						else if (way == 2)
						{
							m = distort2((i - 256), (j - 256), 250)[0];
							n = distort2((i - 256), (j - 256), 250)[1];
						}

						float x = m + 256;
						float y = n + 256;
						int x0, y0;

						x0 = near(x, y)[0];
						y0 = near(x, y)[1];


						if (x >= 510 || y >= 510 || x <= 1 || y <= 1)
						{
							copy_im1.at<Vec3b>(i, j)[0] = 0;
							copy_im1.at<Vec3b>(i, j)[1] = 0;
							copy_im1.at<Vec3b>(i, j)[2] = 0;
						}
						else
						{
							int t;
							int s;
							int q;
							if (method == 1)  //�����
							{
								copy_im1.at<Vec3b>(i, j)[0] = image.at<Vec3b>(x0, y0)[0];
								copy_im1.at<Vec3b>(i, j)[1] = image.at<Vec3b>(x0, y0)[1];
								copy_im1.at<Vec3b>(i, j)[2] = image.at<Vec3b>(x0, y0)[2];
							}
							else if (method == 2)
							{
								q = bilinear(x, y, image)[2];
								s = bilinear(x, y, image)[1];
								t = bilinear(x, y, image)[0];
								copy_im1.at<Vec3b>(i, j)[0] = t;
								copy_im1.at<Vec3b>(i, j)[1] = s;
								copy_im1.at<Vec3b>(i, j)[2] = q;

							}
							else if (method == 3)
							{
								q = bicubic(x, y, image)[2];
								s = bicubic(x, y, image)[1];
								t = bicubic(x, y, image)[0];
								copy_im1.at<Vec3b>(i, j)[0] = t;
								copy_im1.at<Vec3b>(i, j)[1] = s;
								copy_im1.at<Vec3b>(i, j)[2] = q;
							}

						}
					}
				}
			}

			imshow("ԭͼ", image);
			if (mode == 1)
				imshow("Ť��", copy_im);
			if (mode == 2)
				imshow("����", copy_im1);
			waitKey(0);
		}
		
		else if (instruction == 2)
		{
			double K[68][68]; //k����
			double P[68][3];  //P����
			double L[71][71]; //L����
			double Y[71][2];

			char c[4000]; int k = 0;
			char c_[4000]; int k_ = 0;
			int s, d; //ԴͼƬ��Ŀ��ͼƬ���
			cout << "����ԴͼƬ��Ŀ��ͼƬ���(1-9)" << endl;
			cout << "ԴͼƬ:";
			cin >> s;
			cout << "Ŀ��ͼƬ:";
			cin >> d;
			char stxt[13] = "./data/i.txt";
			char dtxt[13] = "./data/i.txt";
			char sjpg[13] = "./data/i.jpg";
			char djpg[13] = "./data/i.jpg";
			for (int i = 0; i < 10; i++)
			{
				if (i == s)
				{
					stxt[7] = i + 48;
					sjpg[7] = i + 48;
				}
			}
			for (int i = 0; i < 10; i++)
			{
				if (i == d)
				{
					dtxt[7] = i + 48;
					djpg[7] = i + 48;
				}
			}
			FILE* fp;
			if ((fp = fopen(dtxt, "r")) == NULL)  //��ȡ1.txt���ļ�
			{
				printf("��ȡ�ļ�ʧ�� \n ");
				exit(1);
			}
			while (!feof(fp))
			{

				c[k] = fgetc(fp);
				k++;
			}

			FILE* fp_;
			if ((fp_ = fopen(stxt, "r")) == NULL)  //��ȡ1.txt���ļ� ԭͼ��
			{
				printf("��ȡ�ļ�ʧ�� \n ");
				exit(1);
			}
			while (!feof(fp_))
			{
				c_[k_] = fgetc(fp_);
				k_++;
			}

			for (int i = 0; i < 68 * 2; i++) // ���Ƶ㸳ֵ
			{
				if (i % 2 == 0)
				{
					char byte[25];
					for (int j = 0; j < 25; j++)
						byte[j] = c[(i / 2) * 50 + j];
					stringstream ss;
					ss << byte;
					ss >> p[i / 2].x;
				}

				else if (i % 2 == 1)
				{
					char byte[25];
					for (int j = 0; j < 25; j++)
						byte[j] = c[(i / 2) * 50 + j + 25];
					stringstream ss;
					ss << byte;
					ss >> p[i / 2].y;
				}
			}

			for (int i = 0; i < 68 * 2; i++) // Ŀ��㸳ֵ
			{
				if (i % 2 == 0)
				{
					char byte[25];
					for (int j = 0; j < 25; j++)
						byte[j] = c_[(i / 2) * 50 + j];
					stringstream ss;
					ss << byte;
					ss >> p_[i / 2].x;
				}

				else if (i % 2 == 1)
				{
					char byte[25];
					for (int j = 0; j < 25; j++)
						byte[j] = c_[(i / 2) * 50 + j + 25];
					stringstream ss;
					ss << byte;
					ss >> p_[i / 2].y;
				}
			}


			for (int i = 0; i < 68; i++)  //k����ֵ
			{
				for (int j = 0; j < 68; j++)
					K[i][j] = U(p[i], p[j]);
			}


			for (int i = 0; i < 68; i++)  //p����ֵ
			{
				for (int j = 0; j < 3; j++)
				{
					if (j == 0)
						P[i][j] = 1;
					else if (j == 1)
						P[i][j] = p[i].x;
					else
						P[i][j] = p[i].y;
				}
			}

			for (int i = 0; i < 71; i++)  //L����ֵ
			{
				for (int j = 0; j < 71; j++)
				{
					if (i < 68 && j < 68)
						L[i][j] = K[i][j];
					else if (i < 68 && j >= 68)
						L[i][j] = P[i][j - 68];
					else if (i >= 68 && j < 68)
						L[i][j] = L[j][i];
					else
						L[i][j] = 0;
				}
			}

			for (int i = 0; i < 71; i++)  //Y����ֵ
			{
				for (int j = 0; j < 2; j++)
				{
					if (i < 68 && j == 0)
						Y[i][j] = p_[i].x;
					else if (i < 68 && j == 1)
						Y[i][j] = p_[i].y;
					else
						Y[i][j] = 0;
				}
			}

			Matrix M_Y(71, 2, Y);
			Matrix M_L(71, 71, L);
			Matrix M = inv(M_L) * M_Y;
			Matrix ans = trv(M);

			Mat tps1 = imread(djpg);
			Mat tps2 = imread(sjpg);
			Mat tps;
			tps2.copyTo(tps);
			cout << "������....�����ĵȴ�"<<endl;
			
			for (int i = 0; i < tps.rows; i++)  //����ͼrgb��ֵ
			{
				for (int j = 0; j < tps.cols; j++)
				{
					//double* posi;
					//posi = f(j, i, ans);
					//f�����ĵ��û�ʹ����������ʽ�����������������
					point point;
					point.x = j; point.y = i;
					double sum1, sum2;
					sum1 = sum2 = 0;
					double result[68];
		
					for (int i = 0; i < 68; i++)
					{
						result[i] = U(p[i], point);
						sum1 += ans[0][i] * result[i];
						sum2 += ans[1][i] * result[i];
					}
					double m = ans[0][68] + ans[0][69] * j + ans[0][70] * i + sum1;
					double n = ans[1][68] + ans[1][69] * j + ans[1][70] * i + sum2;
					
					int x0 = near(m, n)[0];
					if (x0 > tps.cols - 1) x0 = tps.cols - 1;
					else if (x0 < 0) x0 = 0;
					int y0 = near(m, n)[1];
					if (y0 > tps.rows - 1) y0 = tps.rows - 1;
					else if (y0 < 0) y0 = 0;
					
					tps.at<Vec3b>(i, j)[0] = tps2.at<Vec3b>(y0, x0)[0];
					tps.at<Vec3b>(i, j)[1] = tps2.at<Vec3b>(y0, x0)[1];
					tps.at<Vec3b>(i, j)[2] = tps2.at<Vec3b>(y0, x0)[2];
				}
			}
			
			imshow("Ŀ��ͼƬ", tps1);
			imshow("ԴͼƬ", tps2);
			imshow("����ͼƬ", tps);
			waitKey(0);
			
		}

	}
	return 0;
}



