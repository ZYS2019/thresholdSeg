#include "CircleFit.h"



CircleFit::CircleFit()
{
}


CircleFit::~CircleFit()
{
}

Point3f CircleFit::LeastSquareFittingCircle(vector<Point2f> temp_coordinates)//利用opencv solve函数的高斯消元算法（LU分解）解方程组
{
	float x1 = 0;
	float x2 = 0;
	float x3 = 0;
	float y1 = 0;
	float y2 = 0;
	float y3 = 0;
	float x1y1 = 0;
	float x1y2 = 0;
	float x2y1 = 0;
	int num;
	vector<Point2f>::iterator k;
	Point3f tempcircle;
	num = temp_coordinates.size();
	for (k = temp_coordinates.begin(); k != temp_coordinates.end(); k++)
	{

		x1 = x1 + (*k).x;
		x2 = x2 + (*k).x * (*k).x;
		x3 = x3 + (*k).x * (*k).x * (*k).x;
		y1 = y1 + (*k).y;
		y2 = y2 + (*k).y * (*k).y;
		y3 = y3 + (*k).y * (*k).y * (*k).y;
		x1y1 = x1y1 + (*k).x * (*k).y;
		x1y2 = x1y2 + (*k).x * (*k).y * (*k).y;
		x2y1 = x2y1 + (*k).x * (*k).x * (*k).y;
	}
	Mat left_matrix = (Mat_<float>(3, 3) << x2, x1y1, x1, x1y1, y2, y1, x1, y1, num);
	cout << "left_matrix=" << left_matrix << endl;
	Mat right_matrix = (Mat_<float>(3, 1) << -(x3 + x1y2), -(x2y1 + y3), -(x2 + y2));
	cout << "right_matrix=" << right_matrix << endl;
	Mat solution(3, 1, CV_32F);
	solve(left_matrix, right_matrix, solution, CV_LU);
	cout << "solution=" << solution << endl;
	float a, b, c;
	a = solution.at<float>(0);
	b = solution.at<float>(1);
	c = solution.at<float>(2);
	tempcircle.x = -a / 2; //圆心x坐标
	tempcircle.y = -b / 2;//圆心y坐标
	tempcircle.z = sqrt(a*a + b*b - 4 * c) / 2;//圆心半径
	return tempcircle;
}
//Point3f CircleFit::LeastSquareFittingCircle(vector<Point2f> temp_coordinates)//高斯消元法直接求解方程组
//{
//	float x1 = 0;
//	float x2 = 0;
//	float x3 = 0;
//	float y1 = 0;
//	float y2 = 0;
//	float y3 = 0;
//	float x1y1 = 0;
//	float x1y2 = 0;
//	float x2y1 = 0;
//	int num;
//	vector<Point2f>::iterator k;
//	Point3f tempcircle;
//	for (k = temp_coordinates.begin(); k != temp_coordinates.end(); k++)
//	{
//
//		x1 = x1 + (*k).x;
//		x2 = x2 + (*k).x * (*k).x;
//		x3 = x3 + (*k).x * (*k).x * (*k).x;
//		y1 = y1 + (*k).y;
//		y2 = y2 + (*k).y * (*k).y;
//		y3 = y3 + (*k).y * (*k).y * (*k).y;
//		x1y1 = x1y1 + (*k).x * (*k).y;
//		x1y2 = x1y2 + (*k).x * (*k).y * (*k).y;
//		x2y1 = x2y1 + (*k).x * (*k).x * (*k).y;
//	}
//	float C, D, E, G, H, a, b, c;
//	num = temp_coordinates.size();
//	C = num*x2 - x1*x1;
//	D = num*x1y1 - x1*y1;
//	E = num*x3 + num*x1y2 - x1*(x2 + y2);
//	G = num*y2 - y1*y1;
//	H = num*x2y1 + num*y3 - y1*(x2 + y2);
//	a = (H*D - E*G) / (C*G - D*D);
//	b = (H*C - E*D) / (D*D - G*C);
//	c = -(x2 + y2 + a*x1 + b*y1) / num;
//	tempcircle.x = -a / 2; //圆心x坐标
//	tempcircle.y = -b / 2;//圆心y坐标
//	tempcircle.z = sqrt(a*a + b*b - 4 * c) / 2;//圆心半径
//	return tempcircle;
//}