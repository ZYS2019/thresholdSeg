#pragma once
#include<opencv.hpp>

using namespace cv;
using namespace std;
class CircleFit
{
public:
	CircleFit();
	~CircleFit();
	Point3f LeastSquareFittingCircle(vector<Point2f> temp_coordinates);
	//Point3f LeastSquareFittingCircle(vector<Point2f> temp_coordinates);
};

