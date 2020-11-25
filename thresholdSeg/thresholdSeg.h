#pragma once

#include<opencv.hpp>
#include<iostream>
using namespace cv;
using namespace std;

class thresholdSeg
{
public:
	thresholdSeg();
	~thresholdSeg();
	cv::Mat OtsuAlgThreshold(Mat &image);
	void NewOtsuAlgThreshold(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize);

	cv::Mat IterationThreshold(Mat src);
	void NewIterationThreshold(Mat src, cv::Mat &dst, cv::Size wndSize);
	double GetMatAverage(const cv::Mat& src);
	double GetMatStdDev(const cv::Mat& src, double meanValue);
	void Niblack(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize, float k);

	float fastMean(cv::Mat& integral, int x, int y, int window);
	//void Sauvola(cv::Mat& inpImg, cv::Mat& resImg, int window, float k);
	void SauvolaThresh(const cv::Mat& src, cv::Mat& dst, const int k, const cv::Size wndSize);

	void GetMatMaxMin(const cv::Mat& m, int& maxValue, int& minValue);
	void Bernsen(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize);
	void Bernsen(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize, int differMax, int meanMax);
};

