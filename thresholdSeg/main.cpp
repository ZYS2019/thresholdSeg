#include<opencv.hpp>
#include"thresholdSeg.h"
#include"CircleFit.h"

using namespace std;
using namespace cv;


int main()
{
	string srcPath;
	srcPath = "C:\\Users\\SZJ\\Desktop\\阈值分割\\1.jpg";
	
	Mat srcMat = imread(srcPath);
	resize(srcMat, srcMat, Size(2000, 2000));
	Mat polarMat;
	Mat scaleRoIMat; 
	Mat grayMat;
	Mat binaryMat;
	int textNum = 4;
	vector<vector<int>> textRect(textNum, vector<int>(8));
	//step00.圆心提取


	//step0.极坐标变换
	cv::linearPolar(srcMat, polarMat, Point(1474 , 1480), 2000, cv::INTER_LINEAR);
	rotate(polarMat, polarMat, 2);
	//step1.ROI提取
	int x = 967;
	int y = 777;
	int w = 577;
	int h = 399;
	Rect roi(x, y, w, h);
	polarMat(roi).copyTo(scaleRoIMat);
	/*Mat grayMat;
	Mat binaryMat;
	Mat scaleRoIMat = imread("C:\\Users\\SZJ\\Desktop\\阈值分割\\7_roi.jpg");*/
	//step2.灰度转换
	cvtColor(scaleRoIMat, grayMat, COLOR_BGR2GRAY);

	//cv::imwrite("C:\\Users\\SZJ\\Desktop\\阈值分割\\7_roi.jpg", grayMat);
	//step3.阈值分割
	thresholdSeg enThresholdSeg;
	int alg = 3;   //0：大津法  1：迭代求阈值算法 2:Niblack算法  3：Sauvola算法  4：Bernsen算法  5.改进的大津法	double time0 = static_cast<double>(getTickCount());//记录起始时间
	if(alg==0)
	{
		//大津法
		binaryMat = enThresholdSeg.OtsuAlgThreshold(grayMat);
		cv::imwrite("C:\\Users\\SZJ\\Desktop\\阈值分割\\大津\\7.jpg", binaryMat);
	}
	else if (alg == 1)
	{
		//迭代阈值算法
		binaryMat = enThresholdSeg.IterationThreshold(grayMat);
		cv::imwrite("C:\\Users\\SZJ\\Desktop\\阈值分割\\迭代\\3.jpg", binaryMat);
	}
	else if (alg == 2)
	{
		//Niblack算法
		enThresholdSeg.Niblack(grayMat, binaryMat, Size(25,225),-1);
		cv::imwrite("C:\\Users\\SZJ\\Desktop\\阈值分割\\Niblack\\7-1.jpg", binaryMat);
	}
	else if (alg == 3)
	{
		//Sauvola法
		enThresholdSeg.SauvolaThresh(grayMat, binaryMat, 0.1, Size(25, 225));
		cv::imwrite("C:\\Users\\SZJ\\Desktop\\阈值分割\\Sauvola\\1-0.1.jpg", binaryMat);
	}
	else if (alg == 4)
	{
		//Bernsen法
		enThresholdSeg.Bernsen(grayMat, binaryMat, Size(15,15));
	}
	else if (alg == 5)
	{
		//改进的大津法
		enThresholdSeg.NewOtsuAlgThreshold(grayMat, binaryMat, Size(12, 399));
	}
	else if (alg == 6)
	{
		//改进的迭代循环法
		enThresholdSeg.NewIterationThreshold(grayMat, binaryMat, Size(50, 200));
	}
	//time0 = (((double)getTickCount() - time0) / getTickFrequency())*1000;
	//cout << "运行时间为：" << time0 << endl;
	system("pause");
	//imshow("src", srcMat);
	imshow("gray", grayMat);
	imshow("binary", binaryMat);
	waitKey(0);
	return 0;
}
