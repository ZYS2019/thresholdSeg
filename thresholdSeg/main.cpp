#include<opencv.hpp>
#include"thresholdSeg.h"
#include"CircleFit.h"

using namespace std;
using namespace cv;


int main()
{
	string srcPath;
	srcPath = "C:\\Users\\SZJ\\Desktop\\��ֵ�ָ�\\1.jpg";
	
	Mat srcMat = imread(srcPath);
	resize(srcMat, srcMat, Size(2000, 2000));
	Mat polarMat;
	Mat scaleRoIMat; 
	Mat grayMat;
	Mat binaryMat;
	int textNum = 4;
	vector<vector<int>> textRect(textNum, vector<int>(8));
	//step00.Բ����ȡ


	//step0.������任
	cv::linearPolar(srcMat, polarMat, Point(1474 , 1480), 2000, cv::INTER_LINEAR);
	rotate(polarMat, polarMat, 2);
	//step1.ROI��ȡ
	int x = 967;
	int y = 777;
	int w = 577;
	int h = 399;
	Rect roi(x, y, w, h);
	polarMat(roi).copyTo(scaleRoIMat);
	/*Mat grayMat;
	Mat binaryMat;
	Mat scaleRoIMat = imread("C:\\Users\\SZJ\\Desktop\\��ֵ�ָ�\\7_roi.jpg");*/
	//step2.�Ҷ�ת��
	cvtColor(scaleRoIMat, grayMat, COLOR_BGR2GRAY);

	//cv::imwrite("C:\\Users\\SZJ\\Desktop\\��ֵ�ָ�\\7_roi.jpg", grayMat);
	//step3.��ֵ�ָ�
	thresholdSeg enThresholdSeg;
	int alg = 3;   //0�����  1����������ֵ�㷨 2:Niblack�㷨  3��Sauvola�㷨  4��Bernsen�㷨  5.�Ľ��Ĵ��	double time0 = static_cast<double>(getTickCount());//��¼��ʼʱ��
	if(alg==0)
	{
		//���
		binaryMat = enThresholdSeg.OtsuAlgThreshold(grayMat);
		cv::imwrite("C:\\Users\\SZJ\\Desktop\\��ֵ�ָ�\\���\\7.jpg", binaryMat);
	}
	else if (alg == 1)
	{
		//������ֵ�㷨
		binaryMat = enThresholdSeg.IterationThreshold(grayMat);
		cv::imwrite("C:\\Users\\SZJ\\Desktop\\��ֵ�ָ�\\����\\3.jpg", binaryMat);
	}
	else if (alg == 2)
	{
		//Niblack�㷨
		enThresholdSeg.Niblack(grayMat, binaryMat, Size(25,225),-1);
		cv::imwrite("C:\\Users\\SZJ\\Desktop\\��ֵ�ָ�\\Niblack\\7-1.jpg", binaryMat);
	}
	else if (alg == 3)
	{
		//Sauvola��
		enThresholdSeg.SauvolaThresh(grayMat, binaryMat, 0.1, Size(25, 225));
		cv::imwrite("C:\\Users\\SZJ\\Desktop\\��ֵ�ָ�\\Sauvola\\1-0.1.jpg", binaryMat);
	}
	else if (alg == 4)
	{
		//Bernsen��
		enThresholdSeg.Bernsen(grayMat, binaryMat, Size(15,15));
	}
	else if (alg == 5)
	{
		//�Ľ��Ĵ��
		enThresholdSeg.NewOtsuAlgThreshold(grayMat, binaryMat, Size(12, 399));
	}
	else if (alg == 6)
	{
		//�Ľ��ĵ���ѭ����
		enThresholdSeg.NewIterationThreshold(grayMat, binaryMat, Size(50, 200));
	}
	//time0 = (((double)getTickCount() - time0) / getTickFrequency())*1000;
	//cout << "����ʱ��Ϊ��" << time0 << endl;
	system("pause");
	//imshow("src", srcMat);
	imshow("gray", grayMat);
	imshow("binary", binaryMat);
	waitKey(0);
	return 0;
}
