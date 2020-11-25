#include "thresholdSeg.h"



thresholdSeg::thresholdSeg()
{
}


thresholdSeg::~thresholdSeg()
{
}

//------------------------------00.Otsu算法-------------------------
Mat thresholdSeg::OtsuAlgThreshold(Mat &image)
{
	if (image.channels() != 1)
	{
		cout << "Please input Gray-image!" << endl;
	}
	int T = 0; //Otsu算法阈值  
	double varValue = 0; //类间方差中间值保存
	double w0 = 0; //前景像素点数所占比例  
	double w1 = 0; //背景像素点数所占比例  
	double u0 = 0; //前景平均灰度  
	double u1 = 0; //背景平均灰度  
	double Histogram[256] = { 0 }; //灰度直方图，下标是灰度值，保存内容是灰度值对应的像素点总数  
	uchar *data = image.data;

	double totalNum = image.rows*image.cols; //像素总数

	for (int i = 0; i < image.rows; i++)
	{
		for (int j = 0; j < image.cols; j++)
		{
			if (image.at<uchar>(i, j) != 0) Histogram[data[i*image.step + j]]++;
		}
	}
	int minpos, maxpos;
	for (int i = 0; i < 255; i++)
	{
		if (Histogram[i] != 0)
		{
			minpos = i;
			break;
		}
	}
	for (int i = 255; i > 0; i--)
	{
		if (Histogram[i] != 0)
		{
			maxpos = i;
			break;
		}
	}

	for (int i = minpos; i <= maxpos; i++)
	{
		//每次遍历之前初始化各变量  
		w1 = 0;       u1 = 0;       w0 = 0;       u0 = 0;
		//***********背景各分量值计算**************************  
		for (int j = 0; j <= i; j++) //背景部分各值计算  
		{
			w1 += Histogram[j];   //背景部分像素点总数  
			u1 += j*Histogram[j]; //背景部分像素总灰度和  
		}
		if (w1 == 0) //背景部分像素点数为0时退出  
		{
			break;
		}
		u1 = u1 / w1; //背景像素平均灰度  
		w1 = w1 / totalNum; // 背景部分像素点数所占比例
							//***********背景各分量值计算**************************  

							//***********前景各分量值计算**************************  
		for (int k = i + 1; k < 255; k++)
		{
			w0 += Histogram[k];  //前景部分像素点总数  
			u0 += k*Histogram[k]; //前景部分像素总灰度和  
		}
		if (w0 == 0) //前景部分像素点数为0时退出  
		{
			break;
		}
		u0 = u0 / w0; //前景像素平均灰度  
		w0 = w0 / totalNum; // 前景部分像素点数所占比例  
							//***********前景各分量值计算**************************  

							//***********类间方差计算******************************  
		double varValueI = w0*w1*(u1 - u0)*(u1 - u0); //当前类间方差计算  
		if (varValue < varValueI)
		{
			varValue = varValueI;
			T = i;
		}
	}
	Mat dst;
	threshold(image, dst, T, 255, CV_THRESH_OTSU);
	return dst;
}

void thresholdSeg::NewOtsuAlgThreshold(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize)
{
	if (src.channels() != 1)
	{
		cout << "Please input Gray-image!" << endl;
	}
	vector<vector<int>> imageThre;
	for (int y = 0; y < src.rows; y = y + wndSize.height)
	{
		vector<int> xThre;
		for (int x = 0; x < src.cols; x=x + wndSize.width)
		{
			Rect minRoI(x, y, wndSize.width, wndSize.height);
			//矩形框出界处理
			if (minRoI.x + minRoI.width > src.cols)
			{
				minRoI.width = src.cols - minRoI.x;
			}
			if (minRoI.y + minRoI.height > src.rows)
			{
				minRoI.height = src.rows - minRoI.y;
			}

			Mat roiMat = src(minRoI);

			int T = 0; //Otsu算法阈值  
			double varValue = 0; //类间方差中间值保存
			double w0 = 0; //前景像素点数所占比例  
			double w1 = 0; //背景像素点数所占比例  
			double u0 = 0; //前景平均灰度  
			double u1 = 0; //背景平均灰度  
			double Histogram[256] = { 0 }; //灰度直方图，下标是灰度值，保存内容是灰度值对应的像素点总数  
			uchar *data = roiMat.data;

			double totalNum = roiMat.rows*roiMat.cols; //像素总数

			for (int i = 0; i < roiMat.rows; i++)
			{
				for (int j = 0; j < roiMat.cols; j++)
				{
					if (roiMat.at<uchar>(i, j) != 0) Histogram[data[i*roiMat.step + j]]++;
				}
			}
			int minpos, maxpos;
			for (int i = 0; i < 255; i++)
			{
				if (Histogram[i] != 0)
				{
					minpos = i;
					break;
				}
			}
			for (int i = 255; i > 0; i--)
			{
				if (Histogram[i] != 0)
				{
					maxpos = i;
					break;
				}
			}

			for (int i = minpos; i <= maxpos; i++)
			{
				//每次遍历之前初始化各变量  
				w1 = 0;       u1 = 0;       w0 = 0;       u0 = 0;
				//***********背景各分量值计算**************************  
				for (int j = 0; j <= i; j++) //背景部分各值计算  
				{
					w1 += Histogram[j];   //背景部分像素点总数  
					u1 += j*Histogram[j]; //背景部分像素总灰度和  
				}
				if (w1 == 0) //背景部分像素点数为0时退出  
				{
					break;
				}
				u1 = u1 / w1; //背景像素平均灰度  
				w1 = w1 / totalNum; // 背景部分像素点数所占比例
									//***********背景各分量值计算**************************  

									//***********前景各分量值计算**************************  
				for (int k = i + 1; k < 255; k++)
				{
					w0 += Histogram[k];  //前景部分像素点总数  
					u0 += k*Histogram[k]; //前景部分像素总灰度和  
				}
				if (w0 == 0) //前景部分像素点数为0时退出  
				{
					break;
				}
				u0 = u0 / w0; //前景像素平均灰度  
				w0 = w0 / totalNum; // 前景部分像素点数所占比例  
									//***********前景各分量值计算**************************  

									//***********类间方差计算******************************  
				double varValueI = w0*w1*(u1 - u0)*(u1 - u0); //当前类间方差计算  
				if (varValue < varValueI)
				{
					varValue = varValueI;
					T = i;
				}
			}
			xThre.push_back(T);
		}
		imageThre.push_back(xThre);
	}
	dst = cv::Mat::zeros(src.rows, src.cols, CV_8UC1);
	for (int y = 0; y < src.rows; ++y)
	{
		for (int x = 0; x < src.cols; ++x)
		{

			double flagValue = imageThre[y/wndSize.height][x/wndSize.width];
			int value = src.at<uchar>(y, x);
			if (value > flagValue)
			{
				dst.at<uchar>(y, x) = 255;
			}
			else
			{
				dst.at<uchar>(y, x) = 0;
			}
		}
	}
}

//--------------------------------------01.迭代阈值算法------------------------------------------------
Mat thresholdSeg::IterationThreshold(Mat src)
{
	int width = src.cols;
	int height = src.rows;
	int hisData[256] = { 0 };
	for (int j = 0; j < height; j++)
	{
		uchar* data = src.ptr<uchar>(j);
		for (int i = 0; i < width; i++)
			hisData[data[i]]++;
	}

	int T0 = 0;
	for (int i = 0; i < 256; i++)
	{
		T0 += i*hisData[i];
	}
	T0 /= width*height;

	int T1 = 0, T2 = 0;
	int num1 = 0, num2 = 0;
	int T = 0;
	while (1)
	{
		for (int i = 0; i < T0 + 1; i++)
		{
			T1 += i*hisData[i];
			num1 += hisData[i];
		}
		if (num1 == 0)
			continue;
		for (int i = T0 + 1; i < 256; i++)
		{
			T2 += i*hisData[i];
			num2 += hisData[i];
		}
		if (num2 == 0)
			continue;

		T = (T1 / num1 + T2 / num2) / 2;

		if (T == T0)
			break;
		else
			T0 = T;
	}

	Mat dst;
	threshold(src, dst, T, 255, 0);
	return dst;
}
void thresholdSeg::NewIterationThreshold(Mat src, cv::Mat & dst, cv::Size wndSize)
{
	if (src.channels() != 1)
	{
		cout << "Please input Gray-image!" << endl;
	}
	vector<vector<int>> imageThre;
	for (int y = 0; y < src.rows; y = y + wndSize.height)
	{
		vector<int> xThre;
		for (int x = 0; x < src.cols; x = x + wndSize.width)
		{
			Rect minRoI(x, y, wndSize.width, wndSize.height);
			//矩形框出界处理
			if (minRoI.x + minRoI.width > src.cols)
			{
				minRoI.width = src.cols - minRoI.x;
			}
			if (minRoI.y + minRoI.height > src.rows)
			{
				minRoI.height = src.rows - minRoI.y;
			}

			Mat roiMat = src(minRoI);
			int width = roiMat.cols;
			int height = roiMat.rows;
			int hisData[256] = { 0 };
			for (int j = 0; j < height; j++)
			{
				uchar* data = roiMat.ptr<uchar>(j);
				for (int i = 0; i < width; i++)
					hisData[data[i]]++;
			}

			int T0 = 0;
			for (int i = 0; i < 256; i++)
			{
				T0 += i*hisData[i];
			}
			T0 /= width*height;

			int T1 = 0, T2 = 0;
			int num1 = 0, num2 = 0;
			int T = 0;
			while (1)
			{
				for (int i = 0; i < T0 + 1; i++)
				{
					T1 += i*hisData[i];
					num1 += hisData[i];
				}
				if (num1 == 0)
					continue;
				for (int i = T0 + 1; i < 256; i++)
				{
					T2 += i*hisData[i];
					num2 += hisData[i];
				}
				if (num2 == 0)
					continue;

				T = (T1 / num1 + T2 / num2) / 2;

				if (T == T0)
					break;
				else
					T0 = T;
			}
			xThre.push_back(T);
		}
		imageThre.push_back(xThre);
	}
	dst = cv::Mat::zeros(src.rows, src.cols, CV_8UC1);
	for (int y = 0; y < src.rows; ++y)
	{
		for (int x = 0; x < src.cols; ++x)
		{

			double flagValue = imageThre[y / wndSize.height][x / wndSize.width];
			int value = src.at<uchar>(y, x);
			if (value > flagValue)
			{
				dst.at<uchar>(y, x) = 255;
			}
			else
			{
				dst.at<uchar>(y, x) = 0;
			}
		}
	}
}
//---------------------------------02.Niblack算法-------------------------------------
/** @brief 计算单通道灰度图像的平均值

@param src 单通道灰度图
*/
double thresholdSeg::GetMatAverage(const cv::Mat& src)
{
	CV_Assert(src.type() == CV_8UC1);
	double sum = 0.0;
	for (int y = 0; y < src.rows; ++y)
	{
		for (int x = 0; x < src.cols; ++x)
		{
			int value = src.at<uchar>(y, x);
			sum += value;
		}
	}

	return sum / (src.rows * src.cols);
}

/** @brief 计算单通道灰度图像的标准差

@param src 单通道灰度图
*/
double thresholdSeg::GetMatStdDev(const cv::Mat& src, double meanValue)
{
	CV_Assert(src.type() == CV_8UC1);
	double sum = 0.0;
	for (int y = 0; y < src.rows; ++y)
	{
		for (int x = 0; x < src.cols; ++x)
		{
			int value = src.at<uchar>(y, x);
			double var = (value - meanValue)*(value - meanValue);
			sum += var;
		}
	}

	double stdDev = std::sqrt(double(sum) / double(src.rows * src.cols));
	return stdDev;
}

void thresholdSeg::Niblack(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize,float k)
{
	CV_Assert(src.type() == CV_8UC1);
	CV_Assert((wndSize.width % 2 == 1) && (wndSize.height % 2 == 1));
	CV_Assert((wndSize.width <= src.cols) && (wndSize.height <= src.rows));

	cv::Mat flag = cv::Mat::zeros(src.rows, src.cols, CV_64FC1);
	for (int y = 0; y < src.rows; ++y)
	{
		for (int x = 0; x < src.cols; ++x)
		{
			//int value = src.at<uchar>(y, x);
			cv::Point center = cv::Point(x, y);
			cv::Point topLeftPoint = cv::Point(x - wndSize.width / 2, y - wndSize.height / 2);
			cv::Rect wnd = cv::Rect(topLeftPoint.x, topLeftPoint.y, wndSize.width, wndSize.height);
			//出界处理
			if (wnd.x < 0)
			{
				wnd.width = wnd.width + wnd.x;
				wnd.x = 0;
			}if (wnd.y < 0)
			{
				wnd.height = wnd.height + wnd.y;
				wnd.y = 0;
			}if (wnd.x + wnd.width > src.cols)
			{
				wnd.width = src.cols - wnd.x;
			}
			if (wnd.y + wnd.height > src.rows)
			{
				wnd.height = src.rows - wnd.y;
			}
			cv::Mat roiMat = src(wnd);
			double avgValue = GetMatAverage(roiMat);
			double dev = GetMatStdDev(roiMat, avgValue);

			
			//double flagValue = avgValue + k * dev;
			double flagValue = avgValue + k*dev;
			flag.at<double>(y, x) = flagValue;
		}
	}

	dst = cv::Mat::zeros(src.rows, src.cols, CV_8UC1);
	for (int y = 0; y < src.rows; ++y)
	{
		for (int x = 0; x < src.cols; ++x)
		{
			double flagValue = flag.at<double>(y, x);
			int value = src.at<uchar>(y, x);
			if (value > flagValue)
			{
				dst.at<uchar>(y, x) = 255;
			}
			else
			{
				dst.at<uchar>(y, x) = 0;
			}
		}
	}
}

//------------------------------------03.Sauvola算法---------------------------------------------------
//求区域内均值 integral即为积分图
//float thresholdSeg::fastMean(cv::Mat& integral, int x, int y, int window)
//{
//
//	int min_y = std::max(0, y - window / 2);
//	int max_y = std::min(integral.rows - 1, y + window / 2);
//	int min_x = std::max(0, x - window / 2);
//	int max_x = std::min(integral.cols - 1, x + window / 2);
//
//	int topright = integral.at<int>(max_y, max_x);
//	int botleft = integral.at<int>(min_y, min_x);
//	int topleft = integral.at<int>(max_y, min_x);
//	int botright = integral.at<int>(min_y, max_x);
//
//	float res = (float)((topright + botleft - topleft - botright) / (float)((max_y - min_y) *(max_x - min_x)));
//
//	return res;
//}



static int CalcMaxValue(int a, int b)
{
	return (a > b) ? a : b;
}

static double CalcMaxValue(double a, double b)
{
	return (a > b) ? a : b;
}

static int CalcMinValue(int a, int b)
{
	return (a < b) ? a : b;
}

static double CalcMinValue(double a, double b)
{
	return (a < b) ? a : b;
}


/** @brief SauvolaThresh二值算法

此代码不适用与分辨率较大的图像, 此bug准备有空再处理

@param src 单通道灰度图
@param dst 单通道处理后的图
@param k  threshold = mean*(1 + k*((std / 128) - 1))
@param wndSize 处理区域宽高, 一定是奇数

*/
void thresholdSeg::SauvolaThresh(const cv::Mat& src, cv::Mat& dst, const int k, const cv::Size wndSize)
{
	CV_Assert(src.type() == CV_8UC1);
	CV_Assert((wndSize.width % 2 == 1) && (wndSize.height % 2 == 1));
	CV_Assert((wndSize.width <= src.cols) && (wndSize.height <= src.rows));

	dst = cv::Mat::zeros(src.rows, src.cols, CV_8UC1);

	// 产生标志位图像
	unsigned long* integralImg = new unsigned long[src.rows * src.cols];
	unsigned long* integralImgSqrt = new unsigned long[src.rows * src.cols];
	std::memset(integralImg, 0, src.rows *src.cols * sizeof(unsigned long));
	std::memset(integralImgSqrt, 0, src.rows *src.cols * sizeof(unsigned long));

	// 计算直方图和图像值平方的和
	for (int y = 0; y < src.rows; ++y)
	{
		unsigned long sum = 0;
		unsigned long sqrtSum = 0;
		for (int x = 0; x < src.cols; ++x)
		{
			int index = y * src.cols + x;
			sum += src.at<uchar>(y, x);
			sqrtSum += src.at<uchar>(y, x)*src.at<uchar>(y, x);
			if (y == 0)
			{
				integralImg[index] = sum;
				integralImgSqrt[index] = sqrtSum;
			}
			else
			{
				integralImgSqrt[index] = integralImgSqrt[(y - 1)*src.cols + x] + sqrtSum;
				integralImg[index] = integralImg[(y - 1)*src.cols + x] + sum;
			}
		}
	}

	double diff = 0.0;
	double sqDiff = 0.0;
	double diagSum = 0.0;
	double iDiagSum = 0.0;
	double sqDiagSum = 0.0;
	double sqIDiagSum = 0.0;
	for (int x = 0; x < src.cols; ++x)
	{
		for (int y = 0; y < src.rows; ++y)
		{
			int xMin = CalcMaxValue(0, x - wndSize.width / 2);
			int yMin = CalcMaxValue(0, y - wndSize.height / 2);
			int xMax = CalcMinValue(src.cols - 1, x + wndSize.width / 2);
			int yMax = CalcMinValue(src.rows - 1, y + wndSize.height / 2);
			double area = (xMax - xMin + 1)*(yMax - yMin + 1);
			if (area <= 0)
			{
				// blog提供源码是biImage[i * IMAGE_WIDTH + j] = 255;但是i表示的是列, j是行
				dst.at<uchar>(y, x) = 255;
				continue;
			}

			if (xMin == 0 && yMin == 0)
			{
				diff = integralImg[yMax*src.cols + xMax];
				sqDiff = integralImgSqrt[yMax*src.cols + xMax];
			}
			else if (xMin > 0 && yMin == 0)
			{
				diff = integralImg[yMax*src.cols + xMax] - integralImg[yMax*src.cols + xMin - 1];
				sqDiff = integralImgSqrt[yMax * src.cols + xMax] - integralImgSqrt[yMax * src.cols + xMin - 1];
			}
			else if (xMin == 0 && yMin > 0)
			{
				diff = integralImg[yMax * src.cols + xMax] - integralImg[(yMin - 1) * src.cols + xMax];
				sqDiff = integralImgSqrt[yMax * src.cols + xMax] - integralImgSqrt[(yMin - 1) * src.cols + xMax];;
			}
			else
			{
				diagSum = integralImg[yMax * src.cols + xMax] + integralImg[(yMin - 1) * src.cols + xMin - 1];
				iDiagSum = integralImg[(yMin - 1) * src.cols + xMax] + integralImg[yMax * src.cols + xMin - 1];
				diff = diagSum - iDiagSum;
				sqDiagSum = integralImgSqrt[yMax * src.cols + xMax] + integralImgSqrt[(yMin - 1) * src.cols + xMin - 1];
				sqIDiagSum = integralImgSqrt[(yMin - 1) * src.cols + xMax] + integralImgSqrt[yMax * src.cols + xMin - 1];
				sqDiff = sqDiagSum - sqIDiagSum;
			}
			double mean = diff / area;
			double stdValue = sqrt((sqDiff - diff*diff / area) / (area - 1));
			double threshold = mean*(1 + k*((stdValue / 128) - 1));
			if (src.at<uchar>(y, x) < threshold)
			{
				dst.at<uchar>(y, x) = 0;
			}
			else
			{
				dst.at<uchar>(y, x) = 255;
			}

		}
	}

	delete[] integralImg;
	delete[] integralImgSqrt;
}

//void thresholdSeg::Sauvola(cv::Mat& inpImg, cv::Mat& resImg, int window, float k)
//{
//	cv::Mat integral;
//	int nYOffSet = 3;
//	int nXOffSet = 3;
//	cv::integral(inpImg, integral);  //计算积分图像
//	for (int y = 0; y < inpImg.rows; y += nYOffSet)
//	{
//		for (int x = 0; x < inpImg.cols; x += nXOffSet)
//		{
//
//			float fmean = fastMean(integral, x, y, window); float fthreshold = (float)(fmean*(1.0 - k));
//
//			int nNextY = y + nYOffSet;
//			int nNextX = x + nXOffSet;
//			int nCurY = y;
//			while (nCurY < nNextY && nCurY < inpImg.rows)
//			{
//				int nCurX = x;
//				while (nCurX < nNextX && nCurX < inpImg.cols)
//				{
//					uchar val = inpImg.at<uchar>(nCurY, nCurX) < fthreshold;
//					resImg.at<uchar>(nCurY, nCurX) = (val == 0 ? 0 : 255);
//					nCurX++;
//				}
//				nCurY++;
//			}
//
//		}
//	}
//
//	//return resImg;
//}
//-------------------------------04.Bernsen算法-----------------------------
/** @brief 得到矩阵中的最大值与最小值

@param m 单通道CV_8UC1类型矩阵
@param maxValue 最大值
@param minValue 最小值
*/
void thresholdSeg::GetMatMaxMin(const cv::Mat& m, int& maxValue, int& minValue)
{
	CV_Assert(m.type() == CV_8UC1);

	maxValue = INT_MIN;
	minValue = INT_MAX;

	for (int y = 0; y < m.rows; ++y)
	{
		for (int x = 0; x < m.cols; ++x)
		{
			int v = m.at<uchar>(y, x);
			if (v > maxValue) maxValue = v;
			if (v < minValue) minValue = v;
		}
	}
}

void thresholdSeg::Bernsen(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize)
{
	CV_Assert(src.type() == CV_8UC1);
	CV_Assert((wndSize.width % 2 == 1) && (wndSize.height % 2 == 1));
	CV_Assert((wndSize.width <= src.cols) && (wndSize.height <= src.rows));

	cv::Mat meanMat = cv::Mat::zeros(src.rows, src.cols, CV_8UC1);

	for (int y = wndSize.height / 2; y <= src.rows - wndSize.height / 2 - 1; ++y)
	{
		for (int x = wndSize.width / 2; x <= src.cols - wndSize.width / 2 - 1; ++x)
		{
			int value = src.at<uchar>(y, x);
			cv::Point center = cv::Point(x, y);
			cv::Point topLeftPoint = cv::Point(x - wndSize.width / 2, y - wndSize.height / 2);
			cv::Rect wnd = cv::Rect(topLeftPoint.x, topLeftPoint.y, wndSize.width, wndSize.height);
			int maxValue = 0;
			int minValue = 0;
			cv::Mat roiMat = src(wnd);
			GetMatMaxMin(roiMat, maxValue, minValue);
			int meanValue = (maxValue + minValue) / 2.0;
			meanMat.at<uchar>(y, x) = meanValue;
		}
	}

	// 阈值分割
	dst = cv::Mat::zeros(src.rows, src.cols, CV_8UC1);
	for (int y = 0; y < src.rows; ++y)
	{
		for (int x = 0; x < src.cols; ++x)
		{
			int value = src.at<uchar>(y, x);
			int meanValue = meanMat.at<uchar>(y, x);
			if (value > meanValue)
			{
				dst.at<uchar>(y, x) = 255;
			}
			else
			{
				dst.at<uchar>(y, x) = 0;
			}
		}
	}
}

void thresholdSeg::Bernsen(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize, int differMax, int meanMax)
{
	CV_Assert(src.type() == CV_8UC1);
	CV_Assert((wndSize.width % 2 == 1) && (wndSize.height % 2 == 1));
	CV_Assert((wndSize.width <= src.cols) && (wndSize.height <= src.rows));

	// 计算均值矩阵和差异矩阵
	cv::Mat meanMat = cv::Mat::zeros(src.rows, src.cols, CV_8UC1);
	cv::Mat differMat = cv::Mat::zeros(src.rows, src.cols, CV_8UC1);
	for (int y = wndSize.height / 2; y <= src.rows - wndSize.height / 2 - 1; ++y)
	{
		for (int x = wndSize.width / 2; x <= src.cols - wndSize.width / 2 - 1; ++x)
		{
			int value = src.at<uchar>(y, x);
			cv::Point center = cv::Point(x, y);
			cv::Point topLeftPoint = cv::Point(x - wndSize.width / 2, y - wndSize.height / 2);
			cv::Rect wnd = cv::Rect(topLeftPoint.x, topLeftPoint.y, wndSize.width, wndSize.height);
			int maxValue = 0;
			int minValue = 0;
			cv::Mat roiMat = src(wnd);
			GetMatMaxMin(roiMat, maxValue, minValue);
			int meanValue = (maxValue + minValue) / 2.0;
			int differValue = maxValue - minValue;
			meanMat.at<uchar>(y, x) = meanValue;
			differMat.at<uchar>(y, x) = differValue;
		}
	}

	// 赋值
	dst = cv::Mat::zeros(src.rows, src.cols, CV_8UC1);
	for (int y = 0; y < differMat.rows; ++y)
	{
		for (int x = 0; x < differMat.cols; ++x)
		{
			int differValue = differMat.at<uchar>(y, x);
			if (differValue > differMax)
			{
				// blog写的很迷糊, 直说meanValue是阈值
				// 本人认为是边界部分,可以是0,也可以是255
				dst.at<uchar>(y, x) = 255;
			}
			else if (differValue < differMax)
			{
				int meanValue = meanMat.at<uchar>(y, x);
				if (meanValue > meanMax)
				{
					dst.at<uchar>(y, x) = 255;
				}
				else
				{
					dst.at<uchar>(y, x) = 0;
				}
			}
			else
			{
				// TODO
				dst.at<uchar>(y, x) = 0;
			}
		}
	}
}
