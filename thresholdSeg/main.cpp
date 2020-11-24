#include<opencv.hpp>

using namespace std;
using namespace cv;

cv::Mat OtsuAlgThreshold(Mat &image);
cv::Mat IterationThreshold(Mat src);
void Niblack(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize);
void Sauvola(cv::Mat& inpImg, cv::Mat& resImg, int window, float k);
void Bernsen(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize);
void Bernsen(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize, int differMax, int meanMax);
int main()
{
	string srcPath;
	srcPath = "C:\\Users\\SZJ\\Desktop\\������ͼ\\1.jpg";
	int alg = 0;   //0�����  1����������ֵ�㷨 2:Niblack�㷨  3��Sauvola�㷨  4��Bernsen�㷨
	Mat srcMat = imread(srcPath);
	Mat grayMat;
	Mat binaryMat;
	cvtColor(srcMat, grayMat, COLOR_BGR2GRAY);
	
	if(alg==0)
	{
		//���
		binaryMat = OtsuAlgThreshold(grayMat);
	}
	else if (alg == 1)
	{
		//������ֵ�㷨
		binaryMat = IterationThreshold(grayMat);
	}
	else if (alg == 2)
	{
		//Niblack�㷨
		Niblack(srcMat, binaryMat, Size(100, 100));
	}
	else if (alg == 3)
	{
		//Sauvola��
		Sauvola(srcMat, binaryMat, 100, 0.4);
	}
	else if (alg == 4)
	{
		Bernsen(srcMat, binaryMat, Size(100,100));
	}
	
	imshow("src", srcMat);
	imshow("gray", grayMat);
	imshow("binary", binaryMat);
	waitKey(0);
	return 0;
}
//------------------------------00.Otsu�㷨-------------------------
Mat OtsuAlgThreshold(Mat &image)
{
	if (image.channels() != 1)
	{
		cout << "Please input Gray-image!" << endl;
	}
	int T = 0; //Otsu�㷨��ֵ  
	double varValue = 0; //��䷽���м�ֵ����
	double w0 = 0; //ǰ�����ص�����ռ����  
	double w1 = 0; //�������ص�����ռ����  
	double u0 = 0; //ǰ��ƽ���Ҷ�  
	double u1 = 0; //����ƽ���Ҷ�  
	double Histogram[256] = { 0 }; //�Ҷ�ֱ��ͼ���±��ǻҶ�ֵ�����������ǻҶ�ֵ��Ӧ�����ص�����  
	uchar *data = image.data;

	double totalNum = image.rows*image.cols; //��������

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
		//ÿ�α���֮ǰ��ʼ��������  
		w1 = 0;       u1 = 0;       w0 = 0;       u0 = 0;
		//***********����������ֵ����**************************  
		for (int j = 0; j <= i; j++) //�������ָ�ֵ����  
		{
			w1 += Histogram[j];   //�����������ص�����  
			u1 += j*Histogram[j]; //�������������ܻҶȺ�  
		}
		if (w1 == 0) //�����������ص���Ϊ0ʱ�˳�  
		{
			break;
		}
		u1 = u1 / w1; //��������ƽ���Ҷ�  
		w1 = w1 / totalNum; // �����������ص�����ռ����
							//***********����������ֵ����**************************  

							//***********ǰ��������ֵ����**************************  
		for (int k = i + 1; k < 255; k++)
		{
			w0 += Histogram[k];  //ǰ���������ص�����  
			u0 += k*Histogram[k]; //ǰ�����������ܻҶȺ�  
		}
		if (w0 == 0) //ǰ���������ص���Ϊ0ʱ�˳�  
		{
			break;
		}
		u0 = u0 / w0; //ǰ������ƽ���Ҷ�  
		w0 = w0 / totalNum; // ǰ���������ص�����ռ����  
		//***********ǰ��������ֵ����**************************  

		//***********��䷽�����******************************  
		double varValueI = w0*w1*(u1 - u0)*(u1 - u0); //��ǰ��䷽�����  
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


//--------------------------------------01.������ֵ�㷨------------------------------------------------
Mat IterationThreshold(Mat src)
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
//---------------------------------02.Niblack�㷨-------------------------------------
/** @brief ���㵥ͨ���Ҷ�ͼ���ƽ��ֵ

@param src ��ͨ���Ҷ�ͼ
*/
static double GetMatAverage(const cv::Mat& src)
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

/** @brief ���㵥ͨ���Ҷ�ͼ��ı�׼��

@param src ��ͨ���Ҷ�ͼ
*/
static double GetMatStdDev(const cv::Mat& src, double meanValue)
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

void Niblack(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize)
{
	CV_Assert(src.type() == CV_8UC1);
	CV_Assert((wndSize.width % 2 == 1) && (wndSize.height % 2 == 1));
	CV_Assert((wndSize.width <= src.cols) && (wndSize.height <= src.rows));

	cv::Mat flag = cv::Mat::zeros(src.rows, src.cols, CV_64FC1);
	for (int y = wndSize.height / 2; y <= src.rows - wndSize.height / 2 - 1; ++y)
	{
		for (int x = wndSize.width / 2; x <= src.cols - wndSize.width / 2 - 1; ++x)
		{
			int value = src.at<uchar>(y, x);
			cv::Point center = cv::Point(x, y);
			cv::Point topLeftPoint = cv::Point(x - wndSize.width / 2, y - wndSize.height / 2);
			cv::Rect wnd = cv::Rect(topLeftPoint.x, topLeftPoint.y, wndSize.width, wndSize.height);
			cv::Mat roiMat = src(wnd);
			double avgValue = GetMatAverage(roiMat);
			double dev = GetMatStdDev(roiMat, avgValue);

			// ������0.2
			double flagValue = avgValue + 0.2 * dev;
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

//------------------------------------03.Sauvola�㷨---------------------------------------------------
//�������ھ�ֵ integral��Ϊ����ͼ
float fastMean(cv::Mat& integral, int x, int y, int window)
{

	int min_y = std::max(0, y - window / 2);
	int max_y = std::min(integral.rows - 1, y + window / 2);
	int min_x = std::max(0, x - window / 2);
	int max_x = std::min(integral.cols - 1, x + window / 2);

	int topright = integral.at<int>(max_y, max_x);
	int botleft = integral.at<int>(min_y, min_x);
	int topleft = integral.at<int>(max_y, min_x);
	int botright = integral.at<int>(min_y, max_x);

	float res = (float)((topright + botleft - topleft - botright) / (float)((max_y - min_y) *(max_x - min_x)));

	return res;
}


void Sauvola(cv::Mat& inpImg, cv::Mat& resImg, int window, float k)
{
	cv::Mat integral;
	int nYOffSet = 3;
	int nXOffSet = 3;
	cv::integral(inpImg, integral);  //�������ͼ��
	for (int y = 0; y < inpImg.rows; y += nYOffSet)
	{
		for (int x = 0; x < inpImg.cols; x += nXOffSet)
		{

			float fmean = fastMean(integral, x, y, window); float fthreshold = (float)(fmean*(1.0 - k));

			int nNextY = y + nYOffSet;
			int nNextX = x + nXOffSet;
			int nCurY = y;
			while (nCurY < nNextY && nCurY < inpImg.rows)
			{
				int nCurX = x;
				while (nCurX < nNextX && nCurX < inpImg.cols)
				{
					uchar val = inpImg.at<uchar>(nCurY, nCurX) < fthreshold;
					resImg.at<uchar>(nCurY, nCurX) = (val == 0 ? 0 : 255);
					nCurX++;
				}
				nCurY++;
			}

		}
	}

	//return resImg;
}
//-------------------------------04.Bernsen�㷨-----------------------------
/** @brief �õ������е����ֵ����Сֵ

@param m ��ͨ��CV_8UC1���;���
@param maxValue ���ֵ
@param minValue ��Сֵ
*/
static void GetMatMaxMin(const cv::Mat& m, int& maxValue, int& minValue)
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

void Bernsen(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize)
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

	// ��ֵ�ָ�
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

void Bernsen(const cv::Mat & src, cv::Mat & dst, cv::Size wndSize, int differMax, int meanMax)
{
	CV_Assert(src.type() == CV_8UC1);
	CV_Assert((wndSize.width % 2 == 1) && (wndSize.height % 2 == 1));
	CV_Assert((wndSize.width <= src.cols) && (wndSize.height <= src.rows));

	// �����ֵ����Ͳ������
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

	// ��ֵ
	dst = cv::Mat::zeros(src.rows, src.cols, CV_8UC1);
	for (int y = 0; y < differMat.rows; ++y)
	{
		for (int x = 0; x < differMat.cols; ++x)
		{
			int differValue = differMat.at<uchar>(y, x);
			if (differValue > differMax)
			{
				// blogд�ĺ��Ժ�, ֱ˵meanValue����ֵ
				// ������Ϊ�Ǳ߽粿��,������0,Ҳ������255
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