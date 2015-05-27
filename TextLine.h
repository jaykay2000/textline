//Author : Jayant Kumar, jayant@umiacs.umd.edu, UMD college park

#ifndef TextLine_H_
#define TextLine_H_

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

using namespace std;

//Class for text line clustering objects and methods
class TextLine
{
private:

	//properties
	vector <float> centroidX;
	vector <float> centroidY;
	vector <float> ccWidths;
	vector <float> ccHeights;
	vector <float> xTop;
	vector <float> xBot;
	vector <float> yTop;
	vector <float> yBot;
	//thresholds
	float hbandthresh,vbandthresh;
	vector <float> thetamat;
	vector <float> quadmat;
	vector <float> tantheta; 
	vector <int> stquad;
	vector <int> quad;
	vector <int > qcount;
	int TL, BL, TR, BR;
	int **frndmat;
	int **CClabels;
	int **CClabels2;
	int numValidCC;
	int maxCCLabel;
	double orientation;
	float medH,medW;
	vector<vector<int>> mlines;
public : 
	TextLine()
	{

	}
	CvPoint GetCentroidOfCC(CvRect r);
	CvSeq* FilterCC(CvSeq *contours,int hThres, int wThres, int hThresmin, int wThresmin);
	int GetPropertyofCC(CvSeq *contours);
	int** FindLines();
	int GetCountInBands(int orig);
	int ResizeBand();
	int TextLine::ObtainQuad(int orig);
	int StraightOrientation(int orig,double vbandthresh);
	int GetQuadCounts();
	int FindMaxValue(vector<int> &vec);
	int FindMaxIndex(vector<int> &vec,int val);
	int GetMaxIndx(vector<int>& vec,int maxno);
	vector<int> GetQuadIndx(int quadnum);
	double ComputeAvgTheta(int quadno,int sIndex);
	int ComputeFloydAPSP(int **frndmat, int n);
	int FindComponents(int **frndmat, int n);
	int GetNumValidCCs();
	int GetMaxCCLabel(int ** labels,int n);
	double FindOrientation();
	int ObtainTransformedQuad(int orig);
	CvPoint RotateCounterClock(CvPoint r);
	int GetNumLines(){ return maxCCLabel+1;}//Number of textlines
	int CombineLines();
	int LeastSquareFit(vector<int> listcc,vector<float> &param);
	int ExtractLines(int **labels,int lineno, int &n, vector<int> &linecc);
	int FindExtremeCC(vector<int> linecc,int &leftCC,int &rightCC);
	inline int FindMeanCCValue(vector<int> linecc,float &meanCC);
	float FindDistanceToLine(vector<float> param1, int ind);
	//int MergeLines(int label1,int label2);
	int	MakeUniformLabel();
	int MarkLines(int label1,int label2);
	int MergeMarkedLines();
};
#endif 

