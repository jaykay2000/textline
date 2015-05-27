//Author : Jayant Kumar, jayant@umiacs.umd.edu, UMD college park

#include "stdafx.h"

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>

#include "TextLine.h"
#include "graph.h"

//Allocate a 2-D MAT
int** Allocate2DMatrixInt(int Row,int Col)
{
	int **mat = new int*[Row];
	if(mat == NULL)
	{
		cerr<<"Allocation failed";
		exit(0);
	}
	for(int rowIndex = 0;rowIndex<Row;rowIndex++)
	{
		mat[rowIndex] = new int[Col];
	}
	return mat;
}

/*
Obtain centroid using Bounding Box 
*/
CvPoint TextLine::GetCentroidOfCC(CvRect r)
{
	/* get the centroid of the  component */
	CvPoint point = cvPoint(
		r.x + (r.width / 2),
		r.y + (r.height / 2)
	);
	return point;
}

/*
Filter text like Connected-Components 
*/
CvSeq* TextLine::FilterCC(CvSeq *contours,int hThres, int wThres, int hThresmin, int wThresmin)
{
	CvSeq *ptr;
	CvRect r;
	ptr = contours;
	int cnt = 0;
    while (ptr != NULL){
		r = cvBoundingRect(ptr, 1);
		if (r.height > hThres || r.width > wThres || (r.height < hThresmin && r.width < wThresmin)) { //CC is small enough to resemble text
            if (ptr == contours) { //if the first CC itself is not valid, change contours
                contours = ptr->h_next;
                ptr = contours;
                ptr->h_prev = NULL;
            } else { //for second CC onwards just skip the invalid CCs
                ptr->h_prev->h_next = ptr->h_next;
                if (ptr->h_next != NULL) //if its not the last element
                    ptr->h_next->h_prev = ptr->h_prev;
                ptr = ptr->h_next;
            }
        } else {
            ptr = ptr->h_next;
			++cnt;
        }
    }
	numValidCC = cnt;
	return contours;
}

/*
Obtain properties of connected components
*/
int TextLine::GetPropertyofCC(CvSeq *contours)
{
	CvSeq *ptr;
	CvRect r;
	ptr = contours;
	CvPoint point;
	float xT,yT,xB,yB;
    while (ptr != NULL){
		r = cvBoundingRect(ptr, 1);
		point = GetCentroidOfCC(r);
		centroidY.push_back(point.y);
		centroidX.push_back(point.x);
		ccHeights.push_back(r.height);
		ccWidths.push_back(r.width);
		xTop.push_back(point.x - r.width/2);
		xBot.push_back(point.y - r.height/2);
		yTop.push_back(point.x + r.width/2);
		yBot.push_back(point.y + r.height/2);
		ptr = ptr->h_next;
    }
	return 0;
}


inline float Median(vector<float> vec)
{
    sort(vec.begin(),vec.end());
    int indx;
    int size = vec.size();
    float rval,val1,val2;
    if(size%2 == 0) // even
    {
        val1 = vec.at(size/2);
        val2 = vec.at(size/2 + 1);
        rval = (val1 + val2)/2;
    }
    else
    {
        indx = size/2;
        rval = vec.at(indx);
    }
    return rval;
}

//Update TL, BL etc
int TextLine::GetCountInBands(int orig)
{
	int size = centroidX.size();
	int i = 0; 
	float tarcenX,tarcenY;
	float srccenX = centroidX.at(orig);
	float srccenY = centroidY.at(orig);
	TL = 0,BL = 0, TR = 0,BR = 0;
	for(i=0;i < size;i++)
	{
		if(i != orig)
		{
			tarcenX = centroidX.at(i);
            tarcenY = centroidY.at(i);
            if(abs(tarcenY-srccenY)< vbandthresh && abs(tarcenX-srccenX)< hbandthresh)
                if(srccenX - tarcenX > 0 && srccenY-tarcenY > 0)
                    TL = TL+1; // top left
                else if (srccenX - tarcenX > 0 && srccenY-tarcenY < 0) 
                    BL = BL+1; //bottom left 
                else if(srccenX - tarcenX < 0 && srccenY-tarcenY > 0)
                    TR=TR+1; //top right
                else
                    BR = BR+1; //bottom right
		}
	}

	return 0;
}


//Adapt to capture more components
int TextLine::ResizeBand()
{
	//use count and resize
	int flag =  0;
	
	//if skew then diagonally  should have more elements 
	if (((BR + TL < 3) || (BL + TR < 3)) ) 
	{
        vbandthresh = 1.1*vbandthresh;
        flag = 1;
	}
	//right most or left most
    if (((TL +BL < 3) && (TR + BR < 1))  || ((TR + BR < 3) && (TL +BL < 1))) 
	{
        hbandthresh = 1.25*hbandthresh;
        flag = 1;
	}
	
	if ((TL + BL + TR + BR < 4) && flag == 0)
	{
            hbandthresh = 1.3*hbandthresh; // resize
            vbandthresh = 1.1*vbandthresh;
	}
	return 0;
}

double inline rad2deg(double rad)
{
	double const pi = 3.14;
	return (180/pi)*rad;
}


int TextLine::ObtainQuad(int orig)
{

	int size = centroidX.size();
	int i = 0; 
	float tarcenX,tarcenY,tarXTop,tarXBot,tarYBot,tarYTop;
	float srccenX = centroidX.at(orig);
	float srccenY = centroidY.at(orig);
	float srcXTop = xTop.at(orig);
	float srcXBot = xBot.at(orig);
	float srcYTop = yTop.at(orig);
	float srcYBot = yBot.at(orig);
	float distmatX,distmatY;
	double degang,ttheta;
	int pointcnt = 0;
	for(i=0; i < size;i++)
	{
		if(i != orig)
		{
			tarcenX = centroidX.at(i);
            tarcenY = centroidY.at(i);
			tarXTop = xTop.at(i); tarXBot = xBot.at(i);
			//yT = yTop.at(i);
			if((abs(tarcenY-srccenY) < vbandthresh) && ( (abs(tarcenX-srccenX)< hbandthresh) || (abs(srcXTop-tarXBot) < 0.5*hbandthresh) || (abs(srcXBot-tarXTop) < 0.5*hbandthresh) ))
			{
				pointcnt++;
				distmatX = srccenX-tarcenX;
                distmatY = srccenY-tarcenY;
                if (distmatX != 0)
                    ttheta = distmatY/distmatX;
                else
                    ttheta = 100000; //some large number
				
				tantheta[i] = ttheta; // may be positive or negative	
				
				degang = rad2deg(atan(abs(ttheta)));     
				if(ttheta > 0) // 1 or 3 subquadrant  
				{
                        if(degang < 45)     
                             quad[i] = 4;
                        else                                
                             quad[i] = 3;
				}      
                else 
				{
                        if(degang < 45)  
                            quad[i] = 1;
                        else
                            quad[i] = 2;
				}
				
				if(degang < 10)
					stquad[i] = 1;
				else
					stquad[i] = 0;
			}        
				
			
		}
		else
		{
			quad[i] = 0;
			stquad[i]  = 0 ;
		}
                
	}
	return 0;
}

CvPoint inline TextLine::RotateCounterClock(CvPoint r)
{
	/* get the centroid of the  component */
	return cvPoint(r.x*cos(orientation) - r.y*sin(orientation) , r.x*sin(orientation) + r.y*cos(orientation));
}

//Orientation computation after transformation 
int TextLine::ObtainTransformedQuad(int orig)
{
	int size = centroidX.size();
	int i = 0; 
	CvPoint tarpoint,srcpoint;
	float tarcenX,tarcenY,tarXTop,tarXBot,tarYBot,tarYTop;
	//This is origin, so no changes
	float srccenX = centroidX.at(orig);
	float srccenY = centroidY.at(orig);

	srcpoint =  RotateCounterClock(cvPoint( xTop.at(orig),xBot.at(orig)));
	float srcXTop = srcpoint.x;
	float srcXBot = srcpoint.y;
	float distmatX,distmatY;
	double degang,ttheta;
	int pointcnt = 0;
	for(i=0; i < size;i++)
	{
		if(i != orig)
		{
			tarpoint =  RotateCounterClock(cvPoint(centroidX.at(i)-centroidX.at(orig),centroidY.at(i)-centroidY.at(orig))); //Transform centroid
			tarcenX = tarpoint.x + centroidX.at(orig);
            tarcenY = tarpoint.y + centroidY.at(orig);
			//tarpoint =  RotateCounterClock(cvPoint( xTop.at(i),xBot.at(i))); //Transform XTOP, XBOT
			//tarXTop = tarpoint.x; 
			//tarXBot = tarpoint.y;
			if((abs(tarcenY-srccenY) < vbandthresh) && ( (abs(tarcenX-srccenX)< hbandthresh) ))
			{
				pointcnt++;
				distmatX = srccenX-tarcenX;
                distmatY = srccenY-tarcenY;
                if (distmatX != 0)
                    ttheta = distmatY/distmatX;
                else
                    ttheta = 100000; //some large number
				
				tantheta[i] = ttheta; // may be positive or negative	
				
				degang = rad2deg(atan(abs(ttheta)));     
				if(ttheta > 0) // 1 or 3 subquadrant  
				{
                        if(degang < 45)     
                             quad[i] = 4;
                        else                                
                             quad[i] = 3;
				}      
                else 
				{
                        if(degang < 45)  
                            quad[i] = 1;
                        else
                            quad[i] = 2;
				}
				
				if(degang < 10)
					stquad[i] = 1;
				else
					stquad[i] = 0;
			}        
				
			
		}
		else
		{
			quad[i] = 0;
			stquad[i]  = 0 ;
		}
                
	}
	return 0;
}



int SumVec(vector<int> &vec)
{
	int sum = 0,i;
	for(i=0 ; i < vec.size();i++)
		sum = sum +  ((vec.at(i)>0)?1:0);
	return sum;
}

int TextLine::StraightOrientation(int orig,double vbandthresh)
{
		thetamat[orig] = (3.14/180)*0.1; // 0.1 degree
        quadmat.at(orig) = 5;
		int k = 0;
		int i = orig;int p;
		int size = centroidX.size();
		vector<int> stindx1,stindx2,stindx;
		for(k=0;k < size;k++)
		{
			if(k != orig)
			{
				if(stquad.at(k)  == 1)
				{ 
                    frndmat[i][k] = frndmat[k][i] = 1;
				}
			}
		}	 
	return 0; 
}

int TextLine::GetQuadCounts()
{
	int size = quad.size();
	int val = 0,i;
	for( i=0;i < 4;i++)
		qcount[i] = 0;
	
	for( i=0;i<size;i++)
	{
		val = quad[i]-1;
		if(val >= 0)
			qcount[val] = qcount[val] + 1;
	}
	return 0;
}

inline int TextLine::FindMaxValue(vector<int> &vec)
{
	int maxval = -10000;
	for(int i=0;i< vec.size(); i++){
		if (vec[i] > maxval)
			maxval = vec[i];
	}

	return maxval;
}

inline int TextLine::FindMaxIndex(vector<int> &vec, int val)
{
	for(int i=0;i< vec.size(); i++){
		if (vec[i] == val)
			return i;
	}

	return -1;
}

inline vector<int> TextLine::GetQuadIndx(int quadnum)
{
	vector<int> indx;
	for(int k = 0;k < quad.size();k++)
		if(quad.at(k) == quadnum)
			indx.push_back(k);

	return indx;
}

int TextLine::GetNumValidCCs()
{
	return numValidCC;
}

int TextLine::GetMaxIndx(vector<int>& vec,int maxno)
{
	vector<int> maxindx;
	double minstdv = 1000000;
	vector<int> qindx;
	double mean,stdv;
	int quadno,min_indx,k,p;
	for(int k = 0;k < vec.size();k++)
		if(vec.at(k) == maxno)
			maxindx.push_back(k);
	
	if(maxindx.size() < 2)
		return maxindx[0];
	else 
	{
		for( p = 0; p < maxindx.size() ; p++)
		{
			quadno = maxindx[p]+1;
			qindx = GetQuadIndx(quadno);
			mean = 0;
			//mean
			for( k=0;k < qindx.size();k++)
					mean = mean + abs(tantheta[qindx[k]]);
			mean = mean/qindx.size();
			//std dev to decide the correct orientation
			stdv = 0;
			for( k=0;k < qindx.size();k++)
					stdv = stdv + pow(abs(mean - abs(tantheta[qindx[k]])),2); // tantheta may be negative
			stdv = stdv/qindx.size();
			if (minstdv > stdv)
			{
				minstdv = stdv;
				min_indx = p;
			}
		}
		return maxindx[min_indx];
	}
}

inline double TextLine:: ComputeAvgTheta(int quadno,int sIndex)
{
	int size = centroidX.size();
	double avgtheta = 0;
	int cnt =0;
	double SUMx = centroidX[sIndex], SUMy = centroidY[sIndex], SUMxx = centroidX[sIndex]*centroidX[sIndex],SUMxy = centroidX[sIndex]*centroidY[sIndex];
	for (int k =0;k < size;k++)
	{
		if(quad[k] == quadno)
		{
		//fscanf (infile, "%lf %lf", &x[i], &y[i]);
			SUMx = SUMx + centroidX[k];
			SUMy = SUMy + centroidY[k];
			SUMxy = SUMxy + centroidX[k]*centroidY[k];
			SUMxx = SUMxx + centroidX[k]*centroidX[k];
			cnt = cnt + 1;
		}
	}
	int n = cnt + 1; //extra 1 for source 
	double slope = ( SUMx*SUMy - n*SUMxy )/( SUMx*SUMx - n*SUMxx );
	return slope;
}


//Distance in connected points should be less
int TextLine::ComputeFloydAPSP(int **frndmat, int n)
{
  int i,j,k;
   for (i=0; i < n; i++)
   {
       for (j=0; j < n; j++)
        {
            for (k=0; k < n; k++)
            {
				if ( frndmat[i][j] == 1 && frndmat[j][k] == 1)
                {
                    frndmat[i][k] = frndmat[k][i] = 1;
                }
            }
        }
  }
  return 0;
} 

int TextLine::FindComponents(int **frndmat, int n)
{
	vector<int> ccComp ;
	vector<int> visited;
	visited.resize(n);
	int currentlab = -1;
	CClabels = new int*[1];
	CClabels[0] = new int[n];
	//Use BFS for iteratively finding the components in the graph
	Graph g(n);
	for (int i = 0; i < n ; i++)
		for (int k = i+1; k < n ; k++)
		{
			if (frndmat[i][k] == 1)
			{ 
				g.addEdge(i+1,k+1); // node num starting from 1
			}
		}
	for(int k = 0; k < n; k++)
	{
		if (visited[k] == 0)
		{
			ccComp = g.BFS(k+1);
			currentlab = currentlab+1;
			for(int p = 0; p < ccComp.size(); p++)
			{
				visited[ccComp[p]-1] = 1; 
				CClabels[0][ccComp[p]-1] = currentlab;
			}
			visited[k] = 1; // be careful using at
			CClabels[0][k] = currentlab;
		}
	}
  return 0;
}

double TextLine::FindOrientation()
{
	int size = centroidX.size();
	int sIndex = 0, tIndex = 0, p =0;
	int stpointcnt = 0; int pointcnt = 0;
	float medH = Median(ccHeights);
	float medW = Median(ccWidths);
	float vbandthresh1 = 0.6*medH;
	float hbandthresh1 = 3.5*medW;
	thetamat.resize(size);
	quadmat.resize(size);
	qcount.resize(5);
	for (sIndex = 0 ; sIndex < size; sIndex++)
	{
		stquad.clear(); tantheta.clear(); quad.clear(); 
		stquad.resize(size); tantheta.resize(size);quad.resize(size);
		hbandthresh = hbandthresh1;
		vbandthresh = vbandthresh1;
		TL = 0;BL=0;TR=0;BR=0; // reset
		GetCountInBands(sIndex); // Now BL,TL etc should have been updated
		ResizeBand(); // update hbandthres and vbandthresh if required
		ObtainQuad(sIndex);
		pointcnt = SumVec(quad); 
		stpointcnt = SumVec(stquad); //cout<<"\nstpointcnt : "<<stpointcnt;
		if (pointcnt < 2)
		{
			quadmat[sIndex] = -1; //unrealiable
			thetamat[sIndex] = -1;
			continue;
		}
		if (stpointcnt >  2*(pointcnt+ 0.0)/3)
		{
			quadmat[sIndex] = 5; //straight orientation
		}
		else
		{
			GetQuadCounts(); 
			int maxall = FindMaxValue(qcount);	
			//STDDEV for multiple quadnum
			int quadnum = GetMaxIndx(qcount,maxall)+1;
			quadmat[sIndex] = quadnum;
			double avgtheta = ComputeAvgTheta(quadnum,sIndex); // average of the angle in radians
			//double avgttheta = tan(avgtheta); // tan theta of average angle
			thetamat[sIndex] = abs(avgtheta);
		}
	}
	//Now Compute Dominant one (For textlines)
	qcount.clear(); 
	qcount.resize(5);
	for(int i=0;i < size;i++)
	{
		if(quadmat[i] > 0)
			qcount[quadmat[i]-1]++;
	}
	int qmaxval = FindMaxValue(qcount);
	int qmax = FindMaxIndex(qcount,qmaxval);
	double sumtheta = 0;
	for(int i=0;i < size;i++)
	{
		if(quadmat[i] == qmax+1)
			sumtheta = sumtheta + thetamat[i];
	}
	orientation = sumtheta/qmaxval;
	return qmax;
}

int TextLine::LeastSquareFit(vector<int> listcc, vector<float> &param) 
{
  double SUMx, SUMy, SUMxy, SUMxx, SUMres ;
  int i,n;
  n = listcc.size();
  SUMx = 0; SUMy = 0; SUMxy = 0; SUMxx = 0;
  for (i=0; i < n; i++) 
  {
      SUMx = SUMx + centroidX[listcc[i]];
	  SUMy = SUMy + centroidY[listcc[i]];
	  SUMxy = SUMxy + centroidX[listcc[i]]*centroidY[listcc[i]];
	  SUMxx = SUMxx + centroidX[listcc[i]]*centroidX[listcc[i]];
  }
  param[0] = ( ( SUMx*SUMy - n*SUMxy )/( SUMx*SUMx - n*SUMxx )); //slope
  param[1] = (( SUMy-param[0]*SUMx )/n); //y-intercept  
  return 0;
}

//Given labels , it gives the indices of a particular line 
//lineno is the label of line you want, n has number of components in line
int TextLine::ExtractLines(int **labels, int lineno, int &noLineCC, vector<int> &linecc)
{
	int cnt = 0;
	for(int i = 0 ; i < numValidCC;i++)
	{
		if(labels[0][i] == lineno)
		{
			linecc.push_back(i);
			cnt++;
		}
	}
	noLineCC = cnt;
	return 0;
}

inline int FindlabelinCC(int lab, int **cclab,int labno)
{
	int found = 0;
	for(int i =0;i < labno; i++)
	{
		if(cclab[0][i] == lab)
		{
			found = 1;
			break;
		}
	}
	return found;
}

//Find extreme CCs of line
int TextLine::FindExtremeCC(vector<int> linecc,int &leftCC,int &rightCC)
{
	int maxX=-1, minX = 1000000;
	leftCC = -1; rightCC = -1;
	for(int i=0;i<linecc.size();i++)
	{
		if(centroidX[linecc[i]] < minX)
		{
			minX = centroidX[linecc[i]];
			leftCC = linecc[i];
		}
		if(centroidX[linecc[i]] > maxX)
		{
			maxX = centroidX[linecc[i]];
			rightCC = linecc[i];
		}
	}
	return 0;
}

int TextLine::MarkLines(int label1,int label2)
{
	//first check if the labels are already present
	for(int i=0;i < mlines.size();i++)
	{
		vector<int> linelabels = mlines[i];
		for(int j=0; j < linelabels.size();j++)
		{
			if(linelabels[j] == label1 )
			{
				linelabels.push_back(label2);//add another label
				mlines[i] = linelabels; // update the entry
				return 0;
			}
			else if (linelabels[j] == label2 )
			{
				linelabels.push_back(label1);
				mlines[i] = linelabels; //update
				return 0;
			}
		}
	}
	//create a new vector
	vector<int> linelabels ;
	linelabels.push_back(label1);
	linelabels.push_back(label2);
	mlines.push_back(linelabels);
	return 0;
}


//merge two  lines
int TextLine::MergeMarkedLines()
{
	for(int i=0;i < mlines.size(); i++)
	{
		vector<int> linelabels = mlines[i];
		int minlabel = *(min_element(linelabels.begin(),linelabels.end()));
		for(int j=0; j < numValidCC;j++)
		{
			for(int k=0; k < linelabels.size();k++)
				if(CClabels[0][j] == linelabels[k]) //label is present
					CClabels2[0][j] = minlabel;
		}
	}
	return 0;
}

inline int TextLine::FindMeanCCValue(vector<int> linecc,float &meanCC)
{
	float sumvalue = 0;
	for(int i=0;i < linecc.size();i++)
	{
		 sumvalue = sumvalue + centroidX[linecc[i]];
	}
	meanCC = sumvalue/linecc.size();
	return 0;
}
//point to line distance for slope-intercept form
float TextLine::FindDistanceToLine(vector<float> param, int ind)
{
	float num = abs(centroidY[ind] -param[0]*centroidX[ind] - param[1]);
	return num/sqrt(param[0]*param[0]+1);
}

//Copy the matrix
int CopyMat(int **matS,int **matT,int nr,int nc)
{
	for(int rowIndex = 0;rowIndex < nr; rowIndex++)
	{
		for(int colIndex = 0;colIndex < nc; colIndex++)
		{
			matT[rowIndex][colIndex] = matS[rowIndex][colIndex];
		}
	}
	return 0;
}



inline int GetNextHigherLabel(int lab, int **cclab,int labno)
{
	int min = 1000;
	int mindiff = 1000;
	int minlab = -1;
	for(int i =0;i < labno; i++)
	{
		if((cclab[0][i] - lab) > 0 && (cclab[0][i] - lab) < mindiff)
		{
			minlab = cclab[0][i];
			mindiff = (cclab[0][i] - lab);
		}
	}
	return minlab;
}

int	TextLine::MakeUniformLabel()
{
	//check for missing label numbers and make it uniform
	int maxLab = GetMaxCCLabel(CClabels,numValidCC);
	for(int k =0; k <= maxLab;k++)
	{
		
		if (FindlabelinCC(k,CClabels,numValidCC) == 0)
		{
			int lab = GetNextHigherLabel(k,CClabels,numValidCC);
			if(lab!=k)
			{
				for(int i = 0 ; i < numValidCC; i++)
				{
					if(CClabels[0][i] == lab)
						CClabels[0][i] = k;

				}
			}
			maxLab = GetMaxCCLabel(CClabels,numValidCC);
		}
	}
	return 0;
}

//Group broken lines obtained so far using least square fitting...
int TextLine::CombineLines()
{
	int noLineCC1, noLineCC2;
	float meanCC1,meanCC2,dist1,dist2;
	int leftCC1,rightCC1,leftCC2,rightCC2;
	maxCCLabel = GetMaxCCLabel(CClabels,numValidCC);//maxlabel = line-1 (starts with 0)
	CClabels2 = Allocate2DMatrixInt(1,numValidCC);
	CopyMat(CClabels,CClabels2,1,numValidCC);
	for(int i = 0 ; i < maxCCLabel;i++)
	{
		vector<int> linecc1; 
		vector<float> param1; param1.resize(2);
		ExtractLines(CClabels, i, noLineCC1, linecc1);
		if (noLineCC1 < 2)//just the single CC
		{
			//continue;
			param1[0]=0;//slope = 0
			param1[1]=centroidY[linecc1[0]]; //put the only CC's x coord
		}
		else
		{
			LeastSquareFit(linecc1,param1);
		}
		FindExtremeCC(linecc1,leftCC1,rightCC1);
		FindMeanCCValue(linecc1,meanCC1);
		
		for(int j = i+1;j <= maxCCLabel;j++)
		{
			vector<int> linecc2;
			vector<float> param2;param2.resize(2);
			ExtractLines(CClabels, j, noLineCC2, linecc2);
			if (noLineCC2 < 2)
			{
				//continue;
				param2[0]=0;//slope = 0
				param2[1]=centroidY[linecc2[0]]; //put the only CC's x coord
			}
			else
			{
				LeastSquareFit(linecc2,param2);
			}
			FindExtremeCC(linecc2,leftCC2,rightCC2);//should handle single CC
			FindMeanCCValue(linecc2,meanCC2);
			//Find the distance from one extreme to another line and see if its close enough to combine
			if (meanCC1 > meanCC2)//line 2 is left, line 1 is right
			{
				//cout<<"\n dist : "<<centroidX[leftCC1] - centroidX[rightCC2];
				//line2 left, line 1 right
				//if ((centroidX[leftCC1] - centroidX[rightCC2]) > 5*medW)
				//	continue;
				dist1 = FindDistanceToLine(param1,rightCC2);// distance from a point in line 2 to line 1
				dist2 = FindDistanceToLine(param2,leftCC1);
				if(noLineCC1 < 2 && dist2 < medH)
					MarkLines(i,j); 
				else if(noLineCC2 < 2 && dist1 < medH)
					MarkLines(i,j);
				else if ((dist1 < medH || dist2 < medH) && noLineCC1 > 2 && noLineCC2 > 2)
					MarkLines(i,j);//make the label same in CClabels
			}
			else//line 1 is left, line 2 is right
			{
				//cout<<"\n dist : "<<centroidX[leftCC2] - centroidX[rightCC1];
				//if ((centroidX[leftCC2] - centroidX[rightCC1]) > 10*medW) 
				//	continue;
				//Find distance from rightmost of line1 (left) to the right line(line 2)
				dist1 = FindDistanceToLine(param2,rightCC1); // distance from point in line 1 to line 2
				dist2 = FindDistanceToLine(param1,leftCC2); // distance from point in line 2 to line 1
				if(noLineCC2 < 2 && dist2 < medH && dist1 < 2*medH)
					MarkLines(i,j); 
				else if(noLineCC1 < 2 && dist1 < medH && dist2 < 2*medH)
					MarkLines(i,j);
				else if((dist1 < medH || dist2 < medH) && noLineCC1 > 2 && noLineCC2 > 2)
					MarkLines(i,j);//make the label same in CClabels
			}
		}
	}
	MergeMarkedLines();
	CopyMat(CClabels2,CClabels,1,numValidCC);
	MakeUniformLabel();
	maxCCLabel = GetMaxCCLabel(CClabels,numValidCC);//maxlabel = line-1 (starts with 0)
	for(int i=0;i<numValidCC;i++)
		cout<<"\t"<<CClabels[0][i];
	return 0;
}

//compute distance matrix
int** TextLine::FindLines()
{
	int size = centroidX.size();
	FindOrientation(); // find the orientation of whole document using local analysis
	double srcRVal = 0, tarRVal = 0, eRDist = 0,eDistSquare = 0;
	double srcCVal = 0, tarCVal = 0,eCDist = 0;
	float hCC, wCC;
	float vdist,hdist;
	int sIndex = 0, tIndex = 0, p =0,lcnt = 0, rcnt =0;
	int stpointcnt = 0; int pointcnt = 0;
	medH = Median(ccHeights);
	medW = Median(ccWidths);
	float vbandthresh1 = 0.35*medH;
	float hbandthresh1 = 2.5*medW;
	thetamat.resize(size);
	quadmat.resize(size);
	qcount.resize(4);
	frndmat = Allocate2DMatrixInt(size,size);
	for (sIndex = 0 ; sIndex < size; sIndex++)
	{
		stquad.clear(); tantheta.clear(); quad.clear(); 
		stquad.resize(size); tantheta.resize(size);quad.resize(size);
		//srcRVal = centroidX[sIndex]; // select the value for which we have to find the distane
		//srcCVal = centroidY[sIndex]; 
		hbandthresh = hbandthresh1;
		vbandthresh = vbandthresh1;
		TL = 0; BL=0; TR=0; BR=0; // reset
		GetCountInBands(sIndex); // Now BL,TL etc should have been updated
		ResizeBand(); // update hbandthres and vbandthresh if required
		ObtainTransformedQuad(sIndex);
		pointcnt = SumVec(quad); 
		if (pointcnt == 0)
		{
			continue;
		}
		stpointcnt = SumVec(stquad); //cout<<"\nstpointcnt : "<<stpointcnt;
		//fdistmat[sIndex][sIndex] = 0;
		if (stpointcnt >  0.65*(pointcnt+ 0.0))
		{
			StraightOrientation(sIndex,vbandthresh);
			continue;
		}
		GetQuadCounts(); 
		int maxall = FindMaxValue(qcount);	
		thetamat[sIndex] = 0;
		//STDDEV for multiple quadnum
		int quadnum = GetMaxIndx(qcount,maxall) + 1;
		quadmat[sIndex] = quadnum;
		double avgtheta = ComputeAvgTheta(quadnum,sIndex); // average of the angle in radians
		double avgttheta = tan(avgtheta); // tan theta of average angle
		thetamat[sIndex] = avgttheta;
		for(tIndex = 0 ; tIndex < size; tIndex++)
		{	
			if(sIndex != tIndex)
			{
				if(quad[tIndex]  == quadnum)
				{
					frndmat[sIndex][tIndex] = frndmat[tIndex][sIndex] = 1;
				}
			}
		}
		
	}
	//All-pair shortest path
	ComputeFloydAPSP(frndmat, size);
	//Breadth-First-Search for components of graph
	FindComponents(frndmat, size);//Now CClabels has textline information
	//for(int i=0;i<numValidCC;i++)
	//	cout<<"\t"<<CClabels[0][i];
	maxCCLabel = GetMaxCCLabel(CClabels,numValidCC);//maxlabel = line-1 (starts with 0)
	CombineLines(); //Group broken lines if they are aligned and close
	
	return CClabels;
}

int TextLine::GetMaxCCLabel(int ** labels,int n)
{
	int maxcclab = 0;
	for(int k = 0; k < n;k++)
	{
		//cout<<"\t"<<labels[0][k];
		if (maxcclab < labels[0][k])
			maxcclab = labels[0][k];
	}
	return maxcclab; 
}

int _tmain(int argc, _TCHAR* argv[])
{
        IplImage *img2 = cvLoadImage("D:\\experiments\\expr1-images\\series-13\\IMG_20110615_231952.jpg");
        /* always check */
		if( img2 == 0 ) {
			fprintf( stderr, "Cannot load file %s!\n", "image" );
			return 1;
		}
		IplImage *img = cvCreateImage( cvSize( (img2->width)/2, (img2->height)/2 ),img2->depth, img2->nChannels);
		
		//resize image
		cvResize( img2, img);
		/*image properties*/
		int width     = img->width;
		int height    = img->height;
		int depth     = img->depth;
		int nchannels = img->nChannels;
		
		IplImage *grayimg = cvCreateImage( cvSize( width, height ), IPL_DEPTH_8U, 1 );
		
		/* CV_RGB2GRAY: convert RGB image to grayscale */
		cvCvtColor( img, grayimg, CV_RGB2GRAY );
		//IplImage *smoothimg = cvCreateImage( cvSize( grayimg->width, grayimg->height), grayimg->depth, grayimg->nChannels);
		//cvSmooth( grayimg, smoothimg, CV_GAUSSIAN, 3, 3 );
		//cvNamedWindow( "Source", 1 );
        //cvShowImage( "Source", smoothimg );
		//cvWaitKey();
		//grayimg = smoothimg;
		TextLine tLine;
		//Connected Components using Cvcontours
		CvMemStorage *mem; 
        CvSeq *contours, *ptr,*qtr;
		IplImage *cc_color = cvCreateImage(cvGetSize(grayimg), IPL_DEPTH_8U, 3);
		//cvThreshold(grayimg, grayimg, 150, 255, CV_THRESH_BINARY); 
		cvAdaptiveThreshold(grayimg, grayimg,double(255),CV_ADAPTIVE_THRESH_GAUSSIAN_C,1, 105,  5);		//parameter to be subtracted
		
		cvNamedWindow( "Source", 1 );
        cvShowImage( "Source", grayimg );
		//cvWaitKey();
		mem = cvCreateMemStorage(0); 
		cvFindContours(grayimg, mem, &contours, sizeof(CvContour), CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE);
		CvSeq* validCC = tLine.FilterCC(contours, 50, 50,8,8); // Filter text ccs
		CvRect r;
		int numValidCC = tLine.GetNumValidCCs();
		//extract CC properties
		tLine.GetPropertyofCC(validCC);
		int **ccLabels = tLine.FindLines();
		int label = 0,i=0;
		int numlines = tLine.GetNumLines();
		//for (int i=0;i< numValidCC;i++)
		//	cout<<"\t"<<ccLabels[0][i];
		for (label = 0; label < numlines; label++) { 
			
			CvScalar color = CV_RGB( rand()&255, rand()&255, rand()&255 ); 
			qtr = validCC;
			i=0;
			while(qtr != NULL)
			{
				//cout<<"\t"<<ccLabels[0][i];
				if(ccLabels[0][i] == label)
				{
					r = cvBoundingRect(qtr, 1);
					//cvRectangle(cc_color,cvPoint(r.x,r.y),cvPoint(r.x+r.width,r.y+r.height),color);
					cvDrawContours(cc_color, qtr, color, CV_RGB(0,0,0), -1, CV_FILLED, 8, cvPoint(0,0));
				}
				i++;
				qtr = qtr->h_next;
			}
		}

		cvNamedWindow( "Components", 1 );
        cvShowImage( "Components", cc_color );
		//cvSaveImage("D:\\experiments\\expr1-images\\Andreas2\\IMG_20110707_160933_lines.jpg",cc_color);
        cvWaitKey();
        cvDestroyWindow("Image:");
        cvReleaseImage(&img);
		cvDestroyWindow("GrayImage:");
        cvReleaseImage(&grayimg);
        return 0;
}
