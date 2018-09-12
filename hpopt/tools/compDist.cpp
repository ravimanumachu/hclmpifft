/*--------------------------------------------------------*/

/*Copyright 2001-2018 Heterogeneous Computing Laboratory, 
School of Computer Science, University College Dublin

Redistribution and use in source and binary forms, 
with or without modification, are permitted provided 
that the following conditions are met:

1. Redistributions of source code must retain the 
   above copyright notice, this list of conditions 
   and the following disclaimer.

2. Redistributions in binary form must reproduce the 
   above copyright notice, this list of conditions and 
   the following disclaimer in the documentation and/or 
   other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names 
   of its contributors may be used to endorse or promote 
   products derived from this software without specific 
   prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.*/

/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
*/

/*--------------------------------------------------------*/

/*
compDist is a tool which compares heterogeneous workload distribution
with constant and semi functional performance model
*/

/*--------------------------------------------------------*/

#include <algorithm>
#include <vector>
#include <sys/time.h>
#include <time.h>
#include <cstring>

#include "functionclass.h"
#include "wdistmethods.h"
#include "dirhelper.h"

/*--------------------------------------------------------*/

void ShowInParams (char *arg)
{
	cerr << endl << "Heterogeneous vs. Semi-Functional and Constant Performance Models" << endl;
	cerr << "===================================================" << endl;
	cerr << "Usage: " << arg << " <work_size> <point> <num_of_nodes> <versosity> <path to the time functions>." << endl;
	cerr << "===================================================" << endl;
	cerr << "The input parameter \"point\" determines a worksize in the fucntions" << endl;
	cerr <<	"where is used for obtaining relative speeds. The relative speeds are used for" << endl;
	cerr << "Constant Performance Model." << endl;
	cerr << "The input parameter \"num_of_clusters\" determines the number of nodes in the simulated cluster" << endl;
	cerr << "The format of records in time functions:" << endl;
	cerr << "<work_size> <original time> <smooth time>" << endl;
 	cerr << endl;
}

/*--------------------------------------------------------*/

int main(int argc, char** argv)
{
	unsigned int FUNC_GRANURALITY = 0;	// It is the granurality of size in time functions.
	unsigned int numFuncs = 0;
	unsigned int samplingPoint = 0;
	unsigned int workSize;
	FunctionClass **tFunction;				// An array points to origianl time function objects.[numFuncs * numNodes];
	FunctionClass **smoothTFunction;	// An array points to smooth time function objects.[numFuncs * numNodes];
	char **timeFilePathArr;
	char timeDirPath[200];
	int verbosity;
	unsigned int numNodes = 1;
	
	DirOperations *timeDirOpt = new DirOperations ();
	
	if (argc != 6)
	{
		ShowInParams(argv[0]);
   	exit(1);
	}
	
	workSize = atoi(argv[1]);
	samplingPoint = atoi(argv[2]);	
	if (atoi(argv[2]) <= 0)
	{
		cerr << "\"point\" should not be zero or negative" << endl;
		exit(1);
	}
	
	numNodes = atoi(argv[3]);
	if (atoi(argv[3]) <= 0)
	{
		cerr << "\"point\" should not be zero or negative" << endl;
		exit(1);
	}
	
	verbosity = atoi(argv[4]);
	strcpy (timeDirPath, argv[5]);
	
	numFuncs = timeDirOpt->countFiles (timeDirPath);  
  
  if (numFuncs == 0)
  {
  	cerr << "ERROR:There is no file in " << timeDirPath << endl;
  	exit (1);
  }
  
  timeFilePathArr = new char*[numFuncs];
  for (unsigned int i = 0; i < numFuncs; i++)
  {
  	timeFilePathArr[i] = new char[100];
  }
  
  timeDirOpt-> fileNames (timeDirPath, timeFilePathArr);
  
  tFunction = new FunctionClass*[numFuncs * numNodes];
	smoothTFunction = new FunctionClass*[numFuncs * numNodes];
	
	// Creating time functions (original and smooth)
  char tempPath[300];
  for (unsigned int i = 0; i < numFuncs; i++)
  {
  	strcpy (tempPath, timeDirPath);
  	strcat (tempPath, "/");
  	strcat (tempPath, timeFilePathArr[i]);
  	for (unsigned int j = 0; j < numNodes; j++)
  	{
			tFunction[(i * numNodes) + j] = new FunctionClass(tempPath, timeFilePathArr[i], &FUNC_GRANURALITY, 'o');
			smoothTFunction[(i * numNodes) + j] = new FunctionClass(tempPath, timeFilePathArr[i], &FUNC_GRANURALITY, 's');
		}
  }
	
	if (workSize % FUNC_GRANURALITY != 0)
  {
  	cout << "The worksize should be dividable with " << FUNC_GRANURALITY << endl;
  	exit(1);
  }
  
  double execTime_equal_org;
  double execTime_het_org;
  unsigned int parts_opt_org[numFuncs * numNodes];
  
  unsigned int parts_cur[numFuncs * numNodes];
  
  Partitioning *pObj = new Partitioning(numFuncs * numNodes);  
  
  unsigned long maxSizes[numFuncs * numNodes];
  
  // Start point from load equal point
	execTime_equal_org = pObj->timeThreshold (tFunction, numFuncs * numNodes, workSize, parts_opt_org, FUNC_GRANURALITY);
	
	cout << "Time Threshold (Original)= " << execTime_equal_org << endl;
	
	bkToStruct backTrack;
	backTrack.bkToIndex = numeric_limits<unsigned int>::max();
	
	// Initilize memory array
	memorization **memories = new memorization*[numFuncs * numNodes - 1];
	
	for (unsigned int i = 1; i < numFuncs * numNodes - 1; i++)
	{
		memories[i] = new memorization (workSize / FUNC_GRANURALITY + 1, numFuncs * numNodes - i, FUNC_GRANURALITY);
	}
	
	execTime_het_org = execTime_equal_org;
	
	pObj->hetWD (tFunction, numFuncs * numNodes, 0, workSize, parts_cur, &execTime_het_org,
															parts_opt_org, maxSizes, &backTrack, FUNC_GRANURALITY, memories);	
																								
	pObj->printSolutionInOptSolution (tFunction, numFuncs * numNodes, workSize, execTime_het_org, parts_opt_org, 
																			string("Heterogeneous_Original"), verbosity);
    
  for (unsigned int i = 1; i < numFuncs * numNodes - 1; i++)
	{
		delete memories[i];
	}
  
  cout << "--------------------------------------------------------" << endl;
  
	double execTime_equal_smooth;
  double execTime_het_smooth;
  unsigned int parts_opt_smooth[numFuncs * numNodes];
	
	pObj->resetNumCall ();
	// Start point from load equal point
	execTime_equal_smooth = pObj->timeThreshold (smoothTFunction, numFuncs * numNodes, workSize, parts_opt_smooth, FUNC_GRANURALITY);

	cout << "Time Threshold (Smooth) = " << execTime_equal_smooth << endl;

	backTrack.bkToIndex = numeric_limits<unsigned int>::max();

	// Initilize memory array
	for (unsigned int i = 1; i < numFuncs * numNodes - 1; i++)
	{
		memories[i] = new memorization (workSize / FUNC_GRANURALITY + 1, numFuncs * numNodes - i, FUNC_GRANURALITY);
	}

	execTime_het_smooth = execTime_equal_smooth;
	
	pObj->hetWD (smoothTFunction, numFuncs * numNodes, 0, workSize, parts_cur, &execTime_het_smooth,
															parts_opt_smooth, maxSizes, &backTrack, FUNC_GRANURALITY, memories);	
																							
	// This function shows the distribution using original times instead of smooth ones.
	pObj->printSolutionInOptSolution_replace (tFunction, numFuncs * numNodes, workSize, execTime_het_smooth, parts_opt_smooth, 
																			string("Heterogeneous_Smooth"), verbosity);

	for (unsigned int i = 1; i < numFuncs * numNodes - 1; i++)
	{
		delete memories[i];
	}
  
  cout << "--------------------------------------------------------" << endl;
  
	double sumRevTimes = 0;
	
  double execTime_cpm;
  unsigned int parts_cpm[numFuncs * numNodes];
	
	pObj->resetNumCall ();
	
	for (unsigned int i = 0; i < numFuncs * numNodes; i++)
	{
		if (tFunction[i]->getDataBySize (samplingPoint) == -1)
		{
			cerr << "Choose a correct samplingPoint" << endl;
			exit (1);
		}
		sumRevTimes += 1 / tFunction[i]->getDataBySize (samplingPoint);  
	}
	
	bool noSolution = false;
	for (unsigned int i = 0; i < numFuncs * numNodes; i++)
	{
		parts_cpm[i] = (1 / (tFunction[i]->getDataBySize (samplingPoint) * sumRevTimes)) * workSize;
		parts_cpm[i] -= (parts_cpm[i] % FUNC_GRANURALITY);
		
		if (parts_cpm[i] != 0 && tFunction[i]->getDataBySize (parts_cpm[i]) == -1)
		{
			noSolution = true;
			break;
		}
		
		if (parts_cpm[i] > tFunction[i]->getMaxSize())
		{
			noSolution = true;
			break;
		}
	}
	
	// Refinement
	while (1)
	{
		if (noSolution == true)
			break; 
		
		unsigned int sum = 0;
		for (unsigned int i = 0; i < numFuncs * numNodes; i++)
		{
			sum += parts_cpm[i];
		}
	
		if (sum == workSize)
			break;
		
		double minTime = numeric_limits<double>::max();
		unsigned int minIndex = numeric_limits<unsigned int>::max();
		double tempTime;
	
		for (unsigned int i = 0; i < numFuncs * numNodes; i++)
		{
			if (tFunction[i]->getDataBySize(parts_cpm[i] + FUNC_GRANURALITY) != -1 &&
						parts_cpm[i] + FUNC_GRANURALITY <= tFunction[i]->getMaxSize ())
			{
				tempTime = tFunction[i]->getDataBySize (parts_cpm[i] + FUNC_GRANURALITY);
				if (tempTime < minTime)
				{
					minTime = tempTime;
					minIndex = i;
				}
			}
		}
	
		if (minIndex != numeric_limits<unsigned int>::max())
		{
			parts_cpm[minIndex] += FUNC_GRANURALITY; 
		}
		else
		{
			noSolution = true;
		}
	}
	
	// finding time optimal
	execTime_cpm = -1;
	for (unsigned int i = 0; i < numFuncs * numNodes && noSolution == false; i++)
	{		
		if (execTime_cpm < tFunction[i]->getDataBySize (parts_cpm[i]))
			execTime_cpm = tFunction[i]->getDataBySize (parts_cpm[i]);
	}
	
	pObj->printSolutionInOptSolution (tFunction, numFuncs * numNodes, workSize, execTime_cpm, parts_cpm, 
																			string("ConstantPM"), verbosity);
  
	delete[] memories;
	
	for (unsigned int i = 0; i < numFuncs * numNodes; i++)
  {
  	delete tFunction[i];
  	delete smoothTFunction[i];
  }
  delete[] tFunction;
  delete[] smoothTFunction;
}
