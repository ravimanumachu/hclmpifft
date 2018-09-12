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

#include <algorithm>
#include <vector>
#include <sys/time.h>
#include <time.h>
#include <cstring>

#include "functionclass.h"
#include "wdistmethods.h"
#include "dirhelper.h"

/*--------------------------------------------------------*/

unsigned int FUNC_GRANURALITY;	// It is the granurality of size in time functions.

/*--------------------------------------------------------*/

void ShowInParams (char *arg)
{
	cerr << endl << "Usage: " << arg << " <work_size> <dist_method> <versosity> <path to the time functions>.";
 	cerr << endl;
 	cerr << "\tdist_method can be: equal, heuristic, het or hetBasic\n" << endl;
 	cerr << endl; 	
}

/*--------------------------------------------------------*/

int main(int argc, char** argv)
{
	unsigned int numFuncs = 0;
	unsigned int workSize;
	char *partMethod;
	FunctionClass **tFunction;	// An array points to time function objects.[numFuncs];
	char **timeFilePathArr;
	char timeDirPath[200];
	int verbosity;
	
	DirOperations *timeDirOpt = new DirOperations ();
	
	if (argc != 5)
	{
		ShowInParams(argv[0]);
		exit(1);
	}

	workSize = atoi(argv[1]);
	partMethod = argv[2];
	verbosity = atoi(argv[3]);
	strcpy (timeDirPath, argv[4]);
  
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
  	
  tFunction = new FunctionClass*[numFuncs];
  
  // Creating time functions
  char tempPath[300];
  for (unsigned int i = 0; i < numFuncs; i++)
  {
  	strcpy (tempPath, timeDirPath);
  	strcat (tempPath, "/");
  	strcat (tempPath, timeFilePathArr[i]);
  	tFunction[i] = new FunctionClass(tempPath, timeFilePathArr[i], &FUNC_GRANURALITY);
  }
  
  if (workSize % FUNC_GRANURALITY != 0)
  {
  	cout << "The worksize should be dividable with " << FUNC_GRANURALITY << endl;
  	exit(1);
  }
  
  Partitioning *pObj = new Partitioning(numFuncs);
  double executionTime = -1;	// The execution time of the obtained distribution.
  
  unsigned int parts_cur[numFuncs];
  unsigned int parts_opt[numFuncs];
	unsigned long maxSizes[numFuncs];
  
  unsigned int sum = 0;
  
  for (unsigned int i = 0; i < numFuncs; i++)
  {
  	sum += tFunction[i]->getMaxSize();
  }
  
  if (sum < workSize)
  {
  	cerr << "ERROR: Workload exceeds the largest problem size can be solved by current fucntions" << endl;
  	exit (1);
  }
  
  time_t t;
	srand((unsigned)time(&t));
  
  struct timeval start, end;
  
  gettimeofday(&start, NULL);

  if (string(partMethod) == string("equal"))
  {
  	pObj->equalWD (tFunction, numFuncs, workSize, &executionTime, parts_opt, FUNC_GRANURALITY);
  	
  	if (executionTime == -1)
  		cout << "Not found any solution" << endl << endl;
  	else
  		pObj->printSolutionInOptSolution (tFunction, numFuncs, workSize, executionTime, parts_opt, string("Equal"), verbosity);
	}
	else if (string(partMethod) == string("heuristic"))
	{
		pObj->heuristicWD (tFunction, numFuncs, workSize, &executionTime, parts_opt);
		
		if (executionTime == -1)
  		cout << "Not found any solution" << endl << endl;
		else
			pObj->printSolutionInOptSolution (tFunction, numFuncs, workSize, executionTime, parts_opt, string("Heuristic"), verbosity);
	}
	else if (string(partMethod) == string("het"))
	{
		// Start point from load equal point
		cout << endl<< "Finding Time Threshold ..." << endl;
		executionTime = pObj->timeThreshold (tFunction, numFuncs, workSize, parts_opt, FUNC_GRANURALITY);
		
		cout << "Time Threshold= " << executionTime << endl;
		
		bkToStruct backTrack;
		backTrack.bkToIndex = numeric_limits<unsigned int>::max();
		
		// Initilize memory array
		// No record is stored for processors with index = 0 and index = numFuncs - 1
		// At the sake of simplicity, first element in array memories is created for processor 0, but no memory is assigned to it.
		memorization **memories = new memorization*[numFuncs - 1];
		
		for (unsigned int i = 1; i < numFuncs - 1; i++)
		{
			memories[i] = new memorization (workSize / FUNC_GRANURALITY + 1, numFuncs - i, FUNC_GRANURALITY);
		}
		
		pObj->hetWD (tFunction, numFuncs, 0, workSize, parts_cur, &executionTime,
																parts_opt, maxSizes, &backTrack, FUNC_GRANURALITY, memories);
																									
		pObj->printSolutionInOptSolution (tFunction, numFuncs, workSize, executionTime, parts_opt, string("Heterogeneous"), verbosity);
		
		for (unsigned int i = 1; i < numFuncs - 1; i++)
		{
			delete memories[i];
		}
		delete[] memories;
	}
	else if (string(partMethod) == string("hetBasic"))
	{		
		cout << endl<< "Finding The Execution Time for Load-Equal Distribution ..." << endl;
		
		//pObj->equalWD (tFunction, numFuncs, workSize, &executionTime, parts_opt, FUNC_GRANURALITY);
		
		executionTime = pObj->timeThreshold (tFunction, numFuncs, workSize, parts_opt, FUNC_GRANURALITY);
		
		if (executionTime == -1)
		{
			cout << "Not found any solution for Load-Equal Distribution" << endl << endl;
			exit(1);
		}
		
		cout << "Load-Equal Execution Time = " << executionTime << endl;
		unsigned int parts_cur[numFuncs];
		pObj->hetBasicWD (tFunction, numFuncs, 0, workSize, parts_cur, &executionTime, parts_opt);
		
		pObj->printSolutionInOptSolution (tFunction, numFuncs, workSize, executionTime, parts_opt, 
																				string("HeteroBasic"), verbosity);
	}
	else
	{
		cout << "Please select correct value for dist_method" << endl;
		exit(1);
	}
	
	gettimeofday(&end, NULL);
	double tstart = start.tv_sec + start.tv_usec/1000000.;
  double tend = end.tv_sec + end.tv_usec/1000000.;
  double t_sec = (tend - tstart);
  
  // Deconstruction
  for (unsigned int i = 0; i < numFuncs; i++)
  {
  	delete tFunction[i];
  }
  delete[] tFunction;
  
  cout << endl << "Total Execution Time (sec) For Finding The Solution = " << t_sec << endl << endl;
  
  return 0;
}
