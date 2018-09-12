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
#include "hcl_hpopta.h"

/*--------------------------------------------------------*/

unsigned int FUNC_GRANURALITY;	// It is the granurality of size in time functions.

/*--------------------------------------------------------*/


extern int hcl_hpopta(const int verbosity, const unsigned int workSize,
    										const unsigned int numFuncs, const unsigned int* npoints,
										    const unsigned int* psizes, const double* etimes, 
										    unsigned* dOpt, double* eOpt)
{
	if (numFuncs == 0)
  {
  	cerr << "ERROR: p should be greater than 0" << endl;
  	exit (1);
  }
	
	FunctionClass **tFunction;	// An array points to time function objects.[numFuncs];
	tFunction = new FunctionClass*[numFuncs];
	
	// Creating time functions
	unsigned int sum_npoints = 0;
  for (unsigned int i = 0; i < numFuncs; i++)
  {
  	char *dName = new char[10];
  	sprintf(dName, "%u", i);
  	strcat (dName, ".proc");
  	tFunction[i] = new FunctionClass(npoints[i], &psizes[sum_npoints] , &etimes[sum_npoints], dName, &FUNC_GRANURALITY);
  	sum_npoints += npoints[i];
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

  // Finding optimal distribution 
  
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
																								
	int result = pObj->printSolutionInOptSolution (tFunction, numFuncs, workSize, executionTime, parts_opt, 
																											string("Heterogeneous"), verbosity);
																											
	if (result == -1)
	{
		return -1;
	}
	
	*eOpt = executionTime;
	for (unsigned int i = 0; i < numFuncs; i++)
	{
		dOpt[i] = parts_opt[i];
	}
	
	for (unsigned int i = 1; i < numFuncs - 1; i++)
	{
		delete memories[i];
	}
	delete[] memories;
	
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
