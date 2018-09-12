/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
*/

/*--------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <cstddef>
#include <cmath>
#include <stdlib.h>
#include <algorithm>

#include "wdistmethods.h"

/*--------------------------------------------------------*/

unsigned int Partitioning::findMinTime (FunctionClass **tFuncs, unsigned int numDev)
{
	double minTime;
	unsigned int minIndex = 0;
	
	minTime = tFuncs[0]->getNextData ();
	
	for (unsigned int i = 1; i < numDev; i++)
	{
		if (tFuncs[i]->getNextData () < minTime)
		{
			minTime = tFuncs[i]->getNextData ();
			minIndex = i;
		}
	}
	
	return minIndex;
}

/*--------------------------------------------------------*/

double Partitioning::findMaxTime (FunctionClass **tFuncs, unsigned int numDev)
{
	double maxTime;
	
	maxTime = tFuncs[0]->getCurrentData ();
	
	for (unsigned int i = 1; i < numDev; i++)
	{
		if (tFuncs[i]->getCurrentData () > maxTime)
		{
			maxTime = tFuncs[i]->getCurrentData ();
		}
	}
	
	return maxTime;
}

/*--------------------------------------------------------*/

bool Partitioning::sortDesc (const node &lhs, const node &rhs)
{ 
	return lhs.value > rhs.value;
}

/*--------------------------------------------------------*/

bool Partitioning::sortAsc (const node &lhs, const node &rhs)
{ 
	return lhs.value < rhs.value;
}

/*--------------------------------------------------------*/

void Partitioning::sortTimeFuncsWithSpeed (FunctionClass **tFuncs, unsigned int numDev, double timeThreshold)
{
	struct node tempArr[numDev];
	unsigned int cnt;
	double ct;
	unsigned int ranks[numDev];
	unsigned int cRank;
	
	for (unsigned int i = 0; i < numDev; i++)
	{
		tFuncs[i]->setCurrentElm (0);
		ranks[i] = 0;
	}
	
	while (1)
	{
		cnt = 0;
		
		for (unsigned int i = 0; i < numDev; i++)
		{
			ct = tFuncs[i]->getCurrentData ();
			
			if (ct <= timeThreshold)
			{
				tempArr[cnt].index = i;
				tempArr[cnt++].value = ct;						
			}
			
			tFuncs[i]->incCurrentElm ();					
		}
		
		if (cnt == 0)
			break;
		
		if (cnt == 1)
		{
			ranks[tempArr[0].index] += 1;
		}
		else
		{
			sort(tempArr, tempArr + cnt, sortDesc);
		
			cRank = 1;
			ct = tempArr[0].value;
			ranks[tempArr[0].index] += 1;
		
			for (unsigned int i = 1; i < cnt; i++)
			{
				if (tempArr[i].value != ct)
					cRank++;
			
				ranks[tempArr[i].index] += cRank;
			
				ct = tempArr[i].value;
			}
		}
	}
	
	for (unsigned int i = 0; i < numDev; i++)
	{
		tempArr[i].index = i;
		tempArr[i].value = ranks[i];
	}
	
	sort(tempArr, tempArr + numDev, sortDesc);
	
	FunctionClass **functions = new FunctionClass*[numDev];

	for (unsigned int i = 0; i < numDev; i++)
	{
		functions[i] = tFuncs[tempArr[i].index];
	}
	
	for (unsigned int i = 0; i < numDev; i++)
	{
		tFuncs[i] = functions[i];
	}
	
	delete[] functions;
}

/*--------------------------------------------------------*/

void Partitioning::findSolutionOperation (FunctionClass **tFuncs, unsigned int *parts_cur, double *optTime, 
																						unsigned int *parts_opt, memorization **memories, unsigned int totalDev, 
																						int foundIndex, unsigned long *maxSizeArr, bkToStruct *backTrackTo)
{
	// find execution time for current solution
	unsigned long sumSize = parts_cur[totalDev - 1];
	double mxTime = tFuncs[totalDev - 1]->getDataBySize (parts_cur[totalDev - 1]);
	int mxIndex = totalDev - 1;
	
	for (int i = totalDev - 2; i >= 1; i--)
	{
		if (tFuncs[i]->getDataBySize (parts_cur[i]) >= mxTime)
		{
			mxTime = tFuncs[i]->getDataBySize (parts_cur[i]);
			mxIndex = i;
		}
		
		sumSize += parts_cur[i];
		
		if (sumSize != 0 && i < foundIndex)
		{
			memories[i]->storeElm (sumSize, parts_cur[i], mxTime);		// Store in memory
		}
	}
	
	if (tFuncs[0]->getDataBySize (parts_cur[0]) >= mxTime)
	{
		mxTime = tFuncs[0]->getDataBySize (parts_cur[0]);
		mxIndex = 0;
	}
	
	backTrackTo->bkToIndex = mxIndex;
	backTrackTo->bkToTime = mxTime;
	
	// Local variable mxTime has here the maximum time of the solution.
	// Save the latest optimal result
	if (*optTime > mxTime)
	{
		*optTime = mxTime;
		
		for (unsigned int i = 0; i < totalDev; i++)
		{
			parts_opt[i] = parts_cur[i];
		}
		
		// Update maxSizes array
		maxSizeArr[totalDev - 1] = tFuncs[totalDev - 1]->findMaxSizeByData (*optTime);
		for (int i = totalDev - 2; i >= 0; i--)
		{
			maxSizeArr[i] = maxSizeArr[i + 1] + tFuncs[i]->findMaxSizeByData (*optTime);
		}
	}
	else if (*optTime == mxTime)		// Processor optimality
	{
		unsigned int cur_procs = 0;
		unsigned int sol_procs = 0;
		for (unsigned int i = 0; i < totalDev; i++)
		{
			if (parts_opt[i] != 0)
				sol_procs++;
				
			if (parts_cur[i] != 0)
				cur_procs++;
		}
		
		if (cur_procs < sol_procs)
		{
			for (unsigned int i = 0; i < totalDev; i++)
			{
				parts_opt[i] = parts_cur[i];
			}
		}
	}
}

/*--------------------------------------------------------*/

double Partitioning::retrieveSolution (FunctionClass **tFuncs, unsigned int *parts_cur, memorization **memories, 
																				unsigned int totalDev, unsigned int cIndex, int wSize)
{
	unsigned int i;
	double tempT;
	int tempL;
	unsigned int point;
	double mTime = 0;
	
	for (i = cIndex; i < totalDev - 1 && wSize > 0; i++)
	{
		memories[i]->retrieveElm (wSize, &tempT, &tempL, &point);	// retrieve data from memory.
		parts_cur[i] = point;
		wSize -= point;
		
		if (tFuncs[i]->getDataBySize (point) > mTime)
		{
			mTime = tFuncs[i]->getDataBySize (point);
		}
	}
	
	if (wSize > 0)
	{
		parts_cur[totalDev - 1] = wSize;
		
		if (tFuncs[totalDev - 1]->getDataBySize (wSize) > mTime)
		{
			mTime = tFuncs[totalDev - 1]->getDataBySize (wSize);
		}
	}
	else
	{
		for ( ; i < totalDev; i++)
		{
			parts_cur[i] = 0;
		}
	}
	
	return mTime;
}

/*--------------------------------------------------------*/

MEMREADSTATUS Partitioning::readFromMemory (FunctionClass **tFuncs, unsigned int totalDev, unsigned int curDevIndex,
												unsigned int workSize, unsigned int *parts_cur, double *optTime, memorization **memories)
{
	double tempTime; 
	int last;  
	unsigned int point;

	memories[curDevIndex]->retrieveElm (workSize, &tempTime, &last, &point);	// retrieve data from memory.
	
	if (last == _FINALIZED)		// There is no solution
	{
		if (tempTime == _NOT_EVAL)
		{
			return NOT_SOLUTION;
		}
		else	// tempTime != _NOT_EVAL
		{
			// Retrieve finalized solution
			if (tempTime < *optTime)
			{
				parts_cur[curDevIndex] = point;
				retrieveSolution (tFuncs, parts_cur, memories, totalDev, curDevIndex + 1, workSize - point);					
				return SOLUTION;				
			}
				
			return NOT_SOLUTION;	//backtrack
		}
	}
	else if(last != _FINALIZED)
	{
		// Retrieve last stored solution if it is different from the data point determined with varialbe last 
		if (tempTime != _NOT_EVAL && tempTime < *optTime && tFuncs[curDevIndex]->getSize(last) != point)
		{
			parts_cur[curDevIndex] = point;
			retrieveSolution (tFuncs, parts_cur, memories, totalDev, curDevIndex + 1, workSize - point);					
			
			tFuncs[curDevIndex]->setCurrentElm (last);
			
			return SOLUTION_RESUME;			
		}	
		// resume recursive calls from where it was cut.													// Optimization 5
		if (last != _NOT_EVAL)
		{
			tFuncs[curDevIndex]->setCurrentElm (last);
			return RESUME;
		}
	}
	
	return DUMMY;
}
	
/*--------------------------------------------------------*/

Partitioning::Partitioning (unsigned int numDev)
{	
	numCall = 0;
}

/*--------------------------------------------------------*/

Partitioning::~Partitioning ()
{

}

/*--------------------------------------------------------*/

void Partitioning::resetNumCall ()
{
	numCall = 0;
}

/*--------------------------------------------------------*/

int Partitioning::printSolutionInOptSolution (FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize, 
																								double optTime, unsigned int *parts_opt, string pMethod, int verbosity)
{
	unsigned int numUsedProc = 0;
	unsigned long sumSize = 0;
	if (optTime == -1)
	{
		cout << endl << "There is no solution!" << endl;
		return -1;
	}
	
	if (verbosity)
	{
		cout << endl << setw(25) << "Device" << setw(12) << "Size" << setw(12) << "Time" << endl;
	}
	
	for (unsigned int i = 0; i < numDev; i++)
	{
		if (verbosity)
		{
			cout << setw(25) << tFuncs[i]->getDeviceName ();
			cout << setw(12) << parts_opt[i] << setw(12) << tFuncs[i]->getDataBySize (parts_opt[i]) << endl;
		}
		
		sumSize += parts_opt[i];
		
		if (parts_opt[i] != 0)
		{
			numUsedProc++;
		}	
	}
	
	if (sumSize != tWorkSize)
	{
		cerr << "Wrong Solution!" << endl;
		return -1;
	}
	else
	{	
		cout << endl << pMethod << "-> ExecTime= " << optTime << " WSize= " << tWorkSize;
		cout << " TotDevice= " << numDev << " UsedProc= "<< numUsedProc << endl;
		cout << "numCall = " << numCall << endl;
		return 0;
	}
}

/*--------------------------------------------------------*/

void Partitioning::printSolutionInOptSolution_replace (FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize, 
																								double optTime, unsigned int *parts_opt, string pMethod, int verbosity)

{
	unsigned int numUsedProc = 0;
	unsigned long sumSize = 0;
	double maxTime = 0;
	if (optTime == -1)
	{
		cout << endl << "There is no solution!" << endl;
		return;
	}
	
	if (verbosity)
	{
		cout << endl << setw(25) << "Device" << setw(12) << "Size" << setw(12) << "Time" << endl;
	}
	
	for (unsigned int i = 0; i < numDev; i++)
	{		
		if (verbosity)
		{
			cout << setw(25) << tFuncs[i]->getDeviceName ();
			
			if (parts_opt[i] != 0)
			{
				cout << setw(12) << parts_opt[i] << setw(12) << tFuncs[i]->getDataBySize (parts_opt[i]) << endl;
			}
			else
			{
				cout << setw(12) << 0 << setw(12) << 0 << endl;
			}			
		}
		
		sumSize += parts_opt[i];
		
		if (parts_opt[i] != 0)
		{
			numUsedProc++;
		}
		
		if (parts_opt[i] != 0 && maxTime < tFuncs[i]->getDataBySize (parts_opt[i]))
		{
			maxTime = tFuncs[i]->getDataBySize (parts_opt[i]);
		}
	}
	
	if (sumSize != tWorkSize)
	{
		cerr << "Wrong Solution!" << endl;
		exit(1);
	}
	
	cout << endl << pMethod << "-> ExecTime= " << maxTime << " WSize= " << tWorkSize;
	cout << " TotDevice= " << numDev << " UsedProc= "<< numUsedProc << endl;
	cout << "numCall = " << numCall << endl;
}

/*--------------------------------------------------------*/

double Partitioning::refineTimeTh (FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize,double timeTh)
{
	unsigned long maxArray[numDev];
	unsigned long sumSizes = 0;
	unsigned int fIndex = 0;
	double maxTime = 0;
	
	for (unsigned int i = 0; i < numDev; i++)		// initializing maxArray and current index in all time functions.
	{
		int index;
		tFuncs[i]->binSearch (timeTh, 0, tFuncs[i]->funcSize() - 1, 0, &index);
		
		tFuncs[i]->setCurrentElm (index);
		while (tFuncs[i]->getCurrentElm() < (int)tFuncs[i]->funcSize() - 1 &&
						tFuncs[i]->getNextData () <= timeTh)
		{
			tFuncs[i]->incCurrentElm ();
		}
		
		maxArray[i] = tFuncs[i]->findMaxSizeToIndex (tFuncs[i]->getCurrentElm());
		sumSizes += maxArray[i];
	}
	
	while (sumSizes >= tWorkSize)
	{	
		maxTime = 0;
		
		for (unsigned int i = 0; i < numDev; i++)	// Finding the time function with the greatest current time
		{
			if (tFuncs[i]->getCurrentData () > maxTime)
			{
				maxTime = tFuncs[i]->getCurrentData ();
				fIndex = i;
			}
		}
		
		if (maxTime == 0)
			break;
		
		if (tFuncs[fIndex]->getCurrentSize () == maxArray[fIndex])
		{
			unsigned long temp = maxArray[fIndex];
			maxArray[fIndex] = tFuncs[fIndex]->findMaxSizeToIndex (tFuncs[fIndex]->getCurrentElm() - 1);
			
			sumSizes = sumSizes - temp + maxArray[fIndex];
			
			if (sumSizes >= tWorkSize)
			{
				tFuncs[fIndex]->setCurrentElm (tFuncs[fIndex]->getCurrentElm() - 1);
			}
			else
			{
				break;
			}
		}
		else
		{
			tFuncs[fIndex]->setCurrentElm (tFuncs[fIndex]->getCurrentElm() - 1);
		}
	}
	
	maxTime = 0;
	
	for (unsigned int i = 0; i < numDev; i++)	// Finding the final time threshold
	{
		if (tFuncs[i]->getCurrentData () > maxTime)
		{
			maxTime = tFuncs[i]->getCurrentData ();
		}
	}
	
	return maxTime;
}

/*--------------------------------------------------------*/

double Partitioning::timeThreshold (FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize, 
																			unsigned int *parts_opt, int g)
{
	double tTh;
	equalWD (tFuncs, numDev, tWorkSize, &tTh, parts_opt, g);
	
	if (tTh != -1)
	{
		return tTh;
	}
	else
	{
		for (unsigned int i = 0; i < numDev; i++)
		{
			if (tFuncs[i]->getData(tFuncs[i]->funcSize() - 1) > tTh)
			{
				tTh = tFuncs[i]->getData(tFuncs[i]->funcSize() - 1);
			}
		}
		
		return tTh + 0.00001;	// 0.000001 is arbitrary. 
	}
}

/*--------------------------------------------------------*/

// This equalWD finds load-equal distribution p times faster than the other function. However, the 
// found solution by this fucntion is not the best one.

/*void Partitioning::equalWD (FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize, 
																		double *optTime, funcPoints *parts_opt, int g)
{
	if (tWorkSize % g != 0)
	{
		*optTime = -1;
		return;
	}
	
	unsigned int realWorkSize = tWorkSize;
	if(realWorkSize / (double)numDev < g)
	{
		tWorkSize = g * numDev;
	}
	
	int index = -1;
	unsigned int sizePerDevice;
	unsigned int remSize;
	double maxTime = 0;
	
	sizePerDevice = tWorkSize / numDev;
	sizePerDevice = sizePerDevice / g;
	sizePerDevice = sizePerDevice * g;
	
	remSize = tWorkSize - (sizePerDevice * numDev);
	
	if (remSize % g != 0)
	{
		*optTime = -1;
		return;
	}
	
	for (unsigned int i = 0; i < numDev; i++)
	{
		index = tFuncs[i]-> findSize (sizePerDevice);
		if (index == -1)
		{
			cout << "There is no value for work size = " << sizePerDevice << " for " << tFuncs[i]->getDeviceName () << endl;
			*optTime = -1;
			return;
		}
		else
		{					
			parts_opt[i].size = sizePerDevice;
			parts_opt[i].data = tFuncs[i]->getData (index);
		}
	}
	
	// assigning remaining to first processors
	for (unsigned int i = 0; i < remSize / g; i++)
	{
		parts_opt[i].size += g;
		index = tFuncs[i]->findSize (parts_opt[i].size);
		
		if (index == -1)
		{
			cout << "There is no value for work size = " << sizePerDevice << " for " << tFuncs[i]->getDeviceName () << endl;
			*optTime = -1;
			return;
		}
		else
		{
			parts_opt[i].data = tFuncs[i]->getData (index);
		}
	}
	
	if(realWorkSize / (double)numDev < g)
	{
		for (unsigned int i = 0; i < numDev - (realWorkSize / g); i++)
		{
			maxTime = 0;
			for (unsigned int j = 0; j < numDev; j++)
			{
				if (parts_opt[j].data > maxTime)
				{
					maxTime = parts_opt[j].data;
					index = j;
				}
			}
			parts_opt[index].data = 0;
			parts_opt[index].size = 0;			
		}
	}
	
	for (unsigned int i = 0; i < numDev; i++)
	{
		if (parts_opt[i].data > maxTime)
		{
			maxTime = parts_opt[i].data;
		}
	}
	
	*optTime = maxTime;
}*/

/*--------------------------------------------------------*/

void Partitioning::equalWD (FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize, 
																		double *optTime, unsigned int *parts_opt, int g)
{
	if (tWorkSize % g != 0)
	{
		*optTime = -1;
		return;
	}
	
	unsigned int realWorkSize = tWorkSize;
	if(realWorkSize / (double)numDev < g)
	{
		tWorkSize = g * numDev;
	}
	
	unsigned int sizePerDevice;
	unsigned int remSize;
	double maxTime = 0;
	
	sizePerDevice = tWorkSize / numDev;
	sizePerDevice = sizePerDevice / g;
	sizePerDevice = sizePerDevice * g;
	
	remSize = tWorkSize - (sizePerDevice * numDev);
	
	if (remSize % g != 0)
	{
		*optTime = -1;
		return;
	}
	
	for (unsigned int i = 0; i < numDev; i++)
	{
		double eTime = tFuncs[i]-> getDataBySize (sizePerDevice);
		if (eTime == -1)
		{
			*optTime = -1;
			return;
		}
		else
		{					
			parts_opt[i] = sizePerDevice;
		}
	}
	
	// assigning remaining worksize to processors
	for (unsigned int i = 0; i < remSize / g; i++)
	{
		// if adding probelm size to the processor with max time improve total exec time, go ahead with max.
		maxTime = 0;
		int maxIndex = -1;
		for (unsigned int j = 0; j < numDev; j++)
		{
			double tTime = tFuncs[j]->getDataBySize (parts_opt[j]);
			
			if (tTime > maxTime)
			{
				maxTime = tTime;
				maxIndex = j;
			}
		}
		
		double eTime = tFuncs[maxIndex]-> getDataBySize (parts_opt[maxIndex] + g);
		if (eTime != -1)
		{
			if (eTime < maxTime)
			{
				parts_opt[maxIndex] += g;
				continue;
			}
		}
		
		double minTime = numeric_limits<double>::max();
		int minIndex = -1;
		for (unsigned int j = 0; j < numDev; j++)
		{
			double eTime = tFuncs[j]-> getDataBySize (parts_opt[j] + g);
			if (eTime == -1)
			{
				continue;
			}
			else if (eTime < minTime)
			{
				minTime = eTime;
				minIndex = j;
			}
		}
			
		if (minIndex == -1)
		{
			*optTime = -1;
			return;
		}
		else
		{
			parts_opt[minIndex] += g;
		}		
	}
	
	if(realWorkSize / (double)numDev < g)
	{
		for (unsigned int i = 0; i < numDev - (realWorkSize / g); i++)
		{
			maxTime = 0;
			unsigned index = 0;
			for (unsigned int j = 0; j < numDev; j++)
			{
				double tTime = tFuncs[j]->getDataBySize (parts_opt[j]);
				
				if (tTime > maxTime)
				{
					maxTime = tTime;
					index = j;
				}
			}
			parts_opt[index] = 0;			
		}
	}
	
	maxTime = 0;
	
	for (unsigned int i = 0; i < numDev; i++)
	{
		if (tFuncs[i]->getDataBySize (parts_opt[i]) > maxTime)
		{
			maxTime = tFuncs[i]->getDataBySize (parts_opt[i]);
		}
	}
	
	*optTime = maxTime;
}

/*--------------------------------------------------------*/

void Partitioning::balanceWD (FunctionClass **tFuncs, unsigned int totalDev, unsigned int curDevIndex,
																					unsigned int totalWorkSize, unsigned int *parts_cur, double *optTime, 
																					unsigned int *parts_opt, unsigned long *minSizeArr, 
																					unsigned long *maxSizeArr, double delta, bool *backToRoot)
																					// maxSizeArr[i] is the accumulation of maxsize[i] to maxSize[totalDev-1]
																					// minSizeArr[]i is similar to maxSize Arr.
{			
	numCall++; 
	
	// this is a possible solution
	if (curDevIndex == totalDev)
	{
		unsigned long sizes = 0;
		double maxTime = 0;
		
		for (unsigned int i = 0; i < totalDev; i++)
		{
			sizes += parts_cur[i];
			
			if (tFuncs[i]->getDataBySize (parts_cur[i]) > maxTime)
				maxTime = tFuncs[i]->getDataBySize (parts_cur[i]);
		}
		
		// this is a possible solution
		if (*optTime > maxTime && sizes == totalWorkSize)
		{
			*optTime = maxTime;
			
			for (unsigned int i = 0; i < totalDev; i++)
			{
				parts_opt[i] = parts_cur[i];
			}	
			
			// Update minSizes and maxSizes array
			double tMin = tFuncs[0]->getDataBySize (parts_cur[0]) * (1 - delta);
			tFuncs[totalDev - 1]->findMinMaxSizeByData (tMin, *optTime, &minSizeArr[totalDev - 1], &maxSizeArr[totalDev - 1]);
			for (int i = totalDev - 2; i >= 0; i--)
			{
				tFuncs[i]->findMinMaxSizeByData (tMin, *optTime, &minSizeArr[i], &maxSizeArr[i]);
				minSizeArr[i] += minSizeArr[i + 1];
				maxSizeArr[i] += maxSizeArr[i + 1];
			}
		}
		return;
	}
	
	if (curDevIndex == 0)
	{
		unsigned long maxSizeSum = 0;
		unsigned long minSizeSum = 0;
		*optTime = -1;
		*backToRoot = false;
		double lastTime = 0;
		
		for (unsigned int i = 1; i < totalDev; i++)
		{
			minSizeSum += tFuncs[i]->getMinSize ();
			maxSizeSum += tFuncs[i]->getMaxSize ();
		}
						
		for (long i = 0; i < tFuncs[0]->funcSize() && tFuncs[0]->getData(i) <= *optTime; i++)
		{
			if (i != 0)
				lastTime = tFuncs[0]->getDataBySize (parts_cur[0]);
			
			parts_cur[0] = tFuncs[0]->getSize(i);
			
			if (totalWorkSize > parts_cur[0])
			{
				unsigned int remaining = totalWorkSize - parts_cur[0];
				
				if (remaining > 0 && remaining >= minSizeSum && remaining <= maxSizeSum)
				{
					// tMax , tMin: we are looking given T's located at the interval (abs (T - parts_cur[0].time)/parts_cur[0].time) <= delta
					// If this equation is solved, tMax = parts_cur[0].time * (1 + delta) and tMin = parts_cur[0].time * (1 - d)
					double tMax = tFuncs[0]->getDataBySize (parts_cur[0]) * (1 + delta);
					double tMin = tFuncs[0]->getDataBySize (parts_cur[0]) * (1 - delta);
									
					if (tFuncs[0]->getDataBySize (parts_cur[0]) != lastTime)
					{
						tFuncs[totalDev - 1]->findMinMaxSizeByData (tMin, tMax, &minSizeArr[totalDev - 1], &maxSizeArr[totalDev - 1]);
						for (int i = totalDev - 2; i >= 0; i--)
						{
							tFuncs[i]->findMinMaxSizeByData ((tMin < *optTime ? tMin : *optTime), 
																								(tMax < *optTime ? tMax : *optTime), &minSizeArr[i], &maxSizeArr[i]);
							minSizeArr[i] += minSizeArr[i + 1];
							maxSizeArr[i] += maxSizeArr[i + 1];
						}
					}
				
					if (maxSizeArr[curDevIndex + 1] >= remaining && minSizeArr[curDevIndex + 1] <= remaining 
								&& tFuncs[curDevIndex]->getCurrentData() <= *optTime)		// Optimization 1 , 2
					{
						balanceWD (tFuncs, totalDev, 1, totalWorkSize, parts_cur, optTime, parts_opt,
																	minSizeArr, maxSizeArr, delta, backToRoot);
						*backToRoot = false;
					}
				}
			}	
			else if (totalWorkSize == parts_cur[0] && totalDev == 1)
				balanceWD (tFuncs, totalDev, 1, totalWorkSize, parts_cur, optTime, parts_opt, 
															minSizeArr, maxSizeArr, delta, backToRoot);
		}
		return;
	}	
	
	int index;
	
	// There is no solution with current tFuncs[0]->getDataBySize (parts_cur[0]).
	// Backtrack to root and examine the next point.
	if (tFuncs[curDevIndex]->binSearch (tFuncs[0]->getDataBySize (parts_cur[0]), 0, tFuncs[curDevIndex]->funcSize() - 1, 
																				delta, &index) == false)
	{
		*backToRoot = true;
		return;
	}
	
	unsigned int sumSize = 0;
	for (unsigned int k = 0; k < curDevIndex; k++)
		sumSize += parts_cur[k];
	
	for (int j = index; 
				j >= 0 && 
				(abs(tFuncs[0]->getDataBySize (parts_cur[0]) - tFuncs[curDevIndex]->getData(j))) / tFuncs[0]->getDataBySize (parts_cur[0]) <= delta; j--)
	{
		parts_cur[curDevIndex] = tFuncs[curDevIndex]->getSize(j);
		
		if (totalWorkSize == (sumSize + parts_cur[curDevIndex]) && totalDev == curDevIndex + 1)	// it is a solution
		{
			balanceWD (tFuncs, totalDev, curDevIndex + 1, totalWorkSize, parts_cur, optTime, parts_opt,
														minSizeArr, maxSizeArr, delta, backToRoot);
		}
		else if (totalWorkSize != (sumSize + parts_cur[curDevIndex]) && totalDev == curDevIndex + 1)	// it is a solution
		{
			continue;
		}
		else if (totalWorkSize == (sumSize + parts_cur[curDevIndex]) && totalDev != curDevIndex + 1)
		{
			return;
		}
		else if (totalWorkSize > (sumSize + parts_cur[curDevIndex]))
		{
			unsigned int remaining = totalWorkSize - (sumSize + parts_cur[curDevIndex]);
			
			if (tFuncs[curDevIndex]->getCurrentData() <= *optTime && 
						maxSizeArr[curDevIndex + 1] >= remaining && minSizeArr[curDevIndex + 1] <= remaining) // optimization 1, 2
			{
				balanceWD (tFuncs, totalDev, curDevIndex + 1, totalWorkSize, parts_cur, optTime, parts_opt,
															minSizeArr, maxSizeArr, delta, backToRoot);				
				
				if (*backToRoot == true)
					return;
			}
		}
	}
	
	for (unsigned int j = index + 1; 
				j < tFuncs[curDevIndex]->funcSize() && 
				(abs(tFuncs[0]->getDataBySize (parts_cur[0]) - tFuncs[curDevIndex]->getData(j)) / tFuncs[0]->getDataBySize (parts_cur[0])) <= delta; j++)
	{
		parts_cur[curDevIndex] = tFuncs[curDevIndex]->getSize(j);
			
		if (totalWorkSize == (sumSize + parts_cur[curDevIndex]) && totalDev == curDevIndex + 1)	// it is a solution
		{
			balanceWD (tFuncs, totalDev, curDevIndex + 1, totalWorkSize, parts_cur, optTime, parts_opt, 
														minSizeArr, maxSizeArr, delta, backToRoot);
		}
		else if (totalWorkSize != (sumSize + parts_cur[curDevIndex]) && totalDev == curDevIndex + 1)	// it is a solution
		{
			continue;
		}
		else if (totalWorkSize == (sumSize + parts_cur[curDevIndex]) && totalDev != curDevIndex + 1)
		{
			return;
		}
		else if (totalWorkSize > (sumSize + parts_cur[curDevIndex]))
		{
			unsigned int remaining = totalWorkSize - (sumSize + parts_cur[curDevIndex]);
		
			if (tFuncs[curDevIndex]->getCurrentData() <= *optTime &&
						maxSizeArr[curDevIndex + 1] >= remaining && minSizeArr[curDevIndex + 1] <= remaining) // optimization 1, 2
			{
				balanceWD (tFuncs, totalDev, curDevIndex + 1, totalWorkSize, parts_cur, optTime, parts_opt,
																	minSizeArr, maxSizeArr, delta, backToRoot);
				if (*backToRoot == true)
					return;
			}
		}
	}			
}

/*--------------------------------------------------------*/

void Partitioning::heuristicWD(FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize, double *optTime,
																				unsigned int *parts_opt)
{
	//TODO: this method should be modified to store its solution into parts_opt. Now, the founded solution
	// is extracted from the places pointed by currentElm.
	
	unsigned long minIndex = 0;
	unsigned int currentSize = 0;
	
	numCall++;
	
	if (numDev == 1)
	{
		double eTime = tFuncs[0]->getDataBySize (tWorkSize);
		if (eTime == -1)	// Not Found any solution
		{
			*optTime = -1;
			return;
		}
		else
		{
			parts_opt[0] = tWorkSize;
			*optTime = tFuncs[0]->getDataBySize (parts_opt[0]);
			return;
		}				
	}
	
	for (unsigned int i = 0; i < numDev; i++)
	{
		tFuncs[i]->setCurrentElm(-1);	// Pointing no element in time function. It means that there is not any 
																	// selected size to run on device
	}
	
	minIndex = findMinTime (tFuncs, numDev);
	while (currentSize - tFuncs[minIndex]->getCurrentSize () + tFuncs[minIndex]->getNextSize () < tWorkSize)
	{
		// Go ahead with selection of min time
		currentSize -= tFuncs[minIndex]->getCurrentSize ();
		tFuncs[minIndex]->incCurrentElm ();
		currentSize += tFuncs[minIndex]->getCurrentSize ();
		minIndex = findMinTime (tFuncs, numDev);	
	}	
		
	if (currentSize - tFuncs[minIndex]->getCurrentSize () + tFuncs[minIndex]->getNextSize () == tWorkSize)
	{
		// Found a solution
		tFuncs[minIndex]->incCurrentElm();
		
		for (unsigned int i = 0; i < numDev; i++)
		{
			parts_opt[i] = tFuncs[i]->getCurrentSize();
		}
							
		*optTime = findMaxTime (tFuncs, numDev);
		
		return;
	}
	else
	{
		// Removing tFuncs[minIndex] from active time functions.
		unsigned int newWorkSize = tWorkSize - tFuncs[minIndex]->getCurrentSize ();
		if ((numDev - 1) != minIndex)
		{
			FunctionClass *temp = tFuncs[minIndex];
			tFuncs[minIndex] = tFuncs[numDev - 1];
			tFuncs[numDev - 1] = temp;
		}
		
		heuristicWD(tFuncs, numDev - 1, newWorkSize, optTime, parts_opt);	// Solve the partitioning problem for
																																				// remaining workload size and devices
		if (*optTime == -1)
		{
			return;
		}
		else	// comapre eTime with the removed function
		{
			for (unsigned int i = 0; i < numDev; i++)
			{
				parts_opt[i] = tFuncs[i]->getCurrentSize();
			}
							
			*optTime = findMaxTime (tFuncs, numDev);
		
			return;
		}			
	}
}

/*--------------------------------------------------------*/

void Partitioning::hetWD(FunctionClass **tFuncs, unsigned int totalDev, unsigned int curDevIndex,
																				unsigned int workSize, unsigned int *parts_cur, double *optTime, unsigned int *parts_opt,
																				unsigned long *maxSizeArr, bkToStruct *backTrackTo,
																				unsigned int granularity, memorization **memories)	
																				// memories is an array with size of totalNumDev - 1
																				// maxSizeArr[i] is the accumulation of maxsize[i] to maxSize[totalDev-1]		
{
	numCall++; 
	
	double curTime;
	unsigned int curSize;
	
	// Initialisation
	if (curDevIndex == 0)
	{		
		int tempInd;
		int maxIndex = 0;
		
		for (unsigned int i = 0; i < totalDev; i++)
		{
			tFuncs[i]->binSearch (*optTime, 0, tFuncs[i]->funcSize () - 1, 0, &tempInd);
			
			for ( ; tempInd > 0 && tFuncs[i]->getData (tempInd) == *optTime; tempInd--); 	// the largest index with time 
																																										//less than optTime.
			
			if (tempInd > maxIndex)
			{
				maxIndex = tempInd;
			}
		}
		
		// The Maximum number of points should be examined in a funciton.
		cout << "Max Num Points= " << maxIndex + 1 << endl;
		
		/*if (totalDev > 1)
		{
			cout << "Sorting Functions with Time ..." << endl;
			sortTimeFuncsWithSpeed (tFuncs, totalDev, *optTime); //NOTE: This is not an effective optimization.
		}*/
		
		// initilize maxSizeArr
		maxSizeArr[totalDev - 1] = tFuncs[totalDev - 1]->findMaxSizeByData (*optTime);
		for (int i = totalDev - 2; i >= 0; i--)
		{
			maxSizeArr[i] = maxSizeArr[i + 1] + tFuncs[i]->findMaxSizeByData (*optTime);
		}
	}
	
	if (maxSizeArr[curDevIndex] < (workSize))
		return;
	
	// This is a leaf
	if (curDevIndex + 1 == totalDev)
	{
		double eTime = tFuncs[curDevIndex]->getDataBySize (workSize);
		if (eTime == -1)
		{
			//cerr << "There is no value for work size = " << workSize << " for " << tFuncs[curDevIndex]->getDeviceName () << endl;
			return;
		}
		
		// this is not a solution
		if (eTime >= *optTime)
			return;
			
		// This is a solution
		parts_cur[curDevIndex] = workSize;
			
		findSolutionOperation (tFuncs, parts_cur, optTime, parts_opt, memories, totalDev, totalDev - 1, maxSizeArr, backTrackTo);
		
		if (backTrackTo->bkToIndex >=  curDevIndex)		// optimization 3	
		{
			backTrackTo->bkToIndex = numeric_limits<unsigned int>::max();		// backtrack
		}
		
		return;
	}
	
	curTime = 0;
	curSize = 0;
	tFuncs[curDevIndex]->setCurrentElm(-1);		// curDevIndex is not used.	
	
	// Read from memory
	if (curDevIndex > 0 && curDevIndex <= totalDev - 2)	
	{
		MEMREADSTATUS status;
		
		status = readFromMemory(tFuncs, totalDev, curDevIndex, workSize, parts_cur, optTime, memories);
															
		switch (status)
		{
			case NOT_SOLUTION:		// There is no solution or the solution's exec time is greater than current one.
			{
				return;
			}
			case SOLUTION:
			{
				findSolutionOperation (tFuncs, parts_cur, optTime, parts_opt, memories, totalDev, curDevIndex, 
																maxSizeArr, backTrackTo);
			
				if (backTrackTo->bkToIndex >= curDevIndex)	// Optimization 3
				{
					backTrackTo->bkToIndex = numeric_limits<unsigned int>::max();
				}
				return;
			}
			case SOLUTION_RESUME:
			{
				findSolutionOperation (tFuncs, parts_cur, optTime, parts_opt, memories, totalDev, curDevIndex, 
															maxSizeArr, backTrackTo);
		
				curSize = tFuncs[curDevIndex]->getCurrentSize ();
				curTime = tFuncs[curDevIndex]->getCurrentData ();
				
				if (backTrackTo->bkToIndex > curDevIndex)	// Optimization 3
				{
					backTrackTo->bkToIndex = numeric_limits<unsigned int>::max();
				}
				else if (backTrackTo->bkToIndex == curDevIndex)
				{
					backTrackTo->bkToIndex = numeric_limits<unsigned int>::max();
					memories[curDevIndex]->makeFinal (workSize);
					return;
				}
				else
				{
					return;
				}
				
				break;
			}
			case RESUME:
			{
				curSize = tFuncs[curDevIndex]->getCurrentSize ();
				curTime = tFuncs[curDevIndex]->getCurrentData ();
				break;
			}
			case DUMMY:
			{
				break;
			}
		}
	}

	while (curTime < *optTime)
	{
		if (workSize == curSize || 
				 (workSize > curSize && maxSizeArr[curDevIndex + 1] >= (workSize - curSize)))	// Optimization 2 and 1
		{					
			parts_cur[curDevIndex] = curSize;
			
			if (workSize == curSize)
			{
				for (unsigned int i = curDevIndex + 1; i < totalDev; i++)
				{
					parts_cur[i] = 0;
				}
				
				findSolutionOperation (tFuncs, parts_cur, optTime, parts_opt, memories, totalDev, curDevIndex + 1, 
																maxSizeArr, backTrackTo);
			}
			else
			{
				hetWD (tFuncs, totalDev, curDevIndex + 1, workSize - curSize, parts_cur, 
															optTime, parts_opt, maxSizeArr, backTrackTo, granularity, memories);
			}
		
			if (backTrackTo->bkToIndex <  curDevIndex)		//Optimization 3
			{
				if (curDevIndex > 0 && curDevIndex <= totalDev - 2)
				{
					if (tFuncs[curDevIndex]->getCurrentData() == backTrackTo->bkToTime)
					{
						memories[curDevIndex]->makeFinal (workSize);
					}
					else
					{
						memories[curDevIndex]->storeLastIndex(workSize, tFuncs[curDevIndex]->getCurrentElm());
					}
				}
				
				return;		// backtrack
			}
			else if (backTrackTo->bkToIndex ==  curDevIndex)
			{
				backTrackTo->bkToIndex = numeric_limits<unsigned int>::max();
				if (curDevIndex > 0 && curDevIndex <= totalDev - 2)
				{
					memories[curDevIndex]->makeFinal (workSize);
				}
				
				return;		// backtrack
			}
			else
			{
				backTrackTo->bkToIndex = numeric_limits<unsigned int>::max();
			}
		}
		
		// All data points existing at the function have been examined.
		if (tFuncs[curDevIndex]->getCurrentElm() == (int)tFuncs[curDevIndex]->funcSize() - 1)
			break;
		
		tFuncs[curDevIndex]->incCurrentElm();
		curTime = tFuncs[curDevIndex]->getCurrentData();
		curSize = tFuncs[curDevIndex]->getCurrentSize();
	}
	
	if (curDevIndex > 0 && curDevIndex <= totalDev - 2)
	{
		memories[curDevIndex]->makeFinal (workSize);
	}
}

/*--------------------------------------------------------*/

void Partitioning::hetBasicWD(FunctionClass **tFuncs, unsigned int totalDev, unsigned int curDevIndex,
												 										unsigned int workSize, unsigned int *parts_cur, 
												 										double *optTime, unsigned int *parts_opt)									
{			
	double cTime;
	unsigned int cSize;
	
	numCall++;
	
	if (workSize == 0)
	{
		// This is a solution
		// find execution time for current solution
		double maxTime = 0;
		for (unsigned int i = 0; i < curDevIndex; i++)
		{
			if (tFuncs[i]->getDataBySize (parts_cur[i]) > maxTime)
				maxTime = tFuncs[i]->getDataBySize (parts_cur[i]);
		}
		
		for (unsigned int i = curDevIndex; i < totalDev; i++)
		{
			parts_cur[i] = 0;
		}
		
		// Save the latest optimal result
		if (*optTime > maxTime)
		{
			*optTime = maxTime;
			for (unsigned int i = 0; i < totalDev; i++)
			{
				parts_opt[i] = parts_cur[i];
			}
		}
		return;
	}
	
	// Initialisation
	if (curDevIndex == 0)
	{
		for (unsigned int i = 0; i < totalDev; i++)
		{
			tFuncs[i]->setCurrentElm(-1);	// Pointing to first element in every time function.
		}
	}
	
	// This is a leaf
	if (curDevIndex + 1 == totalDev)
	{
		double eTime = tFuncs[curDevIndex]->getDataBySize (workSize);
		if (eTime == -1)
		{
			//cerr << "There is no value for work size = " << workSize << " for " << tFuncs[curDevIndex]->getDeviceName () << endl;
			return;
		}
		
		// this is not a solution
		if (eTime > *optTime)
			return;
			
		// this is a possible solution
		parts_cur[curDevIndex] = workSize;
			
		
		// find execution time for current solution
		double maxTime = 0;
		for (unsigned int i = 0; i < totalDev; i++)
		{
			if (tFuncs[i]->getDataBySize (parts_cur[i]) > maxTime)
				maxTime = tFuncs[i]->getDataBySize (parts_cur[i]);
		}
		
		// Save the latest optimal result
		if (*optTime > maxTime)
		{
			*optTime = maxTime;
			for (unsigned int i = 0; i < totalDev; i++)
			{
				parts_opt[i] = parts_cur[i];
			}
		}
		return;
	}
	
	cTime = tFuncs[curDevIndex]->getCurrentData();
	cSize = tFuncs[curDevIndex]->getCurrentSize();
	//while (cTime <= *optTime)
	while (cTime <= *optTime)
	{		
		if (workSize >= cSize)
		{
			parts_cur[curDevIndex] = cSize;
						
			hetBasicWD (tFuncs, totalDev, curDevIndex + 1, 
																		workSize - cSize, parts_cur, optTime, parts_opt);
		}
		
		// All data points existing at the function have been examined.
		if (tFuncs[curDevIndex]->getCurrentElm() == (int)tFuncs[curDevIndex]->funcSize() - 1)
			break;
		
		tFuncs[curDevIndex]->incCurrentElm();
		cTime = tFuncs[curDevIndex]->getCurrentData();
		cSize = tFuncs[curDevIndex]->getCurrentSize();
	}
	
	tFuncs[curDevIndex]->setCurrentElm(-1);
}

/*--------------------------------------------------------*/
