/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
*/

/*--------------------------------------------------------*/

#ifndef _WDIST_METHODS_H
#define _WDIST_METHODS_H

/*--------------------------------------------------------*/

#include <vector>

#include "functionclass.h"
#include "memorization.h"

using namespace std;

/*--------------------------------------------------------*/

struct bkToStruct
{
	unsigned int bkToIndex;
	double bkToTime;
};

/*--------------------------------------------------------*/

class Partitioning
{
	private:		
		/*--------------------------------------------------------*/
		
		unsigned long numCall;	// count the number of times hetWD is called.
		
		unsigned int findMinTime (FunctionClass **tFuncs, unsigned int numDev);
		
		double findMaxTime (FunctionClass **tFuncs, unsigned int numDev);
		
		struct node
		{
			unsigned int index;
			double value;
		};
		
		// Sort descending
		static bool sortDesc (const node &lhs, const node &rhs);
		
		// Sort descending
		static bool sortAsc (const node &lhs, const node &rhs);

		// This function sort time functions based on the speed so that slower functions will be follwed by
		// faster ones. We need this function for full recursive approach.		
		void sortTimeFuncsWithSpeed (FunctionClass **tFuncs, unsigned int numDev, double worstTime);
		
		// This function is called by hetWD when a solution is found.
		// The input parameter foundIndex determines the index of processor where the solution was found.
		// The input parameter fromMemory is true if we used memory for finding this current solution.
		// If optTime and maxTime of the solution are the same, solution will be the one with lower processors.
		void findSolutionOperation (FunctionClass **tFuncs, unsigned int *parts_cur, double *optTime, unsigned int *parts_opt,
																	memorization **memories, unsigned int totalDev, int foundIndex,
																	unsigned long *maxSizeArr, bkToStruct *backTrackTo);
																
		// This function retrieves solution form memory
		double retrieveSolution (FunctionClass **tFuncs, unsigned int *parts_cur, memorization **memories, unsigned int totalDev, 
															unsigned int cIndex, int wSize);
															
		// This function is called to read the memory and retrieve solution and set the last position in time function to resume 
		// the process. when it outputs true, the caller should retrun.
		MEMREADSTATUS readFromMemory (FunctionClass **tFuncs, unsigned int totalDev, unsigned int curDevIndex,
												unsigned int workSize, unsigned int *parts_cur, double *optTime, memorization **memories);
												
	public:
		
		Partitioning (unsigned int numDev);
		
		~Partitioning ();
		
		void resetNumCall ();
		
		int printSolutionInOptSolution (FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize, 
																			double optTime, unsigned int *parts_opt, string pMethod, int verbosity);
																			
		// This fucntion is similar to printSolutionInOptSolution however it replaces times with those in tFuncs
		void printSolutionInOptSolution_replace (FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize, 
																			double optTime, unsigned int *parts_opt, string pMethod, int verbosity);
		
		// This function scans time functions from right data points to the beginning to find the minimum time threshold. 
		double refineTimeTh (FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize,double timeTh);
		
		double timeThreshold (FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize, 
														unsigned int *parts_opt, int g);
		
		// This function finds equal workload distribution
		// In this method processors do not finish at the same time
		void equalWD (FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize, 
													double *optTime, unsigned int *parts_opt, int g);
		
		// This function finds the best load-balanced workload partitioning
		// In this method all processors finish at the same time
		// This method search all devices to finds points with similar execution times and store them.
		// The solution invloves all processors.
		void balanceWD (FunctionClass **tFuncs, unsigned int totalDev, unsigned int curDevIndex,
															unsigned int totalWorkSize, unsigned int *parts_cur, double *optTime, 
															unsigned int *parts_opt, unsigned long *minSizeArr, 
															unsigned long *maxSizeArr, double delta, bool *backToRoot);
																// maxSizeArr[i] is the accumulation of maxsize[i] to maxSize[totalDev-1]
																// minSizeArr[]i is similar to maxSize Arr.
																
		//This function find an approperiate workload distribution using a
		//recursive method. It selects the point with minimum exection time
		// each step until the summation of all selected points is less than 
		//workload size. Then, it ignore the function with maximum time and 
		//repeat the process for the remaining functions and remaining workload size.
		void heuristicWD(FunctionClass **tFuncs, unsigned int numDev, unsigned int tWorkSize, double *optTime,
														unsigned int *parts_opt);
		
		// This function consider all cominations to find an optimal solution.
		// parts_cur is an array consisting of totalDev elements which saves sizes and times for 
		// each solution from root to one leaf. This array will alanysis in each leaf to find 
		// the workload distribution for that case.
		// optimizations :
		// 1- execution time of each selected point is compared with the execution time of current founded 
		// solution. If it is higher, we should return instead of continuing of work.
		// 2- The maximum possible size is obtained and check if one selection is likely to result in a 
		// solution of not. Sizes can be found at input array maxSize. At first it should be filled by
		// maximum sizes lay between times 1 to equal time.
		// 3- When an solution is found which has lower time than current solution, we should not branch children of 
		// the with maximum time at the solution. So we should bactrack to the ancestor of the processor with maximum time.
		// To this end, we defined input parameter backTrackTo to prevent undue recursive calls. This variable determines
		// until which processor backtracking should be continued. Since there may be more than one
		// processor with maximum time in at the solutin, backTrackTo contains the processor index which
		// is closer to the bactracking tree. If it is -1, there is no backtracking.
		// 4- Storing intermediate results to resue again.
		// 5- Store last examined size for each processor for reuse.
		// This method finds the optimal solution.
		void hetWD(FunctionClass **tFuncs, unsigned int totalDev, unsigned int curDevIndex,
															unsigned int workSize, unsigned int *parts_cur, double *optTime, unsigned int *parts_opt,
															unsigned long *maxSizeArr, bkToStruct *backTrackTo,
															unsigned int granularity, memorization **memories);
															// memories is an array with size of totalNumDev - 1
															// maxSizeArr[i] is the accumulation of maxsize[i] to maxSize[totalDev-1]		
																											
		// This method consider all combinations without any optimization for performance improvement
		void hetBasicWD(FunctionClass **tFuncs, unsigned int totalDev, unsigned int curDevIndex,
						 											unsigned int workSize, unsigned int *parts_cur, 
							 										double *optTime, unsigned int *parts_opt);
};

/*--------------------------------------------------------*/

#endif	// _WDIST_METHODS_H
