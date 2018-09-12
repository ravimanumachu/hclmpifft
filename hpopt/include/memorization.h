/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
*/

/*--------------------------------------------------------*/

#ifndef _MEMORIZATION_H
#define _MEMORIZATION_H

/*--------------------------------------------------------*/

#include "functionclass.h"

/*--------------------------------------------------------*/
#define	_NOT_EVAL			-1				// parameter equals _NOT_EVAL means the node has not been evalulated before.
#define _FINALIZED		-2				// parameter equals _FINALIZED means the node has optimal distribution.

/*--------------------------------------------------------*/

typedef enum {SOLUTION, NOT_SOLUTION, SOLUTION_RESUME, RESUME, DUMMY} MEMREADSTATUS;

/*--------------------------------------------------------*/

// This sturcture stores distribution found on a processor.
struct partitionDetail
{
	double parameter;		// parameter is maximum exucution time of the found distribution in HPOPT
											// parameter is total consumed energy of the found distribution in HEOPT.
											
	int lastEvalIndex;	// The index of last evaluated element. _NOT_EVAL means no element has not been evaluated.
											// _FINALIZED means the node has final result (optimal result) and 
											// all possible answers has been evaluated.
											
	unsigned int part;		
	
	/*--------------------------------------------------------*/
	
	partitionDetail ()
	{
		parameter = _NOT_EVAL;
		lastEvalIndex = _NOT_EVAL;
	}
	
	/*--------------------------------------------------------*/
};

class memorization
{
	private:
		struct partitionDetail *memArr;
		unsigned long numCells;						// Size of memArr.
		unsigned int numDev;							// number of devices.
		unsigned int granularity;
		
		
	public:
		/*--------------------------------------------------------*/
		
		memorization (unsigned long nCells, unsigned int numD, unsigned int gran);
		
		~memorization ();
		
		void makeFinal (unsigned int size);
		
		void storeLastIndex (unsigned int wSize, int value);
		
		int retrieveLastIndex (unsigned int wSize);
		
		// Storig an element at cell with index of wSize. number of elements in array points is 
		void storeElm (unsigned int wSize, unsigned int points, double newParam);
		
		void retrieveElm (unsigned int wSize, double *parameter, int *lastEvalIndex, unsigned int *point);
};

/*--------------------------------------------------------*/

#endif	// _MEMORIZATION_H
