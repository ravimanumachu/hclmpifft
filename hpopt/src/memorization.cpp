/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
*/

/*--------------------------------------------------------*/

#include <iostream>

#include "memorization.h"

/*--------------------------------------------------------*/

using namespace std;

/*--------------------------------------------------------*/
		
memorization::memorization (unsigned long nCells, unsigned int numD, unsigned int gran)
{
	memArr = new partitionDetail[nCells];
	//cout << sizeof (partitionDetail) << endl;
	numCells = nCells;
	numDev = numD;
	granularity = gran;
}	

/*--------------------------------------------------------*/

memorization::~memorization ()
{			
	delete[] memArr;
}

/*--------------------------------------------------------*/

void memorization::makeFinal (unsigned int size)
{
	memArr[size / granularity].lastEvalIndex = _FINALIZED;
}

/*--------------------------------------------------------*/

void memorization::storeLastIndex (unsigned int wSize, int value)
{
	memArr[wSize / granularity].lastEvalIndex = value;
}

/*--------------------------------------------------------*/

int memorization::retrieveLastIndex (unsigned int wSize)
{
	return memArr[wSize / granularity].lastEvalIndex;
}

/*--------------------------------------------------------*/

void memorization::storeElm (unsigned int wSize, unsigned int points, double newParam)
{			
	unsigned int newSize = wSize / granularity;
	
	if (memArr[newSize].parameter == _NOT_EVAL || 
				(memArr[newSize].parameter != _NOT_EVAL && memArr[newSize].parameter > newParam))
	{
		memArr[newSize].parameter = newParam;
		memArr[newSize].part= points;
		
		//cout << "Store: " << numDev << " " << wSize << " " << memArr[wSize / granularity].execTime << endl;
	}
}

/*--------------------------------------------------------*/

void memorization::retrieveElm (unsigned int wSize, double *parameter, int *lastEvalIndex, unsigned int *point)
{						
	unsigned int newSize = wSize / granularity;
	
	if (memArr[newSize].parameter == _NOT_EVAL)
	{
		*parameter = _NOT_EVAL;
		*lastEvalIndex = memArr[newSize].lastEvalIndex;
		return;
	}
	else
	{
		*parameter = memArr[newSize].parameter;
		*lastEvalIndex = memArr[newSize].lastEvalIndex;
		*point = memArr[newSize].part;
		
		//cout << "Retrive: " << numDev << " " << wSize << " " << memArr[wSize / granularity].execTime << endl;
	}
}

/*--------------------------------------------------------*/
