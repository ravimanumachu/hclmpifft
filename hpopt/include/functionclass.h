/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
*/

/*--------------------------------------------------------*/

#ifndef _FUNCTION_CLASS_H
#define _FUNCTION_CLASS_H

/*--------------------------------------------------------*/

#include <iostream>

#include "filehelper.h"

/*--------------------------------------------------------*/

using namespace std;

/*--------------------------------------------------------*/

struct funcPoints		
{
	double data;	// data is execution time in HPOPT
								// data is consumed energy in HEOPT
	unsigned int size;
	
	funcPoints ()
	{
		data = 0;
		size = 0;
	}	
};

/*--------------------------------------------------------*/

class FunctionClass
{
	private:
		int numPoints_sorted;			// number of points in func_sorted
		int numPoints_size;				// number of points in func_sorted
		char *deviceName;
		funcPoints *func_sorted; // an array contains the (time/energy) function of a device sorted by time/energy.
		
		double *func_size;			// an array contain the (time/energy) function sorted by workload size so that
														// func_size[w / granularity] is time/energy for workload w.
														// if there is no point for workload w in the function, func_size[w / granularity]
														// is set to 0.
		
		unsigned int granularity;		// minimum size granularity in functions
		
		int currentElm;		// Point to the current element in the function. 
											// Negative value for currentElm means that there is no selected element.
		
		// The function will be used by qsort () to compare two elements and sort element ascendingly.
		// If two records have the same values of data, they will be sorted according size.
		static bool sortByData(const funcPoints &lhs, const funcPoints &rhs);		
		
	public:
					
		/*--------------------------------------------------------*/
		
		FunctionClass();
		
		~FunctionClass();
		
		// This constructor builts a function from file
		FunctionClass(char *fileName, char *dName, unsigned int *gran);
		
		// This constructor builts a function from array
		FunctionClass(const unsigned int npoints, const unsigned int *sizes, const double *values, char *dname, unsigned int *gran);
		
		// This constructor is used when we are going to compare two workload distibutions.
		FunctionClass(char *fileName, char *dName, unsigned int *gran, char type);
		
		//This function find a value which has same value of its difference with key is delta at most.
		// if key is not found, index is the position it is expected that the key was.
		bool binSearch (double key, int start, int end, double delta, int *index);
		
		unsigned int funcSize ();
		
		void printFunction ();
		
		char* getDeviceName ();
		
		void setCurrentElm (int index);
		
		int getCurrentElm ();
		
		void incCurrentElm ();
		
		double getCurrentData ();
		
		unsigned int getCurrentSize ();
		
		double getNextData ();
		
		unsigned int getNextSize ();
		
		double getData (int index);
		
		unsigned int getSize (int index);
		
		double getDataBySize (unsigned int size);
		
		unsigned int getMaxSize ();
		
		unsigned int getMinSize ();
		
		// This function returns both min and max size. It begins search from the first element in the function. 
		void findMinMaxSizeByData (double data, unsigned long *minSize, unsigned long *maxSize);
		
		// This function returns both min and max. It consider elements between data1 and data2
		void findMinMaxSizeByData (double data1, double data2, unsigned long *minSize, unsigned long *maxSize);
		
		// This function returns max size. It begins search from the first element in speed function and consider all
		// data points with data value less than data. 
		unsigned int findMaxSizeByData (double data);
		
		// This function returns max size. It consider all data points with index in [0 .. index].
		unsigned int findMaxSizeToIndex (int index);
};

/*--------------------------------------------------------*/
#endif	// _FUNCTION_CLASS_H

