/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
*/

/*--------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <numeric>

#include "functionclass.h"
#include "filehelper.h"

/*--------------------------------------------------------*/

using namespace std;

/*--------------------------------------------------------*/
bool FunctionClass::sortByData(const funcPoints &lhs, const funcPoints &rhs)
{ 
	if (lhs.data != rhs.data)
		return lhs.data < rhs.data;
	else
		return lhs.size < rhs.size;
}
			
/*--------------------------------------------------------*/

FunctionClass::FunctionClass()
{

}

/*--------------------------------------------------------*/		

FunctionClass::~FunctionClass()
{
	delete[] func_sorted;
}

/*--------------------------------------------------------*/

FunctionClass::FunctionClass(char *fileName, char *dName, unsigned int *gran)
{
	currentElm = -1;		// No selected element.
	deviceName = dName;
	
	FileOperations fObj (fileName);
	numPoints_sorted = fObj.numPointsCnt();
	func_sorted = new funcPoints[numPoints_sorted];
	if (func_sorted == NULL)
	{
		cerr << "ERROR: Function for " << deviceName << " cannot be built" << endl;
		exit(1);
	}
	
	fObj.goToBegin();
	
	if (*gran == 0)
	{
		double data;
		unsigned int size;
		unsigned int temp_g = 0;
		
		for (int i = 0; i < numPoints_sorted; i++)
		{
			if (!fObj.readNextLine (&size, &data))
			{
				cerr << "There is a problem with reading the function of " << deviceName << endl;
				exit(1);
			}
			
			temp_g = __gcd (size, temp_g);
		}
		
		*gran = temp_g;
	}
	
	granularity = *gran;
	
	fObj.goToBegin();
	
	unsigned int maxSize = 0;
	
	for (int i = 0; i < numPoints_sorted; i++)
	{
		if (!fObj.readNextLine (&func_sorted[i].size, &func_sorted[i].data))
		{
			cerr << "There is a problem with reading the function of " << deviceName << endl;
			exit(1);
		}
		
		if (func_sorted[i].size % granularity != 0)
		{
			cerr << "Size is not dividable by granularity" << endl;
			exit(1);
		}
		
		if (func_sorted[i].size > maxSize)
		{
			maxSize = func_sorted[i].size;
		}
	}
	
	numPoints_size = maxSize / granularity + 1;
	
	func_size = new double[numPoints_size];
	if (func_size == NULL)
	{
		cerr << "ERROR: Function for " << deviceName << " cannot be built" << endl;
		exit(1);
	}
	
	func_size[0] = 0;
	for (int i = 1; i < numPoints_size; i++)
	{
		func_size[i] = -1;
	}
	
	for (int i = 0; i < numPoints_sorted; i++)
	{
		func_size[func_sorted[i].size / granularity] = func_sorted[i].data;
	}	
	
	sort(func_sorted, func_sorted + numPoints_sorted, sortByData);
}	

/*--------------------------------------------------------*/

FunctionClass::FunctionClass(const unsigned int npoints, const unsigned int *sizes, const double *values, 
																char *dName, unsigned int *gran)
{
	currentElm = -1;		// No selected element.
	deviceName = dName;
	//strcpy (deviceName, dName);
	
	numPoints_sorted = npoints;
	func_sorted = new funcPoints[numPoints_sorted];
	if (func_sorted == NULL)
	{
		cerr << "ERROR: Function for " << deviceName << " cannot be built" << endl;
		exit(1);
	}
	
	if (*gran == 0)
	{
		unsigned int temp_g = 0;
		
		for (int i = 0; i < numPoints_sorted; i++)
		{			
			temp_g = __gcd (sizes[i], temp_g);
		}
		
		*gran = temp_g;
	}
	
	granularity = *gran;
		
	unsigned int maxSize = 0;
	
	for (int i = 0; i < numPoints_sorted; i++)
	{
		func_sorted[i].size = sizes[i];
		func_sorted[i].data = values[i];
		
		if (func_sorted[i].size % granularity != 0)
		{
			cerr << "Size is not dividable by granularity" << endl;
			exit(1);
		}
		
		if (func_sorted[i].size > maxSize)
		{
			maxSize = func_sorted[i].size;
		}
	}
	
	numPoints_size = maxSize / granularity + 1;
	
	func_size = new double[numPoints_size];
	if (func_size == NULL)
	{
		cerr << "ERROR: Function for " << deviceName << " cannot be built" << endl;
		exit(1);
	}
	
	func_size[0] = 0;
	for (int i = 1; i < numPoints_size; i++)
	{
		func_size[i] = -1;
	}
	
	for (int i = 0; i < numPoints_sorted; i++)
	{
		func_size[func_sorted[i].size / granularity] = func_sorted[i].data;
	}	
	
	sort(func_sorted, func_sorted + numPoints_sorted, sortByData);
}	

/*--------------------------------------------------------*/

FunctionClass::FunctionClass(char *fileName, char *dName, unsigned int *gran, char type)
{
	currentElm = -1;		// No selected element.
	deviceName = dName;
	
	FileOperations fObj (fileName);
	numPoints_sorted = fObj.numPointsCnt();
	func_sorted = new funcPoints[numPoints_sorted];
	if (func_sorted == NULL)
	{
		cerr << "ERROR: Function for " << deviceName << " cannot be built" << endl;
		exit(1);
	}
	
	fObj.goToBegin();
	
	double data_org, data_smooth;
	
	if (*gran == 0)
	{
		double data;
		unsigned int size;
		unsigned int temp_g = 0;
		
		for (int i = 0; i < numPoints_sorted; i++)
		{
			if (!fObj.readNextLine (&size, &data))
			{
				cerr << "There is a problem with reading the function of " << deviceName << endl;
				exit(1);
			}
			
			temp_g = __gcd (size, temp_g);
		}
		
		*gran = temp_g;
	}
	
	granularity = *gran;
	
	fObj.goToBegin();
	
	unsigned int maxSize = 0;
	
	if (type == 'o') // creating original fucntion
	{
		for (int i = 0; i < numPoints_sorted; i++)
		{
			if (!fObj.readNextLine (&func_sorted[i].size, &func_sorted[i].data, &data_smooth))
			{
				cerr << "There is a problem with reading the function of " << deviceName << endl;
			}
			
			if (func_sorted[i].size % granularity != 0)
			{
				cerr << "Size is not dividable by granularity" << endl;
				exit(1);
			}
			
			if (func_sorted[i].size > maxSize)
			{
				maxSize = func_sorted[i].size;
			}
		}
		
		numPoints_size = maxSize / granularity + 1;
	
		func_size = new double[numPoints_size];
		if (func_size == NULL)
		{
			cerr << "ERROR: Function for " << deviceName << " cannot be built" << endl;
			exit(1);
		}
		
		func_size[0] = 0;
		for (int i = 1; i < numPoints_size; i++)
		{
			func_size[i] = -1;
		}
		
		for (int i = 0; i < numPoints_sorted; i++)
		{
			func_size[func_sorted[i].size / granularity] = func_sorted[i].data;
		}
	}
	else	// creating smooth fucntion
	{
		for (int i = 0; i < numPoints_sorted; i++)
		{
			if (!fObj.readNextLine (&func_sorted[i].size, &data_org, &func_sorted[i].data))
			{
				cerr << "There is a problem with reading the function of " << deviceName << endl;
			}
			
			if (func_sorted[i].size % granularity != 0)
			{
				cerr << "Size is not dividable by granularity" << endl;
				exit(1);
			}
			
			if (func_sorted[i].size > maxSize)
			{
				maxSize = func_sorted[i].size;
			}
		}
		
		numPoints_size = maxSize / granularity + 1;
		func_size = new double[numPoints_size];
		if (func_size == NULL)
		{
			cerr << "ERROR: Function for " << deviceName << " cannot be built" << endl;
			exit(1);
		}
		
		func_size[0] = 0;
		for (int i = 1; i < numPoints_size; i++)
		{
			func_size[i] = -1;
		}
		
		for (int i = 0; i < numPoints_sorted; i++)
		{
			func_size[func_sorted[i].size / granularity] = func_sorted[i].data;
		}
	}
	
	sort(func_sorted, func_sorted + numPoints_sorted, sortByData);
}

/*--------------------------------------------------------*/

bool FunctionClass::binSearch (double key, int start, int end, double delta, int *index)
{	
	// Not found
	if (start > end)
	{
		*index = (end < 0 ? 0 : end);
		return false;
	}
	
	int mid = (start + end) / 2;
	
	//cout << getData (mid) << " " << key << " " << (getData (mid) - key) / denominator << " " << delta << endl;
	
	// Found
	if ((abs(getData (mid) - key) / key) <= delta)
	{
		*index = mid;
		return true;
	}
	
	if (getData (mid) > key)
		return binSearch (key, start, mid -1, delta, index);
	else
		return binSearch (key, mid + 1, end, delta, index);	
}

/*--------------------------------------------------------*/

unsigned int FunctionClass::funcSize ()
{
	return numPoints_sorted;
}

/*--------------------------------------------------------*/

void FunctionClass::printFunction ()
{
	cout << "Function of " << deviceName << endl;
	for (int i = 0; i < numPoints_sorted; i++)
	{
		cout << func_sorted[i].size << " " << func_sorted[i].data << endl;
	}
}

/*--------------------------------------------------------*/

char* FunctionClass::getDeviceName ()
{
	return deviceName;
}

/*--------------------------------------------------------*/

void FunctionClass::setCurrentElm (int index)
{
	currentElm = index;
}

/*--------------------------------------------------------*/

int FunctionClass::getCurrentElm ()
{
	return currentElm;
}

/*--------------------------------------------------------*/

void FunctionClass::incCurrentElm ()
{
	currentElm++;
}

/*--------------------------------------------------------*/

double FunctionClass::getCurrentData ()
{
	if (currentElm == -1)
		return 0;
		
	return func_sorted[currentElm].data;
}

/*--------------------------------------------------------*/

unsigned int FunctionClass::getCurrentSize ()
{
	if (currentElm == -1)
		return 0;
	
	return func_sorted[currentElm].size;
}

/*--------------------------------------------------------*/
		
double FunctionClass::getNextData ()
{				
	if (currentElm >= numPoints_sorted)
		return 0;
	return func_sorted[currentElm + 1].data;
}

/*--------------------------------------------------------*/

unsigned int FunctionClass::getNextSize ()
{
	if (currentElm >= numPoints_sorted)
		return 0;
	return func_sorted[currentElm + 1].size;
}

/*--------------------------------------------------------*/

double FunctionClass::getData (int index)
{				
	if (index == -1)
		return 0;
	
	return func_sorted[index].data;
}

/*--------------------------------------------------------*/

unsigned int FunctionClass::getSize (int index)
{
	if (index == -1)
		return 0;
		
	return func_sorted[index].size;
}

/*--------------------------------------------------------*/

double FunctionClass::getDataBySize (unsigned int size)
{
	if (size % granularity != 0 || size > (numPoints_size - 1) * granularity)
		return -1;
	
	return func_size[size / granularity];
}

/*--------------------------------------------------------*/

unsigned int FunctionClass::getMaxSize ()
{
	return (numPoints_size - 1) * granularity;
}

/*--------------------------------------------------------*/

unsigned int FunctionClass::getMinSize ()
{	
	int i;	
	for (i = 0; i < numPoints_size && (func_size[i] == 0 || func_size[i] == -1); i++);
	
	return i * granularity;
}

/*--------------------------------------------------------*/

void FunctionClass::findMinMaxSizeByData (double data, unsigned long *minSize, unsigned long *maxSize)
{
	*maxSize = 0;
	*minSize = numeric_limits<unsigned int>::max();
	
	for (int i = 0; i < numPoints_sorted && func_sorted[i].data <= data; i++)
	{
		if (*maxSize < func_sorted[i].size)
			*maxSize = func_sorted[i].size;
			
		if (*minSize > func_sorted[i].size)
			*minSize = func_sorted[i].size;
	}
}

/*--------------------------------------------------------*/

void FunctionClass::findMinMaxSizeByData (double data1, double data2, unsigned long *minSize, unsigned long *maxSize)
{
	if (data1 > data2)
	{
		cerr << " ERROR: data1 is supposed not to be larger than data2" << endl;
		
		return;
	}
	
	int i = 0;
	*maxSize = 0;
	*minSize = numeric_limits<unsigned long>::max();
	
	for (i = 0; i < numPoints_sorted && func_sorted[i].data < data1; i++);	// Go to the start point
	
	for ( ; i < numPoints_sorted && func_sorted[i].data <= data2; i++)
	{
		if (*maxSize < func_sorted[i].size)
			*maxSize = func_sorted[i].size;
			
		if (*minSize > func_sorted[i].size)
			*minSize = func_sorted[i].size;
	}
}

/*--------------------------------------------------------*/

unsigned int FunctionClass::findMaxSizeByData (double data)
{
	unsigned int maxSize = 0;
	
	for (int i = 0; i < numPoints_sorted && func_sorted[i].data < data; i++)
	{
		if (maxSize < func_sorted[i].size)
			maxSize = func_sorted[i].size;
	}
	
	return maxSize;
}

/*--------------------------------------------------------*/

unsigned int FunctionClass::findMaxSizeToIndex (int index)
{
	unsigned int maxSize = 0;
	
	for (int i = 0; i < numPoints_sorted && i <= index; i++)
	{
		if (maxSize < func_sorted[i].size)
			maxSize = func_sorted[i].size;
	}
	
	return maxSize;
}

/*--------------------------------------------------------*/
