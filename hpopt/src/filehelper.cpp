/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
*/

/*--------------------------------------------------------*/

#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <sys/types.h>
#include <stdlib.h>

#include "filehelper.h"

/*--------------------------------------------------------*/

using namespace std;

/*--------------------------------------------------------*/

FileOperations::FileOperations (char *fn)
{
	fName = fn;

	ifs.open (fName, ifstream::in);
	if (!ifs.is_open())
	{
		cerr << "ERROR: \"" << fName << "\" cannot be openned" << endl;
		exit(1);
	}
}
		
/*--------------------------------------------------------*/
		
FileOperations::~FileOperations()
{
	ifs.close(); 
}

/*--------------------------------------------------------*/
		
unsigned int FileOperations::numPointsCnt()
{
	string line;
	unsigned int number_of_lines = 0;

  goToBegin();
  
  while(!ifs.eof())
  {
		getline(ifs, line);
		number_of_lines++;
  }
  
  return number_of_lines - 1;
}

/*--------------------------------------------------------*/
		
void FileOperations::goToBegin()
{
	ifs.clear();
	ifs.seekg (0, ios::beg);
}

/*--------------------------------------------------------*/
		
bool FileOperations::readNextLine (unsigned int *s, double *d)
{
	string line;
	
	if (getline(ifs, line))
	{
		istringstream iss(line);			
		if (!(iss >> *s >> *d)) 
		{
			cerr << "ERROR: There is a problem at " << fName << endl;
			exit(1);
		}
		return true;
	}
	else
		return false;
}

/*--------------------------------------------------------*/
		
bool FileOperations::readNextLine (unsigned int *s, double *d1, double *d2)
{
	string line;
	
	if (getline(ifs, line))
	{
		istringstream iss(line);			
		if (!(iss >> *s >> *d1 >> *d2)) 
		{
			cerr << "ERROR: There is a problem at " << fName << endl;
			exit(1);
		}
		return true;
	}
	else
		return false;
}
/*--------------------------------------------------------*/

void FileOperations::removeFileExtension (char *fullName, char *fileName)
{
	unsigned int i = 0;
	unsigned int dotIndex = 0;
	while (fullName[i] != '\0')
	{
		if (fullName[i] == '.')
		{
			dotIndex = i;
		}
		
		i++;
	}	
	
	if (dotIndex == 0)	// fullName is extensionless.
	{
		strcpy (fileName, fullName);
	}
	else
	{
		strncpy (fileName, fullName, dotIndex);
		fileName[dotIndex] = '\0';
	}
}

/*--------------------------------------------------------*/
