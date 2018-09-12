/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
*/

/*--------------------------------------------------------*/

#ifndef _FILE_HELPER_H
#define _FILE_HELPER_H

/*--------------------------------------------------------*/

#include <iostream>
#include <fstream>

/*--------------------------------------------------------*/

using namespace std;

/*--------------------------------------------------------*/

class FileOperations
{
	private:
		char *fName;	// path to file
		ifstream ifs;	
		
	public:
	
		/*--------------------------------------------------------*/
		
		FileOperations (char *fn);
		
		~FileOperations();
		
		unsigned int numPointsCnt();
		
		void goToBegin();
		
		bool readNextLine (unsigned int *s, double *d );
		bool readNextLine (unsigned int *s, double *d1, double *d2 );
		
		static void removeFileExtension (char *fullName, char *fileName);
};

/*--------------------------------------------------------*/
#endif	// _FILE_HELPER_H
