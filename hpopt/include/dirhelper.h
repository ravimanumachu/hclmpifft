/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
*/

/*--------------------------------------------------------*/

#ifndef _DIRECTORY_HELPER_H
#define _DIRECTORY_HELPER_H

/*--------------------------------------------------------*/

class DirOperations
{
	public:
		
		/*--------------------------------------------------------*/
		
		// This function returns the number of all files existing in the path dPath
		int countFiles (char *dPath);
		
		// this function returns the names of all files in the path dPath		
		void fileNames (char *dPath, char **fNames);
};

#endif	// _DIRECTORY_HELPER_H
