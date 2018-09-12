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
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdlib.h>

#include "dirhelper.h"

/*--------------------------------------------------------*/

using namespace std;

/*--------------------------------------------------------*/

int DirOperations::countFiles (char *dPath)
{
	DIR *dp;
	int cnt = 0;
	struct dirent *ep;     
	dp = opendir (dPath);

	if (dp == NULL)
	{
		cerr << "ERROR: Cannot open the directory in " << dPath << endl;
		exit (1);
	}
	else
	{
		while ((ep = readdir (dp)))
		{
			if (ep->d_type == DT_REG)		// If the entry is a regular file
			{
				cnt++;
			}
		}

		(void) closedir (dp);
	}

	return cnt;
}

/*--------------------------------------------------------*/

void DirOperations::fileNames (char *dPath, char **fNames)
{
	DIR *dp;
	int cnt = 0;
	struct dirent *ep;     
	dp = opendir (dPath);
	
	if (dp == NULL)
	{
		cerr << "ERROR: Cannot open the directory in " << dPath << endl;
		exit (1);
	}
	else
	{
		// removing '/' from the end of path if it exists.
		for (cnt = 0; dPath[cnt] != '\0'; cnt++);
		if (dPath[cnt - 1] == '/')
		{
			dPath[cnt - 1] = '\0';
		}
		
		cnt = 0;
		while ((ep = readdir (dp)))
		{
			if (ep->d_type == DT_REG)		// If the entry is a regular file
			{
				strcpy (fNames[cnt], ep->d_name);
				
				cnt++;
			}
		}

		(void) closedir (dp);
	}
}

/*--------------------------------------------------------*/
