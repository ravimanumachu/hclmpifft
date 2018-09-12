/*--------------------------------------------------------*/

/*Copyright 2001-2018 Heterogeneous Computing Laboratory, 
School of Computer Science, University College Dublin

Redistribution and use in source and binary forms, 
with or without modification, are permitted provided 
that the following conditions are met:

1. Redistributions of source code must retain the 
   above copyright notice, this list of conditions 
   and the following disclaimer.

2. Redistributions in binary form must reproduce the 
   above copyright notice, this list of conditions and 
   the following disclaimer in the documentation and/or 
   other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names 
   of its contributors may be used to endorse or promote 
   products derived from this software without specific 
   prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.*/

/*--------------------------------------------------------*/

/*
@file
@author Hamidreza Khaleghzadeh <hamidreza.khaleghzadeh@ucdconnect.ie>
@version 1.0
@The program generate time functions and correspondent speed functions.
*/

/*--------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>

/*--------------------------------------------------------*/

#define _MAX_PATH_SIZE	200

/*--------------------------------------------------------*/

using namespace std;

/*--------------------------------------------------------*/

double fRand(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

int main(int argc, char** argv)
{
	unsigned int numFuncs;		// number of functions should be created
	unsigned int numPoints;		// number of points in each function
	unsigned int stepSize;
	char *fpath;								// the destination to save created functions
	unsigned int size;
	unsigned int i, j;
	double eTime;
	bool isExist;
	struct stat info;
	
	if (argc < 5)
  {
     cerr << "Usage: " << argv[0] << " <num_func> <num_points> <step_size> <path>.";
     cout << endl;
     exit(1);
  }

  numFuncs = atoi(argv[1]);
  numPoints = atoi(argv[2]);
  stepSize = atoi(argv[3]);
  fpath = argv[4];
  
  string timeFiles[numFuncs];	// path to time function files
  string speedFiles[numFuncs];	// path to speed function files
  
  ofstream otfs[numFuncs];
  ofstream osfs[numFuncs];
  
  // initialize random seed
  srand (time(NULL));
  
  // removing '/' from the end of path if it exists.
	for (i = 0; fpath[i] != '\0'; i++);
	if (fpath[i - 1] == '/')
	{
		fpath[i - 1] = '\0';
	}
	
  // create folders
  
  isExist = false;

  stat( (string (fpath)).c_str(), &info );
  if(info.st_mode & S_IFDIR)
	{
		isExist = true;
	}
  
  if (isExist == false)
  {
		if (system(("mkdir -p " + string (fpath)).c_str()) != 0)
		{
			cerr << "ERROR: Folder\"" << string (fpath) << "\" cannot be created" << endl;
			exit(1);
		}
	}
  
  
  
  isExist = false;

  stat( (string (fpath) + string ("/time/")).c_str(), &info );
  if(info.st_mode & S_IFDIR)
	{
		isExist = true;
	}
  
  if (isExist == false)
  {
		if (system(("mkdir -p " + string (fpath) + string ("/time/")).c_str()) != 0)
		{
			cerr << "ERROR: Folder\"" << string (fpath) + string ("/time/") << "\" cannot be created" << endl;
			exit(1);
		}
	}
	
	isExist = false;

  stat( (string (fpath) + string ("/speed/")).c_str(), &info );
  if(info.st_mode & S_IFDIR)
	{
		isExist = true;
	}
	
	if (isExist == false)
	{
		if (system(("mkdir -p " + string (fpath) + string ("/speed/")).c_str()) != 0)
		{
			cerr << "ERROR: Folder\"" << string (fpath) + string ("/speed/") << "\" cannot be created" << endl;
			exit(1);
		}
	}
  
  for (i = 0; i < numFuncs; i++)
  {	
  	size = stepSize;
  	timeFiles[i] = string (fpath) + string ("/time/func") + to_string(i) + string (".ti");
  	speedFiles[i] = string (fpath) + string ("/speed/func") + to_string(i) + string (".sp");
  	
  	otfs[i].open (timeFiles[i]);
  	osfs[i].open (speedFiles[i]);
		if (!otfs[i].is_open())
		{
			cerr << "ERROR: \"" << timeFiles[i] << "\" cannot be openned" << endl;
			exit(1);
		}
		
		if (!osfs[i].is_open())
		{
			cerr << "ERROR: \"" << speedFiles[i] << "\" cannot be openned" << endl;
			exit(1);
		}
		
		for (j = 0; j < numPoints; j++)
		{
			eTime = fRand (0.0001 + 15 * (j / 100), 0.0001 + 15 * (j / 100) + 30);
			otfs[i] << size << " " << eTime << endl;
			osfs[i] << size << " " << size / eTime << endl;
			
			size += stepSize;
		}
		
		otfs[i].close();
		osfs[i].close();
	}
	
	return 0;
}
