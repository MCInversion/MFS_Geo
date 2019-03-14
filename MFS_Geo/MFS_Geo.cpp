// MFS_Geo.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"

int main()
{
	printf("Loading data ... \n");

	std::fstream dataFile;
	dataFile.open("BL-902.dat", std::fstream::in);

	if (!dataFile.is_open()) {
		printf("Unable to open file BL-902.dat\n");
		return 0;
	}
	else {
		printf("BL-902.dat opened successfully\n");


		dataFile.close();
		return 1;
	}
}

