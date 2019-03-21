// MFS_Geo.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"

// global constants
double R = 6371.; // Earth radius [km]
double GM = 389600.5; // [km^3 * s^(-1)]
double uExact = GM / R; // exact surface potential
double qExact = GM / (R * R); // surface acceleration

int N = 902; // dataset size;

void printArray1(std::string name, double *a, int printLim) {
	std::string offset = std::string((name + " = ").length() + 1, ' ');
	std::cout << name << " = " << std::endl;
	for (int i = 0; i < printLim; i++) {
		std::cout << offset << a[i] << std::endl;
	}
	for (int i = 0; i < 3; i++) std::cout << offset << "." << std::endl;
	for (int i = N - printLim; i < N; i++) {
		std::cout << offset << a[i] << std::endl;
	}
	std::cout << std::endl;
}

void printArrayVector3(std::string name, double *vx, double *vy, double *vz, int printLim) {
	std::string offset = std::string((name + " = ").length() + 1, ' ');
	std::cout << name << " = " << std::endl;
	for (int i = 0; i < printLim; i++) {
		std::cout << offset << vx[i] << " " << vy[i] << " " << vz[i] << std::endl;
	}
	for (int i = 0; i < 3; i++) std::cout << offset << "." << std::endl;
	for (int i = N - printLim; i < N; i++) {
		std::cout << offset << vx[i] << " " << vy[i] << " " << vz[i] << std::endl;
	}
	std::cout << std::endl;
}

void loadPointData(double *x, double *y, double *z, double *sx, double *sy, double *sz, double *q, double *d2U) {
	printf("Loading data ... \n");

	std::fstream dataFile;
	dataFile.open("BL-902.dat", std::fstream::in);

	if (!dataFile.is_open()) {
		printf("Unable to open file BL-902.dat\n");
	}
	else {
		printf("BL-902.dat opened successfully\n");

		std::string line;
		int i = 0;

		while (std::getline(dataFile, line)) {
			std::vector<std::string> tokens;			
			std::string s_delimiter = " ";
			size_t pos = 0;
			while (pos < 100) {
				pos = line.find(s_delimiter);
				tokens.push_back(line.substr(0, line.find(s_delimiter)));
				line = line.erase(0, pos + s_delimiter.length());
			}

			double B = std::stod(tokens[0]);
			double L = std::stod(tokens[1]);
			double H = std::stod(tokens[2]);
			double Q = std::stod(tokens[3]);
			double D2U = std::stod(tokens[4]);

			// std::cout << B << " " << L << " " << H << " " << Q << " " << D2U << std::endl;

			sx[i] = 1000 * R * cos(B) * sin(L);
			sy[i] = 1000 * R * cos(B) * sin(L);
			sz[i] = 1000 * R * sin(B);

			x[i] = (1000 * R - H) * cos(B) * sin(L);
			y[i] = (1000 * R - H) * cos(B) * sin(L);
			z[i] = (1000 * R - H) * sin(B);

			q[i] = Q;
			d2U[i] = D2U;

			// std::cout << sx[i] << " " << sy[i] << " " << sz[i] << " " << x[i] << " " << y[i] << " " << z[i] << " " << q[i] << " " << d2U[i] << std::endl << std::endl;
			i++;
		}

		dataFile.close();
	}
}

double norm3(double x, double y, double z) {
	return sqrt(x * x + y * y + z * z);
}

double dist3(double ax, double ay, double az, double bx, double by, double bz) {
	return norm3(
		bx - ax,
		by - ay,
		bz - az
	);
}

void getMatrix(double **M, double *x, double *y, double *z, double *sx, double *sy, double *sz) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {

		}
	}
}

int main() {
	double *x = new double[N];
	double *y = new double[N];
	double *z = new double[N];

	double *sx = new double[N];
	double *sy = new double[N];
	double *sz = new double[N];

	double *q = new double[N];
	double *d2U = new double[N];

	loadPointData(x, y, z, sx, sy, sz, q, d2U);

	printArrayVector3("x", x, y, z, 5);
	printArrayVector3("s", sx, sy, sz, 5);

	printArray1("q", q, 5);
	printArray1("d2U/dn2", d2U, 5);

	double **M = new double*[N];
	for (int i = 0; i < N; i++) M[i] = new double[N];
}

