// MFS_Geo.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"

// global constants
double R = 6371.; // Earth radius [km]
double r_remaining = 11; // [km]
double GM = 389600.5; // [km^3 * s^(-1)]
double uExact = GM / R; // exact surface potential
double qExact = GM / (R * R); // surface acceleration

int N = 902; // dataset size;

void printArray1 (std::string name, double *a, int printLim) {
	std::string offset = std::string((name + " = ").length() + 1, ' ');
	std::cout << name << " = " << std::endl;
	for (int i = 0; i < printLim; i++) {
		std::cout << offset << a[i] << std::endl;
	}
	for (int i = 0; i < 3; i++) std::cout << offset << "  ." << std::endl;
	for (int i = N - printLim; i < N; i++) {
		std::cout << offset << a[i] << std::endl;
	}
	std::cout << std::endl;
}

void printArrayVector3 (std::string name, double *vx, double *vy, double *vz, int printLim) {
	std::string offset = std::string((name + " = ").length() + 1, ' ');
	std::cout << name << " = " << std::endl;
	for (int i = 0; i < printLim; i++) {
		std::cout << offset << vx[i] << " " << vy[i] << " " << vz[i] << std::endl;
	}
	for (int i = 0; i < 3; i++) std::cout << offset << "  ." << std::endl;
	for (int i = N - printLim; i < N; i++) {
		std::cout << offset << vx[i] << " " << vy[i] << " " << vz[i] << std::endl;
	}
	std::cout << std::endl;
}

void printArray2 (std::string name, double **A, int printLim) {
	std::string offset = std::string((name + " = ").length() + 1, ' ');
	std::cout << name << " = " << std::endl;
	for (int i = 0; i < printLim; i++) {
		std::cout << offset;
		for (int j = 0; j < printLim; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << "  ...  ";
		for (int j = N - printLim; j < N; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
	for (int i = 0; i < 3; i++) std::cout << offset << "  ." << std::endl;
	for (int i = N - printLim; i < N; i++) {
		std::cout << offset;
		for (int j = 0; j < printLim; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << "  ...  ";
		for (int j = N - printLim; j < N; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void loadPointData (double *x, double *y, double *z, double *sx, double *sy, double *sz, double *q, double *d2U) {
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

			sx[i] = R * cos(B) * cos(L);
			sy[i] = R * cos(B) * sin(L);
			sz[i] = R * sin(B);

			x[i] = (R + 0.001 * H + r_remaining) * cos(B) * cos(L);
			y[i] = (R + 0.001 * H + r_remaining) * cos(B) * sin(L);
			z[i] = (R + 0.001 * H + r_remaining) * sin(B);

			q[i] = Q;
			d2U[i] = D2U;

			// std::cout << B << " " << L << " " << H << " " << Q << " " << D2U << std::endl;
			// std::cout << sx[i] << " " << sy[i] << " " << sz[i] << " " << x[i] << " " << y[i] << " " << z[i] << " " << q[i] << " " << d2U[i] << std::endl << std::endl;
			i++;
		}

		dataFile.close();
	}
}

void getMatrices (double **dG, double **G, double *x, double *y, double *z, double *sx, double *sy, double *sz) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double dx = x[i] - sx[j];
			double dy = y[i] - sy[j];
			double dz = z[i] - sz[j];

			double d_norm = sqrt(dx * dx + dy * dy + dz * dz);

			double nx = sx[i] / R;
			double ny = sy[i] / R;
			double nz = sz[i] / R;

			double dot = dx * nx + dy * ny + dz * nz;

			dG[i][j] = dot / (4 * M_PI * d_norm * d_norm * d_norm); // dG[i][j]/dn[i]
			G[i][j] = 1 / (4 * M_PI * d_norm);
			if (isnan(dG[i][j])) {
				std::cout << "NAN!: dG[" << i << "][" << j << "] : d_norm = " << d_norm << std::endl;
				std::cout << "dx = " << dx << ", dy = " << dy << ", dz = " << dz << std::endl;
			}
			if (isnan(G[i][j])) {
				std::cout << "NAN!: G[" << i << "][" << j << "] : d_norm = " << d_norm << std::endl;
				std::cout << "dx = " << dx << ", dy = " << dy << ", dz = " << dz << std::endl;
			}
		}
	}
}

double vectorDot (double *a, double *b) {
	double result = 0.;
	for (int i = 0; i < N; i++) {
		result += a[i] * b[i];
	}
	return result;
}

double vectorNorm(double *a) {
	return sqrt(vectorDot(a,a));
}

void mmultVector (double **A, double *x, double *result) {
	for (int i = 0; i < N; i++) {
		result[i] = 0.;
		for (int j = 0; j < N; j++) {
			result[i] += A[i][j] * x[j];
		}
	}
}

void smultVector (double s, double *x, double *result) {
	for (int i = 0; i < N; i++) result[i] = x[i] * s;
}

void subVectors (double *a, double *b, double *result) {
	for (int i = 0; i < N; i++) result[i] = a[i] - b[i];
}

void addVectors (double *a, double *b, double *result) {
	for (int i = 0; i < N; i++) result[i] = a[i] + b[i];
}

void copyVector(double *original, double *target) {
	for (int i = 0; i < N; i++) target[i] = original[i];
}

void Bi_CGSTAB_solve (double **A, double *b, double *x) {
	// ctrl. constants
	int maxIter = 100;
	double tol = 1E-6;

	// iter vectors
	double *x_curr = new double[N];
	double *x_next = new double[N];

	double *r_curr = new double[N];
	double *r_next = new double[N];

	double *rp0 = new double[N];

	double *p_curr = new double[N];
	double *p_next = new double[N];

	double *s = new double[N];

	double *tmp = new double[N];
	double *tmp1 = new double[N];

	// iter scalars
	double alpha, beta, omega;

	// x0 = (1,1,...,1)
	for (int i = 0; i < N; i++) x_curr[i] = 1000.;
	// r0 = b - A x0
	// choose rp0 such that <r0, rp0> != 0
	// p0 = r0
	for (int i = 0; i < N; i++) {
		r_curr[i] = b[i];
		for (int j = 0; j < N; j++) {
			r_curr[i] -= A[i][j] * x_curr[j];
		}
		rp0[i] = r_curr[i] + 100;
		p_curr[i] = r_curr[i];
	}
	std::cout << "==================================================" << std::endl;
	std::cout << "----------- Initializing Bi-CGSTAB Method --------" << std::endl;
	printArray2("systemMatrix", A, 4);
	printArray1("systemRhs", b, 5);
	printArray1("x0", x_curr, 2);
	printArray1("r0", r_curr, 5);

	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "------------ Launching iterations ----------------" << std::endl;
	// begin iterations
	for (int k = 0; k < maxIter; k++) {
		std::cout << "::: iter : " << k << std::endl;
		// alpha[k] = <r[k], rp0> / <Ap[k], rp0>
		mmultVector(A, p_curr, tmp);
		alpha = vectorDot(r_curr, rp0) / vectorDot(tmp, rp0);
		// s[k] = r[k] - alpha[k] * A p[k]
		smultVector(alpha, tmp, tmp);
		subVectors(r_curr, tmp, s);

		if (vectorNorm(s) < tol) {
			// x[k + 1] = x[k] + alpha[k] * p[k]
			smultVector(alpha, p_curr, tmp);
			addVectors(x_curr, tmp, x_next);
			break;
		}

		// omega[k] = <A s[k], s[k]> / <A s[k], A s[k]>
		mmultVector(A, s, tmp);
		omega = vectorDot(tmp, s) / vectorDot(tmp, tmp);
		// x[k + 1] = x[k] + alpha[k] * p[k] + omega[k] * s[k]
		smultVector(alpha, p_curr, tmp);
		addVectors(x_curr, tmp, tmp);
		smultVector(omega, s, tmp1);
		addVectors(tmp, tmp1, x_next);
		printArray1("x", x_next, 5);
		// r[k + 1] = s[k] - omega[k] * A s[k]
		mmultVector(A, s, tmp);
		smultVector(omega, tmp, tmp);
		subVectors(s, tmp, r_next);
		if (vectorNorm(r_next) < tol) {
			break;
		}
		// beta[k] = (alpha[k] / omega[k]) * <r[k + 1], rp0> / <r[k], rp0>
		beta = (alpha / omega) * vectorDot(r_next, rp0) / vectorDot(r_curr, rp0);
		// p[k + 1] = r[k + 1] + beta[k] * (p[k] - omega[k] * A p[k])
		mmultVector(A, p_curr, tmp);
		smultVector(omega, tmp, tmp);
		subVectors(p_curr, tmp, tmp);
		smultVector(beta, tmp, tmp);
		addVectors(r_next, tmp, p_next);
		if (fabs(vectorDot(r_next, rp0)) < tol) {
			// rp0 = r[k + 1]; p[k + 1] = r[k + 1]
			copyVector(r_next, rp0);
			copyVector(r_next, p_next);
		}
		// current = next
		copyVector(x_next, x_curr);
		copyVector(r_next, r_curr);
		copyVector(p_next, p_curr);
		std::cout << "===> finishing iter " << k << std::endl;
	}

	copyVector(x_next, x); // result: x = x_next

	// clean up
	delete[] x_curr; delete[] x_next;
	delete[] r_curr; delete[] r_next;
	delete[] p_curr; delete[] p_next;
	delete[] s; delete[] tmp; delete[] tmp1;
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

	double **dG = new double*[N];
	double **G = new double*[N];
	for (int i = 0; i < N; i++) {
		dG[i] = new double[N];
		G[i] = new double[N];
	}

	getMatrices(dG, G, x, y, z, sx, sy, sz);

	printArray2("dG/dn", dG, 4);
	printArray2("G", G, 4);

	double *alphas = new double[N]; // unknown alpha coeffs
	double *u = new double[N]; // potential solution

	// Bi-CGSTAB solve:
	Bi_CGSTAB_solve(dG, q, alphas);

	// print solution
	printArray1("alphas", alphas, 8);

	// potential solution G . alphas = u
	mmultVector(G, alphas, u);

	printArray1("u", u, 8);

	// clean up
	delete[] x; delete[] y; delete[] z;
	delete[] sx; delete[] sy; delete[] sz;
	delete[] q; delete[] d2U;

	for (int i = 0; i < N; i++) {
		delete[] dG[i]; delete[] G[i];
	}
	delete[] dG; delete[] G;

	return 1;
}

