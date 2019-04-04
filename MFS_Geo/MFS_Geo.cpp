#include "pch.h"

// global constants
double R = 6378000.; // Earth radius [m]
double r_remaining = 500000; // [m]
//double GM = 389600.5; // [km^3 * s^(-1)]
//double uExact = GM / R; // exact surface potential
//double qExact = GM / (R * R); // surface acceleration

// const char* dataFilename = "BL-902.dat";
// const char* dataFilename = "BL-1298.dat";
const char* dataFilename = "BL-3602.dat";
// const char* dataFilename = "BL-8102.dat";
int N = 3602; // dataset size;

void printArray1 (std::string name, double *a, int printLim, bool inRow = true) {
	if (inRow) {
		std::cout << name << " = " << std::endl;
		for (int i = 0; i < printLim; i++) {
			std::cout << " " << a[i];
		}
		std::cout << "  ... ";
		for (int i = N - printLim; i < N; i++) {
			std::cout << " " << a[i];
		}
		std::cout << std::endl;
	}
	else {
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

void loadPointData (
	double *B, double *L, double *H,
	double *x, double *y, double *z,
	double *sx, double *sy, double *sz,
	double *q, double *d2U
) {
	printf("Loading data ... \n");

	std::fstream dataFile;
	dataFile.open(dataFilename, std::fstream::in);

	if (!dataFile.is_open()) {
		printf("Unable to open file %s\n", dataFilename);
	}
	else {
		printf("%s opened successfully\n", dataFilename);

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

			B[i] = std::stod(tokens[0]);
			L[i] = std::stod(tokens[1]);
			H[i] = std::stod(tokens[2]);
			double Q = std::stod(tokens[3]);
			double D2U = std::stod(tokens[4]);

			sx[i] = (R - r_remaining) * cos(B[i] * M_PI / 180) * cos(L[i] * M_PI / 180);
			sy[i] = (R - r_remaining) * cos(B[i] * M_PI / 180) * sin(L[i] * M_PI / 180);
			sz[i] = (R - r_remaining) * sin(B[i] * M_PI / 180);

			x[i] = (R + H[i]) * cos(B[i] * M_PI / 180) * cos(L[i] * M_PI / 180);
			y[i] = (R + H[i]) * cos(B[i] * M_PI / 180) * sin(L[i] * M_PI / 180);
			z[i] = (R + H[i]) * sin(B[i] * M_PI / 180);

			q[i] = 0.00001 * Q;
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

			double nx = sx[i] / (R - r_remaining);
			double ny = sy[i] / (R - r_remaining);
			double nz = sz[i] / (R - r_remaining);

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

double vectorNorm (double *a) {
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
	int maxIter = 1000;
	double tol = 1e-5;

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

		std::cout << "||s|| = " << vectorNorm(s) << std::endl;
		if (vectorNorm(s) < tol) {
			// x[k + 1] = x[k] + alpha[k] * p[k]
			smultVector(alpha, p_curr, tmp);
			addVectors(x_curr, tmp, x_next);
			std::cout << "||s|| < tol = " << tol << ", exiting iterations" << std::endl;
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
		// printArray1("x", x_next, 5);
		// r[k + 1] = s[k] - omega[k] * A s[k]
		mmultVector(A, s, tmp);
		smultVector(omega, tmp, tmp);
		subVectors(s, tmp, r_next);
		std::cout << "||r[k + 1]|| = " << vectorNorm(r_next) << std::endl;
		if (vectorNorm(r_next) < tol) {
			std::cout << "||r[k + 1]|| < tol = " << tol << ", exiting iterations" << std::endl;
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
		std::cout << "|< r[k + 1], rp0 >| = " << fabs(vectorDot(r_next, rp0)) << std::endl;
		if (fabs(vectorDot(r_next, rp0)) < tol) {
			// rp0 = r[k + 1]; p[k + 1] = r[k + 1]
			// std::cout << "|< r[k + 1], rp0 >| < tol = " << tol << ", copying r[k + 1] to rp0 and p[k + 1]" << std::endl;
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


void writeData(double *B, double *L, double *u) {
	std::fstream dataOut;
	dataOut.open("data.dat", std::fstream::out);

	for (int i = 0; i < N; i++) {
		dataOut << B[i] << " " << L[i] << " " << u[i] << std::endl;
	}

	dataOut.close();
}

int main(int argc, char **argv) {
	double *B = new double[N];
	double *L = new double[N];
	double *H = new double[N];

	double *x = new double[N];
	double *y = new double[N];
	double *z = new double[N];

	double *sx = new double[N];
	double *sy = new double[N];
	double *sz = new double[N];

	double *q = new double[N];
	double *d2U = new double[N];

	// MPI Vars:
	int nprocs, myrank;
	int istart, iend, nlocal = 0, nlast = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	auto startLoad = std::chrono::high_resolution_clock::now();
	if (myrank == 0) {		
		// loadPointData(B, L, H, x, y, z, sx, sy, sz, q, d2U); :(
		printf("Loading data ... \n");

		std::fstream dataFile;
		dataFile.open(dataFilename, std::fstream::in);

		if (!dataFile.is_open()) {
			printf("Unable to open file %s\n", dataFilename);
			return 0;
		}
		else {
			printf("%s opened successfully\n", dataFilename);

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

				B[i] = std::stod(tokens[0]);
				L[i] = std::stod(tokens[1]);
				H[i] = std::stod(tokens[2]);
				double Q = std::stod(tokens[3]);
				double D2U = std::stod(tokens[4]);

				sx[i] = (R - r_remaining) * cos(B[i] * M_PI / 180) * cos(L[i] * M_PI / 180);
				sy[i] = (R - r_remaining) * cos(B[i] * M_PI / 180) * sin(L[i] * M_PI / 180);
				sz[i] = (R - r_remaining) * sin(B[i] * M_PI / 180);

				x[i] = (R + H[i]) * cos(B[i] * M_PI / 180) * cos(L[i] * M_PI / 180);
				y[i] = (R + H[i]) * cos(B[i] * M_PI / 180) * sin(L[i] * M_PI / 180);
				z[i] = (R + H[i]) * sin(B[i] * M_PI / 180);

				q[i] = 0.00001 * Q;
				d2U[i] = D2U;

				// std::cout << B << " " << L << " " << H << " " << Q << " " << D2U << std::endl;
				// std::cout << sx[i] << " " << sy[i] << " " << sz[i] << " " << x[i] << " " << y[i] << " " << z[i] << " " << q[i] << " " << d2U[i] << std::endl << std::endl;
				i++;
			}

			dataFile.close();
		}
	}
	auto endLoad = std::chrono::high_resolution_clock::now();

	auto startPrint3 = std::chrono::high_resolution_clock::now();
	if (myrank == 0) printArrayVector3("x", x, y, z, 5);
	auto endPrint3 = std::chrono::high_resolution_clock::now();
	if (myrank == 0) printArrayVector3("s", sx, sy, sz, 5);

	auto startPrint1 = std::chrono::high_resolution_clock::now();
	if (myrank == 0) printArray1("q", q, 5);
	auto endPrint1 = std::chrono::high_resolution_clock::now();
	if (myrank == 0) printArray1("d2U/dn2", d2U, 5);

	// array broadcasts
	MPI_Bcast(q, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(z, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(sx, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(sy, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(sz, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// local/last packet sizes
	nlocal = (N / nprocs) + 1;
	nlast = N - (nprocs - 1) * nlocal;

	// packet indexing
	istart = nlocal * myrank;
	if (myrank == nprocs - 1)
		iend = N - 1;
	else
		iend = istart + nlocal - 1;

	// Calculating system matrices from point data:
	auto startMatrixGen = std::chrono::high_resolution_clock::now();

	double **dGLocal = new double*[nlocal]; // system matrix packet
	double *qLocal = new double[nlocal]; // rhs packet
	for (int i = 0; i < nlocal; i++) {
		dGLocal[i] = new double[N];
	}

	int iGlobal = 0;
	for (int i = 0; i < nlocal; i++) {
		iGlobal = i + istart;

		// sys matrix packet filling
		for (int j = 0; j < N; j++) {
			double dx = x[iGlobal] - sx[j];
			double dy = y[iGlobal] - sy[j];
			double dz = z[iGlobal] - sz[j];

			double d_norm = sqrt(dx * dx + dy * dy + dz * dz);

			double nx = sx[iGlobal] / (R - r_remaining);
			double ny = sy[iGlobal] / (R - r_remaining);
			double nz = sz[iGlobal] / (R - r_remaining);

			double dot = dx * nx + dy * ny + dz * nz;

			dGLocal[i][j] = dot / (4 * M_PI * d_norm * d_norm * d_norm);

			if (isnan(dGLocal[i][j]) && myrank == 0) {
				std::cout << "NAN!: dG[" << i << "][" << j << "] : d_norm = " << d_norm << std::endl;
				std::cout << "dx = " << dx << ", dy = " << dy << ", dz = " << dz << std::endl;

				return 0;
			}
		}

		// rhs packet filling:
		qLocal[i] = q[iGlobal];
		//q[i] = qConst;
	}
	// getMatrices(dG, G, x, y, z, sx, sy, sz); // :(
	auto endMatrixGen = std::chrono::high_resolution_clock::now();

	MPI_Allgather(qLocal, nlocal, MPI_DOUBLE, q, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);


	double *alphas = new double[N]; // unknown alpha coeffs
	double *u = new double[N]; // potential solution

	// Bi-CGSTAB solve:
	// NOTE: This has to be done inside the main() function
	// auto startBi_CGSTAB = std::chrono::high_resolution_clock::now();
	// Bi_CGSTAB_solve(dG, q, alphas); // you had a good life, but I can't use you! :'(
	// auto endBi_CGSTAB = std::chrono::high_resolution_clock::now();

	// =========================================================================================
	// ========== Bi-CGSTAB Implementation inside main() =======================================
	auto startBi_CGSTAB = std::chrono::high_resolution_clock::now();

	if (myrank == 0) {
	std::cout << "==================================================" << std::endl;
	std::cout << "----------- Initializing Bi-CGSTAB Method --------" << std::endl;
	}
	// ctrl. constants
	int maxIter = 100;
	double tol = 1e-5;

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

	// iter scalars
	double alpha, beta, omega;

	// x0 = (1000,1000,...,1000)
	for (int i = 0; i < N; i++) x_curr[i] = 1000.;
	// r0 = b - A x0
	// choose rp0 such that <r0, rp0> != 0
	// p0 = r0

	// local vector packets:
	double *r_curr_local = new double[nlocal];
	double *rp0_local = new double[nlocal];
	double *p_curr_local = new double[nlocal];
	for (int i = 0; i < nlocal; i++) {
		iGlobal = i + istart;
		r_curr_local[i] = q[iGlobal];
		for (int j = 0; j < N; j++) {
			r_curr_local[i] -= dGLocal[i][j] * x_curr[j];
		}
		rp0_local[i] = r_curr_local[i] + 100;
		p_curr_local[i] = r_curr_local[i];
	}

	// gathering vector packets:
	MPI_Allgather(r_curr_local, nlocal, MPI_DOUBLE, r_curr, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(rp0_local, nlocal, MPI_DOUBLE, rp0, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(p_curr_local, nlocal, MPI_DOUBLE, p_curr, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);

	if (myrank == 0) std::cout << "------------ Launching iterations ----------------" << std::endl;
	// begin iterations
	for (int k = 0; k < maxIter; k++) {
		if (myrank == 0) std::cout << "::: iter : " << k << std::endl;

		// alpha[k] = <r[k], rp0> / <Ap[k], rp0>

		double num = 0.; double den = 0.;
		for (int i = 0; i < N; i++) {
			num += r_curr[i] * rp0[i];
		}

		double *tmp_local = new double[nlocal];
		for (int i = 0; i < nlocal; i++) {
			tmp_local[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp_local[i] += dGLocal[i][j] * p_curr[j];
			}
		}
		MPI_Allgather(tmp_local, nlocal, MPI_DOUBLE, tmp, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);

		for (int i = 0; i < N; i++) {
			den += tmp[i] * rp0[i];
		}

		alpha = num / den;

		// s[k] = r[k] - alpha[k] * A p[k]

		for (int i = 0; i < N; i++) {
			s[i] = r_curr[i] - alpha * tmp[i];
		}

		double norm = vectorNorm(s);

		if (myrank == 0) std::cout << "||s|| = " << norm << std::endl;
		if (norm < tol) {
			// x[k + 1] = x[k] + alpha[k] * p[k]

			for (int i = 0; i < N; i++) {
				x_next[i] = x_curr[i] + alpha * p_curr[i];
			}

			if (myrank == 0) std::cout << "||s|| < tol = " << tol << ", exiting iterations" << std::endl;
			delete[] tmp_local;
			break;
		}

		// omega[k] = <A s[k], s[k]> / <A s[k], A s[k]>

		num = 0; den = 0;
		for (int i = 0; i < nlocal; i++) {
			tmp_local[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp_local[i] += dGLocal[i][j] * s[j];
			}
		}
		MPI_Allgather(tmp_local, nlocal, MPI_DOUBLE, tmp, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);

		for (int i = 0; i < N; i++) {
			num += tmp[i] * s[i];
			den += tmp[i] * tmp[i];
		}		
		omega = num / den;

		// x[k + 1] = x[k] + alpha[k] * p[k] + omega[k] * s[k]

		for (int i = 0; i < N; i++) {
			x_next[i] = x_curr[i] + alpha * p_curr[i] + omega * s[i];
		}

		// r[k + 1] = s[k] - omega[k] * A s[k]

		for (int i = 0; i < N; i++) {
			r_next[i] = s[i] - omega * tmp[i];
		}

		norm = vectorNorm(r_next);
		if (myrank == 0) std::cout << "||r[k + 1]|| = " << norm << std::endl;
		if (norm < tol) {
			if (myrank == 0) std::cout << "||r[k + 1]|| < tol = " << tol << ", exiting iterations" << std::endl;
			break;
		}

		// beta[k] = (alpha[k] / omega[k]) * <r[k + 1], rp0> / <r[k], rp0>
		
		num = 0; den = 0;
		for (int i = 0; i < N; i++) {
			num += r_next[i] * rp0[i];
			den += r_curr[i] * rp0[i];
		}

		beta = (alpha / omega) * num / den;

		// p[k + 1] = r[k + 1] + beta[k] * (p[k] - omega[k] * A p[k])

		iGlobal = 0;
		for (int i = 0; i < nlocal; i++) {
			iGlobal = i + istart;
			tmp_local[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp_local[i] += dGLocal[i][j] * p_curr[j];
			}
			// p_next[iGlobal] = r_next[iGlobal] + beta * (p_curr[iGlobal] - omega * tmp[i]);
		}
		MPI_Allgather(tmp_local, nlocal, MPI_DOUBLE, tmp, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);
		delete[] tmp_local;

		for (int i = 0; i < N; i++) {
			p_next[i] = r_next[i] + beta * (p_curr[i] - omega * tmp[i]);
		}

		norm = fabs(vectorDot(r_next, rp0));
		if (myrank == 0) std::cout << "|< r[k + 1], rp0 >| = " << norm << std::endl;
		if (norm < tol) {
			// rp0 = r[k + 1]; p[k + 1] = r[k + 1]
			
			for (int i = 0; i < N; i++) {
				rp0[i] = r_next[i]; p_next[i] = r_next[i];
			}
		}
		// current = next

		for (int i = 0; i < N; i++) {
			x_curr[i] = x_next[i];
			r_curr[i] = r_next[i];
			p_curr[i] = p_next[i];
		}

		if (myrank == 0) std::cout << "===> finishing iter " << k << std::endl;
	}

	// result: x = x_next
	for (int i = 0; i < N; i++) alphas[i] = x_next[i];

	// clean up
	delete[] r_curr_local; delete[] rp0_local; delete[] p_curr_local;

	delete[] x_curr; delete[] x_next;
	delete[] r_curr; delete[] r_next;
	delete[] p_curr; delete[] p_next;
	delete[] s;

	auto endBi_CGSTAB = std::chrono::high_resolution_clock::now();
	// ========================= Solution done ===============================================

	// print solution
	if (myrank == 0) printArray1("alphas", alphas, 8);

	// potential solution G . alphas = u
	auto startMMult = std::chrono::high_resolution_clock::now();	
	for (int i = 0; i < nlocal; i++) {
		u[i] = 0;
		for (int j = 0; j < N; j++) {
			double dx = x[i] - sx[j];
			double dy = y[i] - sy[j];
			double dz = z[i] - sz[j];

			double d = sqrt(dx * dx + dy * dy + dz * dz);
			double G = 1 / (4 * M_PI * d);

			u[i] += G * alphas[j];
		}
	}
	auto endMMult = std::chrono::high_resolution_clock::now();

	if (myrank == 0) {
		printArray1("u", u, 8);
		
		std::fstream dataOut;
		dataOut.open("data.dat", std::fstream::out);
		if (!dataOut.is_open()) {
			std::cout << "unable to open output file data.dat!" << std::endl;
			return 0;
		}

		for (int i = 0; i < N; i++) {
			dataOut << B[i] << " " << L[i] << " " << u[i] << std::endl;
		}

		dataOut.close();

		// ================================ Summary ============================================

		printf("\n============================================================================\n");
		printf("======================= Program summary =====================================\n");
		std::cout << "data size = " << N << std::endl;
		std::chrono::duration<double> elapsedTotal = (endMMult - startLoad);
		std::cout << "total runtime :    " << elapsedTotal.count() << " s" << std::endl;
		std::cout << "--------------------------------------------------------------------------" << std::endl;
		std::chrono::duration<double> elapsedLoad = (endLoad - startLoad);
		std::cout << "loading data file :    " << elapsedLoad.count() << " s" << std::endl;
		std::chrono::duration<double> elapsedPrint3 = (endPrint3 - startPrint3);
		std::cout << "printing vector 3 array :    " << elapsedPrint3.count() << " s" << std::endl;
		std::chrono::duration<double> elapsedPrint1 = (endPrint1 - startPrint1);
		std::cout << "printing vector 1 array :    " << elapsedPrint3.count() << " s" << std::endl;
		std::cout << "..........................................................................." << std::endl;
		std::chrono::duration<double> elapsedMatrixGen = (endMatrixGen - startMatrixGen);
		std::cout << "generating system matrix :    " << elapsedMatrixGen.count() << " s" << std::endl;
		std::chrono::duration<double> elapsedMMult = (endMMult - startMMult);
		std::cout << "matrix * vector multiplication :    " << elapsedMMult.count() << " s" << std::endl;
		std::cout << ".... Bi-CGSTAB: .........................................................." << std::endl;
		std::chrono::duration<double> elapsedBi_CGSTAB = (endBi_CGSTAB - startBi_CGSTAB);
		std::cout << "Bi-CGSTAB solution :    " << elapsedBi_CGSTAB.count() << " s" << std::endl;
		std::cout << "--------------------------------------------------------------------------" << std::endl;
	}

	// clean up
	delete[] x; delete[] y; delete[] z;
	delete[] sx; delete[] sy; delete[] sz;
	delete[] q; delete[] d2U;

	for (int i = 0; i < N; i++) {
		delete[] dGLocal[i];
	}
	delete[] dGLocal;

	MPI_Finalize();

	return 1;
}

