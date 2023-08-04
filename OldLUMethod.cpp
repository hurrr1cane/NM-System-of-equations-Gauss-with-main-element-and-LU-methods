//double* CSystemOfLinearEquations::solveUsingLUFactorisation() {
//
//	double** U = new double* [countOfUnknown], ** L = new double* [countOfUnknown];
//	double* coefficients = new double[(int)(((float)countOfUnknown - 1) / 2) * countOfUnknown];
//	double* yVector = new double[countOfUnknown];
//
//	for (int i = 0; i < countOfUnknown; i++) {
//		U[i] = new double[countOfUnknown];
//		L[i] = new double[countOfUnknown];
//	}
//	//copying the original array to U array
//	for (int i = 0; i < countOfUnknown; i++) {
//		for (int j = 0; j < countOfUnknown; j++) {
//			U[i][j] = matrixOfCoefficients2[i][j];
//			L[i][j] = 0;
//		}
//		yVector[i] = 0;
//	}
//
//	std::cout << "\nLU method:";
//	printMatrix(U, countOfUnknown, "Original");
//
//	double m;
//	//making an U matrix from original by getting a triangular from it
//	int counter = 0;
//	for (int i = 0; i < countOfUnknown; i++) {
//		for (int j = i + 1; j < countOfUnknown; j++) {
//			//find m
//			m = -(U[j][i] / U[i][i]);
//			coefficients[counter++] = m; //saving m for making an L - matrix
//			for (int k = 0; k < countOfUnknown; k++) {
//				U[j][k] = U[j][k] + m * U[i][k];
//			}
//		}
//	}
//	printMatrix(U, countOfUnknown, "U");
//	//making an L - matrix
//	for (int i = 0; i < countOfUnknown; i++) {
//		L[i][i] = 1;
//	}
//	counter = 0;
//	int counter2 = 0;
//	for (int i = 0; i < countOfUnknown; i++) {
//		for (int j = 0; j < counter2; j++) {
//			L[i][j] = -coefficients[counter++];
//		}
//		counter2++;
//	}
//
//	printMatrix(L, countOfUnknown, "L");
//	//finding y- vector
//	for (int i = 0; i < countOfUnknown; i++) {
//		double sum = 0;
//		for (int j = 0; j < countOfUnknown; j++) {
//			sum += -L[i][j] * yVector[j];
//		}
//		yVector[i] = (vectorOfFreeMembers[i] + sum) / L[i][i];
//	}
//	printVector(yVector, countOfUnknown, "Y-vector");
//
//	//finding x-vector
//	for (int i = countOfUnknown - 1; i >= 0; i--) {
//		double sum = 0;
//		for (int j = 0; j < countOfUnknown; j++) {
//			sum += -U[i][j] * vectorOfSolutions2[j];
//		}
//		vectorOfSolutions2[i] = (yVector[i] + sum) / U[i][i];
//	}
//
//	for (int i = 0; i < countOfUnknown; i++) {
//		delete[] U[i];
//		delete[] L[i];
//	}
//	delete[] yVector;
//	delete[] U;
//	delete[] L;
//	return vectorOfSolutions2;
//
//}