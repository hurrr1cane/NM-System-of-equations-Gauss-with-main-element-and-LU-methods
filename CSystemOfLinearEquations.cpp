#include "CSystemOfLinearEquations.h"
#include <iomanip>


CSystemOfLinearEquations::CSystemOfLinearEquations(int countOfUnknown) {
	this->countOfUnknown = countOfUnknown;
	matrixOfCoefficients = new double*[countOfUnknown];
	for (int i = 0; i < countOfUnknown; i++) {
		matrixOfCoefficients[i] = new double[countOfUnknown + 1];
	}
	matrixOfCoefficients2 = new double* [countOfUnknown];
	for (int i = 0; i < countOfUnknown; i++) {
		matrixOfCoefficients2[i] = new double[countOfUnknown + 1];
	}
	vectorOfSolutions = new double[countOfUnknown];
	vectorOfSolutions2 = new double[countOfUnknown];
	vectorOfFreeMembers = new double[countOfUnknown];
	for (int i = 0; i < countOfUnknown; i++) {
		vectorOfSolutions[i] = 0;
		vectorOfSolutions2[i] = 0;
		vectorOfFreeMembers[i] = matrixOfCoefficients[i][countOfUnknown];
	}
}

double* CSystemOfLinearEquations::solveUsingGaussWithChoosingMainElement() {
	std::cout << "Gauss method:\nMaking triangular matrix...\n";
	for (int i = 0; i < countOfUnknown; i++) {
		//this->findMaximumInColumnAndMakeItFirst(i);
		std::cout << "Step " << i + 1 << ":\n" << *this;
		for (int j = i + 1; j < countOfUnknown; j++) {
			//find m
			double m = -(matrixOfCoefficients[j][i] / matrixOfCoefficients[i][i]);
			for (int k = 0; k < countOfUnknown + 1; k++) {
				matrixOfCoefficients[j][k] = matrixOfCoefficients[j][k] + m * matrixOfCoefficients[i][k];
			}
		}
	}
	std::cout << "\nUpper-triangle matrix:\n" << *this << "\n";
	for (int i = countOfUnknown - 1; i >= 0; i--) {
		double sum = 0;
		for (int j = 0; j < countOfUnknown; j++) {
			sum += -matrixOfCoefficients[i][j] * vectorOfSolutions[j];
		}
		vectorOfSolutions[i] = (matrixOfCoefficients[i][countOfUnknown] + sum) / matrixOfCoefficients[i][i];
	}
	return vectorOfSolutions;
}



double* CSystemOfLinearEquations::solveUsingLUFactorisation() {

	double** U = new double* [countOfUnknown], ** L = new double* [countOfUnknown];
	double* coefficients = new double[(int)(((float)countOfUnknown - 1) / 2) * countOfUnknown];
	double* yVector = new double[countOfUnknown];

	for (int i = 0; i < countOfUnknown; i++) {
		U[i] = new double[countOfUnknown];
		L[i] = new double[countOfUnknown];
		yVector[i] = 0;
	}

	for (int i = 0; i < countOfUnknown; i++) {
		for (int j = 0; j < countOfUnknown; j++) {
			U[i][j] = 0;
			L[i][j] = 0;
		}
	}

	std::cout << "\nLU method:";
	printMatrix(matrixOfCoefficients2, countOfUnknown, "Original");

	//first steps in defining L and U matrixes
	for (int i = 0; i < countOfUnknown; i++) {
		L[i][0] = matrixOfCoefficients2[i][0];
		U[i][i] = 1;
		U[0][i] = matrixOfCoefficients2[0][i] / matrixOfCoefficients2[0][0];
	}

	//Doing all the rest
	for (int line = 1; line < countOfUnknown; line++) {
		
		//Finding l-elements
		for (int column = 1; column <= line; column++) {

			//sum of l u
			double sum = 0;

			for (int k = 0; k < column; k++) {
				sum -= L[line][k] * U[k][column];
			}

			L[line][column] = matrixOfCoefficients2[line][column] + sum;
		}

		//Finding u-elements
		for (int column = line+1; column < countOfUnknown; column++) {
			
			//sum of l u
			double sum = 0;
			for (int k = 0; k < line; k++) {
				sum -= L[line][k] * U[k][column];
			}

			U[line][column] = (1 / L[line][line]) * (matrixOfCoefficients2[line][column] + sum);
		}
	}

	printMatrix(L, countOfUnknown, "L");
	printMatrix(U, countOfUnknown, "U");

	//Finding Y-vector
	for (int i = 0; i < countOfUnknown; i++) {
		double sum = 0;
		for (int j = 0; j < countOfUnknown; j++) {
			sum += -L[i][j] * yVector[j];
		}
		yVector[i] = (vectorOfFreeMembers[i] + sum) / L[i][i];
	}

	printVector(yVector, countOfUnknown, "Y-vector");

	//finding x-vector
	for (int i = countOfUnknown - 1; i >= 0; i--) {
		double sum = 0;
		for (int j = 0; j < countOfUnknown; j++) {
			sum += -U[i][j] * vectorOfSolutions2[j];
		}
		vectorOfSolutions2[i] = (yVector[i] + sum) / U[i][i];
	}

	for (int i = 0; i < countOfUnknown; i++) {
		delete[] U[i];
		delete[] L[i];
	}
	delete[] yVector;
	delete[] U;
	delete[] L;
	return vectorOfSolutions2;

}

void CSystemOfLinearEquations::findMaximumInColumnAndMakeItFirst(int column) {
	double max = fabs(matrixOfCoefficients[column][column]);
	int savedIndex = column;
	for (int i = column; i < countOfUnknown; i++) {
		if (fabs(matrixOfCoefficients[i][column]) > max) {
			max = fabs(matrixOfCoefficients[i][column]);
			savedIndex = i;
		}
	}
	
	//Saving one row
	double* savedRow = matrixOfCoefficients[savedIndex];
	//Switching two rows
	matrixOfCoefficients[savedIndex] = matrixOfCoefficients[column];
	matrixOfCoefficients[column] = savedRow;
}

std::istream& operator>>(std::istream& in, CSystemOfLinearEquations object) {
	for (int i = 0; i < object.countOfUnknown; i++) {
		std::cout << "Enter coefficients of " << i + 1 << " line: ";
		for (int j = 0; j < object.countOfUnknown + 1; j++) {
			in >> object.matrixOfCoefficients[i][j];
			object.matrixOfCoefficients2[i][j] = object.matrixOfCoefficients[i][j];
			object.vectorOfFreeMembers[i] = object.matrixOfCoefficients[i][object.countOfUnknown];
		}
	}
	return in;
}

std::ostream& operator<<(std::ostream& out, CSystemOfLinearEquations object) {
	out << std::setprecision(4);
	for (int i = 0; i < object.countOfUnknown; i++) {
		for (int j = 0; j < object.countOfUnknown + 1; j++) {
			if (j == 0) {
				out << std::setw(8) << object.matrixOfCoefficients[i][j] << std::resetiosflags(std::ios_base::adjustfield) << std::resetiosflags(std::ios_base::showpos) << "x" << j + 1;
			}
			else if (j > 0 && j != object.countOfUnknown) {
				out << std::showpos << std::setw(8) << object.matrixOfCoefficients[i][j] << std::resetiosflags(std::ios_base::adjustfield) << std::resetiosflags(std::ios_base::showpos) << "x" << j + 1;
			}
			else if (j = object.countOfUnknown) {
				out << "=" << object.matrixOfCoefficients[i][j];
			}
		}
		out << std::endl;
	}
	return out;
}

void printMatrix(double** Matrix, int a, const char* string) {
	std::cout << "\n" << string << "-matrix:\n";
	for (int i = 0; i < a; i++) {
		for (int j = 0; j < a; j++) {
			std::cout << std::setw(6) << std::setprecision(3) << Matrix[i][j] << " ";
		}
		std::cout << std::resetiosflags(std::ios_base::adjustfield) << std::resetiosflags(std::ios_base::showpos) << std::endl;
	}
}

void printVector(double* Vector, int a, const char* string) {
	std::cout << "\n" << string << ":\n";
	for (int i = 0; i < a; i++) {
		std::cout << Vector[i] << " ";
	}
	std::cout << std::resetiosflags(std::ios_base::adjustfield) << std::resetiosflags(std::ios_base::showpos) << std::endl;
}

double Vyznachnyk(CSystemOfLinearEquations& object) {
	double** Matrix = object.matrixOfCoefficients;
	return (Matrix[0][0] * Matrix[1][1] * Matrix[2][2]\
		+ Matrix[0][1] * Matrix[1][2] * Matrix[2][0]\
		+ Matrix[0][2] * Matrix[1][0] * Matrix[2][1]\
		- (Matrix[0][2] * Matrix[1][1] * Matrix[2][0]\
			+ Matrix[0][1] * Matrix[1][0] * Matrix[2][2]\
			+ Matrix[0][0] * Matrix[1][2] * Matrix[2][1]));
}