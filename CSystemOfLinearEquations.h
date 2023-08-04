#pragma once
#include <iostream>

class CSystemOfLinearEquations {
	int countOfUnknown;
	double** matrixOfCoefficients;
	double* vectorOfSolutions;

	double** matrixOfCoefficients2;
	double* vectorOfSolutions2;
	double* vectorOfFreeMembers;
public:
	CSystemOfLinearEquations(int countOfUnknown);
	friend std::istream& operator>> (std::istream& in, CSystemOfLinearEquations object);
	friend std::ostream& operator<< (std::ostream& out, CSystemOfLinearEquations object);
	double* solveUsingGaussWithChoosingMainElement();
	double* solveUsingLUFactorisation();
	friend double Vyznachnyk(CSystemOfLinearEquations& object);

private:
	void findMaximumInColumnAndMakeItFirst(int column);
};

void printVector(double* Vector, int a, const char* string);
void printMatrix(double** Matrix, int a, const char* string);
