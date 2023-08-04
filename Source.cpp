#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "CSystemOfLinearEquations.h"

int main(void) {
	int rang = 3;

	std::cout << "Do you want to read data from file? Yes - 1, No - 0: ";
	int checker = 0;
	std::cin >> checker;
	std::ifstream file;
	if (checker) {
		std::string FileName;
		std::cout << "Enter file path and name: ";
		std::cin >> FileName;
		file.open(FileName);
		file >> rang;
	}
	else {
		std::cout << "Enter the rang of the matrix: ";
		std::cin >> rang;		
	}

	CSystemOfLinearEquations MySystem1(rang);
	if (checker) {
		file >> MySystem1;
	}
	else {
		std::cin >> MySystem1;
	}
	
	std::cout << "\nEntered system of equations:\n" << MySystem1 << "\n";
	
	double det = 1;
	if (rang == 3) {
		det = Vyznachnyk(MySystem1);
		std::cout << "The determinant is: " << det << std::endl;
	}
	if (det != 0 || rang != 3) {

		double* SolutionsByGauss = MySystem1.solveUsingGaussWithChoosingMainElement();

		for (int i = 0; i < rang; i++) {
			std::cout << "x" << i + 1 << " = " << SolutionsByGauss[i] << ", ";
		}
		std::cout << "\n";

		double* SolutionsByLU = MySystem1.solveUsingLUFactorisation();

		std::cout << std::endl;
		for (int i = 0; i < rang; i++) {
			std::cout << "x" << i + 1 << " = " << std::setprecision(4) << SolutionsByLU[i] << ", ";
		}
		std::cout << "\n\n\n";
	}
	else if (det == 0) {
		std::cout << "System of equations cannot be solved";
		return 1;
	}
	return 0;
}