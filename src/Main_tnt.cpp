#include <iostream>
#include "tnt.h"

int main()
{

        std::cout << "Hello" << std::endl; 
		TNT::Array1D<double> A (3,1);
		TNT::Array1D<double> B (3,1);


		A+=B;

		std::cout << A;
        return 0;
}
