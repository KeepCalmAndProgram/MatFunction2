#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <complex>

class Display {

protected:

	virtual void display(std::vector<int>& myVector)
	{
		for (int i = 0; i < myVector.size(); i++) {
			std::cout << "arry[" << i << "] = " << myVector[i] << std::endl;
		}
	}

	double add(double* sum, double* difference, double* sumr, double* sumu, double* pdelta, double* x1r, double* x2r, double* x1u, double* x2u) {

		if (*pdelta > 0) {
			*sum = *x1r + *x2r;
		}
		else if (pdelta < 0) {
			*sumr = *x1r + *x2r;
			*sumu = *x1u + *x2u;
			*sum = *sumr + *sumu;
		}
		return *sum;
	}
	double subtract(double* difference, double* differencer, double* differenceu, double* pdelta, double* x1r, double* x2r, double* x1u, double* x2u) {

		if (*pdelta > 0) {
			*difference = *x1r - *x2r;
		}
		else if (pdelta < 0) {
			*differencer = *x1r - *x2r;
			*differenceu = *x1u - *x2u;
			*difference = *differencer + *differenceu;
		}
		return *difference;
	}
	double multiplication(double* sum_multiplication, double* sumr, double* sumu, double* pdelta, double* x1r, double* x2r, double* x1u, double* x2u) {

		if (*pdelta > 0) {
			*sum_multiplication = ((*x1r) * (*x2r));
		}
		else if (pdelta < 0) {
			*sumr = ((*x1r) * (*x2r));
			*sumu = ((*x1u) * (*x2u));
			*sum_multiplication = *sumr + *sumu;
		}
		return *sum_multiplication;
	}
};

class MyClass : public Display {

public:

	MyClass(std::vector<int>& myVector, double err) {
		myVector.push_back(0);
		myVector.push_back(2);
		myVector.push_back(-3);

		err = 0.000001;
	}

	using Display::display;
	using Display::add;
	using Display::subtract;
	using Display::multiplication;

	double calculateDelte(double* delta, std::vector<int>& myVector);
	double sqrtNewton(double* delta, double* pdelta, double err);
	double sqrtHeron(std::vector<int>& myVector, double* pdelta);

	void addZespolS(double sr, double su);
	void addZespolR(double rr, double ru);
	void calculateTheRoot(double* delta, double* pdelta, std::vector<int>& myVector, double x1r, double x2r, double x1u, double x2u, double sr, double su, double rr, double ru);
	void displayResult(double* delta, double* pdelta, std::complex<double>& myVectorOutput, std::vector<int>& myVector, double x1r, double x2r, double x1u, double x2u);

};

double MyClass::calculateDelte(double* delta, std::vector<int>& myVector) {
	std::cout << "\n" << myVector.at(0) << "xx + " << myVector.at(1) << "x + " << myVector.at(2) << " = 0" << std::endl;
	*delta = ((myVector.at(1) * myVector.at(1)) - (4 * myVector.at(0) * myVector.at(2)));
	std::cout << "\n" << *delta << " = " << myVector.at(1) << " - 4*" << myVector.at(0) << "*" << myVector.at(2) << std::endl;

	return *delta;
}
double MyClass::sqrtNewton(double* delta, double* pdelta, double err) {
	*pdelta = *delta / 2;
	while ((*pdelta - *delta / *pdelta) > err) //0.000001
	{
		*pdelta = (*pdelta + *delta / *pdelta) / 2;
		if ((*pdelta) * (*pdelta) == *delta) break;
	}
	std::cout << "Sqrt Newtona " << *pdelta << std::endl;

	return *pdelta;
}
double MyClass::sqrtHeron(std::vector<int>& myVector, double* pdelta) {
	double p;

	p = ((myVector.at(0) + myVector.at(1) + myVector.at(2)) / 2);
	p *= -1;
	*pdelta = sqrt((p * (p - myVector.at(0)) * (p - myVector.at(1)) * (p - myVector.at(2))));
	std::cout << "Sqrt Herona " << *pdelta << std::endl;

	return *pdelta;
}
void MyClass::addZespolS(double sr, double su) {
	su += sr;
	std::cout << "\nsu = " << su << std::endl;
}
void MyClass::addZespolR(double rr, double ru) {
	ru += rr;
	std::cout << "ru = " << ru << std::endl;
}

void MyClass::calculateTheRoot(double* delta, double* pdelta, std::vector<int>& myVector, double x1r, double x2r, double x1u, double x2u, double sr, double su, double rr, double ru) {

	if (*delta >= 0) {
		*pdelta = sqrt(*pdelta);
		x1r = ((-myVector.at(1) - *pdelta) / (2.0 * myVector.at(0))); //delta +
		x2r = ((-myVector.at(1) + *pdelta) / (2.0 * myVector.at(0)));

	}
	else if (*delta < 0) {
		//FOR COMPLEX NUMBERS
		*pdelta = sqrt((sqrt(myVector.at(0) * myVector.at(0) + myVector.at(1) * myVector.at(1) + myVector.at(0)) / 2.0) + (sqrt(myVector.at(0) * myVector.at(0) + myVector.at(1) * myVector.at(1) - myVector.at(0)) / 2.0));
		x1r = ((-myVector.at(1)) / (2.0 * myVector.at(0)));   //delta -
		x2r = x1r;
		x1u = ((-*pdelta) / (2.0 * myVector.at(0)));
		x2u = -x1u;
		x1r += x1u;
		x2r += x2u;
	}

	sr = x1r + x2r;
	su = x1u + x2u;
	rr = x1r - x2r;
	ru = x1u - x2u;
	addZespolS(sr, su);
	addZespolR(rr, ru);

}
void MyClass::displayResult(double* delta, double* pdelta, std::complex<double>& myVectorOutput, std::vector<int>& myVector, double x1r, double x2r, double x1u, double x2u) {

	if (delta > 0) {
		std::cout << "          _____" << std::endl;
		std::cout << "pdelta = Vdelta" << std::endl;
		std::cout << "x1r = (" << myVector.at(1) << "-" << *pdelta << ")/2*" << myVector.at(0) << std::endl;
		std::cout << "x1r = (" << myVector.at(0) << "+" << *pdelta << ")/2*" << myVector.at(0) << std::endl;
	}
	else if (delta < 0) {
		std::cout << "          _______" << std::endl;
		std::cout << "pdelta = V|delta|" << std::endl;

		std::cout << "x1r = -" << myVector.at(1) << "/2*" << myVector.at(0) << std::endl;
		std::cout << "x2r = " << x1r << std::endl;
		std::cout << "x1u = " << *pdelta << "/2*" << myVector.at(0) << std::endl;
		std::cout << "x2u = " << x1u << std::endl;
		std::cout << x1r << "+" << x1u << std::endl;
		std::cout << x2r << "+" << x2u << std::endl;
	}
	else if (delta == 0) {
		std::cout << "x1r = " << myVector.at(1) << "/2*" << myVector.at(0) << std::endl;
	}
	
	if (myVectorOutput.real(myVector.at(0)) != NULL) {
		std::cout << " a) REAL : " << myVectorOutput.real(myVector.at(0)) << std::endl;

		if (myVectorOutput.imag(myVector.at(0)) != NULL) {
			if (myVectorOutput.imag(myVector.at(0)) > 0)
				std::cout << "\t + " << myVectorOutput.imag(myVector.at(0)) << "i\n" << std::endl;
			else if (myVectorOutput.imag(myVector.at(0)) < 0)
				std::cout << "\t + (" << myVectorOutput.imag(myVector.at(0)) << ") * i\n" << std::endl;
		}
	}
	else {
		if (myVectorOutput.imag(myVector.at(0)) != NULL) {
			std::cout << " a) IMAG : " << myVectorOutput.real(myVector.at(0)) << std::endl;

			if (myVectorOutput.imag(myVector.at(0)) != NULL) {
				if (myVectorOutput.imag(myVector.at(0)) > 0)
					std::cout << "\t + " << myVectorOutput.imag(myVector.at(0)) << "i\n" << std::endl;
				else if (myVectorOutput.imag(myVector.at(0)) < 0)
					std::cout << "\t + (" << myVectorOutput.imag(myVector.at(0)) << ") * i\n" << std::endl;
			}
		}
	}

	if (myVectorOutput.real(myVector.at(1)) != NULL) {
		std::cout << " b) REAL : " << myVectorOutput.real(myVector.at(1)) << std::endl;

		if (myVectorOutput.imag(myVector.at(1)) != NULL) {
			if (myVectorOutput.imag(myVector.at(1)) > 0)
				std::cout << "\t + " << myVectorOutput.imag(myVector.at(1)) << "i\n" << std::endl;
			else if (myVectorOutput.imag(myVector.at(1)) < 0)
				std::cout << "\t + (" << myVectorOutput.imag(myVector.at(1)) << ") * i\n" << std::endl;
		}
	}
	else {
		if (myVectorOutput.imag(myVector.at(1)) != NULL) {
			std::cout << " b) IMAG : " << myVectorOutput.real(myVector.at(1)) << std::endl;

			if (myVectorOutput.imag(myVector.at(1)) != NULL) {
				if (myVectorOutput.imag(myVector.at(1)) > 0)
					std::cout << "\t + " << myVectorOutput.imag(myVector.at(1)) << "i\n" << std::endl;
				else if (myVectorOutput.imag(myVector.at(1)) < 0)
					std::cout << "\t + (" << myVectorOutput.imag(myVector.at(1)) << ") * i\n" << std::endl;
			}
		}
	}

	if (myVectorOutput.real(myVector.at(2)) != NULL) {
		std::cout << " c) REAL : " << myVectorOutput.real(myVector.at(2)) << std::endl;

		if (myVectorOutput.imag(myVector.at(2)) != NULL) {
			if (myVectorOutput.imag(myVector.at(2)) > 0)
				std::cout << "\t + " << myVectorOutput.imag(myVector.at(2)) << "i\n" << std::endl;
			else if (myVectorOutput.imag(myVector.at(2)) < 0)
				std::cout << "\t + (" << myVectorOutput.imag(myVector.at(2)) << ") * i\n" << std::endl;
		}
	}
	else {
		if (myVectorOutput.imag(myVector.at(2)) != NULL) {
			std::cout << " c) IMAG : " << myVectorOutput.real(myVector.at(2)) << std::endl;

			if (myVectorOutput.imag(myVector.at(2)) != NULL) {
				if (myVectorOutput.imag(myVector.at(2)) > 0)
					std::cout << "\t + " << myVectorOutput.imag(myVector.at(2)) << "i\n" << std::endl;
				else if (myVectorOutput.imag(myVector.at(2)) < 0)
					std::cout << "\t + (" << myVectorOutput.imag(myVector.at(2)) << ") * i\n" << std::endl;
			}
		}
	}

	std::vector<double>myVectorSecond;

	myVectorSecond.push_back(1.0);
	myVectorSecond.push_back(1.0);
	for (int i = 0; i < myVectorSecond.size(); i++)
	{
		myVectorSecond.at(0) *= myVectorOutput.real(myVector.at(i));
		myVectorSecond.at(1) *= myVectorOutput.imag(myVector.at(i));
	}
	std::cout << "Complex numbers : " << myVectorSecond.front() << " -- " << myVectorSecond.back() << std::endl;
}

int main()
{
	std::vector<int> myVector;
	std::complex<double> myVectorOutput;

	double sr, su, rr, ru;
	sr = su = rr = ru = 0;
	double delta, pdelta, x1r, err, x2r, x1u, x2u, x;
	delta = pdelta = err = x1r = x2r = x1u = x2u = x = 0;
	double sum, difference, sum_multiplication, sumr, sumu, differencer, differenceu;
	sum = difference = sumr = sumu = differencer = differenceu = 0;

	MyClass* obj = new MyClass(myVector, err);


	myVector.reserve(10);
	myVector.shrink_to_fit();

	std::cout << " capasity vector : " << myVector.capacity() << std::endl;


	if (myVector.at(0) == 0 && myVector.at(1) == 0 && myVector.at(2) == 0) {
		std::cout << "Wrong data !!!" << std::endl;
		return 0;
	}
	if (myVector.at(0) == 0) {
		if (myVector.at(1) != 0) {
			x1r = -(double)myVector.at(2) / (double)myVector.at(1);
			std::cout << x1r << " = " << myVector.at(2) << "/" << myVector.at(1) << std::endl;
			
			if (x1r > 0) {
				x1r = sqrt(x1r);
				std::cout << "Element = " << x1r << std::endl;
			}
			else if (x1r < 0) {
				x1r = sqrt(-1 * x1r);
				std::cout << "Element = " << x1r << "i" << std::endl;
			}
			return 0;
		}
		else if (myVector.at(1) == 0 && myVector.at(2) != 0) {
			std::cout << "Contradictory equation" << std::endl;
			return 0;
		}
		else if (myVector.at(1) == 0 && myVector.at(2) == 0) {
			std::cout << "Identity equation" << std::endl;
			return 0;
		}
	}

	obj->display(myVector);
	obj->calculateDelte(&delta, myVector);

	if (err == 0) {
		pdelta = sqrt(delta);
		obj->calculateTheRoot(&delta, &pdelta, myVector, x1r, x2r, x1u, x2u, sr, su, rr, ru);
		obj->displayResult(&delta, &pdelta, myVectorOutput, myVector, x1r, x2r, x1u, x2u);
	}
	else if (err > 0.1) {
		obj->sqrtNewton(&delta, &pdelta, err);
		obj->calculateTheRoot(&delta, &pdelta, myVector, x1r, x2r, x1u, x2u, sr, su, rr, ru);
		obj->displayResult(&delta, &pdelta, myVectorOutput, myVector, x1r, x2r, x1u, x2u);
	}
	else if (err > 0 && err < 0.1) {
		obj->sqrtHeron(myVector, &pdelta);
		obj->calculateTheRoot(&delta, &pdelta, myVector, x1r, x2r, x1u, x2u, sr, su, rr, ru);
		obj->displayResult(&delta, &pdelta, myVectorOutput, myVector, x1r, x2r, x1u, x2u);
	}
	obj->add(&sum, &difference, &sumr, &sumu, &pdelta, &x1r, &x2r, &x1u, &x2u);
	obj->subtract(&difference, &differencer, &differenceu, &pdelta, &x1r, &x2r, &x1u, &x2u);
	obj->multiplication(&sum_multiplication, &sumr, &sumu, &pdelta, &x1r, &x2r, &x1u, &x2u);

	myVector.clear();

	return 0;
}