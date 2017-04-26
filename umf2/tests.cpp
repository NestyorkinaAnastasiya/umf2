/*tests.cpp*/
#include "tests.h"

namespace tests
{
	double Tests::Lambda(int formNumber, double x, double y)
	{
		if (test == 1) return 1;
		if (test == 2) return 1;
		if (test == 3) return 1;
		return 1;
	}
	double Tests::Sigma(int formNumber, double x, double y)
	{
		if (test == 1) return 1;
		if (test == 2) return 1;
		if (test == 3) return 1;
		return 1;
	}
	double Tests::Ug(int formNumber, double x, double y, double t)
	{
		if (test == 1) return 1;
		if (test == 2) return 2 * t;
		if (test == 3) return 4 * t * t;
		if (test == 4) return x*x*x;
		if (test == 5) return pow(x+y,2);

	}
	double Tests::Betta(int formNumber, double x, double y)
	{
		if (test == 1)
			switch (formNumber)
			{
				//левое ребро
			case 0: return 0;
				//правое ребро
			case 1: return 0;
				//нижнее ребро 
			case 2: return 0;
				//верхнее ребро
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//левое ребро
			case 0: return 1;
				//правое ребро
			case 1: return 1;
				//нижнее ребро 
			case 2: return 1;
				//верхнее ребро
			case 3: return 1;
			}

		if (test == 3)
			switch (formNumber)
			{
				//левое ребро
			case 0: return 1;
				//правое ребро
			case 1: return 1;
				//нижнее ребро 
			case 2: return 1;
				//верхнее ребро
			case 3: return 1;
			}
	}
	double Tests::Ubetta(int formNumber, double x, double y)
	{
		if (test == 1)
			switch (formNumber)
			{
				//левое ребро
			case 0: return 0;
				//правое ребро
			case 1: return 0;
				//нижнее ребро 
			case 2: return 0;
				//верхнее ребро
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//левое ребро
			case 0: return x + y - 1;
				//правое ребро
			case 1: return x + y + 1;
				//нижнее ребро 
			case 2: return x + y - 1;
				//верхнее ребро
			case 3: return x + y + 1;
			}
		if (test == 3)
			switch (formNumber)
			{
				//левое ребро
			case 0: return 0;
				//правое ребро
			case 1: return 2*exp(x + y);
				//нижнее ребро 
			case 2: return 0;
				//верхнее ребро
			case 3: return 2 * exp(x + y);
			}
	}
	double Tests::Tetta(int formNumber, double x, double y)
	{
		if (test == 1)
			switch (formNumber)
			{
			//левое ребро
			case 0: return 0;
			//правое ребро
			case 1: return 0;
			//нижнее ребро 
			case 2: return 0;
			//верхнее ребро
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//левое ребро
			case 0: return -1;
				//правое ребро
			case 1: return 1;
				//нижнее ребро 
			case 2: return -1;
				//верхнее ребро
			case 3: return 1;
			}

		if (test == 3)
			switch (formNumber)
			{
				//левое ребро
			case 0: return -exp(x+y);
				//правое ребро
			case 1: return exp(x + y);
				//нижнее ребро 
			case 2: return -exp(x + y);
				//верхнее ребро
			case 3: return exp(x + y);
			}
	}
	double Tests::Fi(double u)
	{
		if (test == 1) return 0;
		if (test == 2) return 2;
		if (test == 3) return 4*sqrt(u);
		if (test == 4) return -6*pow(u,1/3);
		if (test == 5) return 4;
		return 0;
	}
	Tests::Tests()
	{
		ifstream fo;
		fo.open("number_of_test.txt");
		fo >> test;
	}
}