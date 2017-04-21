/*tests.cpp*/
#include "tests.h"

namespace tests
{
	double Tests::Lambda(int formNumber, double x, double y)
	{
		if (test == 1) return 1;
		if (test == 2)
			switch (formNumber)
			{
			case 0: return 1;
			}
		if (test == 3)
			switch (formNumber)
			{
			case 0: return 1;
			}
		return 0.0;
	}
	double Tests::Sigma(int formNumber, double x, double y)
	{
		if (test == 1) return 1;
		if (test == 2)
			switch (formNumber)
			{
			case 0: return x + y;
			}
		if (test == 3)
			switch (formNumber)
			{
			case 0: return x + y;
			}
		return 0.0;
	}
	double Tests::Ug(int formNumber, double x, double y, double t)
	{
		if (test == 1) return 1;

		if (test == 2)
			switch (formNumber)
			{
			case 0: return x+y;
			}

		if (test == 3)
			switch (formNumber)
			{
			case 0:  exp(x + y);
			}
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
	double Tests::Fi(double x, double y)
	{
		if (test == 1) return 0;
		if (test == 2) return x*x + y*y + 2*x*y;
		if (test == 3) return (x + y - 2) *  exp(x + y);
		return 0.0;
	}
	Tests::Tests()
	{
		ifstream fo;
		fo.open("number_of_test.txt");
		fo >> test;
	}
}