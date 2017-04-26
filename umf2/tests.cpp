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
				//����� �����
			case 0: return 0;
				//������ �����
			case 1: return 0;
				//������ ����� 
			case 2: return 0;
				//������� �����
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//����� �����
			case 0: return 1;
				//������ �����
			case 1: return 1;
				//������ ����� 
			case 2: return 1;
				//������� �����
			case 3: return 1;
			}

		if (test == 3)
			switch (formNumber)
			{
				//����� �����
			case 0: return 1;
				//������ �����
			case 1: return 1;
				//������ ����� 
			case 2: return 1;
				//������� �����
			case 3: return 1;
			}
	}
	double Tests::Ubetta(int formNumber, double x, double y)
	{
		if (test == 1)
			switch (formNumber)
			{
				//����� �����
			case 0: return 0;
				//������ �����
			case 1: return 0;
				//������ ����� 
			case 2: return 0;
				//������� �����
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//����� �����
			case 0: return x + y - 1;
				//������ �����
			case 1: return x + y + 1;
				//������ ����� 
			case 2: return x + y - 1;
				//������� �����
			case 3: return x + y + 1;
			}
		if (test == 3)
			switch (formNumber)
			{
				//����� �����
			case 0: return 0;
				//������ �����
			case 1: return 2*exp(x + y);
				//������ ����� 
			case 2: return 0;
				//������� �����
			case 3: return 2 * exp(x + y);
			}
	}
	double Tests::Tetta(int formNumber, double x, double y)
	{
		if (test == 1)
			switch (formNumber)
			{
			//����� �����
			case 0: return 0;
			//������ �����
			case 1: return 0;
			//������ ����� 
			case 2: return 0;
			//������� �����
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//����� �����
			case 0: return -1;
				//������ �����
			case 1: return 1;
				//������ ����� 
			case 2: return -1;
				//������� �����
			case 3: return 1;
			}

		if (test == 3)
			switch (formNumber)
			{
				//����� �����
			case 0: return -exp(x+y);
				//������ �����
			case 1: return exp(x + y);
				//������ ����� 
			case 2: return -exp(x + y);
				//������� �����
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