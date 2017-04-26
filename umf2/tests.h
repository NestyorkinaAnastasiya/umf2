/*tests.h*/
#include "integration.h"

namespace tests
{
	struct Tests
	{
		int test;
		double Lambda(int formNumber, double x, double y);
		double Sigma(int formNumber, double x, double y);
		double Ug(int formNumber, double x, double y, double t);
		double Betta(int formNumber, double x, double y);
		double Ubetta(int formNumber, double x, double y);
		double Tetta(int formNumber, double x, double y);
		double Fi(double u);
		Tests();
	};
}
