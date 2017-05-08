/*tests.h*/
#include "integration.h"

namespace tests
{
	struct Tests
	{
		int test;
		double Lambda(double x, double y, double t);
		double Sigma(double x, double y, double t);
		double Ug(double x, double y, double t);
		double Betta(int formNumber, double x, double y);
		double Ubetta(int formNumber, double x, double y);
		double Tetta(int formNumber, double x, double y);
		double Fi(double u, double x, double y, double t);
		double dFdq(double u, double x, double y, double t);
		Tests();
	};
}
