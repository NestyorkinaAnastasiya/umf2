/*slae.h*/
#pragma once
#include "tests.h"
using namespace matrix;
using namespace basis;
using namespace integration;
using namespace tests;

namespace slae
{
	class SLAE : private Basis, private GaussIntegration
	{
		//����������� ������
		int n;
		//������������ ���������� �������� � ��������
		int maxiter = 10000;
		//������������ ���������� �������� �� �������
		int maxtime = 8;
		//��� �� �������
		double ht = 0.1;
		//�������� ������� ����
		const double eps = 1e-14;
		//�������� ����������
		double w;
		//�����
		Grid grid;
		//��������� �������� �������
		Tests tests;
		//���������� �������
		Matrix A;
		Matrix globalM;
		//��������� �������
		//������� ��������
		array<array<double, 9>, 9> G;
		//������� �����
		array<array<double, 9>, 9> M;
		//��������� ������ ������ �����
		array <double, 9> locF;
		//���������� ������ ������ �����
		vector <double> F;
		//������ ������������� ������� �� ���������� ��������
		vector <double> u_;
		//������ ������������� �������
		vector <double> u;
		//����� ������� ������ �����
		double normF;

		//������ ��������� ������ ��������
		void CalculateG(int elementNumber);
		//������ ��������� ������ ����
		void CalculateM(int elementNumber);
		//������ ��������� ������ ������
		void CalculateLocalF(int elementNumber);
		//������� ���������� �������� � ����������
		void AddElementToGlobalMatrix(Matrix &B, int i, int j, double element);
		//������ ��������� ������(��������) � ���������� � ����������
		void CalculateLocals(int elementNumber);

		//������ ������� ����� ��� ������� �������� 
		array<double, 3> g;
		//���������� ������ ����� ��� 1��� �������� �������
		void Calculate_g(int formNumber, int orientation, int elNumber);
		//���������� 1��� �������� ������� ��� ������ ����
		void CalculateBoundaries1ForNode(int node, double gi, double weight);
		//���� ������� �������� �������
		void CalculateBoundaries1(int number);
		//���� ������� �������� �������
		void CalculateBoundaries2(int number);
		//���� �������� �������� �������
		void CalculateBoundaries3(int number);

		//���������� ������� � �������������
		vector <double> L;
		vector <double> D;
		vector <double> U;

		//������� �������� �������
		double t;
		//������ �������
		vector <double> r;
		//������ ������
		vector <double> z;

		//���������� ����� �������
		double Norm(const vector<double>& x);
		//��������� ������������ ��������
		double Scalar(const vector<double>& x, const vector<double>& y);

		//��������� ���� �� i-�� �������� �� �������
		void GenerateSLAE();
		//LU-������������
		void LU();
		//��������������� ������� ��� ��������
		void LYF(const vector<double>& C, vector<double>& yl);
		void UXY(const vector<double>& C, vector<double>& yu);
		double Rel_Discrepancy();
		//�������� ��� � LU-�������������
		void LULOS();

		double StopIteration();

	public:
		SLAE();

		void TSolve();

		~SLAE() {};
	};
}
