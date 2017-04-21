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
		//������������ ���������� ��������
		const int maxiter = 10000;
		//�������� ������� ����
		const double eps = 1e-14;
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
		//��� �� �������
		double ht;
		//������ �������
		vector <double> r;
		//������ ������
		vector <double> z;

		//���������� ����� �������
		double Norm(const vector<double>& x);
		//��������� ������������ ��������
		double Scalar(const vector<double>& x, const vector<double>& y);

		void GenerateSLAE();
		//LU-������������
		void LU();
		//��������������� ������� ��� ��������
		void MultiplyUx(vector<double> a, vector<double>& result);
		void LFx(vector<double> b, vector<double>& result);
		void LTFx(vector<double> b, vector<double>& result);
		void UFx(vector<double> b, vector<double>& result);
		void UTFx(vector<double> b, vector<double>& result);
		double IterMSG(const vector<double>& Az);
	public:
		SLAE();
		void TSolve();
		void LYF(const vector<double>& C, vector<double>& yl);
		void UXY(const vector<double>& C, vector<double>& yu);
		void LULOS(FILE * fd);
		double Rel_Discrepancy();
		void LOS(FILE * Out);
		//�������� ��� � LU-�������������
		void LU_MSG(FILE * fo);
		~SLAE() {};
	};
}
