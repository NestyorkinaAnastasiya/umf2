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
		//Размерность задачи
		int n;
		//Максимальное количество итераций
		const int maxiter = 10000;
		//Точность решения СЛАУ
		const double eps = 1e-14;
		//Сетка
		Grid grid;
		//Хранилище тестовых функций
		Tests tests;
		//Глобальная матрица
		Matrix A;
		Matrix globalM;
		//Локальные матрицы
		//Матрица жёсткости
		array<array<double, 9>, 9> G;
		//Матрица массы
		array<array<double, 9>, 9> M;
		//Локальный вектор правой части
		array <double, 9> locF;
		//Глобальный вектор правой части
		vector <double> F;
		//Вектор приближенного решения
		vector <double> u;
		//Норма вектора правой части
		double normF;

		//Сборка локальных матриц жёсткости
		void CalculateG(int elementNumber);
		//Сборка локальных матриц масс
		void CalculateM(int elementNumber);
		//Сборка локальных правых частей
		void CalculateLocalF(int elementNumber);
		//Добавка локального элемента в глобальный
		void AddElementToGlobalMatrix(Matrix &B, int i, int j, double element);
		//Сборка локальных матриц(векторов) и добавление в глобальные
		void CalculateLocals(int elementNumber);

		//Вектор праввой части для первого краевого 
		array<double, 3> g;
		//Нахождение правой части для 1ого краевого условия
		void Calculate_g(int formNumber, int orientation, int elNumber);
		//Вычисление 1ого краевого условия для одного узла
		void CalculateBoundaries1ForNode(int node, double gi, double weight);
		//Учёт первого краевого условия
		void CalculateBoundaries1(int number);
		//Учёт второго краевого условия
		void CalculateBoundaries2(int number);
		//Учёт третьего краевого условия
		void CalculateBoundaries3(int number);

		//Компоненты матрицы с факторизацией
		vector <double> L;
		vector <double> D;
		vector <double> U;

		//Текущее значение времени
		double t;
		//шаг по времени
		double ht;
		//Вектор невязки
		vector <double> r;
		//Вектор спуска
		vector <double> z;

		//Вычисление нормы вектора
		double Norm(const vector<double>& x);
		//Скалярное произведение векторов
		double Scalar(const vector<double>& x, const vector<double>& y);

		void GenerateSLAE();
		//LU-факторизация
		void LU();
		//Вспомогательные функции для решателя
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
		//Решатель МСГ с LU-факторизацией
		void LU_MSG(FILE * fo);
		~SLAE() {};
	};
}
