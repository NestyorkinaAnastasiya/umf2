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
		//Максимальное количество итераций в решателе
		int maxiter = 10000;
		//Максимальное количество итераций по времени
		int maxtime = 8;
		//шаг по времени
		double ht = 0.1;
		//Точность решения СЛАУ
		const double eps = 1e-14;
		//параметр релаксации
		double w;
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
		//Вектор приближенного решения на предыдущей итерации
		vector <double> u_;
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
		//Вектор невязки
		vector <double> r;
		//Вектор спуска
		vector <double> z;

		//Вычисление нормы вектора
		double Norm(const vector<double>& x);
		//Скалярное произведение векторов
		double Scalar(const vector<double>& x, const vector<double>& y);

		//Генерация СЛАУ на i-ой итерации по времени
		void GenerateSLAE();
		//LU-факторизация
		void LU();
		//Вспомогательные функции для решателя
		void LYF(const vector<double>& C, vector<double>& yl);
		void UXY(const vector<double>& C, vector<double>& yu);
		double Rel_Discrepancy();
		//Решатель ЛОС с LU-факторизацией
		void LULOS();

		double StopIteration();

	public:
		SLAE();

		void TSolve();

		~SLAE() {};
	};
}
