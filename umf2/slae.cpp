/*slae.cpp*/
#include "slae.h"

namespace slae
{
	SLAE::SLAE()
	{
		//Построение сетки
		grid.DoPartition();
		//Размерность задачи соответствует общему числу базисных функций
		n = (2 * grid.nx - 1)*(2 * grid.ny - 1);

		F.resize(n);
		u.resize(n);
		u_.resize(n);
		r.resize(n);
		z.resize(n);
		maxtime = 8;
		w = 1;
		//Генерация портрета матрицы и её инициализация
		A.CreatePortret(n, grid);
		globalM.CreatePortret(n, grid);

	}

	//Сборка локальных матриц жёсткости
	void SLAE::CalculateG(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double g1, g2, ksi, etta, x_, y_, lambda,
			x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x,
			y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y,
			hx = x3 - x1, hy = y3 - y1,
			hx2 = hx * hx, hy2 = hy * hy,
			jacobian = hx * hy / 4.0;

		for (int i = 0; i < 9; i++)
		{	
			for (int j = i; j < 9; j++)
			{	
				g1 = 0; g2 = 0;
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi*hx; y_ = y1 + etta*hy;
					lambda = tests.Lambda(element.numberOfArea, x_, y_);

					g1 += gaussWeights[k] * dphiksi[i](ksi, etta) * dphiksi[j](ksi, etta) * lambda;
					g2 += gaussWeights[k] * dphietta[i](ksi, etta) * dphietta[j](ksi, etta) * lambda;
				}
				G[i][j] = g1 * jacobian / hx2 + g2 * jacobian / hy2;
			}
		}
		//матрица симметричная, заполняем нижний треугольник
		for (int i = 1; i < 9; i++)
			for (int j = 0; j < i; j++)
				G[i][j] = G[j][i];
	}

	//Сборка локальных матриц масс
	void SLAE::CalculateM(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double  g, ksi, etta, sigma, x_, y_, 
			x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x,
			y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y,
			hx = x3 - x1, hy = y3 - y1,
			jacobian = hx * hy / 4.0;

		for (int i = 0; i < 9; i++)
		{	
			for (int j = i; j < 9; j++)
			{	
				g = 0;
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi*hx; y_ = y1 + etta*hy;
					sigma = tests.Sigma(element.numberOfArea, x_, y_);

					g +=  sigma *  gaussWeights[k] * phi[i](ksi, etta) * phi[j](ksi, etta);
				}
				M[i][j] = g * jacobian / ht;
			}
		}
		//матрица симметричная, заполняем нижний треугольник
		for (int i = 1; i < 9; i++)
			for (int j = 0; j < i; j++)
				M[i][j] = M[j][i];
	}

	//Сборка локальных правых частей
	void SLAE::CalculateLocalF(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double fi, ksi, etta, x_, y_, ui,
			x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x,
			y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y,
			hx = x3 - x1, hy = y3 - y1,
			jacobian = hx * hy / 4.0;
		//интегрирование(Гаусс 5)
		for (int i = 0; i < 9; i++)
		{	
			locF[i] = 0;
			
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi*hx; y_ = y1 + etta*hy;

					fi = 0;
					//разложение по базису правой части
					for (int l = 0; l < 9; l++)
					{ 
						ui = 0;
						//разложение по базису предыдущего решения
						for (int k = 0; k < 9; k++)
							ui += u[grid.GetGlobalFuncNumber(elementNumber, k)] * phi[k](ksi, etta);
						fi += tests.Fi(ui) * phi[l](ksi, etta);
					}
					locF[i] += fi * gaussWeights[k] * phi[i](ksi, etta);
				}
			locF[i] *= jacobian ;
		}
	}

	//Добавка локального элемента в глобальный
	void SLAE::AddElementToGlobalMatrix(Matrix &B, int i, int j, double element)
	{
		int id;
		bool flag;

		if (i == j)
			B.di[i] += element;
		else
		{
			if (i < j)
			{
				flag = false;
				for (id = B.ig[j]; !flag && id < B.ig[j + 1]; id++)
					if (B.jg[id] == i) flag = true;

				if (flag) B.ggu[id - 1] += element;
			}
			else
			{
				flag = false;
				for (id = B.ig[i]; !flag && id < B.ig[i + 1]; id++)
					if (B.jg[id] == j) flag = true;

				if (flag) B.ggl[id - 1] += element;
			}
		}
	}

	//Сборка локальных матриц(векторов) и добавление в глобальные
	void SLAE::CalculateLocals(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		int ki, kj;

		//вычисление локальных матриц
		CalculateG(elementNumber);
		CalculateM(elementNumber);
		CalculateLocalF(elementNumber);

		for (int i = 0; i < 9; i++)
		{
			ki = element.dof[i];
			for (int j = 0; j < 9; j++)
			{
				kj = element.dof[j];
				//добавка в глобальную матрицу А
				AddElementToGlobalMatrix(A, ki, kj, G[i][j] + M[i][j]);
				//добавка в глобальную матрицу масс
				AddElementToGlobalMatrix(globalM, ki, kj, M[i][j]);
			}
			//добавка в глобальную правую часть
			F[ki] += locF[i];
		}
	}

	//Генерация СЛАУ
	void SLAE::GenerateSLAE()
	{
		for (int i = 0; i < n; i++)
			A.di[i] = F[i] = globalM.di[i] = 0;
		for (int i = 0; i < A.ggl.size(); i++)
			A.ggl[i] = A.ggu[i] = globalM.ggl[i] = globalM.ggu[i] = 0;

		//Высчитывание локальных матриц(векторов) и добавление в глобальные
		for (int i = 0; i < grid.elements.size(); i++)
			CalculateLocals(i);

		vector<double> b;
		b.resize(n);

		//добавка от матрицы масс по времени
		globalM.MultiplyAx(u, b);
		for (int i = 0; i < F.size(); i++)
			F[i] += b[i];

		//Учёт краевых условий
		for (int i = 0; i < grid.ku[0].size(); i++)
			CalculateBoundaries1(i);

		for (int i = 0; i < grid.ku[1].size(); i++)
			CalculateBoundaries2(i);

		for (int i = 0; i < grid.ku[2].size(); i++)
			CalculateBoundaries3(i);

		normF = Norm(F);
	}

	//Нахождение правой части для 1ого краевого условия
	void SLAE::Calculate_g(int formNumber, int orientation, int elNumber)
	{
		Element element = grid.elements[elNumber];

		switch (orientation)
		{
			//левое ребро
		case 0:
		{
			double x = grid.nodes[element.nodes[0]].x;
			double y1 = grid.nodes[element.nodes[0]].y;
			double y3 = grid.nodes[element.nodes[2]].y;
			double y2 = (y1 + y3) / 2;

			g[0] = tests.Ug(formNumber, x, y1,t);
			g[1] = tests.Ug(formNumber, x, y2, t);
			g[2] = tests.Ug(formNumber, x, y3, t);
		}
		break;
		//правое ребро
		case 1:
		{
			double x = grid.nodes[element.nodes[1]].x;
			double y1 = grid.nodes[element.nodes[1]].y;
			double y3 = grid.nodes[element.nodes[3]].y;
			double y2 = (y1 + y3) / 2;

			g[0] = tests.Ug(formNumber, x, y1, t);
			g[1] = tests.Ug(formNumber, x, y2, t);
			g[2] = tests.Ug(formNumber, x, y3, t);
		}
		break;
		//нижнее ребро
		case 2:
		{
			double y = grid.nodes[element.nodes[0]].y;
			double x1 = grid.nodes[element.nodes[0]].x;
			double x3 = grid.nodes[element.nodes[1]].x;
			double x2 = (x1+x3)/2;
			g[0] = tests.Ug(formNumber, x1, y, t);
			g[1] = tests.Ug(formNumber, x2, y, t);
			g[2] = tests.Ug(formNumber, x3, y, t);
		}
		break;
		//верхнее ребро
		case 3:
		{
			double y = grid.nodes[element.nodes[2]].y;
			double x1 = grid.nodes[element.nodes[2]].x;
			double x3 = grid.nodes[element.nodes[3]].x;
			double x2 = (x1 + x3) / 2;
			g[0] = tests.Ug(formNumber, x1, y, t);
			g[1] = tests.Ug(formNumber, x2, y, t);
			g[2] = tests.Ug(formNumber, x3, y, t);
		}
		break;
		default:; break;
		}
	}

	//Вычисление 1ого краевого условия для одного узла
	void SLAE::CalculateBoundaries1ForNode(int node, double gi, double weight)
	{
		int id;
		F[node] = gi;
		A.di[node] = weight;

		for (int j = 0; j < n; j++)
			if (node < j)
			{
				bool flag = false;
				for (id = A.ig[j]; !flag && id <= A.ig[j + 1] - 1; id++)
					if (A.jg[id] == node) flag = true;
				if (flag) A.ggu[id - 1] = 0.0;
			}
			else
			{
				bool flag = false;
				for (id = A.ig[node]; !flag && id <= A.ig[node + 1] - 1; id++)
					if (A.jg[id] == j) flag = true;
				if (flag) A.ggl[id - 1] = 0.0;
			}
	}

	//Учёт первого краевого условия
	void SLAE::CalculateBoundaries1(int number)
	{
		Element element = grid.elements[grid.ku[0][number].elem];

		if (grid.ku[0][number].edges[0] == 1)
		{
			int indexes[3] = { element.dof[0], element.dof[3], element.dof[6] };
			Calculate_g(grid.ku[0][number].formNumber[0], 0, grid.ku[0][number].elem);
			for (int i = 0; i < 3; i++)
				CalculateBoundaries1ForNode(indexes[i], g[i], 1);

		}
		if (grid.ku[0][number].edges[1] == 1)
		{
			int indexes[3] = { element.dof[2], element.dof[5], element.dof[8] };
			Calculate_g(grid.ku[0][number].formNumber[1], 1, grid.ku[0][number].elem);
			for (int i = 0; i < 3; i++)
				CalculateBoundaries1ForNode(indexes[i], g[i], 1);
		}
		if (grid.ku[0][number].edges[2] == 1)
		{
			int indexes[3] = { element.dof[0], element.dof[1], element.dof[2] };
			Calculate_g(grid.ku[0][number].formNumber[2], 2, grid.ku[0][number].elem);
			for (int i = 0; i < 3; i++)
				CalculateBoundaries1ForNode(indexes[i], g[i], 1);
		}
		if (grid.ku[0][number].edges[3] == 1)
		{
			int indexes[3] = { element.dof[6], element.dof[7], element.dof[8] };
			Calculate_g(grid.ku[0][number].formNumber[3], 3, grid.ku[0][number].elem);
			for (int i = 0; i < 3; i++)
				CalculateBoundaries1ForNode(indexes[i], g[i], 1);
		}
	}

	//Учёт второго краевого условия
	void SLAE::CalculateBoundaries2(int number)
	{
		Element element = grid.elements[grid.ku[1][number].elem];
		double g, jacobian;
		double x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x;
		double y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y;
		double hx = x3 - x1;
		double hy = y3 - y1;
		//Левое ребро
		if (grid.ku[1][number].edges[0] == 1)
		{
			//Глобальные номера базисных функций
			array<int, 3> indexes = { element.dof[0], element.dof[3], element.dof[6] };
			//Локальные номера
			array<int, 3> id = { 0, 3, 6 };

			jacobian = hy / 2;

			for (int i = 0; i < 3; i++)
			{
				g = 0;
				for (int k = 0; k < 5; k++)
				{
					double etta = 0.5 + 0.5 * gaussPoints1[k];
					double y_ = y1 + etta*hy;
					double tetta_ = tests.Tetta(grid.ku[1][number].formNumber[0], grid.nodes[element.dof[id[i]]].x, y_);
					g += gaussWeights1[k] * phi[id[i]](0, etta) * tetta_;
				}
				//Добавляем к глобальному вектору правой части
				F[indexes[i]] += g * jacobian;
			}
		}
		//Правое ребро
		if (grid.ku[1][number].edges[1] == 1)
		{
			int indexes[3] = { element.dof[2], element.dof[5], element.dof[8] };
			int id[3] = { 2, 5, 8 };

			jacobian = hy / 2;

			for (int i = 0; i < 3; i++)
			{
				g = 0;
				for (int k = 0; k < 5; k++)
				{
					double etta = 0.5 + 0.5 * gaussPoints1[k];
					double y_ = y1 + etta*hy;
					double tetta_ = tests.Tetta(grid.ku[1][number].formNumber[1], grid.nodes[element.dof[id[i]]].x, y_);
					g += gaussWeights1[k] * phi[id[i]](1, etta) * tetta_;
				}
				F[indexes[i]] += g * jacobian;
			}
		}
		//Нижнее ребро
		if (grid.ku[1][number].edges[2] == 1)
		{
			int indexes[3] = { element.dof[0], element.dof[1], element.dof[2] };
			int id[3] = { 0, 1, 2 };

			jacobian = hx / 2;

			for (int i = 0; i < 3; i++)
			{
				g = 0;
				for (int k = 0; k < 5; k++)
				{
					double ksi = 0.5 + 0.5 * gaussPoints1[k];
					double x_ = x1 + ksi*hx;
					double tetta_ = tests.Tetta(grid.ku[1][number].formNumber[2], x_, grid.nodes[element.dof[id[i]]].y);
					g += gaussWeights1[k] * phi[id[i]](ksi, 0) * tetta_;
				}
				F[indexes[i]] += g * jacobian;
			}
		}
		//Верхнее ребро
		if (grid.ku[1][number].edges[3] == 1)
		{
			int indexes[3] = { element.dof[6], element.dof[7], element.dof[8] };
			int id[3] = { 6, 7, 8 };

			jacobian = hx / 2;

			for (int i = 0; i < 3; i++)
			{
				g = 0;
				for (int k = 0; k < 5; k++)
				{
					double ksi = 0.5 + 0.5 * gaussPoints1[k];
					double x_ = x1 + ksi*hx;
					double tetta_ = tests.Tetta(grid.ku[1][number].formNumber[3], x_, grid.nodes[element.dof[id[i]]].y);
					g += gaussWeights1[k] * phi[id[i]](ksi, 1) * tetta_;
				}
				F[indexes[i]] += g * jacobian;
			}
		}
	}

	//Учёт третьего краевого условия
	void SLAE::CalculateBoundaries3(int number)
	{
		Element element = grid.elements[grid.ku[2][number].elem];
		double g, jacobian, u_betta, betta;
		double x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x;
		double y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y;
		double hx = x3 - x1;
		double hy = y3 - y1;

		//Левое ребро
		if (grid.ku[2][number].edges[0] == 1)
		{
			int indexes[3] = { element.dof[0], element.dof[3], element.dof[6] };
			int id[3] = { 0, 3, 6 };

			jacobian = hy / 2;

			for (int i = 0; i < 3; i++)
			{
				g = 0;
				for (int k = 0; k < 5; k++)
				{
					double etta = 0.5 + 0.5 * gaussPoints1[k];
					double y_ = y1 + etta*hy;
					betta = tests.Betta(grid.ku[2][number].formNumber[0], grid.nodes[element.dof[id[i]]].x, y_);
					u_betta = tests.Ubetta(grid.ku[2][number].formNumber[0], grid.nodes[element.dof[id[i]]].x, y_);
					g += gaussWeights1[k] * phi[id[i]](0, etta)* u_betta * betta;
				}
				F[indexes[i]] += g * jacobian;
			}

			for (int i = 0; i < 9; i++)
			{
				for (int j = 0; j < 9; j++)
				{
					g = 0;
					for (int k = 0; k < 5; k++)
					{
						double etta = 0.5 + 0.5 * gaussPoints1[k];
						double y_ = y1 + etta*hy;
						betta = tests.Betta(grid.ku[2][number].formNumber[0], grid.nodes[element.dof[id[i]]].x, y_);
						g += gaussWeights1[k] * phi[i](0, etta) * phi[j](0, etta);
					}
					AddElementToGlobalMatrix(A,element.dof[i], element.dof[j], g * jacobian * betta);
				}
			}
		}
		//Правое ребро
		if (grid.ku[2][number].edges[1] == 1)
		{

			int indexes[3] = { element.dof[2], element.dof[5], element.dof[8] };
			int id[3] = { 2, 5, 8 };

			jacobian = hy / 2;

			for (int i = 0; i < 3; i++)
			{
				for (int k = 0; k < 5; k++)
				{
					double etta = 0.5 + 0.5 * gaussPoints1[k];
					double y_ = y1 + etta*hy;
					betta = tests.Betta(grid.ku[2][number].formNumber[1], grid.nodes[element.dof[id[i]]].x, y_);
					u_betta = tests.Ubetta(grid.ku[2][number].formNumber[1], grid.nodes[element.dof[id[i]]].x, y_);
					g += gaussWeights1[k] * phi[id[i]](1, etta) * u_betta * betta;
				}
				F[indexes[i]] += g * jacobian;
			}

			for (int i = 0; i < 9; i++)
			{
				for (int j = 0; j < 9; j++)
				{
					g = 0;
					for (int k = 0; k < 5; k++)
					{
						double etta = 0.5 + 0.5 * gaussPoints1[k];
						double y_ = y1 + etta*hy;
						betta = tests.Betta(grid.ku[2][number].formNumber[1], grid.nodes[element.dof[id[i]]].x, y_);
						g += gaussWeights1[k] * phi[i](1, etta) * phi[j](1, etta) * betta;
					}
					AddElementToGlobalMatrix(A,element.dof[i], element.dof[j], g * jacobian);
				}
			}
		}
		//Нижнее ребро
		if (grid.ku[2][number].edges[2] == 1)
		{
			int indexes[3] = { element.dof[0], element.dof[1], element.dof[2] };
			int id[3] = { 0, 1, 2 };

			jacobian = hx / 2;

			for (int i = 0; i < 3; i++)
			{
				g = 0;
				for (int k = 0; k < 5; k++)
				{
					double ksi = 0.5 + 0.5 * gaussPoints1[k];
					double x_ = x1 + ksi*hx;
					betta = tests.Betta(grid.ku[2][number].formNumber[2], x_, grid.nodes[element.dof[id[i]]].y);
					u_betta = tests.Ubetta(grid.ku[2][number].formNumber[2], x_, grid.nodes[element.dof[id[i]]].y);
					g += gaussWeights1[k] * phi[id[i]](ksi, 0) * u_betta * betta;
				}
				F[indexes[i]] += g * jacobian;
			}

			for (int i = 0; i < 9; i++)
			{
				for (int j = 0; j < 9; j++)
				{
					g = 0;
					for (int k = 0; k < 5; k++)
					{
						double ksi = 0.5 + 0.5 * gaussPoints1[k];
						double x_ = x1 + ksi*hx;
						betta = tests.Betta(grid.ku[2][number].formNumber[2], x_, grid.nodes[element.dof[id[i]]].y);
						g += gaussWeights1[k] * phi[i](ksi, 0) * phi[j](ksi, 0)* betta;
					}
					AddElementToGlobalMatrix(A, element.dof[i], element.dof[j], g * jacobian);
				}
			}
		}
		//Верхнее ребро
		if (grid.ku[2][number].edges[3] == 1)
		{
			int indexes[3] = { element.dof[6], element.dof[7], element.dof[8] };
			int id[3] = { 6, 7, 8 };

			jacobian = hx / 2;

			for (int i = 0; i < 3; i++)
			{
				g = 0;
				for (int k = 0; k < 5; k++)
				{
					double ksi = 0.5 + 0.5 * gaussPoints1[k];
					double x_ = x1 + ksi*hx;
					betta = tests.Betta(grid.ku[2][number].formNumber[3], x_, grid.nodes[element.dof[id[i]]].y);
					u_betta = tests.Ubetta(grid.ku[2][number].formNumber[3], x_, grid.nodes[element.dof[id[i]]].y);
					g += gaussWeights1[k] * phi[id[i]](ksi, 1) * u_betta * betta;
				}
				F[indexes[i]] += g * jacobian;
			}

			for (int i = 0; i < 9; i++)
			{
				for (int j = 0; j < 9; j++)
				{
					g = 0;
					for (int k = 0; k < 5; k++)
					{
						double ksi = 0.5 + 0.5 * gaussPoints1[k];
						double x_ = x1 + ksi*hx;
						betta = tests.Betta(grid.ku[2][number].formNumber[3], x_, grid.nodes[element.dof[id[i]]].y);
						g += gaussWeights1[k] * phi[i](ksi, 1) * phi[j](ksi, 1) * betta;
					}
					AddElementToGlobalMatrix(A, element.dof[i], element.dof[j], g * jacobian);
				}
			}
		}
	}

	//Вычисление нормы вектора
	double SLAE::Norm(const vector<double> &x)
	{
		double norm = 0;
		int size = x.size();

		for (int i = 0; i < size; i++)
			norm += x[i] * x[i];

		return sqrt(norm);
	}

	//Скалярное произведение векторов
	double SLAE::Scalar(const vector<double> &x, const vector<double> &y)
	{
		double sum = 0;
		int size = x.size();
		for (int i = 0; i < size; i++)
			sum += x[i] * y[i];

		return sum;
	}

	//LU-факторизация
	void SLAE::LU()
	{
		int i, i0, j0, iend, num, ki, kj, jend;
		double suml, sumu, sumdg;

		L.resize(A.ggl.size());
		L = A.ggl;
		U.resize(A.ggu.size());
		U = A.ggu;
		D.resize(A.di.size());
		D = A.di;

		for (i = 0; i < n; i++)
		{
			i0 = A.ig[i];
			iend = A.ig[i + 1];

			for (num = i0, sumdg = 0; num < iend; num++)
			{
				j0 = A.ig[A.jg[num]];
				jend = A.ig[A.jg[num] + 1];
				ki = i0;
				kj = j0;
				//для num учитываются все предыдущие элементы
				for (suml = 0, sumu = 0, ki = i0; ki < num; ki++)
				{
					for (int m = kj; m < jend; m++)
						//ищем соответствующие ненулевые элементы для умножения
						if (A.jg[ki] == A.jg[m])
						{
							suml += L[ki] * U[m];
							sumu += L[m] * U[ki];
						}
				}
				L[num] -= suml;
				U[num] = (U[num] - sumu) / D[A.jg[num]];
				//умножаются симметричные элементы
				sumdg += L[num] * U[num];
			}
			D[i] -= sumdg;
		}
	}
	
	void SLAE::LYF(const vector <double> &C, vector <double> &yl)
	{
		int i, i0, iend; //i0-адрес начала строки, iend-адрес конца строки
		double sum;
		for (i = 0; i < n; i++)
		{
			i0 = A.ig[i]; iend = A.ig[i + 1];

			yl[i] = C[i];

			for (i0, sum = 0; i0 < iend; i0++)
				yl[i] -= yl[A.jg[i0]] * L[i0];
			yl[i] /= D[i];
		}
	}

	void SLAE::UXY(const vector <double> &C, vector <double> &yu)
	{
		int i, i0, iend;

		for (i = 0; i < n; i++)
			yu[i] = 0.0;

		for (i = n - 1; i >= 0; i--)//проход по столбцам с конца
		{
			yu[i] += C[i];

			i0 = A.ig[i]; iend = A.ig[i + 1]; iend--;

			for (; iend >= i0; iend--)//идём по столбцу с конца
				yu[A.jg[iend]] -= yu[i] * U[iend];
		}
	}

	double SLAE::Rel_Discrepancy()
	{
		double dis1, dis2;
		dis1 = Scalar(r, r);
		dis2 = Scalar(F, F);
		double dis = dis1 / dis2;
		return sqrt(dis);
	}

	void SLAE::LULOS()
	{
		double a, b, pp, dis, rr;
		int i,k;
		vector <double> Ax(n), C(n), p(n);

		LU();
		//Ax0
		A.MultiplyAx(u, Ax);
		//f-Ax0
		for (i = 0; i < n; i++)
			r[i] = F[i] - Ax[i];

		//r0=L^(-1)(f-Ax0)
		LYF(r, r);

		//z0=U^(-1)r0->r0=Uz0
		UXY(r, z);

		//p0=L^(-1)Az0
		A.MultiplyAx(z, Ax);//Az0
		LYF(Ax, p);

		rr = Scalar(r, r);
		dis = Scalar(r, r) / rr;
		dis = sqrt(dis);
		k = 0;

		for (k = 1; dis > eps && k <= maxiter; k++)
		{
			//Аk
			pp = Scalar(p, p);
			a = Scalar(p, r) / pp;

			//Xk, Rk
			for (i = 0; i < n; i++)
			{
				u[i] = u[i] + a*z[i];
				r[i] = r[i] - a*p[i];
			}

			//UY=rk->Y=U^(-1)rk
			UXY(r, C);
			//AU^(-1)rk=Ax
			A.MultiplyAx(C, Ax);
			//L^(-1)AU^(-1)rk=Y2->L^(-1)B=Y2->LY2=B->Y2=L^(-1)AU^(-1)rk
			LYF(Ax, Ax);

			//bk
			b = -Scalar(p, Ax) / pp;

			//zk=U^(-1)rk+bkz[k-1]
			//pk
			for (i = 0; i < n; i++)
			{
				z[i] = C[i] + b * z[i];
				p[i] = Ax[i] + b * p[i];
			}
			dis = Scalar(r, r) / rr;
			dis = sqrt(dis);
		}

		dis = Rel_Discrepancy();
	}

	double SLAE::StopIteration()
	{
		vector <double> b(n), z(n);
		//генерируем СЛАУ с текущим решением
		GenerateSLAE();
		//умножаем матрицу на текущее решение
		A.MultiplyAx(u, z);
		//и находим абсолютную невязку
		for (int i = 0; i < n; i++)
			b[i] = z[i] - F[i];

		return Norm(b) / Norm(F);
	}
	
	void SLAE::TSolve()
	{
		FILE *fo;
		fopen_s(&fo, "result.txt", "w");
		int k;
		double exit_condition;
		t = 1;

		vector <double> b(n), z(n);
		u = { 0,1,8,0,1,8,0,1,8 };

		//for (auto &i : u) i = 1;
		//цикл по времени
		while (t < maxtime)
		{
			printf("Calculating......([t] = %f)\n", t);
			//цикл по нелинейности
			for (k = 1; k < maxiter; k++)
			{
				//сохраняем предыдущее решение
				u_ = u;
				//находим новое решение
				GenerateSLAE();
				//умножаем матрицу на текущее решение
				A.MultiplyAx(u, z);
				//и находим абсолютную невязку
				for (int i = 0; i < n; i++)
					b[i] = z[i] - F[i];
				LULOS();
				//проверяем условие останова 
				//+ генерация СЛАУ с текущим решением
				exit_condition = StopIteration();
				//если оно не выполнено, то, используя параметр
				//релаксации, находим новое решение
				if (exit_condition > eps)
				{
					for (int i = 0; i < n; i++)
						u[i] = w*u[i] + (1 - w)*u_[i];
				}
				//иначе выходим из цикла
				else break;
			}

			for (int i = 0; i < n; i++)
				fprintf(fo, "%.14lf\n", u[i]);
			fprintf(fo, "%f\t", t);
			fprintf(fo, "%d\t", k);
			fprintf(fo, "%e\t", exit_condition);
			fprintf(fo, "\n\n");

			printf("Vector:");
			for (int i = 0; i < n; i++)
				printf("\t\t%.14lf\n", u[i]);
			printf("%d\t", k);
			printf("%e\t", exit_condition);
			printf("\n\n");

			t += ht;
		}

		fclose(fo);
	}

}