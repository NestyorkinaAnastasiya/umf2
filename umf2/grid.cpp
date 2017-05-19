/*grid.cpp*/
#include "grid.h"

namespace grid
{
	void Area::Input(FILE *fo)
	{
		fscanf_s(fo, "%d", &leftX);
		fscanf_s(fo, "%d", &rightX);
		fscanf_s(fo, "%d", &lowY);
		fscanf_s(fo, "%d", &upY);

		fscanf_s(fo, "%d", &ku[0]);
		fscanf_s(fo, "%d", &ku[1]);
		fscanf_s(fo, "%d", &ku[2]);
		fscanf_s(fo, "%d", &ku[3]);

		fscanf_s(fo, "%d", &kuForm[0]);
		fscanf_s(fo, "%d", &kuForm[1]);
		fscanf_s(fo, "%d", &kuForm[2]);
		fscanf_s(fo, "%d", &kuForm[3]);
	}

	AreasLines::AreasLines()
	{
		int n;
		double tmp;
		FILE *fo;
		fopen_s(&fo,"AreasLines.txt", "r");
		fscanf_s(fo, "%d", &n);
		x.reserve(n);
		for (int i = 0; i < n; i++)
		{
			fscanf_s(fo, "%lf", &tmp);
			x.push_back(tmp);
		}

		fscanf_s(fo, "%d", &n);
		y.reserve(n);
		for (int j = 0; j < n; j++)
		{
			fscanf_s(fo, "%lf", &tmp);
			y.push_back(tmp);
		}
		fclose(fo);
	}

	//Ïîèñê ãëîáàëüíîãî íîìåðà á.ô. dofsNumber â ýëåìåíòå
	bool Element::SearchDof(int dofsNumber)
	{
		for (int i = 0; i < 9; i++)
			if (dofsNumber == dof[i]) return true;

		return false;
	}

	Grid::Grid()
	{
		int n;
		FILE *fo;
		fopen_s(&fo,"Areas.txt", "r");
		fscanf_s(fo, "%d", &n);
		areas.resize(n);
		for (int i = 0; i < n; i++)
			areas[i].Input(fo);
		fclose(fo);
	}

	Grid::~Grid() {}

	//Ãåíåðàöèÿ êîîðäèíàòû ñ ó÷¸òîì ðàçáèåíèÿ íà âñåõ èíòåðâàëàõ
	void Grid::PartitionÑoordinate(vector <double> &x, vector <double> areasLines,
		vector <double> coefficient, vector <int> nIntervals)
	{
		int count;
		//äëèíà èíòåðâàëà
		double l;
		//øàã
		double h;
		//÷èñëî èíòåðâàëîâ
		int nLines = areasLines.size();

		count = 0;
		for (int i = 0; i < nLines - 1; i++)
		{
			x.push_back(areasLines[i]);
			count++;

			//äëèíà èíòåðâàëà
			l = abs(areasLines[i + 1] - areasLines[i]);

			//ðàññ÷èòûâàåì ïåðâûé øàã
			//ðàâíîìåðíàÿ
			if (abs(1.0 - coefficient[i]) < 1E-14)
				h = l / nIntervals[i];
			else //ñãóùàåì âïðàâî
				if (coefficient[i] < 1)
				{
					h = l * (1.0 - coefficient[i]);
					h /= 1.0 - pow(coefficient[i], nIntervals[i]);
				}
				else //ñãóùàåì âëåâî
				{
					h = l * (coefficient[i] - 1.0);
					h /= coefficient[i] * (pow(coefficient[i], nIntervals[i] - 1.0) - 1.0);
					h += coefficient[i] - 1.0;
				}

			//ïîëó÷àåì ñåòêó âíóòðè èíòåðâàëà
			for (int j = 1; j < nIntervals[i]; j++)
			{
				if (j != 1) h *= coefficient[i];
				x.push_back(x[count - 1] + h);
				count++;
			}
		}
		x.push_back(areasLines[nLines - 1]);
	}

	//Äîáàâëåíèå óçëà
	void Grid::PushNode(double x, double y)
	{
		Point p;
		p.y = y;
		p.x = x;
		//êëàä¸ì óçåë â êîíåö ìàññèâà óçëîâ
		nodes.push_back(p);
	}

	//Ïîñòðîåíèå ñåòêè
	void Grid::BuildGrid()
	{
		//äëÿ iîãî èíòåðâàëà è êàæäîé êîîðäèíàòû
		//êîýôôèöèåíò ðàçðÿäêè
		vector <double> xCoefficient;
		vector <double> yCoefficient;
		//÷èñëî ïîäèíòåðâàëîâ
		vector <int> xIntervals;
		vector <int> yIntervals;

		//êîëè÷åñòâî êîîðäèíàòíûõ ëèíèé ïî x è y
		int xLines = areasLines.x.size();
		int yLines = areasLines.y.size();

		//ãåîìåòðè÷åñêèå ëèíèè ðàçáèåíèÿ ïî õ è y
		vector <double> xi;
		vector <double> yj;

		//âðåìåííàÿ ïåðåìåííàÿ äëÿ ñ÷èòûâàíèÿ
		int tmp1;
		//âðåìåííàÿ ïåðåìåííàÿ äëÿ ñ÷èòûâàíèÿ
		double tmp2;

		//÷èñëî óçëîâ è ÷èñëî ýëåìåíòîâ
		int nNodes, nElements;

		xIntervals.reserve(xLines - 1); xCoefficient.reserve(xLines - 1);
		yIntervals.reserve(yLines - 1); yCoefficient.reserve(yLines - 1);

		FILE *fo;
		fopen_s(&fo,"Intervals.txt", "r");
		//ââîä êîëè÷åñòâà èíòåðâàëîâ ïî õ è ó è êîýôôèöèåíòîâ ðàçðÿäêè
		for (int i = 0; i < xLines - 1; i++)
		{
			fscanf_s(fo, "%d", &tmp1);
			xIntervals.push_back(tmp1);
			fscanf_s(fo, "%lf", &tmp2);
			xCoefficient.push_back(tmp2);
		}

		for (int j = 0; j < yLines - 1; j++)
		{
			fscanf_s(fo, "%d", &tmp1);
			yIntervals.push_back(tmp1);
			fscanf_s(fo, "%lf", &tmp2);
			yCoefficient.push_back(tmp2);
		}
		fclose(fo);
		nx = 0;
		for (int i = 0; i < xLines - 1; i++)
			//îáùåå êîëè÷åñòâî èíòåðâàëîâ ïî õ
			nx += xIntervals[i];
		nx++;

		ny = 0;
		for (int j = 0; j < yLines - 1; j++)
			//îáùåå êîëè÷åñòâî èíòåðâàëîâ ïî ó
			ny += yIntervals[j];
		ny++;

		xi.reserve(nx); yj.reserve(ny);

		//ïîñòðîåíèå ñåòîê ïî õ è ó
		PartitionÑoordinate(xi, areasLines.x, xCoefficient, xIntervals);
		PartitionÑoordinate(yj, areasLines.y, yCoefficient, yIntervals);

		xIntervals.clear(); xCoefficient.clear();
		yIntervals.clear(); yCoefficient.clear();


		nElements = (nx - 1) * (ny - 1);
		nNodes = xi.size() * yj.size();

		elements.reserve(nElements);
		nodes.reserve(nNodes);

		//çàïîëíÿåì ñïèñîê óçëîâ
		for (int j = 0; j < yj.size(); j++)
			for (int i = 0; i < xi.size(); i++)
				PushNode(xi[i], yj[j]);
	}

	//Ïîëó÷åíèå ãëîáàëüíîãî íîìåðà
	int Grid::GetGlobalNumber(int elementNumber, int localNumber)
	{
		//÷èñëî èíòåðâàëîâ ïî x
		int nxint = nx - 1;
		//íîìåð ãîðèçîíòàëüíîé ëèíèè, êîòîðàÿ ÿâëÿåòñÿ íèæíèì ðåáðîì ýëåìåíòà
		int nGorLine = elementNumber / nxint;
		//íà÷àëüíûé íîìåð óçíà íà ðåáðå
		int nodeOnLine = nGorLine * (nxint + 1);
		//îòñòóï îò íà÷àëüíîãî óçëà, ÷òîáû ïîëó÷èòü òåêóùèé íîìåð
		//ëåâîãî íèæíåãî óçëà ýëåìåíòà;
		int offsetX = elementNumber % nxint;
		//íîìåð íèæíåãî ëåâîãî óçëà
		int nodeElem = nodeOnLine + offsetX;


		if (localNumber == 0)
			return nodeElem;

		if (localNumber == 1)
			return nodeElem + 1;

		if (localNumber == 2)
			return nodeElem + (nxint + 1);

		if (localNumber == 3)
			return nodeElem + (nxint + 1) + 1;

		return -1;

	}

	//Ïîëó÷åíèå ãëîáàëüíûõ áàçèñíûõ íîìåðîâ
	int Grid::GetGlobalFuncNumber(int elNum, int localFuncNum)
	{
		//÷èñëî èíòåðâàëîâ ïî x
		int nxint = nx - 1;
		//íîìåð ãîðèçîíòàëüíîé ëèíèè, êîòîðàÿ ÿâëÿåòñÿ íèæíèì ðåáðîì ýëåìåíòà
		int nGorLine = elNum / nxint;
		//Êîëè÷åñòâî óçëîâ íà îäíîé ãîðèçîíòàëüíîé ëèíèè
		int nGorNodes = 2 * nx - 1;
		//îòñòóï îò íà÷àëüíîãî óçëà, ÷òîáû ïîëó÷èòü òåêóùèé íîìåð
		//ëåâîãî íèæíåãî óçëà ýëåìåíòà;
		int offsetX = elNum % nxint;

		//ãëîáàëüíûé íîìåð íèæíåãî ëåâîãî óçëà
		int res = 2 * nGorLine * nGorNodes + 2 * offsetX;
		//â íóìåðàöèè îò åäèíèöû 1-3
		if (localFuncNum <= 2)
			return res + localFuncNum;
		else //... 4-6
			if (localFuncNum <= 5)
				return res + nGorNodes + (localFuncNum - 3);
			else //... 7-9
				return res + 2 * nGorNodes + (localFuncNum - 6);
	}

	//Íàõîæäåíèå íîìåðà îáëàñòè, â êîòîðîé ëåæèò çàäàííàÿ òî÷êà
	int Grid::FindArea(double x, double y)
	{
		bool xInterval, yInterval;
		int size = areas.size();
		for (int i = 0; i < size; i++)
		{
			xInterval = x > areasLines.x[areas[i].leftX] && x < areasLines.x[areas[i].rightX];
			yInterval = y > areasLines.y[areas[i].lowY] && y < areasLines.y[areas[i].upY];

			//åñëè òî÷êà ïîïàäàåò â ióþ îáëàñòü, âîçâðàùàåì íîìåð ýòîé îáëàñòè
			if (xInterval && yInterval) return i;
		}
	}

	//Íàõîæäåíèå ñîñåäíèõ êîíå÷íûõ ýëåìåíòîâ
	void Grid::FindNeighbors(int elementNumber)
	{
		int count = 0;
		//êîëè÷åñòâî êîíå÷íûõ ýëåìåíòîâ
		int nElements = elements.size();
		Element element;
		element = elements[elementNumber];

		bool tmp[4] = { false, false, false, false };

		for (int i = 0; count < 4 && i < nElements; i++)
		{
			if (i != elementNumber) 
			{
				//ñîñåä ïî ëåâîìó ðåáðó
				if (nodes[element.nodes[0]] == nodes[elements[i].nodes[1]] && nodes[element.nodes[2]] == nodes[elements[i].nodes[3]])
				{
					element.neighbors[0] = i;
					tmp[0] = true;
					count++;
				}
				//ñîñåä ïî ïðàâîìó ðåáðó
				if (nodes[element.nodes[1]] == nodes[elements[i].nodes[0]] && nodes[element.nodes[3]] == nodes[elements[i].nodes[2]])
				{
					element.neighbors[1] = i;
					tmp[1] = true;
					count++;
				}
				//ñîñåä ïî íèæíåìó ðåáðó
				if (nodes[element.nodes[0]] == nodes[elements[i].nodes[2]] && nodes[element.nodes[1]] == nodes[elements[i].nodes[3]])
				{
					element.neighbors[2] = i;
					tmp[2] = true;
					count++;
				}
				//ñîñåä ïî âåðõíåìó ðåáðó
				if (nodes[element.nodes[2]] == nodes[elements[i].nodes[0]] && nodes[element.nodes[3]] == nodes[elements[i].nodes[1]])
				{
					element.neighbors[3] = i;
					tmp[3] = true;
					count++;
				}
			}
		}

		for (int i = 0; i < 4; i++)
			//îòñóòñòâèå ñîñåäà íà ñîîòâåòñòâóþùåé ñòîðîíå
			if (tmp[i] == false) element.neighbors[i] = -1;

		elements[elementNumber] = element;
	}

	//Âû÷èñëåíèå êîíå÷íûõ ýëåìåíòîâ
	void Grid::ComputeElements()
	{
		Element tmpEl;
		double x, y;

		//âû÷èñëÿåì ãëîáàëüíûå íîìåðà óçëîâ è áàçèñíûõ ôóíêöèé
		//êàæäîãî êý è çàïîëíÿåì ñïèñîê êý
		for (int i = 0; i < elements.capacity(); i++)
		{
			for (int j = 0; j < 4; j++)
				tmpEl.nodes[j] = GetGlobalNumber(i, j);

			for (int j = 0; j < 9; j++)
				tmpEl.dof[j] = GetGlobalFuncNumber(i, j);

			elements.push_back(tmpEl);
		}

		//íàõîäèì íîìåð ïîäîáëàñòè äëÿ êàæäîãî êý
		for (int i = 0; i < elements.size(); i++)
		{

			//áåð¸ì öåíòðàëüíûé óçåë
			x = (nodes[elements[i].nodes[0]].x + nodes[elements[i].nodes[1]].x) / 2;
			y = (nodes[elements[i].nodes[0]].y + nodes[elements[i].nodes[2]].y) / 2;

			elements[i].numberOfArea = FindArea(x, y);
		}

		//íàõîäèì ñîñåäíèå êý
		for (int i = 0; i < elements.size(); i++)
			FindNeighbors(i);
	}

	//Ôîðìèðîâàíèå ìàññèâîâ êðàåâûõ óñëîâèé
	void Grid::FormKU()
	{
		int size = elements.size();
		BoundaryCondition tmp;

		//äëÿ êàæäîãî ýëåìåíòà íàõîäèì,êàêèå êó íà åãî ãðàíèöàõ
		for (int i = 0; i < size; i++)
		{
			bool left[3] = { false, false, false }, right[3] = { false, false, false },
				low[3] = { false, false, false }, up[3] = { false, false, false };
			if (elements[i].neighbors[0] == -1)
			{
				if (areas[elements[i].numberOfArea].ku[0] == 1) left[0] = true;
				else
					if (areas[elements[i].numberOfArea].ku[0] == 2) left[1] = true;
					else
						if (areas[elements[i].numberOfArea].ku[0] == 3) left[2] = true;
			}
			if (elements[i].neighbors[1] == -1)
			{
				if (areas[elements[i].numberOfArea].ku[1] == 1) right[0] = true;
				else
					if (areas[elements[i].numberOfArea].ku[1] == 2) right[1] = true;
					else
						if (areas[elements[i].numberOfArea].ku[1] == 3) right[2] = true;
			}
			if (elements[i].neighbors[2] == -1)
			{
				if (areas[elements[i].numberOfArea].ku[2] == 1) low[0] = true;
				else
					if (areas[elements[i].numberOfArea].ku[2] == 2) low[1] = true;
					else
						if (areas[elements[i].numberOfArea].ku[2] == 3) low[2] = true;
			}
			if (elements[i].neighbors[3] == -1)
			{
				if (areas[elements[i].numberOfArea].ku[3] == 1) up[0] = true;
				else
					if (areas[elements[i].numberOfArea].ku[3] == 2) up[1] = true;
					else
						if (areas[elements[i].numberOfArea].ku[3] == 3) up[2] = true;
			}

			for (int j = 0; j < 3; j++)
				if (left[j] || right[j] || low[j] || up[j])
				{
					tmp.elem = i;
					if (left[j])
					{
						tmp.edges[0] = 1;
						tmp.formNumber[0] = areas[elements[i].numberOfArea].kuForm[0];
					}
					else
					{
						tmp.edges[0] = 0;
						tmp.formNumber[0] = -1;
					}
					if (right[j])
					{
						tmp.edges[1] = 1;
						tmp.formNumber[1] = areas[elements[i].numberOfArea].kuForm[1];
					}
					else
					{
						tmp.edges[1] = 0;
						tmp.formNumber[1] = -1;
					}
					if (low[j])
					{
						tmp.edges[2] = 1;
						tmp.formNumber[2] = areas[elements[i].numberOfArea].kuForm[2];
					}
					else
					{
						tmp.edges[2] = 0;
						tmp.formNumber[2] = -1;
					}
					if (up[j])
					{
						tmp.edges[3] = 1;
						tmp.formNumber[3] = areas[elements[i].numberOfArea].kuForm[3];
					}
					else
					{
						tmp.edges[3] = 0;
						tmp.formNumber[3] = -1;
					}

					ku[j].push_back(tmp);
				}
		}
	}

	//Îáðàçîâàíèå ìíîæåñòâà êý
	void Grid::DoPartition()
	{
		//Ïîñòðîåíèå ñåòêè
		BuildGrid();
		//Âû÷èñëåíèå êîíå÷íûõ ýëåìåíòîâ
		ComputeElements();
		//Ôîðìèðîâàíèå ìàññèâîâ êðàåâûõ óñëîâèé
		FormKU();
	}

	//Ôîðìèðîâàíèå ñïèñêà ýëåìåíòîâ, ñîäåðæàùèõ ãëîáàëüíûé íîìåð á.ô.
	//ðàâíûé dofsNumber
	void Grid::SearchElements(int dofsNumber, vector <int> &elList)
	{
		int count;
		int size = elements.size();
		elList.reserve(4);

		count = 0;
		for (int i = 0; i < size && count < 4; i++)
		{
			if (elements[i].SearchDof(dofsNumber))
			{
				count++;
				elList.push_back(i);
			}
		}
	}
}
