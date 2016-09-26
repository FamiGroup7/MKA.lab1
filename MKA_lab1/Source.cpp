#include <stdio.h>
#include <conio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "Point.cpp"
#include "Dividing.cpp"
#include "nvtr.cpp"
using namespace std;

int Nx, Ny;
int xTotalCells = 0, yTotalCells = 0;
Point**lines;
Dividing*xDivide;
Dividing*yDivide;
ofstream cells("cells.txt");
ofstream out("out.txt");
vector<Point> xy;
vector<nvtr> KE;

void output() {
	ofstream xyFile("xy.txt");
	ofstream nvtrFile("nvtr.txt");

	xyFile << xy.size() << endl;
	for (size_t i = 0; i < xy.size(); i++)
	{
		xyFile << xy[i].x << " " << xy[i].y << endl;
	}

	nvtrFile << KE.size() << endl;
	for (size_t i = 0; i < KE.size(); i++)
	{
		for (size_t k = 0; k < 4; k++)
		{
			nvtrFile << KE[i].uzel[k] << " ";
		}
		nvtrFile << KE[i].numberField << endl;
	}
}
void inputLines() {
	ifstream wFile("W.txt");
	wFile >> Nx >> Ny;
	xDivide = new Dividing[Nx - 1];
	yDivide = new Dividing[Ny - 1];

	lines = new Point*[Nx];
	for (size_t i = 0; i < Nx; i++)
	{
		lines[i] = new Point[Ny];
	}
	for (size_t j = 0; j < Ny; j++)
	{
		for (size_t i = 0; i < Nx; i++)
		{
			wFile >> lines[i][j].x >> lines[i][j].y;
		}
	}
	wFile.close();

	ifstream divFile("Dividing.txt");
	for (size_t i = 0; i < Nx - 1; i++)
	{
		divFile >> xDivide[i].nIntervals >> xDivide[i].koefRazr;
		xTotalCells += xDivide[i].nIntervals;
	}
	xTotalCells++;
	for (size_t i = 0; i < Ny - 1; i++)
	{
		divFile >> yDivide[i].nIntervals >> yDivide[i].koefRazr;
		yTotalCells += yDivide[i].nIntervals;
	}
	yTotalCells++;

}

double findH(double length, double alfa, int n) {
	if (alfa < 0)
		alfa = -1 / alfa;
	return length*(1 - alfa) / (1 - pow(alfa, n));
}
double relation(double alfa, int n) {
	if (alfa < 0)
		alfa = -1 / alfa;
	return (1 - alfa) / (1 - pow(alfa, n));
}
double relation(double alfa, double alfa1) {
	if (alfa < 0)
		alfa = -1 / alfa;
	return alfa1*alfa;
}

bool Inside(Point source,  vector<Point> poligons)
{
	int size = poligons.size();

	Point p, last = poligons.back();
	int intersect_counter = int();
	for (int i = 0; i<size; ++i)
	{
		Point p = poligons[i];
		if ((p.y >= source.y) ^ (last.y >= source.y))
		{
			float t = (source.y - last.y) / (p.y - last.y);
			float x = last.x + t * (p.x - last.x);
			if (x>=source.x)
				++intersect_counter;
		}
		last = p;
	}
	return (intersect_counter & 1);
}

double length(Point p1, Point p2)
{
	return sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
}
//проверка принадлежности точки внутренней границе
bool belongToLine(Point source, vector<Point> poligons)
{
	Point currentPoint = poligons[0];
	for (size_t i = 1; i < poligons.size(); i++)
	{
		if (length(currentPoint, poligons[i]) == length(source, currentPoint) + length(source, poligons[i]))
			return true;
		currentPoint = poligons[i];
	}
	return false;
}
//индекс точки в массиве xy
int indexXY(Point source)
{
	for (size_t i = 0; i < xy.size(); i++)
	{
		if (source.x == xy[i].x && source.y == xy[i].y)
			return i;
	}
	return -1;
}

void generateArrayOfCells() {
	int countY = 0, countX = 0;
	Point**arrayCells = new Point*[xTotalCells];
	for (size_t i = 0; i < xTotalCells; i++)
	{
		arrayCells[i] = new Point[yTotalCells];
	}

	for (size_t j = 0; j < Ny; j++)
	{
		int startY = countY;
		int startX;
		countX = 0;
		for (size_t i = 0; i < Nx - 1; i++)
		{
			startX = countX;
			double lengthX = lines[i + 1][j].x - lines[i][j].x;
			double lengthY = lines[i + 1][j].y - lines[i][j].y;
			arrayCells[countX][countY] = lines[i][j];
			countX++;
			if(xDivide[i].koefRazr==1)
				for (size_t indX = 1; indX < xDivide[i].nIntervals; indX++, countX++)
				{
					arrayCells[countX][countY].x = lines[i][j].x + lengthX *indX / xDivide[i].nIntervals;
					arrayCells[countX][countY].y = lines[i][j].y + lengthY *indX / xDivide[i].nIntervals;
				}
			else
			{
				double otnosh = relation(xDivide[i].koefRazr, xDivide[i].nIntervals);
				double alfa = relation(xDivide[i].koefRazr, 1);
				for (size_t indX = 1; indX < xDivide[i].nIntervals; indX++, countX++)
				{
					arrayCells[countX][countY].x = lines[i][j].x + otnosh*lengthX;
					arrayCells[countX][countY].y = lines[i][j].y + otnosh*lengthY;
					otnosh *= alfa;
				}
			}
		}
		arrayCells[xTotalCells - 1][startY] = lines[Nx - 1][j];

		countY = startY + yDivide[j].nIntervals;
		
	}

	countX = countY = 0;
	for (countX = 0; countX < xTotalCells; countX++)
	{
		countY = 0;
		int startY = countY;
		for (size_t j = 0; j < Ny - 1; j++)
		{
			startY = countY;
			countY += 1;
			//double lengthX = lines[countX][j + 1].x - lines[countX][j].x;
			//double lengthY = lines[countX][j + 1].y - lines[countX][j].y;
			double lengthX = arrayCells[countX][startY + yDivide[j].nIntervals].x - arrayCells[countX][startY].x;
			double lengthY = arrayCells[countX][startY + yDivide[j].nIntervals].y - arrayCells[countX][startY].y;

			if (yDivide[j].koefRazr == 1)
				for (size_t indY = 1; indY < yDivide[j].nIntervals; indY++, countY++)
				{
					arrayCells[countX][countY].x = arrayCells[countX][startY].x + lengthX *indY / yDivide[j].nIntervals;
					arrayCells[countX][countY].y = arrayCells[countX][startY].y + lengthY *indY / yDivide[j].nIntervals;
				}
			else
			{
				double otnosh = relation(yDivide[j].koefRazr, yDivide[j].nIntervals);
				double alfa = relation(yDivide[j].koefRazr, 1);
				for (size_t indY = 1; indY < yDivide[j].nIntervals; indY++, countY++)
				{
					arrayCells[countX][countY].x = arrayCells[countX][startY].x + otnosh*lengthX;
					arrayCells[countX][countY].y = arrayCells[countX][startY].y + otnosh*lengthY;
					otnosh *= alfa;
				}
			}
		}
	}

	for (int j = yTotalCells-1; j >=0; j--)
	{
		for (int i = 0; i < xTotalCells; i++)
		{
			cells << setw(10) << arrayCells[i][j].x << " " << setw(4) << arrayCells[i][j].y;
		}
		cells << endl;
	}
	
	ifstream ignoredFieldFile("ignoredFields.txt");
	vector<Point> ignoredField;
	int nFields;
	ignoredFieldFile >> nFields;
	for (size_t i = 0; i < nFields; i++)
	{
		double x, y;
		ignoredFieldFile >> x >> y;
		ignoredField.push_back(Point(x, y));
	}
	for (size_t j = 0; j < yTotalCells; j++)
	{
		for (size_t i = 0; i < xTotalCells; i++)
		{
			if (!Inside(arrayCells[i][j], ignoredField)) {
				xy.push_back(Point(arrayCells[i][j]));
			}
			else
				if (belongToLine(arrayCells[i][j], ignoredField)) {
					xy.push_back(Point(arrayCells[i][j]));
				}
			/*	else
				cout << setw(8) << arrayCells[i][j].x << setw(8) << arrayCells[i][j].y << endl;*/
		}
	}

	for (size_t j = 0; j < yTotalCells - 1; j++)
	{
		for (size_t i = 0; i < xTotalCells - 1; i++)
		{
			nvtr temp = nvtr();
			if ((temp.uzel[0] = indexXY(arrayCells[i][j])) == -1)
				continue;
			if ((temp.uzel[1] = indexXY(arrayCells[i+1][j])) == -1)
				continue;
			if ((temp.uzel[2] = indexXY(arrayCells[i][j+1])) == -1)
				continue;
			if ((temp.uzel[3] = indexXY(arrayCells[i+1][j+1])) == -1)
				continue;
			temp.numberField = 0;
			KE.push_back(temp);
		}
	}

}
void main() {
	inputLines();
	generateArrayOfCells();
	output();
}