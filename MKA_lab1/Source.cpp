#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <conio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <set>
#include <map>
#include "Point.cpp"
#include "Dividing.cpp"
#include "nvtr.cpp"
#include "glut.h"
#define isDisplayingFunctions false
#define isDisplayingEdges false
#define isDisplayingAreas false
#define isDisplayingQuadrangles false
#define isDisplayingIsolines true

using namespace std;

int Nx, Ny;
int xTotalCells = 0, yTotalCells = 0;
Point**lines;
Dividing*xDivide;
Dividing*yDivide;
ofstream cells("cells.txt");
vector<Point> xy;
vector<nvtr> KE;
vector<Point> ignoredField;
vector<Point> allBounds;

int*ig, *jg, *igEdge, *jgEdge;
double *ggl, *ggu, *di, *b, *q;

const int COUNT_OF_ISOLINES = 10;
const int COUNT_OF_COLOR_AREAS = 5;
const int COUNT_OF_X_INNER_QUADRES = 20;
double isolinesValues[COUNT_OF_ISOLINES];
struct Map {
	double value;
	unsigned short red, green, blue;
};
Map rainbow[COUNT_OF_COLOR_AREAS];

struct Triangle {
	int nvtr[3];
	Triangle(int uzel1, int uzel2, int uzel3) {
		nvtr[0] = uzel1; nvtr[1] = uzel2; nvtr[2] = uzel3;
	}
	Triangle() {}
};

struct TriangleGeneral {
	Point nvtr[3];
	double solutions[3];
};

struct Quadrangle {
	Point nvtr[4];
	double solutions[4];
	Quadrangle(){}
};

GLint Width, Height;
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
		for (size_t k = 0; k < 9; k++)
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
		if (xDivide[i].koefRazr == -1)
			xDivide[i].koefRazr = 1;
		xTotalCells += xDivide[i].nIntervals;
	}
	xTotalCells++;
	for (size_t i = 0; i < Ny - 1; i++)
	{
		divFile >> yDivide[i].nIntervals >> yDivide[i].koefRazr;
		if (yDivide[i].koefRazr == -1)
			yDivide[i].koefRazr = 1;
		yTotalCells += yDivide[i].nIntervals;
	}
	yTotalCells++;

}

double BasicFunc1d(int num, double ksi) {
	switch (num)
	{
	case 0: return 2.0 * (ksi - 0.5) * (ksi - 1.0);
	case 1: return -4.0 * ksi * (ksi - 1.0);
	case 2: return 2.0 * ksi * (ksi - 0.5);
	default:
		cerr << "Error in Basic Function" << endl;
		system("pause");
		exit(1);
	}
}

double BasicFunc2d(int i, double ksi, double eta) {
	return BasicFunc1d(i % 3, ksi)*BasicFunc1d(i / 3, eta);
}

double difBasicFunc1d(int num, double ksi) {
	switch (num)
	{
	case 0: return 4.0 * ksi - 3;
	case 1: return 4.0 - 8.0 * ksi;
	case 2: return 4.0 * ksi - 1;
	default:
		cerr << "Error in Basic Function" << endl;
		system("pause");
		exit(1);
	}
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
double reverseKoef(double alfa, double alfa1) {
	if (alfa < 0)
		alfa = -1 / alfa;
	return alfa1*alfa;
}

bool Inside(Point source,  vector<Point> poligons)
{
	int size = poligons.size();
	if (size < 3)return false;

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

bool BelongToLine(Point source, Point start, Point end) {
	if (fabs(length(start, end) - (length(source, start) + length(source, end))) < 1e-10)
		return true;
	else return false;
}

bool BelongToPolygon(Point source, vector<Point> poligons)
{
	if (poligons.size() < 2)return false;
	Point currentPoint = poligons[0];
	for (size_t i = 1; i < poligons.size(); i++)
	{
		if (BelongToLine(source, currentPoint, poligons[i]))
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
				double alfa = reverseKoef(xDivide[i].koefRazr, 1);
				double hX_pred = 0;
				double hY_pred = 0;
				for (size_t indX = 1; indX < xDivide[i].nIntervals; indX++, countX++)
				{
					hX_pred += otnosh*lengthX;
					hY_pred += otnosh*lengthY;
					arrayCells[countX][countY].x = lines[i][j].x + hX_pred;
					arrayCells[countX][countY].y = lines[i][j].y + hY_pred;
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
				double alfa = reverseKoef(yDivide[j].koefRazr, 1);
				double hX_pred = 0;
				double hY_pred = 0;
				for (size_t indY = 1; indY < yDivide[j].nIntervals; indY++, countY++)
				{
					hX_pred += otnosh*lengthX;
					hY_pred += otnosh*lengthY;
					arrayCells[countX][countY].x = arrayCells[countX][startY].x + hX_pred;
					arrayCells[countX][countY].y = arrayCells[countX][startY].y + hY_pred;
					otnosh *= alfa;
				}
			}
		}
	}

	for (int j = yTotalCells-1; j >=0; j--)
	{
		for (int i = 0; i < xTotalCells; i++)
		{
			cells << setw(15) << arrayCells[i][j].x << " " << setw(8) << arrayCells[i][j].y;
		}
		cells << endl;
	}
	
	ifstream ignoredFieldFile("ignoredFields.txt");
	int nFields;
	ignoredFieldFile >> nFields;
	for (size_t i = 0; i < nFields; i++)
	{
		double x, y;
		ignoredFieldFile >> x >> y;
		ignoredField.push_back(Point(x, y));
	}
	for (size_t j = 0; j < yTotalCells - 1; j++)
	{
		//*************пробежались по Х вдоль основной линии
		for (size_t i = 0; i < xTotalCells - 1; i++)
		{
			if (!Inside(arrayCells[i][j], ignoredField)) {
				xy.push_back(arrayCells[i][j]);
			}
			else
				if (BelongToPolygon(arrayCells[i][j], ignoredField)) {
					xy.push_back(arrayCells[i][j]);
				}

			//добавляем доп. точку для квадратичных базисных функций
			Point dopTochka((arrayCells[i][j].x + arrayCells[i + 1][j].x) / 2, 
				(arrayCells[i][j].y + arrayCells[i + 1][j].y) / 2);
			if (!Inside(dopTochka, ignoredField)) {
				xy.push_back(dopTochka);
			}
			else
				if (BelongToPolygon(dopTochka, ignoredField)) {
					xy.push_back(dopTochka);
				}
		}
		//добавляем последнюю по Х точку
		if (!Inside(arrayCells[xTotalCells-1][j], ignoredField)) {
			xy.push_back(arrayCells[xTotalCells - 1][j]);
		}
		else
			if (BelongToPolygon(arrayCells[xTotalCells - 1][j], ignoredField)) {
				xy.push_back(arrayCells[xTotalCells - 1][j]);
			}

		//*************закончили пробег по Х
		//*************пробежались по Х вдоль средней линии
		for (size_t i = 0; i < xTotalCells - 1; i++)
		{
			Point mainPoint((arrayCells[i][j].x + arrayCells[i][j + 1].x) / 2,
				(arrayCells[i][j].y + arrayCells[i][j + 1].y) / 2);
			if (!Inside(mainPoint, ignoredField)) {
				xy.push_back(mainPoint);
			}
			else
				if (BelongToPolygon(mainPoint, ignoredField)) {
					xy.push_back(mainPoint);
				}

			//добавляем доп. точку для квадратичных базисных функций
			Point nextPoint((arrayCells[i + 1][j].x + arrayCells[i + 1][j + 1].x) / 2, 
				(arrayCells[i + 1][j].y + arrayCells[i + 1][j + 1].y) / 2);
			Point dopTochka((mainPoint.x + nextPoint.x) / 2, (mainPoint.y + nextPoint.y) / 2);
			if (!Inside(dopTochka, ignoredField)) {
				xy.push_back(dopTochka);
			}
			else
				if (BelongToPolygon(dopTochka, ignoredField)) {
					xy.push_back(dopTochka);
				}
		}
		//добавляем последнюю по Х точку
		Point dopTochka((arrayCells[xTotalCells - 1][j].x + arrayCells[xTotalCells - 1][j + 1].x) / 2, 
			(arrayCells[xTotalCells - 1][j].y + arrayCells[xTotalCells - 1][j + 1].y) / 2);
		if (!Inside(dopTochka, ignoredField)) {
			xy.push_back(dopTochka);
		}
		else
			if (BelongToPolygon(dopTochka, ignoredField)) {
				xy.push_back(dopTochka);
			}

		//*************закончили пробег по Х

	}

	//*************пробежались по Х вдоль последней основной горизонтальной линии
	for (size_t i = 0; i < xTotalCells - 1; i++)
	{
		if (!Inside(arrayCells[i][yTotalCells - 1], ignoredField)) {
			xy.push_back(arrayCells[i][yTotalCells - 1]);
		}
		else
			if (BelongToPolygon(arrayCells[i][yTotalCells - 1], ignoredField)) {
				xy.push_back(arrayCells[i][yTotalCells - 1]);
			}

		//добавляем доп. точку для квадратичных базисных функций
		Point dopTochka((arrayCells[i][yTotalCells - 1].x + arrayCells[i + 1][yTotalCells - 1].x) / 2,
			(arrayCells[i][yTotalCells - 1].y + arrayCells[i + 1][yTotalCells - 1].y) / 2);
		if (!Inside(dopTochka, ignoredField)) {
			xy.push_back(dopTochka);
		}
		else
			if (BelongToPolygon(dopTochka, ignoredField)) {
				xy.push_back(dopTochka);
			}
	}
	//добавляем последнюю по Х точку
	if (!Inside(arrayCells[xTotalCells - 1][yTotalCells - 1], ignoredField)) {
		xy.push_back(arrayCells[xTotalCells - 1][yTotalCells - 1]);
	}
	else
		if (BelongToPolygon(arrayCells[xTotalCells - 1][yTotalCells - 1], ignoredField)) {
			xy.push_back(arrayCells[xTotalCells - 1][yTotalCells - 1]);
		}

	//*************закончили пробег по Х

	//формируем КЭ
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

			temp.uzel[4] = temp.uzel[0] + 1;
			temp.uzel[8] = temp.uzel[2] + 1;

			temp.uzel[5] = indexXY(Point((arrayCells[i][j].x + arrayCells[i][j + 1].x) / 2, 
				(arrayCells[i][j].y + arrayCells[i][j + 1].y) / 2));
			temp.uzel[6] = temp.uzel[5] + 1;
			temp.uzel[7] = temp.uzel[6] + 1;
			//temp.uzel[6] = indexXY(Point((arrayCells[i][j].x + arrayCells[i + 1][j + 1].x) / 2,
			//	(arrayCells[i][j].y + arrayCells[i + 1][j + 1].y) / 2));
			//temp.uzel[7] = indexXY(Point((arrayCells[i + 1][j].x + arrayCells[i + 1][j + 1].x) / 2,
			//	(arrayCells[i + 1][j].y + arrayCells[i + 1][j + 1].y) / 2));
			temp.numberField = 0;
			KE.push_back(temp);
		}
	}

	//Очистка памяти
	for (size_t i = 0; i < xTotalCells; i++)
	{
		delete arrayCells[i];
	}
	delete arrayCells;

}

void PrintGlobalMatrix() {
	ofstream A("A.txt");
	for (size_t i = 0; i < xy.size(); i++)
	{
		for (size_t j = 0; j < xy.size(); j++)
		{
			if (i == j) {
				A << setw(18) << di[i];
				continue;
			}
			else {
				int I, J;
				bool flag;
				if (i > j) {
					I = i; J = j;
					flag = true;
				}
				else {
					I = j; J = i;
					flag = false;
				}
				for (size_t l = ig[I]; l < ig[I+1]; l++)
				{
					if (jg[l] == J) {
						flag ? A << setw(18) << ggl[l] : A << setw(18) << ggu[l];
					}
				}
			}
		}
		A << endl;
	}
	A << "B:" << endl;
	for (size_t i = 0; i < xy.size(); i++)
	{
		A << setw(18) << b[i];
	}
	A.close();
}

void drawText(const char *text, int length, int x, int y)
{
	glMatrixMode(GL_PROJECTION);
	double *matrix = new double[16];
	glGetDoublev(GL_PROJECTION_MATRIX, matrix);
	glLoadIdentity();
	glOrtho(0, Width, 0, Height, -5, 5);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPushMatrix();
	glLoadIdentity();
	glRasterPos2i(x, y);
	for (int i = 0; i<length; i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, (int)text[i]);//GLUT_BITMAP_9_BY_15
	}
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(matrix);
	glMatrixMode(GL_MODELVIEW);

}

void GetColorForPoint(int indPoint) {
	int i;
	for (i = 0; i < COUNT_OF_COLOR_AREAS - 1 && q[indPoint] >= rainbow[i + 1].value; i++);

		//glColor3ub((rainbow[i].red+rainbow[i+1].red)/2,
		//	(rainbow[i].green + rainbow[i + 1].green) / 2, 
		//	(rainbow[i].blue + rainbow[i + 1].blue) / 2);
	
	double basFuncOfValue = (q[indPoint] - rainbow[i].value) / (rainbow[i + 1].value - rainbow[i].value);

	double red = rainbow[i].red + basFuncOfValue*(rainbow[i + 1].red - rainbow[i].red);
	double green = rainbow[i].green + basFuncOfValue*(rainbow[i + 1].green - rainbow[i].green);
	double blue = rainbow[i].blue + basFuncOfValue*(rainbow[i + 1].blue - rainbow[i].blue);

	glColor3ub(rainbow[i].red + basFuncOfValue*(rainbow[i + 1].red - rainbow[i].red), 
		rainbow[i].green + basFuncOfValue*(rainbow[i + 1].green - rainbow[i].green),
		rainbow[i].blue + basFuncOfValue*(rainbow[i + 1].blue - rainbow[i].blue));
}

bool CheckIsolineOnTriangleEdge(double isoline, Point v0, Point v1, double q0, double q1, double E_x, double E_y, double tempX, double tempY) {
	double minQ, maxQ;
	if (q0 < q1) {
		minQ = q0; maxQ = q1;
	}
	else {
		maxQ = q0; minQ = q1;
	}
	if (isoline >= minQ && isoline <= maxQ) {
		//double basiqFunc = (isoline - q[indMinQ]) / (q[indMaxQ] - q[indMinQ]);
		//double x = xy[indMinQ].x + basiqFunc*fabs(xy[indMinQ].x - xy[indMaxQ].x);
		//double y = xy[indMinQ].y + basiqFunc*fabs(xy[indMinQ].y - xy[indMaxQ].y);
		//glVertex2f((E_x + xy[indMinQ].x + basiqFunc*fabs(xy[indMinQ].x- xy[indMaxQ].x)) * tempX, (E_y + xy[indMinQ].y + basiqFunc*fabs(xy[indMinQ].y - xy[indMaxQ].y)) * tempY);
		glVertex2f((E_x + (v0.x + v1.x) / 2) * tempX, (E_y + (v0.y + v1.y) / 2) * tempY);
		return true;
	}
	return false;
}

void DrawIsolineInTriangle(TriangleGeneral triangle, double E_x, double E_y, double tempX, double tempY) {
	glColor3f(0, 0, 0);
	for (size_t i = 0; i < COUNT_OF_ISOLINES; i++)
	{
		glBegin(GL_LINES);
		CheckIsolineOnTriangleEdge(isolinesValues[i], triangle.nvtr[0], triangle.nvtr[1], triangle.solutions[0], triangle.solutions[1], E_x, E_y, tempX, tempY);
		CheckIsolineOnTriangleEdge(isolinesValues[i], triangle.nvtr[1], triangle.nvtr[2], triangle.solutions[1], triangle.solutions[2], E_x, E_y, tempX, tempY);
		CheckIsolineOnTriangleEdge(isolinesValues[i], triangle.nvtr[0], triangle.nvtr[2], triangle.solutions[0], triangle.solutions[2], E_x, E_y, tempX, tempY);
		glEnd();
	}
}
void DrawFieldOnTriangle(Triangle triangle, double E_x, double E_y, double tempX, double tempY) {
	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i < 3; i++)
	{
		GetColorForPoint(triangle.nvtr[i]);
		glVertex2f((E_x + xy[triangle.nvtr[i]].x) * tempX, (E_y + xy[triangle.nvtr[i]].y) * tempY);
	}
	glEnd();

}

void DrawFieldOnQuad(int v0, int v1, int v2, int v3, double E_x, double E_y, double tempX, double tempY) {
	Triangle triangle;
	triangle.nvtr[0] = v0; triangle.nvtr[1] = v1; triangle.nvtr[2] = v2;
	DrawFieldOnTriangle(triangle, E_x, E_y, tempX, tempY);
	triangle.nvtr[0] = v1; triangle.nvtr[1] = v2; triangle.nvtr[2] = v3;
	DrawFieldOnTriangle(triangle, E_x, E_y, tempX, tempY);

}

Point CalcQuadranglesCoordinates(double ksi, double eta, double beta[6], double x0, double y0) {
	Point point;
	point.x = x0 + beta[1] * ksi + beta[0] * eta + beta[4] * ksi*eta;
	point.y = y0 + beta[3] * ksi + beta[2] * eta + beta[5] * ksi*eta;
	return point;
}
double CalcSolutionInPointOfSquare(double ksi, double eta, int ielem) {
	double result = 0;
	int relation[] = {
		0,4,1,5,6,7,2,8,3
	};
	for (size_t i = 0; i < 9; i++)
	{
		result += BasicFunc2d(i, ksi, eta)*q[KE[ielem].uzel[relation[i]]];
	}
	return result;
}

vector<Quadrangle> GenerateInnerQuadres(int ielem) {
	double beta[6] = {
		xy[KE[ielem].uzel[2]].x - xy[KE[ielem].uzel[0]].x,
		xy[KE[ielem].uzel[1]].x - xy[KE[ielem].uzel[0]].x,
		xy[KE[ielem].uzel[2]].y - xy[KE[ielem].uzel[0]].y,
		xy[KE[ielem].uzel[1]].y - xy[KE[ielem].uzel[0]].y,
		xy[KE[ielem].uzel[0]].x - xy[KE[ielem].uzel[1]].x - xy[KE[ielem].uzel[2]].x + xy[KE[ielem].uzel[3]].x,
		xy[KE[ielem].uzel[0]].y - xy[KE[ielem].uzel[1]].y - xy[KE[ielem].uzel[2]].y + xy[KE[ielem].uzel[3]].y
	};
	double hSquare = 1. / COUNT_OF_X_INNER_QUADRES;
	vector<Quadrangle> quadrangles;
	Quadrangle quad;
	for (size_t j = 0; j < COUNT_OF_X_INNER_QUADRES; j++)
	{
		double eta = hSquare*j;
		for (size_t i = 0; i < COUNT_OF_X_INNER_QUADRES; i++)
		{
			double ksi = hSquare*i;
			quad.nvtr[0] = CalcQuadranglesCoordinates(ksi, eta, beta, xy[KE[ielem].uzel[0]].x, xy[KE[ielem].uzel[0]].y);
			quad.nvtr[1] = CalcQuadranglesCoordinates(ksi + hSquare, eta, beta, xy[KE[ielem].uzel[0]].x, xy[KE[ielem].uzel[0]].y);
			quad.nvtr[2] = CalcQuadranglesCoordinates(ksi, eta + hSquare, beta, xy[KE[ielem].uzel[0]].x, xy[KE[ielem].uzel[0]].y);
			quad.nvtr[3] = CalcQuadranglesCoordinates(ksi + hSquare, eta + hSquare, beta, xy[KE[ielem].uzel[0]].x, xy[KE[ielem].uzel[0]].y);
			quad.solutions[0] = CalcSolutionInPointOfSquare(ksi, eta, ielem);
			quad.solutions[1] = CalcSolutionInPointOfSquare(ksi + hSquare, eta, ielem);
			quad.solutions[2] = CalcSolutionInPointOfSquare(ksi, eta + hSquare, ielem);
			quad.solutions[3] = CalcSolutionInPointOfSquare(ksi + hSquare, eta + hSquare, ielem);
			quadrangles.push_back(quad);
		}
	}
	return quadrangles;
}

void Display(void) // функция вывода
{
	glClearColor(1, 1, 1, 1); // очистка буфера
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3ub(0, 0, 0);
	//glShadeModel(GL_FLAT); // отключение интерполяции
	
	double h_x = lines[Nx - 1][0].x - lines[0][0].x;
	double h_y = lines[0][Ny - 1].y - lines[0][0].y;
	double E_x = -lines[0][0].x/*+0.03*h_x*/;
	double E_y = -lines[0][0].y/* + 0.03*h_y*/;
	double tempX,tempY;
	tempX = (Width - 50) / h_x;
	tempY = (Height - 50) / h_y;

	//закрашивание КЭ
	for (size_t iKE = 0; iKE < KE.size(); iKE++)
	{
		//деление по диагонали снизу-вверх слева направо
		DrawFieldOnQuad(KE[iKE].uzel[0], KE[iKE].uzel[4], KE[iKE].uzel[5], KE[iKE].uzel[6], E_x, E_y, tempX, tempY);
		DrawFieldOnQuad(KE[iKE].uzel[4], KE[iKE].uzel[1], KE[iKE].uzel[6], KE[iKE].uzel[7], E_x, E_y, tempX, tempY);
		DrawFieldOnQuad(KE[iKE].uzel[5], KE[iKE].uzel[6], KE[iKE].uzel[2], KE[iKE].uzel[8], E_x, E_y, tempX, tempY);
		DrawFieldOnQuad(KE[iKE].uzel[6], KE[iKE].uzel[7], KE[iKE].uzel[8], KE[iKE].uzel[3], E_x, E_y, tempX, tempY);

		//рисуем изолинии
		if (isDisplayingIsolines) {
			vector<Quadrangle> quadrangles = GenerateInnerQuadres(iKE);
			for (size_t iQuad = 0; iQuad < quadrangles.size(); iQuad++)
			{
				TriangleGeneral triangle;
				triangle.nvtr[0] = quadrangles[iQuad].nvtr[0]; triangle.nvtr[1] = quadrangles[iQuad].nvtr[1]; triangle.nvtr[2] = quadrangles[iQuad].nvtr[2];
				triangle.solutions[0] = quadrangles[iQuad].solutions[0]; triangle.solutions[1] = quadrangles[iQuad].solutions[1]; triangle.solutions[2] = quadrangles[iQuad].solutions[2];
				DrawIsolineInTriangle(triangle, E_x, E_y, tempX, tempY);

				triangle.nvtr[0] = quadrangles[iQuad].nvtr[3]; triangle.nvtr[1] = quadrangles[iQuad].nvtr[1]; triangle.nvtr[2] = quadrangles[iQuad].nvtr[2];
				triangle.solutions[0] = quadrangles[iQuad].solutions[3]; triangle.solutions[1] = quadrangles[iQuad].solutions[1]; triangle.solutions[2] = quadrangles[iQuad].solutions[2];
				DrawIsolineInTriangle(triangle, E_x, E_y, tempX, tempY);

				////отрисовка мелких четырехугольников
				if (isDisplayingQuadrangles) {
					glColor3f(0, 0, 0);
					glBegin(GL_LINE_LOOP);
					glVertex2f((E_x + quadrangles[iQuad].nvtr[0].x) * tempX, (E_y + quadrangles[iQuad].nvtr[0].y) * tempY);
					glVertex2f((E_x + quadrangles[iQuad].nvtr[1].x) * tempX, (E_y + quadrangles[iQuad].nvtr[1].y) * tempY);
					glVertex2f((E_x + quadrangles[iQuad].nvtr[3].x) * tempX, (E_y + quadrangles[iQuad].nvtr[3].y) * tempY);
					glVertex2f((E_x + quadrangles[iQuad].nvtr[2].x) * tempX, (E_y + quadrangles[iQuad].nvtr[2].y) * tempY);
					glEnd();
				}
			}
		}

		//отрисовка КЭ
		glColor3ub(100, 100, 100);
		glLineWidth(1);
		glBegin(GL_LINE_LOOP);
		glVertex2f((E_x + xy[KE[iKE].uzel[0]].x) * tempX, (E_y + xy[KE[iKE].uzel[0]].y) * tempY);
		glVertex2f((E_x + xy[KE[iKE].uzel[1]].x) * tempX, (E_y + xy[KE[iKE].uzel[1]].y) * tempY);
		glVertex2f((E_x + xy[KE[iKE].uzel[3]].x) * tempX, (E_y + xy[KE[iKE].uzel[3]].y) * tempY);
		glVertex2f((E_x + xy[KE[iKE].uzel[2]].x) * tempX, (E_y + xy[KE[iKE].uzel[2]].y) * tempY);
		glEnd();
	}

	//отрисовка линий подобластей
	glLineWidth(2);
	glColor3ub(255, 0, 0);
	for (int i = 0; i < Nx && isDisplayingAreas; i++)
	{
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < Ny; j++)
		{
			glVertex2f((E_x + lines[i][j].x) * tempX, (E_y + lines[i][j].y) * tempY);
		}
		glEnd();
	}
	for (int j = 0; j < Ny && isDisplayingAreas; j++)
	{
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < Nx; i++)
		{
			glVertex2f((E_x + lines[i][j].x) * tempX, (E_y + lines[i][j].y) * tempY);
		}
		glEnd();
	}

	//нумерация функций
	glColor3ub(0, 0, 0);
	for (int i = 0; i < xy.size() && isDisplayingFunctions; i++)
	{
					char*text = new char[3];
					//sprintf_s(text, 2, "%d", i);
					sprintf(text, "%d", i);
					drawText(text, strlen(text), (E_x + xy[i].x) * tempX, (E_y + xy[i].y) * tempY);
	}


	//нумерация ребер
	glColor3ub(0, 255, 0);
	for (size_t i = 0; i < xy.size() && isDisplayingEdges; i++)
	{
		for (size_t j = igEdge[i]; j < igEdge[i + 1]; j++)
		{
			char*text = new char[3];
			//sprintf_s(text, 2, "%d", i);
			sprintf(text, "%d", j);
			drawText(text, strlen(text), (E_x + (xy[i].x + xy[jgEdge[j]].x) / 2) * tempX, 
				(E_y + (xy[i].y + xy[jgEdge[j]].y) / 2) * tempY);
		}
	}
	glFinish();
}
void Reshape(GLint w, GLint h) // При изменении размеров окна
{
	Width = w;
	Height = h;
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, w, 0, h, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void GeneratePortrait() {
	set<size_t> * portrait = new set<size_t>[xy.size()];
	for (size_t k = 0; k < KE.size(); k++)
	{
		for (size_t i = 0; i < 9; i++)
		{
			size_t a = KE[k].uzel[i];
			for (size_t j = 0; j < i; j++)
			{
				size_t b = KE[k].uzel[j];
				if (b > a)
					portrait[b].insert(a);
				else
					portrait[a].insert(b);
			}
		}
	}

	ig = new int[xy.size() + 1];
	di = new double[xy.size()];
	b = new double[xy.size()];
	q = new double[xy.size()];
	ig[0] = ig[1] = 0;
	for (size_t i = 0; i < xy.size(); i++)
	{
		di[i] = b[i] = q[i] = 0;
		ig[i + 1] = ig[i] + portrait[i].size();
	}
	jg = new int[ig[xy.size()]];
	ggl = new double[ig[xy.size()]];
	ggu = new double[ig[xy.size()]];

	for (size_t i = 0; i < ig[xy.size()]; i++)
	{
		ggl[i] = 0;
		ggu[i] = 0;
	}

	size_t tmp = 0;
	for (size_t i = 0; i < xy.size(); i++)
	{
		for (set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
		{
			jg[tmp] = *j;
			tmp++;
		}
		portrait[i].clear();
	}
	delete[] portrait;

	int i, j;
	ofstream igOut("ig.txt");
	ofstream jgOut("jg.txt");
	for (i = 0; i <= xy.size(); i++)
	{
		igOut << ig[i] << " ";
	}
	for (j = 0; j < ig[xy.size()]; j++)
	{
		jgOut << jg[j] << " ";
	}
}

bool IsEdgeExist(int ielem, int ind1, int ind2) {
	if (abs(ind2 - ind1) == 2 ||
		(abs(ind2 - ind1) == 1 && ind2 != 2))return true;
	else return false;
}
void GenerateEdges()
{
	int KEcount = KE.size();

	int *list[2], *listbeg;//список
	list[0] = new int[xy.size()*xy.size()];
	list[1] = new int[xy.size()*xy.size()];
	listbeg = new int[xy.size()];
	int listsize = -1;//количество элементов в списке, а также в jg

	for (int i = 0; i < xy.size(); i++)	listbeg[i] = -1;

	for (int ielem = 0; ielem < KEcount; ielem++)//проходим по всем КЭ
	{
		if (ielem == 9)
			ielem = ielem;
		for (int i = 0; i < 4; i++)//перебираем базисные функции
		{
			int kk = KE[ielem].uzel[i];//kk - глобальный номер текущей рассматриваемой базисной функции
									   //int kk = 2 * (ielem / (Nx - 1))*(2 * Nx - 1) + 2 * (ielem % (Nx - 1));
			for (int j = i + 1; j < 4; j++)//перебираем следующие за ней функции на элементе
			{
				int ind1 = kk;
				int ind2 = KE[ielem].uzel[j];
				if (!IsEdgeExist(ielem, i, j))
				continue;
				if (ind2 < ind1)//вносится связь большего с меньшим
				{
					ind1 = ind2;
					ind2 = kk;
				}
				int iaddr = listbeg[ind2];
				if (iaddr < 0)//список был пуст
				{
					//создание списка
					listsize++;
					listbeg[ind2] = listsize;
					list[0][listsize] = ind1;
					list[1][listsize] = -1;
				}
				else//список не был пуст
				{
					//ищем в списке ind1
					while (list[0][iaddr] < ind1 && list[1][iaddr] >= 0)		iaddr = list[1][iaddr];
					if (list[0][iaddr] > ind1)
					{
						//если найденный там элемент имеет больший номер, 
						//то нужно добавить перед ним, чтобы список был упорядоченным
						listsize++;
						list[0][listsize] = list[0][iaddr];	//перекладываем вперед найденный элемент
						list[1][listsize] = list[1][iaddr];
						list[0][iaddr] = ind1;			//на его место ложим новый
						list[1][iaddr] = listsize;
					}
					else
						if (list[0][iaddr] < ind1)
						{
							//не нашли, а список закончился
							//добавляем в конец списка
							listsize++;
							list[1][iaddr] = listsize;
							list[0][listsize] = ind1;
							list[1][listsize] = -1;	//указываем, что это последний элемент списка
						}
				}
			}
		}
	}
	int gk = listsize + 1;
	jgEdge = new int[gk];
	igEdge = new int[xy.size() + 1];
	igEdge[0] = 0;
	for (int i = 0; i < xy.size(); i++)
	{
		igEdge[i + 1] = igEdge[i];			//igEdge[i+1] - номер ячейки массива jgEdge, куда надо поместить следующий элемент
		int iaddr = listbeg[i];		//iaddr - индекс начала списка
		while (iaddr >= 0)			//просматриваем список
		{
			jgEdge[igEdge[i + 1]] = list[0][iaddr];
			igEdge[i + 1]++;
			iaddr = list[1][iaddr]; //переходим к следующему элементу
		}
	}
	delete list[0]; delete list[1]; delete listbeg;


	ofstream igEdgeOut("igEdge.txt");
	ofstream jgEdgeOut("jgEdge.txt");
	for (int i = 0; i <= xy.size(); i++)
	{
		igEdgeOut << igEdge[i] << " ";
	}
	for (int j = 0; j < igEdge[xy.size()]; j++)
	{
		jgEdgeOut << jgEdge[j] << " ";
	}
}
void genProfile()
{
	int KEcount = KE.size();

	int *list[2], *listbeg;//список
	list[0] = new int[xy.size()*xy.size()];
	list[1] = new int[xy.size()*xy.size()];
	listbeg = new int[xy.size()];
	int listsize = -1;//количество элементов в списке, а также в jg

	for (int i = 0; i < xy.size(); i++)	listbeg[i] = -1;

	for (int ielem = 0; ielem < KEcount; ielem++)//проходим по всем КЭ
	{
		for (int i = 0; i < 9; i++)//перебираем базисные функции
		{
			int kk = KE[ielem].uzel[i];//kk - глобальный номер текущей рассматриваемой базисной функции
									//int kk = 2 * (ielem / (Nx - 1))*(2 * Nx - 1) + 2 * (ielem % (Nx - 1));
			for (int j = i + 1; j < 9; j++)//перебираем следующие за ней функции на элементе
			{
				int ind1 = kk;
				int ind2 = KE[ielem].uzel[j];
				if (ind2 < ind1)//вносится связь большего с меньшим
				{
					ind1 = ind2;
					ind2 = kk;
				}
				int iaddr = listbeg[ind2];
				if (iaddr < 0)//список был пуст
				{
					//создание списка
					listsize++;
					listbeg[ind2] = listsize;
					list[0][listsize] = ind1;
					list[1][listsize] = -1;
				}
				else//список не был пуст
				{
					//ищем в списке ind1
					while (list[0][iaddr] < ind1 && list[1][iaddr] >= 0)		iaddr = list[1][iaddr];
					if (list[0][iaddr] > ind1)
					{
						//если найденный там элемент имеет больший номер, 
						//то нужно добавить перед ним, чтобы список был упорядоченным
						listsize++;
						list[0][listsize] = list[0][iaddr];	//перекладываем вперед найденный элемент
						list[1][listsize] = list[1][iaddr];
						list[0][iaddr] = ind1;			//на его место ложим новый
						list[1][iaddr] = listsize;
					}
					else
						if (list[0][iaddr] < ind1)
						{
							//не нашли, а список закончился
							//добавляем в конец списка
							listsize++;
							list[1][iaddr] = listsize;
							list[0][listsize] = ind1;
							list[1][listsize] = -1;	//указываем, что это последний элемент списка
						}
				}
			}
		}
	}
	int gk = listsize + 1;
	ig = new int[xy.size() + 1];
	jg = new int[gk];
	ig[0] = 0;
	for (int i = 0; i < xy.size(); i++)
	{
		ig[i + 1] = ig[i];			//ig[i+1] - номер ячейки массива jg, куда надо поместить следующий элемент
		int iaddr = listbeg[i];		//iaddr - индекс начала списка
		while (iaddr >= 0)			//просматриваем список
		{
			jg[ig[i + 1]] = list[0][iaddr];
			ig[i + 1]++;
			iaddr = list[1][iaddr]; //переходим к следующему элементу
		}
	}
	delete list[0]; delete list[1]; delete listbeg;


	ofstream igOut("ig.txt");
	ofstream jgOut("jg.txt");
	for (int i = 0; i <= xy.size(); i++)
	{
		igOut << ig[i] << " ";
	}
	for (int j = 0; j < ig[xy.size()]; j++)
	{
		jgOut << jg[j] << " ";
	}

	ggl = new double[ig[xy.size()]];
	ggl = new double[gk];
	ggu = new double[gk];
	di = new double[xy.size()];
	b = new double[xy.size()];
	q = new double[xy.size()];

	for (int i = 0; i < gk; i++)	ggl[i] = ggu[i] = 0;
	for (int i = 0; i < xy.size(); i++) di[i] = b[i] = q[i] = 0;
}

double Lambda(int numberField) {
	switch (numberField) {
	case 0: return 1;
	default:
		cerr << "Error in Lambda" << endl;
		system("pause");
		exit(1);
	}
}

double Gamma(int numberField) {
	switch (numberField) {
	case 0: return 1;
	default:
		cerr << "Error in Gamma" << endl;
		system("pause");
		exit(1);
	}
}

double Func(Point p) {
	//return p.x;
	//return p.y;
	//return 1 + p.x + p.y;
	//return - 2 * p.x + p.x * p.x - p.y;
	return (p.x - 15)*(p.x - 15) + (p.y - 10)*(p.y - 10) - 4;
	//return p.x*p.x + p.y*p.y - p.x*p.y - 4;
	//return p.x + p.y + p.x*p.y + p.x*p.x + p.y*p.y - 3;
	//return p.x*p.x*p.x + p.y*p.y*p.y - 6 * p.x - 6 * p.y;
	//return -0.02*exp(0.1*(p.x + p.y)) + exp(0.1*(p.x + p.y));
}

double AnaliticSolution(Point p) {
	//return p.x;
	//return p.y;
	//return 1 + p.x + p.y;
	//return 2 - 2 * p.x + p.x * p.x - p.y;
	return (p.x-15)*(p.x - 15) + (p.y - 10)*(p.y - 10);
	//return p.x*p.x + p.y*p.y - p.x*p.y;
	//return 1 + p.x + p.y + p.x*p.y + p.x*p.x + p.y*p.y;
	//return p.x*p.x*p.x + p.y*p.y*p.y;
	//return exp(0.1*(p.x + p.y));
}

void CreateLocalMatrix(int ielem, double integrPoints[], double tauKoefs[], int countExtraPoints, double A[9][9], double localB[9]) {
	double alfa0, alfa1, alfa2, signAlfa0;
	alfa0 = (xy[KE[ielem].uzel[1]].x - xy[KE[ielem].uzel[0]].x)*(xy[KE[ielem].uzel[2]].y - xy[KE[ielem].uzel[0]].y)
		- (xy[KE[ielem].uzel[1]].y - xy[KE[ielem].uzel[0]].y)*(xy[KE[ielem].uzel[2]].x - xy[KE[ielem].uzel[0]].x);
	alfa1= (xy[KE[ielem].uzel[1]].x - xy[KE[ielem].uzel[0]].x)*(xy[KE[ielem].uzel[3]].y - xy[KE[ielem].uzel[2]].y)
		- (xy[KE[ielem].uzel[1]].y - xy[KE[ielem].uzel[0]].y)*(xy[KE[ielem].uzel[3]].x - xy[KE[ielem].uzel[2]].x);
	alfa2 = (xy[KE[ielem].uzel[2]].y - xy[KE[ielem].uzel[0]].y)*(xy[KE[ielem].uzel[3]].x - xy[KE[ielem].uzel[1]].x)
		- (xy[KE[ielem].uzel[2]].x - xy[KE[ielem].uzel[0]].x)*(xy[KE[ielem].uzel[3]].y - xy[KE[ielem].uzel[1]].y);
	double beta[6] = {
		xy[KE[ielem].uzel[2]].x - xy[KE[ielem].uzel[0]].x,
		xy[KE[ielem].uzel[1]].x - xy[KE[ielem].uzel[0]].x,
		xy[KE[ielem].uzel[2]].y - xy[KE[ielem].uzel[0]].y,
		xy[KE[ielem].uzel[1]].y - xy[KE[ielem].uzel[0]].y,
		xy[KE[ielem].uzel[0]].x - xy[KE[ielem].uzel[1]].x - xy[KE[ielem].uzel[2]].x + xy[KE[ielem].uzel[3]].x,
		xy[KE[ielem].uzel[0]].y - xy[KE[ielem].uzel[1]].y - xy[KE[ielem].uzel[2]].y + xy[KE[ielem].uzel[3]].y
	};
	if (alfa0 > 0)signAlfa0 = 1;
	else signAlfa0 = -1;

	for (size_t i = 0; i < 9; i++)
	{
		for (size_t j = 0; j <= i; j++)
		{
			double resultG = 0;
			double resultM = 0;
			for (size_t k = 0; k < countExtraPoints; k++)
			{
				for (size_t l = 0; l < countExtraPoints; l++)
				{
					//double ksi = (beta[2] * (integrPoints[k] - xy[KE[ielem].uzel[0]].x) - beta[0] * (integrPoints[l] - xy[KE[ielem].uzel[0]].y)) /
					//	(beta[1] * beta[2] - beta[0] * beta[3]);
					//double eta = (beta[1] * (integrPoints[l] - xy[KE[ielem].uzel[0]].y) - beta[3] * (integrPoints[k] - xy[KE[ielem].uzel[0]].x)) /
					//		(beta[1] * beta[2] - beta[0] * beta[3]);

					double J = alfa0 + alfa1*integrPoints[k] + alfa2*integrPoints[l];
					double dfi_I_dksi = difBasicFunc1d(i % 3, integrPoints[k])*BasicFunc1d(i / 3, integrPoints[l]);
					double dfi_I_deta = BasicFunc1d(i % 3, integrPoints[k])*difBasicFunc1d(i / 3, integrPoints[l]);
					double dfi_J_dksi = difBasicFunc1d(j % 3, integrPoints[k])*BasicFunc1d(j / 3, integrPoints[l]);
					double dfi_J_deta = BasicFunc1d(j % 3, integrPoints[k])*difBasicFunc1d(j / 3, integrPoints[l]);

					resultG +=
						((dfi_I_dksi*(beta[5] * integrPoints[k] + beta[2]) - dfi_I_deta*(beta[5] * integrPoints[l] + beta[3])) *
							(dfi_J_dksi*(beta[5] * integrPoints[k] + beta[2]) - dfi_J_deta*(beta[5] * integrPoints[l] + beta[3])) +
							(dfi_I_deta*(beta[4] * integrPoints[l] + beta[1]) - dfi_I_dksi*(beta[4] * integrPoints[k] + beta[0])) *
							(dfi_J_deta*(beta[4] * integrPoints[l] + beta[1]) - dfi_J_dksi*(beta[4] * integrPoints[k] + beta[0])))*
						tauKoefs[k] * tauKoefs[l] / J ;

					resultM += BasicFunc1d(i % 3, integrPoints[k])*BasicFunc1d(i / 3, integrPoints[l]) *
						BasicFunc1d(j % 3, integrPoints[k])*BasicFunc1d(j / 3, integrPoints[l])*
						tauKoefs[k] * tauKoefs[l] * J ;
				}
			}
			A[i][j] += resultG / 4.*signAlfa0*Lambda(KE[ielem].numberField) + resultM / 4.*signAlfa0*Gamma(KE[ielem].numberField);
		}

		//right part
		double resultRightPart = 0;
		for (size_t k = 0; k < countExtraPoints; k++)
		{
			for (size_t l = 0; l < countExtraPoints; l++)
			{
				double J = alfa0 + alfa1*integrPoints[k] + alfa2*integrPoints[l];

				resultRightPart += BasicFunc1d(i % 3, integrPoints[k])*BasicFunc1d(i / 3, integrPoints[l]) *
					Func(Point(
						xy[KE[ielem].uzel[0]].x + beta[1] * integrPoints[k] + beta[0] * integrPoints[l] + beta[4] * integrPoints[k] * integrPoints[l], 
						xy[KE[ielem].uzel[0]].y + beta[3] * integrPoints[k] + beta[2] * integrPoints[l] + beta[5] * integrPoints[k] * integrPoints[l]))*
					tauKoefs[k] * tauKoefs[l] * J;
			}
		}
		localB[i] += resultRightPart / 4. *signAlfa0;
	}
}

void AddToMatrix(int posI, int posJ, double el)
{
	int tmp;
	if (posI == posJ)
	{
		di[posI] += el;
		return;
	}
	else
	{
		if (posI < posJ)
		{
			return;
			tmp = posI;
			posI = posJ;
			posJ = tmp;
		}
		for (tmp = ig[posI]; tmp < ig[posI + 1]; tmp++)
		{
			if (jg[tmp] == posJ)
			{
				ggl[tmp] += el;
				return;
			}
		}
	}
}

void Addition(int ielem, double A[9][9], double localB[9]) {
	int i, j;
	int relation[] = {
		0,4,1,5,6,7,2,8,3
	};
	for (i = 0; i < 9; i++)
	{
		b[KE[ielem].uzel[relation[i]]] += localB[i];
		for (j = 0; j <= i; j++)
		{
			AddToMatrix(KE[ielem].uzel[relation[i]], KE[ielem].uzel[relation[j]], A[i][j]);
		}
	}
}

struct BoundType {
	int b[3];
	BoundType(int b1new,int b2new,int b3new){
		b[0] = b1new; b[1] = b2new; b[2] = b3new;
	}
	BoundType(const BoundType &newBoundType) {
		b[0] = newBoundType.b[0]; b[1] = newBoundType.b[1]; b[2] = newBoundType.b[2];
	}
};

vector<vector<BoundType>> GetBoundOfGlobalField() {
	vector<vector<BoundType>> mapBounds;
	mapBounds.resize(3);


	ifstream ku1bounds("edge1bounds.txt");
	int nKU1;
	ku1bounds >> nKU1;
	Point **boundsOfEdge1 = new Point*[nKU1];
	for (size_t i = 0; i < nKU1; i++)
	{
		boundsOfEdge1[i] = new Point[2];
		ku1bounds >> boundsOfEdge1[i][0].x >> boundsOfEdge1[i][0].y >> boundsOfEdge1[i][1].x >> boundsOfEdge1[i][1].y;
	}

	for (int i = 0; i < Ny; i++)
	{
		allBounds.push_back(lines[0][i]);
	}
	for (int i = 0; i < Nx && ignoredField.size()==0; i++) {
		allBounds.push_back(lines[i][Ny - 1]);
	}
	for (int i = 0; i < ignoredField.size(); i++)
	{
		allBounds.push_back(ignoredField[i]);
	}
	for (int i = Ny - 1; i > 0; i--)
	{
		allBounds.push_back(lines[Nx - 1][i]);
	}
	for (int i = Nx-1; i >= 0; i--)
	{
		allBounds.push_back(lines[i][0]);
	}
	for (size_t i = 0; i < xy.size(); i++)
	{
		if (i == 26)
			i = i;
		for (size_t j = igEdge[i]; j < igEdge[i+1]; j++)
		{
			Point dopPoint((xy[jgEdge[j]].x + xy[i].x) / 2, (xy[jgEdge[j]].y + xy[i].y) / 2);

			if (BelongToPolygon(Point(xy[i].x, xy[i].y), allBounds) &&
				BelongToPolygon(dopPoint, allBounds) &&
				BelongToPolygon(Point(xy[jgEdge[j]].x, xy[jgEdge[j]].y), allBounds)) {
				int k = indexXY(dopPoint);
				if (k == -1) {
					cerr << "Error in GetBoundOfGlobal" << endl;
					system("pause");
					exit(1);
				}
				bool flagIsEdge1 = false;
				for (size_t indEdge1 = 0; indEdge1 < nKU1; indEdge1++)
				{
					if (BelongToLine(Point(xy[i].x, xy[i].y), boundsOfEdge1[indEdge1][0], boundsOfEdge1[indEdge1][1]) &&
						BelongToLine(Point(xy[jgEdge[j]].x, xy[jgEdge[j]].y), boundsOfEdge1[indEdge1][0], boundsOfEdge1[indEdge1][1]))
					{
						flagIsEdge1 = true;
						break;
					}
				}
				BoundType edge(i, jgEdge[j], k);
				if (flagIsEdge1) {
					mapBounds[0].push_back(edge);
					//mapBounds[0].push_back(edge);
					cout << i << " " << jgEdge[j] << " " << k << endl;
				}
				else {
					mapBounds[1].push_back(edge);
				}
			}
		}
	}
	return mapBounds;
}

void Edge1_sim(vector<BoundType> boundsEdge1) {
	ofstream ku1("ku1.txt");
	for (size_t iBound = 0; iBound < boundsEdge1.size(); iBound++)
	{
		for (int ind = 0; ind < 3; ind++)//!!!!!!!!
		{
			int k = boundsEdge1[iBound].b[ind];
			di[k] = 1;
			for (int m = ig[k]; m < ig[k + 1]; m++)
			{
				b[jg[m]] -= ggl[m] * AnaliticSolution(xy[k]);
				ggl[m] = 0;
			}
			b[k] = AnaliticSolution(xy[k]);
			for (int l = 0; l < xy.size(); l++)
			{
				for (int m = ig[l]; m < ig[l + 1]; m++)
				{
					if (k == jg[m])
					{
						b[l] -= b[k] * ggl[m];
						ggl[m] = 0;
					}
				}
			}
			ku1 << k << '\t' << b[k] << endl;
		}
	}
}

void Edge1_not_sim(vector<BoundType> boundsEdge1) {
	ofstream ku1("ku1.txt");
	for (size_t iBound = 0; iBound < boundsEdge1.size(); iBound++)
	{
		for (int ind = 0; ind < 3; ind++)//!!!!!!!!
		{
			int k = boundsEdge1[iBound].b[ind];
			di[k] = 1;
			b[k] = AnaliticSolution(xy[k]);
			for (int m = ig[k]; m < ig[k + 1]; m++)
			{
				ggl[m] = 0;
			}
			for (int l = 0; l < xy.size(); l++)
			{
				for (int m = ig[l]; m < ig[l + 1]; m++)
				{
					if (k == jg[m])
					{
						ggu[m] = 0;
					}
				}
			}
			ku1 << k << '\t' << b[k] << endl;
		}
	}
}

void Edge2(vector<BoundType> boundsEdge2) {
	double h = 0.001;
	int relation[] = { 0,2,1 };
	int loc1dA[][3] = {
		{4, 2, -1},
		{ 2, 16, 2 },
		{ -1, 2, 4 }
	};
	for (size_t iBound = 0; iBound < boundsEdge2.size(); iBound++)
	{
		double tetta[3];
		double dh;
		int koef;
		if (xy[boundsEdge2[iBound].b[0]].x == xy[boundsEdge2[iBound].b[1]].x){
			if (Inside(Point(xy[boundsEdge2[iBound].b[2]].x - h, xy[boundsEdge2[iBound].b[2]].y), allBounds)) {
				koef = 1;
			}
			else {
				koef = -1;
			}
			dh = fabs(xy[boundsEdge2[iBound].b[1]].y - xy[boundsEdge2[iBound].b[0]].y);
			for (size_t indPos = 0; indPos < 3; indPos++)
			{
				tetta[relation[indPos]] = koef*(AnaliticSolution(Point(xy[boundsEdge2[iBound].b[indPos]].x + h, xy[boundsEdge2[iBound].b[indPos]].y)) -
					AnaliticSolution(Point(xy[boundsEdge2[iBound].b[indPos]].x - h, xy[boundsEdge2[iBound].b[indPos]].y))) / (2 * h);
			}
		}
		else{
			if (Inside(Point(xy[boundsEdge2[iBound].b[2]].x, xy[boundsEdge2[iBound].b[2]].y - h), allBounds)) {
				koef = 1;
			}
			else {
				koef = -1;
			}
			dh = fabs(xy[boundsEdge2[iBound].b[1]].x - xy[boundsEdge2[iBound].b[0]].x);
			for (size_t indPos = 0; indPos < 3; indPos++)
			{
				tetta[relation[indPos]] = koef*(AnaliticSolution(Point(xy[boundsEdge2[iBound].b[indPos]].x, xy[boundsEdge2[iBound].b[indPos]].y + h)) -
					AnaliticSolution(Point(xy[boundsEdge2[iBound].b[indPos]].x, xy[boundsEdge2[iBound].b[indPos]].y - h))) / (2 * h);
			}
		}
		for (size_t i = 0; i < 3; i++)
		{
			double value = 0;
			for (size_t j = 0; j < 3; j++)
			{
				value += loc1dA[i][j] * tetta[j];
			}
			b[boundsEdge2[iBound].b[relation[i]]] += value*dh / 30;
		}
	}
}

void GenerateGlobalMatrix() {
	int ielem, i, j;
	double tKoef = sqrt(3. / 5.);
	double integrationPoints[3] = {
		0.5,
		(1 + tKoef) / 2,
		(1 - tKoef) / 2
	};
	double tauKoefs[3] = {
		8. / 9.,
		5. / 9.,
		5. / 9.
	};
	double A[9][9];
	double localB[9];
	for (ielem = 0; ielem < KE.size(); ielem++)
	{
		for (size_t i = 0; i < 9; i++)
		{
			localB[i] = 0;
			for (size_t j = 0; j < 9; j++)
			{
				A[i][j] = 0;
			}
		}
		CreateLocalMatrix(ielem, integrationPoints, tauKoefs, 3, A, localB);
		Addition(ielem, A, localB);
	}
	PrintGlobalMatrix();
	vector<vector<BoundType>> bounds = GetBoundOfGlobalField();
	Edge2(bounds[0]);
	Edge1_sim(bounds[1]);
	for (i = 0; i < ig[xy.size()]; i++)
	{
		ggu[i] = ggl[i];
	}
	//Edge1_not_sim(bounds[0]);
	//PrintGlobalMatrix();
}

//	Умножение матрицы на вектор
void MultMatrixOnVector(double *in, double *out)
{
	int i, j;
	double *out1;
	out1 = new double[xy.size()];
	for (i = 0; i<xy.size(); i++)
	{
		out1[i] = di[i] * in[i];
		for (j = ig[i]; j<ig[i + 1]; j++)
		{
			out1[i] += ggl[j] * in[jg[j]];
			out1[jg[j]] += ggu[j] * in[i];
		}
	}
	for (i = 0; i<xy.size(); i++)
		out[i] = out1[i];
	delete[] out1;
}

double ScalarMult(double *v1, double *v2)
{
	int i;
	double result;
	result = 0;
	for (i = 0; i < xy.size(); i++)
	{
		result += v1[i] * v2[i];
	}
	return result;
}

void LOS()
{
	int maxiter = 10000, i;
	double alfa, alfachisl, alfaznam, beta, betachisl, betaznam, checkE, epsMSG = 1e-16, startNeviazka;
	double*r = new double[xy.size()];
	double*s = new double[xy.size()];
	double*z = new double[xy.size()];
	double*p = new double[xy.size()];
	double*rout = new double[xy.size()];
	for (i = 0; i < xy.size(); i++)
	{
		s[i] = rout[i] = r[i] = q[i] = z[i] = p[i] = 0;
	}
	MultMatrixOnVector(q, r);
	for (i = 0; i < xy.size(); i++)
	{
		r[i] = b[i] - r[i];
		z[i] = r[i];
	}
	MultMatrixOnVector(z, p);
	checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	//startNeviazka = checkE = ScalarMult(r, r);
	for (int iter = 0; iter < maxiter && checkE >= epsMSG; iter++)
	{
		alfachisl = ScalarMult(p, r);
		alfaznam = ScalarMult(p, p);
		alfa = alfachisl / alfaznam;
		for (i = 0; i < xy.size(); i++)
		{
			q[i] = q[i] + alfa*z[i];
			r[i] = r[i] - alfa*p[i];
		}
		MultMatrixOnVector(r, rout);
		betachisl = ScalarMult(p, rout);
		betaznam = ScalarMult(p, p);
		beta = -betachisl / betaznam;
		for (i = 0; i < xy.size(); i++)
		{
			z[i] = r[i] + beta*z[i];
			p[i] = rout[i] + beta*p[i];
		}
		checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	}
}

void MSG()
{
	int maxiter = 10000, i;
	double alfa, alfachisl, alfaznam, beta, betachisl, betaznam, checkE, epsMSG = 1e-16;
	double*r = new double[xy.size()];
	double*s = new double[xy.size()];
	double*z = new double[xy.size()];
	double*rout = new double[xy.size()];
	for (i = 0; i < xy.size(); i++)
	{
		s[i] = rout[i] = r[i] = q[i] = z[i] = 0;
	}
	MultMatrixOnVector(q, r);
	for (i = 0; i < xy.size(); i++)
	{
		r[i] = b[i] - r[i];
		z[i] = r[i];
	}
	checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	for (int iter = 0; iter < maxiter && checkE >= epsMSG; iter++)
	{
		alfachisl = ScalarMult(r, r);
		MultMatrixOnVector(z, s);
		alfaznam = ScalarMult(s, z);
		alfa = alfachisl / alfaznam;
		for (i = 0; i < xy.size(); i++)
		{
			q[i] = q[i] + alfa*z[i];
			rout[i] = r[i] - alfa*s[i];
		}
		betachisl = ScalarMult(rout, rout);
		betaznam = ScalarMult(r, r);
		beta = betachisl / betaznam;
		for (i = 0; i < xy.size(); i++)
		{
			z[i] = rout[i] + beta*z[i];
			r[i] = rout[i];
		}
		checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	}
}

double CalcPogreshnost(double *qVect, int dimension) {
	double sumPogr = 0;
	double sumU = 0;
	double correctSolution;
	for (size_t i = 0; i < dimension; i++)
	{
		correctSolution = AnaliticSolution(xy[i]);
		sumPogr += (qVect[i] - correctSolution)*(qVect[i] - correctSolution);
		sumU += correctSolution*correctSolution;
	}
	return sqrt(sumPogr / sumU);
}

double CalcPogreshnost(double *qVect, int dimension, string fileName) {
	ofstream output(fileName);
	double sumPogr = 0;
	double sumU = 0;
	double correctSolution;
	for (size_t i = 0; i < dimension; i++)
	{
		correctSolution = AnaliticSolution(xy[i]);
		output << setw(20) << qVect[i] << setw(20) << correctSolution << setw(20) << qVect[i] - correctSolution << setw(5) << i << endl;
		sumPogr += (qVect[i] - correctSolution)*(qVect[i] - correctSolution);
		sumU += correctSolution*correctSolution;
	}
	correctSolution = sqrt(sumPogr / sumU);
	output << correctSolution;
	return correctSolution;
}

void PrepareSolutionToDisplay() {
	double minValue, maxValue;
	minValue = maxValue = q[0];
	for (size_t i = 1; i < xy.size(); i++)
	{
		if (q[i] < minValue) minValue = q[i];
		else if (q[i]>maxValue) maxValue = q[i];
	}
	//изолинии
	double hValues = (maxValue - minValue) / (COUNT_OF_ISOLINES + 1);
	double predVal = minValue;
	for (size_t i = 0; i < COUNT_OF_ISOLINES; i++)
	{
		isolinesValues[i] = predVal + hValues;
		predVal = isolinesValues[i];
	}
	//цветовые подобласти
	hValues = (maxValue - minValue) / (COUNT_OF_COLOR_AREAS - 1);
	for (size_t i = 0; i < COUNT_OF_COLOR_AREAS; i++)
	{
		rainbow[i].value = minValue + i*hValues;
	}
	////радуга
	if (COUNT_OF_COLOR_AREAS == 7) {
		rainbow[0].red = 255; rainbow[0].green = 0; rainbow[0].blue = 0;	//каждый
		rainbow[1].red = 255; rainbow[1].green = 128; rainbow[1].blue = 0;	//охотник
		rainbow[2].red = 255; rainbow[2].green = 255; rainbow[2].blue = 0;	//желает
		rainbow[3].red = 0; rainbow[3].green = 255; rainbow[3].blue = 0;	//знать
		rainbow[4].red = 0; rainbow[4].green = 255; rainbow[4].blue = 255;	//где
		rainbow[5].red = 0; rainbow[5].green = 0; rainbow[5].blue = 255;	//сидит
		rainbow[6].red = 128; rainbow[6].green = 0; rainbow[6].blue = 128;	//фазан
	}

	if (COUNT_OF_COLOR_AREAS == 5) {
		rainbow[0].red = 255; rainbow[0].green = 0; rainbow[0].blue = 0;	//каждый
		rainbow[1].red = 255; rainbow[1].green = 255; rainbow[1].blue = 0;	//желает
		rainbow[2].red = 0; rainbow[2].green = 255; rainbow[2].blue = 0;	//знать
		rainbow[3].red = 0; rainbow[3].green = 0; rainbow[3].blue = 255;	//сидит
		rainbow[4].red = 255; rainbow[4].green = 0; rainbow[4].blue = 255;	//фазан
	}

}

int main(int argc, char *argv[])
{
	setlocale(LC_ALL, "rus");
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);

	inputLines();
	generateArrayOfCells();
	output();
	genProfile();
	GenerateEdges();
	GenerateGlobalMatrix();
	if (true) {
		LOS();
	}
	else MSG();
	CalcPogreshnost(q, xy.size(), "output.txt");
	PrepareSolutionToDisplay();

	Width = 1000;
	Height = 600;
	glutInitWindowSize(Width, Height);
	glutCreateWindow("GLUT");
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutIdleFunc(Display);
	glutMainLoop();
	return 1;
}