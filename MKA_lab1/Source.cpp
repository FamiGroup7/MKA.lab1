#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <conio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "Point.cpp"
#include "Dividing.cpp"
#include "nvtr.cpp"
#include "glut.h"

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
double reverseKoef(double alfa, double alfa1) {
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


void drawText(const char *text, int length, int x, int y)
{
	glMatrixMode(GL_PROJECTION);
	double *matrix = new double[16];
	glGetDoublev(GL_PROJECTION_MATRIX, matrix);
	glLoadIdentity();
	glOrtho(0, 600, 0, 600, -5, 5);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPushMatrix();
	glLoadIdentity();
	glRasterPos2i(x, y);
	for (int i = 0; i<length; i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, (int)text[i]);//GLUT_BITMAP_9_BY_15
	}
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(matrix);
	glMatrixMode(GL_MODELVIEW);

}

void Display(void) // функция вывода
{
	int i, j, k, t, colorSreda;
	glClearColor(1, 1, 1, 1); // очистка буфера
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3ub(0, 0, 0);
	glShadeModel(GL_FLAT); // отключение интерполяции
	double h_x = lines[Nx - 1][0].x - lines[0][0].x;
	double h_y = lines[0][Ny - 1].y - lines[0][0].y;
	double E_x = -lines[0][0].x+0.03*h_x;
	double E_y = -lines[0][0].y + 0.03*h_y;
	double tempY;
	if (h_x < h_y)tempY = 550 / h_y;
	else tempY = 550 / h_x;

	glColor3ub(0, 0, 0);
	glLineWidth(1);
	for (int keIndex = 0; keIndex < KE.size(); keIndex++)
	{
		glBegin(GL_LINE_LOOP);
		glVertex2f((E_x + xy[KE[keIndex].uzel[0]].x) * tempY, (E_y + xy[KE[keIndex].uzel[0]].y) * tempY);
		glVertex2f((E_x + xy[KE[keIndex].uzel[1]].x) * tempY, (E_y + xy[KE[keIndex].uzel[1]].y) * tempY);
		glVertex2f((E_x + xy[KE[keIndex].uzel[3]].x) * tempY, (E_y + xy[KE[keIndex].uzel[3]].y) * tempY);
		glVertex2f((E_x + xy[KE[keIndex].uzel[2]].x) * tempY, (E_y + xy[KE[keIndex].uzel[2]].y) * tempY);
		glEnd();
	}

	glLineWidth(2);
	glColor3ub(255, 0, 0);
	for (i = 0; i < Nx; i++)
	{
		glBegin(GL_LINE_STRIP);
		for (j = 0; j < Ny; j++)
		{
			glVertex2f((E_x + lines[i][j].x) * tempY, (E_y + lines[i][j].y) * tempY);
		}
		glEnd();
	}
	for (j = 0; j < Ny; j++)
	{
		glBegin(GL_LINE_STRIP);
		for (i = 0; i < Nx; i++)
		{
			glVertex2f((E_x + lines[i][j].x) * tempY, (E_y + lines[i][j].y) * tempY);
		}
		glEnd();
	}

	glColor3ub(0, 0, 0);
	for (i = 0; i < xy.size(); i++)
	{
					char*text = new char[3];
					//sprintf_s(text, 2, "%d", i);
					sprintf(text, "%d", i);
					drawText(text, strlen(text), (E_x + xy[i].x) * tempY, (E_y + xy[i].y) * tempY);
	}
	//glColor3ub(255, 0, 0);
	//glPointSize(2);
	//xy.size();
	//glBegin(GL_POINTS);
	//for (i = 0; i < nX; i++)
	//{
	//	for (j = 0; j < nY; j++)
	//	{
	//		glVertex2f((E_x + xNet[i])*tempY, (E_y + yNet[j])*tempY);
	//	}
	//}
	//glEnd();

	//glPointSize(7);
	//glBegin(GL_POINTS);
	//if (koordSourceX >= leftX && koordSourceX <= rightX && koordSourceY >= leftY && koordSourceY <= rightY)
	//	glVertex2f((E_x + koordSourceX)*tempY, (E_y + koordSourceY)*tempY);
	//glEnd();

	//if (DrowText)
	//{
		//glColor3ub(0, 0, 255);
		//for (j = 0; j < nY - 1; j++)
		//{
		//	for (i = 0; i < nX; i++)
		//	{
		//		if (matrixNode[i][j] != 'Y')
		//		{
		//			char*text = new char[3];
		//			sprintf(text, "%d", NumberNode(xNet[i], yNet[j]));
		//			drawText(text, strlen(text), (E_x + xNet[i]) * tempY, (E_y + yNet[j]) * tempY);
		//		}
		//	}
	//	}
	//	for (i = 0; i < nX; i++)
	//	{
	//		if (matrixNode[i][nY - 1] != 'Y')
	//		{
	//			char*text = new char[3];
	//			sprintf(text, "%d", NumberNode(xNet[i], yNet[nY - 1]));
	//			drawText(text, strlen(text), (E_x + xNet[i]) * tempY, (E_y + yNet[nY - 1] - (yNet[nY - 1] - yNet[nY - 2]) / 20) * tempY);
	//		}
	//	}
	//}
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

int main(int argc, char *argv[])
{
	setlocale(LC_ALL, "rus");
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);

	inputLines();
	generateArrayOfCells();
	output();

	Width = 600;
	Height = 600;
	glutInitWindowSize(Width, Height);
	glutCreateWindow("GLUT");
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutIdleFunc(Display);
	glutMainLoop();
	return 1;
}