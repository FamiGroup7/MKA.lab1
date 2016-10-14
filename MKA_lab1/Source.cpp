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

int*ig, *jg;
double *ggl, *ggu, *di, *b, *q;

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
//�������� �������������� ����� ���������� �������
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
//������ ����� � ������� xy
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
	for (size_t j = 0; j < yTotalCells - 1; j++)
	{
		//*************����������� �� � ����� �������� �����
		for (size_t i = 0; i < xTotalCells - 1; i++)
		{
			if (!Inside(arrayCells[i][j], ignoredField)) {
				xy.push_back(arrayCells[i][j]);
			}
			else
				if (belongToLine(arrayCells[i][j], ignoredField)) {
					xy.push_back(arrayCells[i][j]);
				}

			//��������� ���. ����� ��� ������������ �������� �������
			Point dopTochka((arrayCells[i][j].x + arrayCells[i + 1][j].x) / 2, 
				(arrayCells[i][j].y + arrayCells[i + 1][j].y) / 2);
			if (!Inside(dopTochka, ignoredField)) {
				xy.push_back(dopTochka);
			}
			else
				if (belongToLine(dopTochka, ignoredField)) {
					xy.push_back(dopTochka);
				}
		}
		//��������� ��������� �� � �����
		if (!Inside(arrayCells[xTotalCells-1][j], ignoredField)) {
			xy.push_back(arrayCells[xTotalCells - 1][j]);
		}
		else
			if (belongToLine(arrayCells[xTotalCells - 1][j], ignoredField)) {
				xy.push_back(arrayCells[xTotalCells - 1][j]);
			}

		//*************��������� ������ �� �
		//*************����������� �� � ����� ������� �����
		for (size_t i = 0; i < xTotalCells - 1; i++)
		{
			Point mainPoint((arrayCells[i][j].x + arrayCells[i][j + 1].x) / 2, 
				(arrayCells[i][j].y + arrayCells[i][j + 1].y) / 2);
			if (!Inside(mainPoint, ignoredField)) {
				xy.push_back(mainPoint);
			}
			else
				if (belongToLine(mainPoint, ignoredField)) {
					xy.push_back(mainPoint);
				}

			//��������� ���. ����� ��� ������������ �������� �������
			Point nextPoint((arrayCells[i + 1][j].x + arrayCells[i + 1][j + 1].x) / 2, 
				(arrayCells[i + 1][j].y + arrayCells[i + 1][j + 1].y) / 2);
			Point dopTochka((mainPoint.x + nextPoint.x) / 2, (mainPoint.y + nextPoint.y) / 2);
			if (!Inside(dopTochka, ignoredField)) {
				xy.push_back(dopTochka);
			}
			else
				if (belongToLine(dopTochka, ignoredField)) {
					xy.push_back(dopTochka);
				}
		}
		//��������� ��������� �� � �����
		Point dopTochka((arrayCells[xTotalCells - 1][j].x + arrayCells[xTotalCells - 1][j + 1].x) / 2, 
			(arrayCells[xTotalCells - 1][j].y + arrayCells[xTotalCells - 1][j + 1].y) / 2);
		if (!Inside(dopTochka, ignoredField)) {
			xy.push_back(dopTochka);
		}
		else
			if (belongToLine(dopTochka, ignoredField)) {
				xy.push_back(dopTochka);
			}

		//*************��������� ������ �� �

	}

	//*************����������� �� � ����� ��������� �������� �������������� �����
	for (size_t i = 0; i < xTotalCells - 1; i++)
	{
		if (!Inside(arrayCells[i][yTotalCells - 1], ignoredField)) {
			xy.push_back(arrayCells[i][yTotalCells - 1]);
		}
		else
			if (belongToLine(arrayCells[i][yTotalCells - 1], ignoredField)) {
				xy.push_back(arrayCells[i][yTotalCells - 1]);
			}

		//��������� ���. ����� ��� ������������ �������� �������
		Point dopTochka((arrayCells[i][yTotalCells - 1].x + arrayCells[i + 1][yTotalCells - 1].x) / 2,
			(arrayCells[i][yTotalCells - 1].y + arrayCells[i + 1][yTotalCells - 1].y) / 2);
		if (!Inside(dopTochka, ignoredField)) {
			xy.push_back(dopTochka);
		}
		else
			if (belongToLine(dopTochka, ignoredField)) {
				xy.push_back(dopTochka);
			}
	}
	//��������� ��������� �� � �����
	if (!Inside(arrayCells[xTotalCells - 1][yTotalCells - 1], ignoredField)) {
		xy.push_back(arrayCells[xTotalCells - 1][yTotalCells - 1]);
	}
	else
		if (belongToLine(arrayCells[xTotalCells - 1][yTotalCells - 1], ignoredField)) {
			xy.push_back(arrayCells[xTotalCells - 1][yTotalCells - 1]);
		}

	//*************��������� ������ �� �

	//��������� ��
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

	//������� ������
	for (size_t i = 0; i < xTotalCells; i++)
	{
		delete arrayCells[i];
	}
	delete arrayCells;

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

void Display(void) // ������� ������
{
	int i, j, k, t, colorSreda;
	glClearColor(1, 1, 1, 1); // ������� ������
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3ub(0, 0, 0);
	glShadeModel(GL_FLAT); // ���������� ������������
	double h_x = lines[Nx - 1][0].x - lines[0][0].x;
	double h_y = lines[0][Ny - 1].y - lines[0][0].y;
	double E_x = -lines[0][0].x/*+0.03*h_x*/;
	double E_y = -lines[0][0].y/* + 0.03*h_y*/;
	double tempX,tempY;
	tempX = (Width - 50) / h_x;
	tempY = (Height - 50) / h_y;

	glColor3ub(0, 0, 0);
	glLineWidth(1);
	for (int keIndex = 0; keIndex < KE.size(); keIndex++)
	{
		glBegin(GL_LINE_LOOP);
		glVertex2f((E_x + xy[KE[keIndex].uzel[0]].x) * tempX, (E_y + xy[KE[keIndex].uzel[0]].y) * tempY);
		glVertex2f((E_x + xy[KE[keIndex].uzel[1]].x) * tempX, (E_y + xy[KE[keIndex].uzel[1]].y) * tempY);
		glVertex2f((E_x + xy[KE[keIndex].uzel[3]].x) * tempX, (E_y + xy[KE[keIndex].uzel[3]].y) * tempY);
		glVertex2f((E_x + xy[KE[keIndex].uzel[2]].x) * tempX, (E_y + xy[KE[keIndex].uzel[2]].y) * tempY);
		glEnd();
	}

	glLineWidth(2);
	glColor3ub(255, 0, 0);
	for (i = 0; i < Nx; i++)
	{
		glBegin(GL_LINE_STRIP);
		for (j = 0; j < Ny; j++)
		{
			glVertex2f((E_x + lines[i][j].x) * tempX, (E_y + lines[i][j].y) * tempY);
		}
		glEnd();
	}
	for (j = 0; j < Ny; j++)
	{
		glBegin(GL_LINE_STRIP);
		for (i = 0; i < Nx; i++)
		{
			glVertex2f((E_x + lines[i][j].x) * tempX, (E_y + lines[i][j].y) * tempY);
		}
		glEnd();
	}

	glColor3ub(0, 0, 0);
	for (i = 0; i < xy.size(); i++)
	{
					char*text = new char[3];
					//sprintf_s(text, 2, "%d", i);
					sprintf(text, "%d", i);
					drawText(text, strlen(text), (E_x + xy[i].x) * tempX, (E_y + xy[i].y) * tempY);
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
void Reshape(GLint w, GLint h) // ��� ��������� �������� ����
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

void genProfile()
{
	int KEcount = KE.size();
	ig = new int[xy.size() + 1];

	int *list[2], *listbeg;//������
	list[0] = new int[xy.size()*xy.size()];
	list[1] = new int[xy.size()*xy.size()];
	listbeg = new int[xy.size()];
	int listsize = -1;//���������� ��������� � ������, � ����� � jg

	for (int i = 0; i < xy.size(); i++)	listbeg[i] = -1;

	for (int ielem = 0; ielem < KEcount; ielem++)//�������� �� ���� ��
	{
		for (int i = 0; i < 9; i++)//���������� �������� �������
		{
			int kk = KE[ielem].uzel[i];//kk - ���������� ����� ������� ��������������� �������� �������
									//int kk = 2 * (ielem / (Nx - 1))*(2 * Nx - 1) + 2 * (ielem % (Nx - 1));
			for (int j = i + 1; j < 9; j++)//���������� ��������� �� ��� ������� �� ��������
			{
				int ind1 = kk;
				int ind2 = KE[ielem].uzel[j];
				if (ind2 < ind1)//�������� ����� �������� � �������
				{
					ind1 = ind2;
					ind2 = kk;
				}
				int iaddr = listbeg[ind2];
				if (iaddr < 0)//������ ��� ����
				{
					//�������� ������
					listsize++;
					listbeg[ind2] = listsize;
					list[0][listsize] = ind1;
					list[1][listsize] = -1;
				}
				else//������ �� ��� ����
				{
					//���� � ������ ind1
					while (list[0][iaddr] < ind1 && list[1][iaddr] >= 0)		iaddr = list[1][iaddr];
					if (list[0][iaddr] > ind1)
					{
						//���� ��������� ��� ������� ����� ������� �����, 
						//�� ����� �������� ����� ���, ����� ������ ��� �������������
						listsize++;
						list[0][listsize] = list[0][iaddr];	//������������� ������ ��������� �������
						list[1][listsize] = list[1][iaddr];
						list[0][iaddr] = ind1;			//�� ��� ����� ����� �����
						list[1][iaddr] = listsize;
					}
					else
						if (list[0][iaddr] < ind1)
						{
							//�� �����, � ������ ����������
							//��������� � ����� ������
							listsize++;
							list[1][iaddr] = listsize;
							list[0][listsize] = ind1;
							list[1][listsize] = -1;	//���������, ��� ��� ��������� ������� ������
						}
				}
			}
		}
	}
	int gk = listsize + 1;
	jg = new int[gk];
	ig[0] = 0;
	for (int i = 0; i < xy.size(); i++)
	{
		ig[i + 1] = ig[i];			//ig[i+1] - ����� ������ ������� jg, ���� ���� ��������� ��������� �������
		int iaddr = listbeg[i];		//iaddr - ������ ������ ������
		while (iaddr >= 0)			//������������� ������
		{
			jg[ig[i + 1]] = list[0][iaddr];
			ig[i + 1]++;
			iaddr = list[1][iaddr]; //��������� � ���������� ��������
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
	return p.x + p.y;
}

double AnaliticSolve(Point p) {
	return p.x + p.y;
}

double BasicFunc1d(int num, double ksi) {
	switch (num)
	{
	case 0: return 2 * (ksi - 0.5)*(ksi - 1);
	case 1: return -4 * ksi * (ksi - 1);
	case 2: return 2 * ksi * (ksi - 0.5); 
	default:
		cerr << "Error in Basic Function" << endl;
		system("pause");
		exit(1);
	}
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
					double J = alfa0 + alfa1*integrPoints[k] + alfa2*integrPoints[l];
					double dfi_I_dksi = difBasicFunc1d(i % 3, integrPoints[k])*BasicFunc1d(i / 3, integrPoints[l]);
					double dfi_I_deta = BasicFunc1d(i % 3, integrPoints[k])*difBasicFunc1d(i / 3, integrPoints[l]);
					double dfi_J_dksi = difBasicFunc1d(j % 3, integrPoints[k])*BasicFunc1d(j / 3, integrPoints[l]);
					double dfi_J_deta = BasicFunc1d(j % 3, integrPoints[k])*difBasicFunc1d(j / 3, integrPoints[l]);
					resultG +=
						((dfi_I_dksi*(beta[5] * integrPoints[k] + beta[2]) - dfi_I_deta*(beta[5] * integrPoints[l] * beta[3])) *
							(dfi_J_dksi*(beta[5] * integrPoints[k] + beta[2]) - dfi_J_deta*(beta[5] * integrPoints[l] * beta[3])) +
							(dfi_I_deta*(beta[4] * integrPoints[l] + beta[1]) - dfi_I_dksi*(beta[4] * integrPoints[k] + beta[0])) *
							(dfi_J_deta*(beta[4] * integrPoints[l] + beta[1]) - dfi_J_dksi*(beta[4] * integrPoints[k] + beta[0])))*
						tauKoefs[k] * tauKoefs[l] / J / 4;

					resultM += BasicFunc1d(i % 3, integrPoints[k])*BasicFunc1d(i / 3, integrPoints[l]) *
						BasicFunc1d(j % 3, integrPoints[k])*BasicFunc1d(j / 3, integrPoints[l])*
						tauKoefs[k] * tauKoefs[l] * J / 4;
				}
			}
			A[i][j] += resultG*signAlfa0*Lambda(KE[ielem].numberField) + resultM*signAlfa0*Gamma(KE[ielem].numberField);
		}

		//right part
		double resultRightPart = 0;
		for (size_t k = 0; k < countExtraPoints; k++)
		{
			for (size_t l = 0; l < countExtraPoints; l++)
			{
				double J = alfa0 + alfa1*integrPoints[k] + alfa2*integrPoints[l];

				resultRightPart += BasicFunc1d(i % 3, integrPoints[k])*BasicFunc1d(i / 3, integrPoints[l]) *
					Func(Point(xy[KE[ielem].uzel[0]].x + beta[1] * integrPoints[k] + beta[0] * integrPoints[l] + beta[4] * integrPoints[k] * integrPoints[l], 
						xy[KE[ielem].uzel[0]].y + beta[3] * integrPoints[k] + beta[2] * integrPoints[l] + beta[5] * integrPoints[k] * integrPoints[l]))*
					tauKoefs[k] * tauKoefs[l] * J / 4;
			}
		}
		localB[i] += resultRightPart*signAlfa0;
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
	int i, j, k, l, posI, posJ;
	double koefI, koefJ;
	for (i = 0; i < 9; i++)
	{
		b[KE[ielem].uzel[i]] += localB[i];
		for (j = 0; j <= i; j++)
		{
			AddToMatrix(KE[ielem].uzel[i], KE[ielem].uzel[j], A[i][j]);
		}
	}
}

void GenerateGlobalMatrix() {
	int ielem, i, j;
	double tKoef = sqrt(3. / 5.);
	double integrationPoints[3] = {
		(1 - tKoef) / 2,
		0.5,
		(1 + tKoef) / 2
	};
	double tauKoefs[3] = {
		5. / 9.,
		8. / 9.,
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

	//Edge2(1, 1, 1, 1);
	//Edge3(1, 1, 1, 1);
	//Edge1_sim(1, 1, 1, 1);
	//for (i = 0; i < ig[kolvoRegularNode]; i++)
	//{
	//	ggu[i] = ggl[i];
	//}
	//Edge1_not_sim(1, 1, 1, 1);
}

//	��������� ������� �� ������
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

void runLOS()
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

int main(int argc, char *argv[])
{
	setlocale(LC_ALL, "rus");
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);

	inputLines();
	generateArrayOfCells();
	output();
	//GeneratePortrait();
	//genProfile();
	//GenerateGlobalMatrix();
	//runLOS();
	//ofstream output("output.txt");
	//double sumPogr = 0;
	//for (size_t i = 0; i < xy.size(); i++)
	//{
	//	output << setw(20) << q[i] << setw(20) << AnaliticSolve(xy[i]) << setw(20) << q[i] - AnaliticSolve(xy[i]) << endl;
	//	sumPogr += (q[i] - AnaliticSolve(xy[i]))*(q[i] - AnaliticSolve(xy[i]));
	//}
	//output << sqrt(sumPogr);

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