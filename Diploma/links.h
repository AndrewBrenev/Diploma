#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include<fstream>
#include <windows.h>
#include <cmath>
#include <vector>
#include <algorithm>

# pragma once;

#define IRON		1
#define OIL			2
#define GROUND		3

using namespace std;

typedef double real;


struct Point{
	//координаты, глобальный номер вершины
	real x, y, z;
	int id;
	//конструктор
	Point(){
		x = 0.0;
		y = 0.0;
		z = 0.0;
		id = 0;
	};
	Point(real p1, real p2, real p3, int k){
		x = p1;
		y = p2;
		z = p3;
		id = k;
	};
	Point(real p1, real p2, real p3){
		x = p1;
		y = p2;
		z = p3;
	};
};

struct Circle{
	//Координаты центра окружности
	Point O;
	//real x, y, z;
	//Внешний и внутренний радиусы
	real R, r;
	Circle(){};
	Circle(real x0, real y0, real z0, real R0, real r0){
		O.x = x0;
		O.y = y0;
		O.z = z0;
		R = R0;
		r = r0;
	};
	Circle(Point O, real R0, real r0){
		O.x = O.x;
		O.y = O.y;
		O.z = O.z;
		R = R0;
		r = r0;
	};
};

struct Plane{
	//координаты, глобальный номер вершины
	real A, B, C, D;
	//конструктор
	Plane(){
		A = 0.0;
		B = 0.0;
		C = 0.0;
		D = 0;
	};
	Plane(real a, real b, real c, real d){
		A = a;
		B = b;
		C = c;
		D = d;
	};

};
struct Edge {
	int a, b;

	//1 - соседи по горизонтали 2 - соседи по вертикали
	int axis;
	Edge() {};
	Edge(int i, int j, int k) { a = i; b = j; axis = k; };
};
struct Material{
	//параметры материала
	real mu, sigma;
	// номер материала
	int id;
	//конструктор
	Material(){};
};

struct NVTR{
	//массив глобальных номаров вершин КЭ
	vector <int> n;
	//Глобальный номер КЭ
	int id;
	//Номер материала
	int material;
	NVTR(){
		n.resize(8);
		n[0] = 0; n[1] = 0; n[2] = 0; n[3] = 0;
		n[4] = 0; n[5] = 0; n[6] = 0; n[7] = 0;
		id = 0; material = 0;
	};
	NVTR(int p0, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int i, int ma){
		n.resize(8);
		n[0] = p0; n[1] = p1; n[2] = p2; n[3] = p3;
		n[4] = p4; n[5] = p5; n[6] = p6; n[7] = p7;
		id = i; material = ma;
	};
};

struct nvtr_point
{
public:
	int p1, p2, p3, p4;
	int material;
	nvtr_point(int a, int b, int c, int d, int e){
		p1 = a;
		p2 = b;
		p3 = c;
		p4 = d;
		material = e;
	};
};

struct T_Point{
	// Номер добавляемой вершины
	int point_id;
	// A -  номер вершины имеющей большее значение , B -  меньшее значение
	int A, B,C,D;
	real A_value, B_value,C_value,D_value;
	// 1 - по Х, 2 - по  Z, 3 - по Y 
	int way;
	T_Point() {
		A = -1, B = -1, C = -1, D = -1;
		A_value = -1; B_value = -1; C_value = -1; D_value = -1;
		way = -1;
		point_id = -1;
	};
	T_Point(int p, int move) {
		A = -1, B = -1, C = -1, D = -1;
		A_value = -1; B_value = -1; C_value = -1; D_value = -1;
		way = move;
		point_id = p;
	};
};

struct Face_Point {
	// Номер добавляемой вершины
	int point_id;
	// соседи вершины
	int up,down;
	int left, right;
	Face_Point() {};
	Face_Point(int t) {
		point_id = t;
		up = -1; down = -1;
		left = -1; right = -1;
	};
};

struct T_Matrix{
	vector <int> ig;
	vector <int> jg;
	vector <real> gg;
	// Подпрограмма вывода в файл
	void writeMatrixInfo(){
		FILE *b = fopen("ig.txt", "w");
		//ofstream fout("ig1.txt"); 
		
		for (int i = 0; i < ig.size(); i++)
			fprintf(b, "%.1f ", (float)ig[i]);
		fclose(b);

		b = fopen("jg.txt", "w");
		for (int i = 0; i < jg.size(); i++)
			fprintf(b, "%d ", jg[i]);
		fclose(b);

		b = fopen("gg.txt", "w");
		for (int i = 0; i < gg.size(); i++)
			fprintf(b, "%lf ", gg[i]);
		fclose(b);
	}
};

struct Matrix_Elem{ int id; real val; };