#include "vect.h"


//Класс для сортировки массива точек
class sort_coord_vector
{
public:
	bool operator() (Point i, Point j)
	{
		return (i.id<j.id);
	}
};

class sort_coord_edges
{
public:
	bool operator() (Edge i, Edge j)
	{
		return (i.a < j.a);
	}
};

class sort_coord_t
{
public:
	bool operator() (T_Point i, T_Point j)
	{
		return (i.point_id < j.point_id);
	}
};
class sort_coord_faces
{
public:
	bool operator() (Face_Point i, Face_Point j)
	{
		return (i.point_id < j.point_id);
	}
};


class sort_coord_column
{
public:
	bool operator() (Matrix_Elem i, Matrix_Elem j)
	{
		return (i.id<j.id);
	}
};

class net_generator{
public:
	void buildNet();

private:
	//Класс для сортировки
	sort_coord_vector sort_coord_vect;
	sort_coord_column sort_coord_col;
	sort_coord_edges sort_coord_ed;
	sort_coord_faces sort_faces;
	sort_coord_t sort_t_points;

	vector <Point> coord;	//вектор 3-мерных точек
	vector <Point> tmp;		//Дополнительный вектор
	vector <Circle> circle_centers; //вектор центров сечений
	vector <Circle> path;	//вектор 3-мерных конечных элементов
	vector <NVTR> nvtr;		//вектор 3-мерных конечных элементов 
	vector <nvtr_point> nv; //двумерные элементы в сечении
	vector <Material> materials; //Вектор используемых материалов
	vector <vect> normals; //Вектор нормалей
	vector <Edge> edges; // Вектор рёбер
	//Кол-во узлов, и конечных элементов в сетке
	int n_uzlov, n_el;
	//количество контрольных точек, кол-во материалов
	int n_path, n_material;

	//Границы области и шаг разбиения
	real X0, Xn, Xh;
	real Y0, Yn, Yh;
	real Z0, Zn, Zh;

	// массив координат "траншеи"
	// 0,1 - X; 2,3 - Y; 4,5 - Z
	int *trench_size;

	// Число разбиений по стороне квадрата, кол-во внутренних окружностей, кол-во разбиений при повороте, кол-во разбиений в траншее
	int n, l,m,p;
	int coor_on_layer = 0, el_on_layer=0;
	//Итерация  и начальный индекс
	int iter = 0, start_ind, start_elem;

	// Коэфициент при поворотое
	real stretch_coeff;

	bool initial_zone, end_zone;
	// Режим построения сетки 1 - только труба 0 - с внешней сеткой
	bool tube_only;
	
	//Подпрограммы решения геометрич. задач
	// поиск пересечения точки прямой, заданной двумя точками и окружностью
	void circle_point(real *res_x, real *res_y, real x1, real y1, real x2, real y2, real x0, real y0, real r);
	// Нахождения точки вращения при повороте
	Point findPointOfRotation(Point A, Point B,Plane P, int inorm);
	// Вычисление определителся матрица
	real countMatrixDeterminant(real a[][3]);
	// Нахождение окружности между двумя КТ
	Circle findCircleOnLoyer(Circle Begin, Circle End, real step);
	// Проверка на перехлёст по y
	int testExternalGrid(real more, real less);

	// Генерация файлов с сеткой
	void writeIntoFile();
	//Подпрограммы вывода информации в файлы
	void writeMaterialInfo();
	void writeElemInfo();
	void writeCooordInfo();
	void writeKuzlovInfo();

	void write_special();

	//Определение размера траншеи
	void findTrench();
	
	// Считывание основных параметров сетки
	void input();
	//Подпрограммы Считывания информации о сетке
	void readMaterialInfo();
	void readElemInfo();
	void readCooordInfo();
	void readCirclePath();

	//Генерация внешней сетки
	void generatOuterNet();

	//Построениет траншеи вокруг трубы
	void generatTrenchGround(int iters);

	//перемещения внешней сетки
	void moveSection(Point A);
	//поворот слоя в направлении заданного вектора
	void rotateSection(int j);

	// Координаты трубы в сечении слоя
	void calculate2DLayer(Circle c, int norma);

	//Нахождение координат трубы в сечении
	void coordTubeOnly(Circle c);
	//Задание конечных элементов на слое
	void nvtrTubeOnly();

	// Алгоритм сборки трёхмерной сетки
	void combine3DNet();
	void compyteTubeFE(int k);
	// задание материала каждому конкретному КЭ
	void setMaterial();

	// Задание граничных условий
	void setBoundaryClause();

	//---------------------------------------------------------------------------------------
	
	T_Matrix t_matrix;
	vector <T_Point> t_points;		//Вектор глобальных номеров терминальных вершин

	vector <Face_Point > face_points;
	vector <T_Point> points_with_overlap; 

	//нахождение соседей для каждой точки
	void findFacePointNeighbour();
	//построение матрицы для торцов
	void facetPointsBuilder();
	//Вычисление значения Т-вершины
	void calculateTVertexValue(T_Point &a);
	//рекурсивное заполнение солонки Т-матрицы
	void recColomBuild(int id,real val, vector <Matrix_Elem> &column, bool horizontal);
	// Построение Т-матрицы
	void buildT_Matrix();
	// перенумерование точек
	void renumbePoints();
	// Проверка на перехлёст
	bool checkTheOverlap();
	// Исправлкение перехлёста
	void fixTheOverlap();

	//Обработка терминальных узлов
	void setPointsInfo();
	// Заполнение Т-матрицы
	void setTMatrix();

	int findPoint(Point A);
};