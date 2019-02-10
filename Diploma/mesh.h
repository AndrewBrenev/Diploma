#include "vect.h"


//����� ��� ���������� ������� �����
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
	//����� ��� ����������
	sort_coord_vector sort_coord_vect;
	sort_coord_column sort_coord_col;
	sort_coord_edges sort_coord_ed;
	sort_coord_faces sort_faces;
	sort_coord_t sort_t_points;

	vector <Point> coord;	//������ 3-������ �����
	vector <Point> tmp;		//�������������� ������
	vector <Circle> circle_centers; //������ ������� �������
	vector <Circle> path;	//������ 3-������ �������� ���������
	vector <NVTR> nvtr;		//������ 3-������ �������� ��������� 
	vector <nvtr_point> nv; //��������� �������� � �������
	vector <Material> materials; //������ ������������ ����������
	vector <vect> normals; //������ ��������
	vector <Edge> edges; // ������ ����
	//���-�� �����, � �������� ��������� � �����
	int n_uzlov, n_el;
	//���������� ����������� �����, ���-�� ����������
	int n_path, n_material;

	//������� ������� � ��� ���������
	real X0, Xn, Xh;
	real Y0, Yn, Yh;
	real Z0, Zn, Zh;

	// ������ ��������� "�������"
	// 0,1 - X; 2,3 - Y; 4,5 - Z
	int *trench_size;

	// ����� ��������� �� ������� ��������, ���-�� ���������� �����������, ���-�� ��������� ��� ��������, ���-�� ��������� � �������
	int n, l,m,p;
	int coor_on_layer = 0, el_on_layer=0;
	//��������  � ��������� ������
	int iter = 0, start_ind, start_elem;

	// ���������� ��� ���������
	real stretch_coeff;

	bool initial_zone, end_zone;
	// ����� ���������� ����� 1 - ������ ����� 0 - � ������� ������
	bool tube_only;
	
	//������������ ������� ���������. �����
	// ����� ����������� ����� ������, �������� ����� ������� � �����������
	void circle_point(real *res_x, real *res_y, real x1, real y1, real x2, real y2, real x0, real y0, real r);
	// ���������� ����� �������� ��� ��������
	Point findPointOfRotation(Point A, Point B,Plane P, int inorm);
	// ���������� ������������� �������
	real countMatrixDeterminant(real a[][3]);
	// ���������� ���������� ����� ����� ��
	Circle findCircleOnLoyer(Circle Begin, Circle End, real step);
	// �������� �� �������� �� y
	int testExternalGrid(real more, real less);

	// ��������� ������ � ������
	void writeIntoFile();
	//������������ ������ ���������� � �����
	void writeMaterialInfo();
	void writeElemInfo();
	void writeCooordInfo();
	void writeKuzlovInfo();

	void write_special();

	//����������� ������� �������
	void findTrench();
	
	// ���������� �������� ���������� �����
	void input();
	//������������ ���������� ���������� � �����
	void readMaterialInfo();
	void readElemInfo();
	void readCooordInfo();
	void readCirclePath();

	//��������� ������� �����
	void generatOuterNet();

	//����������� ������� ������ �����
	void generatTrenchGround(int iters);

	//����������� ������� �����
	void moveSection(Point A);
	//������� ���� � ����������� ��������� �������
	void rotateSection(int j);

	// ���������� ����� � ������� ����
	void calculate2DLayer(Circle c, int norma);

	//���������� ��������� ����� � �������
	void coordTubeOnly(Circle c);
	//������� �������� ��������� �� ����
	void nvtrTubeOnly();

	// �������� ������ ��������� �����
	void combine3DNet();
	void compyteTubeFE(int k);
	// ������� ��������� ������� ����������� ��
	void setMaterial();

	// ������� ��������� �������
	void setBoundaryClause();

	//---------------------------------------------------------------------------------------
	
	T_Matrix t_matrix;
	vector <T_Point> t_points;		//������ ���������� ������� ������������ ������

	vector <Face_Point > face_points;
	vector <T_Point> points_with_overlap; 

	//���������� ������� ��� ������ �����
	void findFacePointNeighbour();
	//���������� ������� ��� ������
	void facetPointsBuilder();
	//���������� �������� �-�������
	void calculateTVertexValue(T_Point &a);
	//����������� ���������� ������� �-�������
	void recColomBuild(int id,real val, vector <Matrix_Elem> &column, bool horizontal);
	// ���������� �-�������
	void buildT_Matrix();
	// ��������������� �����
	void renumbePoints();
	// �������� �� ��������
	bool checkTheOverlap();
	// ������������ ���������
	void fixTheOverlap();

	//��������� ������������ �����
	void setPointsInfo();
	// ���������� �-�������
	void setTMatrix();

	int findPoint(Point A);
};