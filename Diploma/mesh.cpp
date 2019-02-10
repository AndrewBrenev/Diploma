#include "mesh.h"


//Общий вызов функции
void net_generator::buildNet(){
	input(); 
	combine3DNet();
	if (!tube_only)	buildT_Matrix();
	writeIntoFile();
	//write_special();
	
}

//---------------------------------------------------------------------------------------
// поиск  точки пересечения прямой, заданной двумя точками, и окружностью
void net_generator::circle_point(real *res_x, real *res_y, real x1, real y1, real x2, real y2, real x0, real y0, real r){
	if (x1 == x2){
		*res_x = x1;
		*res_y = y0 + r;
	}
	if (y1 == y2){
		*res_x = x0 + r;
		*res_y = y1;
	}
	else{
		real dx = x1 - x2;
		real dy = y1 - y2;
		real l = sqrt(dx*dx + dy*dy);
		real c = dx / l;
		real s = dy / l;
		real x = x0+c*r;
		real y = y0+s*r;
		*res_x = x;
		*res_y = y;
	}
}

// Вычисление определителся матрицы 3х3
real net_generator::countMatrixDeterminant(real a[][3]){
	return  a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[0][1] * a[1][2] * a[2][0] -
		a[2][0] * a[1][1] * a[0][2] - a[2][1] * a[1][2] * a[0][0] - a[1][0] * a[0][1] * a[2][2];
}

// Нахождения точки вращения при повороте
Point net_generator::findPointOfRotation(Point Ap, Point B,Plane rotatePlane, int i){
	Point O;

	real alfa = acos((normals[i].x*normals[i + 1].x + normals[i].y*normals[i + 1].y + normals[i].z*normals[i + 1].z) /
		(sqrt(normals[i].x*normals[i].x + normals[i].y*normals[i].y + normals[i].z*normals[i].z)*
		sqrt(normals[i + 1].x*normals[i + 1].x + normals[i + 1].y*normals[i + 1].y + normals[i + 1].z*normals[i + 1].z)));
	real gamma = M_PI_2 - alfa / 2.0;
	real la = sqrt((Ap.x - path[i + 1].O.x)*(Ap.x - path[i + 1].O.x) + (Ap.y - path[i + 1].O.y)*(Ap.y - path[i + 1].O.y) + (Ap.z - path[i + 1].O.z)*(Ap.z - path[i + 1].O.z));

	real k1 = -normals[i].x*Ap.x - normals[i].y*Ap.y - normals[i].z*Ap.z;
	real k2 = -normals[i + 1].x*B.x - normals[i + 1].y*B.y - normals[i + 1].z*B.z;
	
	// Построим матрицу
	real A[3][3], b[3];
	A[0][0] = normals[i].x; A[0][1] = normals[i].y; A[0][2] = normals[i].z; b[0] = -k1;
	A[1][0] = normals[i + 1].x; A[1][1] = normals[i + 1].y; A[1][2] = normals[i+1].z; b[1] = -k2;
	A[2][0] = rotatePlane.A; A[2][1] = rotatePlane.B; A[2][2] = rotatePlane.C; b[2] = -rotatePlane.D;
	// Решим матрицу методом Крамера
	real det = countMatrixDeterminant(A);

	A[0][0] = -k1;	A[1][0] =  -k2;	A[2][0] =  -rotatePlane.D;
	real detX = countMatrixDeterminant(A);

	A[0][0] = normals[i].x; A[0][1] =  -k1;
	A[1][0] = normals[i + 1].x; A[1][1] =  -k2;
	A[2][0] = rotatePlane.A; A[2][1] =  -rotatePlane.D;
	real detY = countMatrixDeterminant(A);

	 A[0][1] = normals[i].y; A[0][2] =  -k1;
	 A[1][1] = normals[i + 1].y; A[1][2] =  -k2;
	 A[2][1] = rotatePlane.B; A[2][2] =  -rotatePlane.D;
	real detZ = countMatrixDeterminant(A);

	O.x = detX / det; O.y = detY / det; O.z = detZ / det;

	/*
	//Найдём квадрат радиуса поворота
	real R = la*tan(gamma);

	 k1 = normals[i + 1].x*B.x + normals[i + 1].y*B.y + normals[i + 1].z*B.z;
	 k2 = normals[i].x*Ap.x + normals[i].y*Ap.y + normals[i].z*Ap.z;
	real a1 = -normals[i + 1].y / normals[i + 1].x, a2 = -normals[i + 1].z / normals[i + 1].x;
	real b1 = a1 * normals[i].x, b2 = a2 * normals[i].x;
	k1 = k1 / normals[i + 1].x;
	k2 = k2 - k1*normals[i].x;
	real c = -(b2 + normals[i].z) / (b1 + normals[i].y);
	k2 = k2 / (b1 + normals[i].y);
	real k3 = k1 + a1*k2;
	real a3 = a1*c + a2;

	real e = a3*a3 + c*c + 1;
	real d = a3*(k3 - B.x) + c*(k2 - B.y) - B.z;
	real f = (k3 - B.x)*(k3 - B.x) + (k2 - B.y)*(k2 - B.y) + B.z*B.z - R*R;

	real D = d*d - e*f;
	if (D < 1e-5) D = 0;
	Point o1;
	//И таки находим точку, вокруг которой будем поворачивать трубу
	o1.z = (-d + sqrt(D)) / e;
	o1.y = c*O.z + k2;
	o1.x = a3*O.z + k3;
	

	return o1;
	*/
	return O;
} 

// Определение границ траншеи
void net_generator::findTrench(){

	trench_size = new int[6];
	int i;	real A;
	vector<real> minX, maxX;	minX.resize(n_path);	maxX.resize(n_path);
	vector<real> minY, maxY;	minY.resize(n_path);	maxY.resize(n_path);
	vector<real> minZ, maxZ;	minZ.resize(n_path);	maxZ.resize(n_path);

	real res_min_x, res_max_x;
	real res_min_y, res_max_y;
	real res_min_z, res_max_z;

	// Для кажой контрольной точки определим её границы
	for (i = 0; i < n_path; i++){
		A = 1.15*path[i].R;
		minX[i] = path[i].O.x - A;	maxX[i] = path[i].O.x + A;
		minY[i] = path[i].O.y;	maxY[i] = path[i].O.y;
		minZ[i] = path[i].O.z - A;	maxZ[i] = path[i].O.z + A;
	}

	vector<real>::iterator res_min_X, res_max_X;
	vector<real>::iterator res_min_Y, res_max_Y;
	vector<real>::iterator res_min_Z, res_max_Z;

	// По кажой координате найдём максимальное и минимальное значение
	res_min_X = min_element(minX.begin(), minX.end());
	if (res_min_X._Ptr[0] < X0) res_min_x = X0; else res_min_x = res_min_X._Ptr[0];
	res_max_X = max_element(maxX.begin(), maxX.end());
	if (res_max_X._Ptr[0] >= Xn) res_max_x = Xn; else res_max_x = res_max_X._Ptr[0];

	res_min_Y = min_element(minY.begin(), minY.end());
	if (res_min_Y._Ptr[0] < Y0) res_min_y = Y0; else res_min_y = res_min_Y._Ptr[0];
	res_max_Y = max_element(maxY.begin(), maxY.end());
	if (res_max_Y._Ptr[0] >= Yn) res_max_y = Yn; else res_max_y = res_max_Y._Ptr[0];

	res_min_Z = min_element(minZ.begin(), minZ.end());
	if (res_min_Z._Ptr[0] < Z0) res_min_z = Z0; else res_min_z = res_min_Z._Ptr[0];
	res_max_Z = max_element(maxZ.begin(), maxZ.end());
	if (res_max_Z._Ptr[0] >= Zn) res_max_z = Zn; else res_max_z = res_max_Z._Ptr[0];

	trench_size[0] = (res_min_x - X0) / Xh;
	trench_size[1] = (res_max_x - X0) / Xh;
	if (res_max_x - X0 != trench_size[1] * Xh) trench_size[1]++;

	trench_size[2] = (res_min_y - Y0) / Yh;
	trench_size[3] = (res_max_y - Y0) / Yh;
	if (res_max_y - Y0 != trench_size[3] * Yh) trench_size[3]++;

	trench_size[4] = (res_min_z - Z0) / Zh;
	trench_size[5] = (res_max_z - Z0) / Zh;
	if (res_max_z - Z0 != trench_size[5] * Zh) trench_size[5]++;
	
}

Circle net_generator::findCircleOnLoyer(Circle Begin, Circle End, real step){
	Circle Res;
	Res.O.y = step;
	real ratio = (step - Begin.O.y) / (End.O.y - Begin.O.y);
	Res.O.x = ratio * (End.O.x - Begin.O.x) + Begin.O.x;
	Res.O.z = ratio * (End.O.z - Begin.O.z) + Begin.O.z;
	Res.R = ratio * ((End.O.z + End.R) - (Begin.O.z + Begin.R)) + (Begin.O.z + Begin.R) - Res.O.z;
	Res.r = ratio * ((End.O.z + End.r) - (Begin.O.z + Begin.r)) + (Begin.O.z + Begin.r) - Res.O.z;
	return Res;
}

int net_generator::testExternalGrid(real less,real more ){
	bool flag = true;
	real y_st= Y0;
	int i;
	for ( i = 0; flag && i < trench_size[3] - trench_size[2]; i++)
		if (Y0 + i*Yh > less && Y0 + i*Yh < more )
			flag = false;
	i--;
	if (flag) return -1; else return i;

}

//---------------------------------------------------------------------------------------

//Считывание входных данных
void net_generator::input(){
	//считывание координат
	FILE* file = fopen("input.txt", "r");
	fscanf(file, "X : ( %lf ; %lf ) h = %lf\n", &X0, &Xn, &Xh);
	fscanf(file, "Y : ( %lf ; %lf ) h = %lf\n", &Y0, &Yn, &Yh);
	fscanf(file, "Z : ( %lf ; %lf ) h = %lf\n", &Z0, &Zn, &Zh);
	fscanf(file, "n = %d l = %d p = %d m = %d k = %lf\n", &n, &l, &p, &m, &stretch_coeff);
	fscanf(file, "%d\n", &tube_only);
	readCirclePath();

	if (path[0].O.y > Y0) initial_zone = true;
	else initial_zone = false;

	if (path[n_path - 1 ].O.y < Yn) end_zone = true;
	else end_zone = false;
}
void net_generator::readCirclePath(){
	//считывание координат
	FILE* file = fopen("Path.txt", "r");
	fscanf(file, "%d\n", &n_path);
	path.resize(n_path);
	for (int i = 0; i < n_path; i++)
		fscanf(file, " ( %lf ; %lf ; %lf ) R = %lf r = %lf\n", &path[i].O.x, &path[i].O.y, &path[i].O.z, &path[i].R, &path[i].r);

	//Создание вектора нормалей
	normals.resize(n_path - 1);
	for (int i = 1; i < n_path; i++){
		bool Z, X,Y;
		if (path[i].O.x - path[i - 1].O.x >= 0) X = true; else X = false;
		if (path[i].O.z - path[i - 1].O.z >= 0) Z = true; else Z = false;
		if (path[i].O.y - path[i - 1].O.y > 0) Y = true; else Y = false;
		vect temp(path[i].O.x - path[i - 1].O.x,
			path[i].O.y - path[i - 1].O.y,
			path[i].O.z - path[i - 1].O.z, Z, X, Y);
		normals[i-1]=temp;
	}
		
}

//---------------------------------------------------------------------------------------
// Подпрораммы вывода информации в бинарные файлы 
void net_generator::writeIntoFile(){
	writeMaterialInfo();
	writeElemInfo();
	writeCooordInfo();
	writeKuzlovInfo();
}
void net_generator::writeElemInfo(){

	//vector<int> tmp; tmp.resize(n_el*14);
	size_t t;
	int *tmp=new int[n_el*14];
	for (int i = 0,k=0; k < n_el;k++, i=i+14)
		for (int j = 0; j < 8; j++)
		tmp[i+j] = nvtr[k].n[j];
	
	//прежде следует вычленить информацию
	FILE *fout;
	fopen_s(&fout, "nver.dat", "wb");
	t=fwrite(tmp, sizeof(int), n_el*14, fout);
	fclose(fout);
}
void net_generator::writeMaterialInfo(){

	int *tmp = new int[n_el];
	size_t t;
	for (int i = 0; i < n_el; i++)	tmp[i ] = nvtr[i].material;
	FILE *fout;
	fopen_s(&fout, "nvkat.dat", "wb");
    t=fwrite(tmp, sizeof(int), n_el, fout);
	fclose(fout);

}
void net_generator::writeCooordInfo(){

	double *tmp = new double[3 * n_uzlov];
	size_t t;
	for (int i = 0,k=0; k < n_uzlov; k++, i=i+3){
		tmp[i] = coord[k].x;
		tmp[i + 1] = coord[k].y;
		tmp[i + 2] = coord[k].z;
	}
	FILE *fout;
	fopen_s(&fout, "xyz.dat", "wb");
	t=fwrite(tmp, sizeof(double), 3*n_uzlov, fout);
	fclose(fout);

}
void net_generator::writeKuzlovInfo(){

	FILE *b = fopen("inftry.dat", "w");
	fprintf(b, " ISLAU=       0 INDKU1=       0 INDFPO=       0\nKUZLOV=%8d   KPAR=%8d    KT1=       0   KTR2=       0   KTR3=       0\nKISRS1=       0 KISRS2=       0 KISRS3=       0   KBRS=       0\n   KT7=       0   KT10=       0   KTR4=       0  KTSIM=       0\n   KT6=       0\n", n_uzlov, n_el);
	fclose(b);

}

void net_generator::write_special(){
	int fst = 0;
		//краевые
		FILE *ba = fopen("output/l1.txt", "w");
		for (int i = 0; i < coord.size(); i++)
			if (coord[i].x == X0 || coord[i].x == Xn ||
				coord[i].y == Y0 || coord[i].y == Yn ||
				coord[i].z == Z0 || coord[i].z == Zn){
					fprintf(ba, "%d\n", coord[i].id);
					fst++;
			}
		fclose(ba);
		//осн инфо
		FILE *b = fopen("output/info.txt", "w");
		fprintf(b, "KUZLOV = %d KPAR = %d KT1 = %d\n", n_uzlov, n_el,fst);
		fclose(b);
		//координаты
		FILE *a = fopen("output/nodes.txt", "w");
		for (int i = 0; i < coord.size(); i++)
			fprintf(a, "%lf %lf %lf\n", coord[i].x, coord[i].y, coord[i].z);
		fclose(a);
		//элементы
		FILE *af = fopen("output/elements.txt", "w");
		for (int i = 0; i < nvtr.size(); i++)
			fprintf(af, "%d %d %d %d %d %d %d %d\n", nvtr[i].n[0], nvtr[i].n[1], nvtr[i].n[2], nvtr[i].n[3], 
														nvtr[i].n[4], nvtr[i].n[5], nvtr[i].n[6], nvtr[i].n[7]);
		fclose(af);
		// материалы
		FILE *m = fopen("output/materials.txt", "w");
		for (int i = 0; i < nvtr.size(); i++)
			fprintf(m, "%d\n", nvtr[i].material);
		fclose(m);

		//Т-матрица
		b = fopen("output/t_ig.txt", "w");
		for (int i = 0; i < t_matrix.ig.size(); i++)
			fprintf(b, "%d ", t_matrix.ig[i]);
		fclose(b);

		b = fopen("output/t_jg.txt", "w");
		for (int i = 0; i < t_matrix.jg.size(); i++)
			fprintf(b, "%d ", t_matrix.jg[i]);
		fclose(b);

		b = fopen("output/t_gg.txt", "w");
		for (int i = 0; i < t_matrix.gg.size(); i++)
			fprintf(b, "%lf ", t_matrix.gg[i]);
		fclose(b);
}

//---------------------------------------------------------------------------------------
// Подпрограммы считывания информации о сетки из бинарных файлов
void net_generator::readMaterialInfo(){
	int *tmp1 = new int[n_el];
	size_t t;
	FILE *fout1;
	fopen_s(&fout1, "nvkat.dat", "rb");
	t = fread(tmp1, sizeof(int), n_el, fout1);
	fclose(fout1);
}
void net_generator::readCooordInfo(){
	real *tmp1 = new real[3 * n_uzlov];
	size_t t1;
	FILE *fout1;
	fopen_s(&fout1, "xyz.dat", "rb");
	t1 = fread(tmp1, sizeof(double), 3 * n_uzlov, fout1);
	int k = feof(fout1);
	fclose(fout1);
}
void net_generator::readElemInfo(){
	FILE *fout1;	size_t t;
	fopen_s(&fout1, "nver.dat", "rb");
	int  *tmp1 = new int[n_el * 14];
	t = fread(tmp1, sizeof(int), n_el * 14, fout1);
}

//---------------------------------------------------------------------------------------
//Построение 3-х мерной сетки
void net_generator::combine3DNet(){

	// Определим размер траншеи
	findTrench();
	
	// номер нормы с которой работаем
	int i_n;
	// Шаг по прямой
	int j = 0;
	// флаг окончания движения по данной прямой
	bool end = false;
	for (int i = 0; i < n_path - 1; i++){

		end = false;
		//Вычисляем шаг на данном отрезке
		real dx = (path[i + 1].O.x - path[i].O.x) / n;
		real dy = (path[i + 1].O.y - path[i].O.y) / n;
		real dz = (path[i + 1].O.z - path[i].O.z) / n;

		while (!end){

			// Вычисляется точка на прямой
			real nx = j*dx;
			real ny = j*dy;
			real nz = j*dz;
			vect tmp(nx, ny, nz, normals[i].clockwiseZ, normals[i].clockwiseX, normals[i].clockwiseY);

			if (normals[i].length() - tmp.length() <= stretch_coeff * path[i + 1].R && i != n_path - 2)
				// Если попали в поворот
			{
				
				real alfa = acos((normals[i].x*normals[i + 1].x + normals[i].y*normals[i + 1].y + normals[i].z*normals[i + 1].z) /
					(sqrt(normals[i].x*normals[i].x + normals[i].y*normals[i].y + normals[i].z*normals[i].z)*
					sqrt(normals[i + 1].x*normals[i + 1].x + normals[i + 1].y*normals[i + 1].y + normals[i + 1].z*normals[i + 1].z)));

				//Находим точки на прямых
				real Ax, Ay, Az;
				real Bx, By, Bz;

				Ax = path[i + 1].O.x - (normals[i].x * stretch_coeff * path[i + 1].R) / normals[i].length();
				Ay = path[i + 1].O.y - (normals[i].y * stretch_coeff * path[i + 1].R) / normals[i].length();
				Az = path[i + 1].O.z - (normals[i].z * stretch_coeff * path[i + 1].R) / normals[i].length();

				Bx = path[i + 1].O.x + (normals[i + 1].x * stretch_coeff * path[i + 1].R) / normals[i + 1].length();
				By = path[i + 1].O.y + (normals[i + 1].y * stretch_coeff * path[i + 1].R) / normals[i + 1].length();
				Bz = path[i + 1].O.z + (normals[i + 1].z * stretch_coeff * path[i + 1].R) / normals[i + 1].length();
				
				Point A(Ax, Ay, Az), B(Bx, By, Bz),O;

				// Выисляем плоскость, в которой они лежат
				real Nx, Ny, Nz,D;
				real a21, a22, a23;
				real a31, a32, a33;
				a21 = Ax - path[i + 1].O.x; a22 = Ay - path[i + 1].O.y; a23 = Az - path[i + 1].O.z;
				a31 = Bx - path[i + 1].O.x; a32 = By - path[i + 1].O.y; a33 = Bz - path[i + 1].O.z;

				// Вектор нормали к плосткости, в которой происходит поворот
				Nx = a22*a33 - a32*a23;
				Ny = a23*a31 - a21*a33;
				Nz = a21*a32 - a22*a31;
				D = -path[i + 1].O.x*Nx - path[i + 1].O.y*Ny - path[i + 1].O.z*Nz;

				Plane rotatePlane(Nx, Ny, Nz, D);
				
				//Находим точку поворота
				O = findPointOfRotation(A, B, rotatePlane, i);
				
				//Нормируем ветор нормали для поворота точек вокруг произвольного единичного вектора
				real N_length = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
				Nx = - Nx / N_length;
				Ny = - Ny / N_length;
				Nz = -Nz / N_length;

				//Вычисляем шаг по углу
				real alfa_step = alfa / m;

				Circle Res_old, Res;
				Res_old.O.x = path[i].O.x; Res_old.O.y = path[i].O.y; Res_old.O.z = path[i].O.z;
				Res_old.R = path[i].R; Res_old.r = path[i].r;
				
				//Нашли крайние точки, теперь реализуем сам поворот
				for (int section = 0; section <= m; section++){

					real aop = section*alfa_step*180.0 / M_PI;
					real sn = sin(section*alfa_step);
					real cs = cos(section*alfa_step);

					// Повернули центр окружности
					Res.O.x = (Ax - O.x) *(cs + (1 - cs)*Nx*Nx) + (Ay - O.y)*((1 - cs)*Nx*Ny - sn*Nz) + (Az - O.z)*((1 - cs)*Nx*Nz + sn*Ny) + O.x;
					Res.O.y = (Ax - O.x) *((1 - cs)*Ny*Nx + sn*Nz) + (Ay - O.y)*(cs + (1 - cs)*Ny*Ny) + (Az - O.z)*((1 - cs)*Ny*Nz - sn*Nx) + O.y;
					Res.O.z = (Ax - O.x) *((1 - cs)*Nz*Nx - sn*Ny) + (Ay - O.y)*((1 - cs)*Nz*Ny + sn*Nx) + (Az - O.z)*(cs + (1 - cs)*Nz*Nz) + O.z;
					Res.R = path[i + 1].R; Res.r = path[i + 1].r;

					// Определили направление поворота
					bool Z, X,Y(true);
					if (Res.O.x - Res_old.O.x >= 0) X = true; else X = false;
					if (Res.O.z - Res_old.O.z >= 0) Z = true; else Z = false;
					//Добавляем вектор, в направлении которого движимся
					vect norml(Res.O.x - Res_old.O.x, Res.O.y - Res_old.O.y, Res.O.z - Res_old.O.z, Z, X, Y);
					normals.push_back(norml);

					{
						// Если попали в стыковку с внешней сеткой
						int n_centers = circle_centers.size() - 1;
						int extra_loyer = testExternalGrid(circle_centers[n_centers].O.y, Res.O.y);
						if (extra_loyer != -1)
						{
							Circle prev_circle = circle_centers[n_centers];
							Circle oneMore = findCircleOnLoyer(prev_circle, Res, Y0 + extra_loyer* Yh);
							// Если мы на пути следования трубы
							calculate2DLayer(oneMore, normals.size() - 1);
							iter++;
						}
					}
					//Строим слой
					calculate2DLayer(Res, normals.size() - 1);

					Res_old = Res;
					iter++;
				}
					// Определяем с какого в следующем
		
				dx = (path[i + 2].O.x - path[i + 1].O.x) / n;
				dy = (path[i + 2].O.y - path[i + 1].O.y) / n;
				dz = (path[i + 2].O.z - path[i + 1].O.z) / n;
					
				real n_ny = abs(path[i + 2].O.y - path[i + 1].O.y) / n;
				bool j_f = true;

				for (int w = 0; w < n && j_f; w++){
					// Вычисляется точка на прямой
					nx = w*dx;
					ny = w*dy;
					nz = w*dz;
					vect tmp1(nx, ny, nz, normals[i].clockwiseZ, normals[i].clockwiseX, normals[i].clockwiseY);
					if (tmp1.length() >= stretch_coeff * path[i + 1].R ){
							j = w;
							j_f = false;
						}
					}
					end = true;
				}
				else
					// Если движемся по прямой
				{
					//Если мы вошли в радиус поворота
					if (i == n_path - 2 && normals[i].length() - tmp.length() < stretch_coeff * path[i + 1].R) 
						end = true;
					else{
						
						// находим нужную точку
						Circle Res(nx + path[i].O.x, ny + path[i].O.y, nz + path[i].O.z, path[i].R, path[i].r);
						if (iter){
							// Если попали в стыковку с внешней сеткой
							int n_centers = circle_centers.size() - 1;
							
							int extra_loyer = testExternalGrid(circle_centers[n_centers].O.y, Res.O.y);
							if (extra_loyer != -1)
							{
								Circle prev_circle= circle_centers[n_centers];
								Circle oneMore = findCircleOnLoyer(prev_circle, Res, Y0 + extra_loyer* Yh);
								// Если мы на пути следования трубы
								calculate2DLayer(oneMore, i);
								iter++;
							}
							
						}
						// Если мы на пути следования трубы
						calculate2DLayer(Res, i);

						// На первой итерации находим трубу
						if (iter == 0)
						{
							// Узнали сколько точек на одном слое
							coor_on_layer = coord.size();
							// Строим сету трубы на слое
							nvtrTubeOnly();
							el_on_layer = nv.size();

						}
						iter++;
					}					
					j++;
				}
			}
		}

	{
		// Если попали в стыковку с внешней сеткой
		int n_centers = circle_centers.size() - 1;
		int extra_loyer = testExternalGrid(circle_centers[n_centers].O.y, path[n_path - 1].O.y);
		if (extra_loyer != -1)
		{
			Circle prev_circle = circle_centers[n_centers];
			Circle oneMore = findCircleOnLoyer(prev_circle, path[n_path - 1], Y0 + extra_loyer* Yh);
			// Если мы на пути следования трубы
			calculate2DLayer(oneMore, normals.size() - 1);
			iter++;
		}
					}
	//Отдельно обработаем последнюю точку
	calculate2DLayer(path[n_path - 1], n_path - 2);
	iter++;
	//считаем все остальные точки
	compyteTubeFE(iter);
	start_ind = coord.size();
	start_elem = nvtr.size();
	if (!tube_only){
		//Высчитываем траншею
		generatTrenchGround(iter);

		//Считаем внешнюю сетку
		generatOuterNet();

		//задаём материалы
		setMaterial();
	}
	//Сортируем точки
	std::sort(coord.begin(), coord.end(), sort_coord_vect);

	n_uzlov = coord.size();
	n_el = nvtr.size();

	
}

// Построение сетки на слое
void net_generator::calculate2DLayer(Circle c, int norma){
	//Вычисляем
	coordTubeOnly(c);
	Point A(c.O.x, c.O.y, c.O.z, c.R);
	circle_centers.push_back(c);
	rotateSection(norma);
	moveSection(A);
	for (int i = 0; i < tmp.size(); i++)  coord.push_back(tmp[i]); 

	//Если труба не выходит на границу, а находится в среде
	if(c.O.y == path[0].O.y && initial_zone || c.O.y == path[n_path - 1].O.y && end_zone)
		for (int i = 0; i < tmp.size(); i++) {
			Face_Point A(tmp[i].id);
			face_points.push_back(A);
		}
	tmp.clear();
};

// построение КЭ, зная чило слоёв
void net_generator::compyteTubeFE(int k){
	for (int q = 1; q < k; q++)
		for (int j = 0; j < el_on_layer; j++)
		{
		NVTR tmp_FE(
			nv[j].p1 + (q - 1)*coor_on_layer + 1,
			nv[j].p2 + (q - 1)*coor_on_layer + 1,
			nv[j].p3 + (q - 1)*coor_on_layer + 1,
			nv[j].p4 + (q - 1)*coor_on_layer + 1,
			nv[j].p1 + q *coor_on_layer + 1,
			nv[j].p2 + q *coor_on_layer + 1,
			nv[j].p3 + q * coor_on_layer + 1,
			nv[j].p4 + q *coor_on_layer + 1,
			(q - 1)*coor_on_layer + j + 1,
			nv[j].material
			);
		nvtr.push_back(tmp_FE);
		// Добавим рёбра
		if (q == 1 && initial_zone || q == k - 1 && end_zone) {
			int c = 0;
			if (q == k - 1 && end_zone) c = 4;
			Edge A(tmp_FE.n[0 + c] - 1, tmp_FE.n[1 + c] - 1, 2);
			Edge B(tmp_FE.n[2 + c] - 1, tmp_FE.n[3 + c] - 1, 2);
			Edge C(tmp_FE.n[0 + c] - 1, tmp_FE.n[2 + c] - 1, 1);
			Edge D(tmp_FE.n[1 + c] - 1, tmp_FE.n[3 + c] - 1, 1);
			edges.push_back(A); edges.push_back(B); edges.push_back(C); edges.push_back(D);

		}
		}
}

//Нахождение координат трубы в сечении
void net_generator::coordTubeOnly(Circle c){

	Point Temp(0, 0, 0, 0);
	int i, j;
	real a = c.R*1.1;
	real b = c.r*0.6;
	real a_step = 2 * a / (int)n;
	real b_step = 2 * b / (int)n;
	real p = a / a_step;

	real step = (c.R - c.r) / (int)l;

	//Точки на внутреннем и внешнем квадратах
	for (int k = 0; k < n; k++)
	{
		//Верхняя сторона
		// Точка с внутреннего квадрата
		Temp.x = b - k*b_step;
		Temp.z = b;
		Temp.id = (l + 1)*(k+1) + k + iter*coor_on_layer;
		tmp.push_back(Temp);

		//Левая сторона
		// Точка с внутреннего квадрата
		Temp.x = - b;
		Temp.z =  b - k*b_step;
		Temp.id = (l + 2)*(n + k)  + l + 1 + iter*coor_on_layer;
		tmp.push_back(Temp);

		//Нижняя сторона
		// Точка с внутреннего квадрата
		Temp.x = - b + k*b_step;
		Temp.z = - b;
		Temp.id = (k + 2 * n) * (l + 2) + l + 1 + +iter*coor_on_layer;
		tmp.push_back(Temp);

		//Правая сторона
		// Точка с внутреннего квадрата
		Temp.x =  b;
		Temp.z = - b + k*b_step;
		Temp.id = (k + 3 * n) * (l + 2) + l + 1 + +iter*coor_on_layer;
		tmp.push_back(Temp);

	}
	//Точки на радиусах
	for (int k = 0; k <= p; k++)
	{
		// Если попали ровно в середину стороны
		if (k == p && (a - p*a_step <= 2))
			for (i = 0; i <= l; i++)
			{
			//Вертикальная верх
			Temp.x = 0;
			Temp.z =  c.R - i*step;
			Temp.id = (l + 2)*k + i  + iter*coor_on_layer;
			tmp.push_back(Temp);
			//Горизонтальная лево
			Temp.x = - c.R + i*step;
			Temp.z = 0;
			Temp.id = (l + 2) * n + (l + 2)*k + i + +iter*coor_on_layer;
			tmp.push_back(Temp);
			//Вертикальная низ
			Temp.x = 0;
			Temp.z = - c.R + i*step;
			Temp.id = 2 * n *(l + 2) + (l + 2)*k + i + +iter*coor_on_layer;
			tmp.push_back(Temp);
			//Горизонтальная право
			Temp.x =  c.R - i*step;
			Temp.z = 0;
			Temp.id = 3 * n*(l + 2) + (l + 2)*k + i + +iter*coor_on_layer;
			tmp.push_back(Temp);
			}
		else{

			real x1, y1, x2, y2;
			//Точка с внешнего квадрата
			x1 = a - k*a_step;
			y1 = a;

			// Точка с внутреннего квадрата
			x2 = b - k*b_step;
			y2 = b;

			real r_x, r_z, R_x, R_z;
			circle_point(&R_x, &R_z, x1, y1, x2, y2, 0.0, 0.0, c.R);
			circle_point(&r_x, &r_z, x1, y1, x2, y2, 0.0, 0.0, c.r);

			real x_step = abs(R_x - r_x) / l;
			real z_step = abs(R_z - r_z) / l;

			//концентрические окружности
			for (i = 0; i <= l; i++){

				//Верх право
				Temp.x = R_x - i*x_step;
				Temp.z = R_z - i*z_step;
				Temp.id = (l + 2)*k + i + iter*coor_on_layer;
				tmp.push_back(Temp);

				real dx = Temp.x;
				real dy = Temp.z;
				Point Temp1(0, 0, 0, 0);

				//Верх лево
				if (k != 0){
					Temp1.x = -dx;
					Temp1.z = Temp.z;
					Temp1.id = (l + 2)*n - (l + 2)*k + i + iter*coor_on_layer;
					tmp.push_back(Temp1);
				}

				//низ лево
				Temp1.x = -dx;
				Temp1.z = -dy;
				Temp1.id = (l + 2)*n * 2 + (l + 2)*k + i  + iter*coor_on_layer;
				tmp.push_back(Temp1);

				if (k != 0){
					//низ право
					Temp1.x = Temp.x;
					Temp1.z = -dy;
					Temp1.id = (l + 2)*n * 3 - (l + 2)*k + i + iter*coor_on_layer;
					tmp.push_back(Temp1);
				}

				//Право верх
				if (k != 0){
					Temp1.x = dy;
					Temp1.z = dx;
					Temp1.id = (l + 2)* n * 4 + i - (l + 2)*k + iter*coor_on_layer;
					tmp.push_back(Temp1);
				}

				//Право низ
				Temp1.x = dy;
				Temp1.z = -dx;
				Temp1.id = (l + 2)*n * 3 + i + (l + 2)*k + iter*coor_on_layer;
				tmp.push_back(Temp1);

				//Лево верх
				Temp1.x = -dy;
				Temp1.z = dx;
				Temp1.id = (l + 2)*n + i + (l + 2)*k + iter*coor_on_layer;
				tmp.push_back(Temp1);

				if (k != 0){
					//Лево низ
					Temp1.x = -dy;
					Temp1.z = -dx;
					Temp1.id = (l + 2)*n * 2 + i  - (l + 2)*k + iter*coor_on_layer;
					tmp.push_back(Temp1);
				}

			}
		}
	}
	if (iter == 0)  start_ind = tmp.size() - 1;
	//Добавляем внутренюю сетку
	for (i = 1; i < n; i++)
		for (j = 1; j < n; j++){
		Temp.x = b - j*b_step;
		Temp.z = b - i*b_step;
		Temp.id = start_ind + (n - 1)*(i - 1) + j + iter*coor_on_layer;
		tmp.push_back(Temp);
		}

}
//Задание двумерных КЭ трубы
void net_generator::nvtrTubeOnly(){
	//Займёмся сеткой
	int i, j, material = 0;
	int a, b, c, d;
	//Сетка по окружностям
	for (int k = 0; k < 4 * n; k++)
		for (i = 0; i <= l; i++)
		{
			if (i == l ) material = OIL; else material = IRON;
		//Склейка конца с началом
		if (k == 4 * n - 1){
			a = (l + 2)*k + i;
			b = (l + 2)*k + i + 1;
			c = i;
			d = i + 1;
			nvtr_point A(c, a, d, b, material);
			nv.push_back(A);
		}
		else {
			a = (l + 2)*k + i;
			b = (l + 2)*k + i + 1;
			c = (l + 2)*(k + 1) + i;
			d = (l + 2)*(k + 1) + i + 1;
			if (k / n == 0 || k / n == 2)
			{
				nvtr_point A(a, b, c, d, material); nv.push_back(A);
			}
			else
			{
				nvtr_point A(b, d, a, c, material); nv.push_back(A);
			}
			
		}
		}
	//Верхняя полоса внутренней сетки
	for (i = 0; i < n - 1; i++){
		a = (l + 2)*i + l + 1;
		b = (l + 2)*(i + 1) + l + 1;
		c = start_ind + i;
		d = start_ind + i + 1;
		nvtr_point A(a, c, b, d, OIL);
		nv.push_back(A);
	}
	//Верхняя левая в квадрате
	a = (l + 2)*n - 1;
	b = (l + 2)*n + l + 1;
	c = start_ind + n - 1;
	d = (l + 2)*(n + 1) + l + 1;
	nvtr_point A1(a, c, b, d, OIL);
	nv.push_back(A1);
	//Нижняя левая
	a = (l + 2)*(2 * n) - 1;
	b = coor_on_layer - 1;
	c = (l + 2)*(2 * n) + l + 1;
	d = (l + 2)*(2 * n + 2) - 1;
	nvtr_point A2(b, d,a,c, OIL);
	nv.push_back(A2);
	//Нижняя правая
	a = (l + 2)*(3 * n) - 1;
	b = coor_on_layer - n + 1;
	c = (l + 2)*(3 * n + 1) - 1;
	d = (l + 2)*(3 * n + 2) - 1;
	nvtr_point A3(d,c,b,a, OIL);
	nv.push_back(A3);

	for (i = 1; i < n - 1; i++)
	{
		// по левой стенке
		a = (l + 2)*(n + i) + l + 1;
		b = start_ind + i*(n - 1);
		c = (l + 2)*(n + i + 1) + l + 1;
		d = start_ind + (i + 1)*(n - 1);
		nvtr_point A4(b,d,a,c, OIL);
		nv.push_back(A4);
		//по низу
		a = (l + 2)*(2 * n + 1 + i) - 1;
		b = coor_on_layer - i;
		c = (l + 2)*(2 * n + 2 + i) - 1;
		d = coor_on_layer - i - 1;
		nvtr_point A5(d,c,b,a, OIL);
		nv.push_back(A5);
		// По правой стенке
		a = start_ind + (n - 1)*(i - 1) + 1;
		b = start_ind - (i - 1)*(l + 2);
		c = start_ind + (n - 1)*i + 1;
		d = start_ind - (i)*(l + 2);
		nvtr_point A6(b,d,a,c, OIL);
		nv.push_back(A6);
	}
	//Внутренние КЭ
	for (i = 1; i < n - 1; i++)
		for (j = 1; j < n - 1; j++){
		a = start_ind + (n - 1)*(i - 1) + j;
		b = start_ind + (n - 1)*(i - 1) + j + 1;
		c = start_ind + (n - 1)*(i)+j;
		d = start_ind + (n - 1)*i + j + 1;
		nvtr_point A(a, c,b, d, OIL);
		nv.push_back(A);
		}
}

//---------------------------------------------------------------------------------------
// Поворот заданных точек на ветор
void net_generator::rotateSection(int j){
	Point normal(0, 1, 0, 0);
	for (int i = 0; i < tmp.size(); i++)
		tmp[i] = normals[j].rotatePoint(tmp[i], normal);
}
// Сдвиг вектора точек на вектор
void net_generator::moveSection(Point A){
	for (int i = 0; i < tmp.size(); i++){
		tmp[i].x += A.x;
		tmp[i].y += A.y;
		tmp[i].z += A.z;
	}
}

//---------------------------------------------------------------------------------------
//Построение внешней сетки
void net_generator::generatOuterNet(){

	Point T; NVTR N;
	vector <NVTR> n_tmp;
	
	int a;
	int id = coord.size(), start_id = id;
	int nX = (Xn - X0) / Xh;
	int nY = (Yn - Y0) / Yh;
	int nZ = (Zn - Z0) / Zh;
	int el_layer = (nZ + 1)*(nX + 1);
	int t_el = 0;
	// Движение по OY
	for (int j = 0; j <= nY; j++)
		// Движени по OZ
		for (int k = 0; k <= nZ; k++)
			//Движение по OX
			for (int i = 0; i <= nX; i++)
			{
				T.id = id;
				T.x = X0 + i*Xh ;
				T.y = Y0 + j*Yh;
				T.z = Z0 + k*Zh;
				tmp.push_back(T);
				id++;

				if (i != nX && j != nY && k != nZ)
				{
					NVTR N(
						start_id + (k + 1) *(nX + 1) + (i + 1) + j * el_layer + 1,
						start_id + k * (nX + 1) + (i + 1) + j * el_layer + 1,
						start_id + (k + 1) * (nX + 1) + i + j * el_layer + 1,
						start_id + k*(nX + 1) + i + j * el_layer + 1,
						start_id + (k + 1) *(nX + 1) + (i + 1) + (j+1) * el_layer + 1,
						start_id + k * (nX + 1) + (i + 1) + (j + 1) * el_layer + 1,
						start_id + (k + 1) * (nX + 1) + i + (j + 1) * el_layer + 1,
						start_id + k*(nX + 1) + i + (j + 1) * el_layer + 1,
						t_el,
						GROUND
						);
					n_tmp.push_back(N);
					t_el++;
				}
			}
	int n_del = start_id + trench_size[2] * el_layer + trench_size[4] * (nX+1) + trench_size[0];
	bool flag = true;

	//Вырежем ненужные

	for (int i = 0; i < tmp.size(); i++){
		//Если данная точка лежит на границе, либо в области 
		if ((tmp[i].id >= n_del)
			&& (tmp[i].x >= X0 + trench_size[0] * Xh) && (tmp[i].x <= X0 + trench_size[1] * Xh)
			&& (tmp[i].y >= Y0 + trench_size[2] * Yh) && (tmp[i].y <= Y0 + trench_size[3] * Yh)
			&& (tmp[i].z >= Z0 + trench_size[4] * Zh) && (tmp[i].z <= Z0 + trench_size[5] * Zh))
		{
			// Флаг для удаления КЭ
			flag = true;
		
			//Если данный элемент не граничный
			if ((tmp[i].x != X0 + trench_size[1] * Xh) && (tmp[i].z != Z0 + trench_size[5] * Zh) && (tmp[i].y != Y0 + trench_size[3] * Yh))
			{
				// Для всех конечных элементов будет лишь один КЭ с этой вершиной
				for (int j = 0; j < n_tmp.size() && flag; j++)
					if (n_tmp[j].n[3] == tmp[i].id + 1)
					{
					//Удоляем данный КЭ из списка КЭ
					n_tmp.erase(n_tmp.begin() + j);
					flag = false;
					}
			}
				bool flag_double = true;
				int numer;
				// Находим вершину с той же координатой
				for (numer = 0; numer < coord.size() && flag_double; numer++)
					if ((tmp[i].x == coord[numer].x) &&
						(tmp[i].y == coord[numer].y) &&
						(tmp[i].z == coord[numer].z))
						flag_double = false;
				numer--;
				// Заменяем в списке КЭ вершину на данную
				if (!flag_double){
					for (int p = 0; p < n_tmp.size(); p++)
						for (int t = 0; t < 8; t++)
							if (n_tmp[p].n[t] == tmp[i].id + 1)
								n_tmp[p].n[t] = coord[numer].id + 1;

					tmp[i].id = -1;
				}
		}
	}
	// Меняем индексы на новые остальным вершинам и КЭ
	int id_ad = 0;
	for (int i = 0; i < tmp.size(); i++)
		if (tmp[i].id != -1 && tmp[i].id > n_del)
		{
			for (int p = 0; p < n_tmp.size(); p++)
				for (int t = 0; t < 8; t++)
					if (n_tmp[p].n[t] == tmp[i].id + 1)
						n_tmp[p].n[t] = n_del + id_ad + 1;
			tmp[i].id = n_del + id_ad;
			id_ad++;
		}
	//Добавляем нужные
	for (int i = 0; i < tmp.size(); i++){
		if (tmp[i].id != -1)
		{
			
		T_Point overlap_point;
		// Добавляем X-овые
		if (
			(tmp[i].x > X0 + trench_size[0] * Xh) && (tmp[i].x < X0 + trench_size[1] * Xh)
			&& (tmp[i].y >= Y0 + trench_size[2] * Yh) && (tmp[i].y <= Y0 + trench_size[3] * Yh)
			&& ((tmp[i].z == Z0 + trench_size[4] * Zh) || (tmp[i].z == Z0 + trench_size[5] * Zh))){

			overlap_point.point_id = tmp[i].id;
			overlap_point.way = 1;
			points_with_overlap.push_back(overlap_point);
		}
		// Добавляем Z-овые
		if (
			((tmp[i].x == X0 + trench_size[0] * Xh) || (tmp[i].x == X0 + trench_size[1] * Xh))
			&& (tmp[i].y >= Y0 + trench_size[2] * Yh) && (tmp[i].y <= Y0 + trench_size[3] * Yh)
			&& (tmp[i].z > Z0 + trench_size[4] * Zh) && (tmp[i].z < Z0 + trench_size[5] * Zh)){

			overlap_point.point_id = tmp[i].id;
			overlap_point.way = 2;
			points_with_overlap.push_back(overlap_point);
		}
		//// Добавляем Y-овые
		//if (
		//	((tmp[i].x == X0 + trench_size[0] * Xh) || (tmp[i].x == X0 + trench_size[1] * Xh))
		//	&& (tmp[i].y > Y0 + trench_size[2] * Yh) && (tmp[i].y < Y0 + trench_size[3] * Yh)
		//	&& ((tmp[i].z == Z0 + trench_size[4] * Zh) || (tmp[i].z == Z0 + trench_size[5] * Zh))){

		//	overlap_point.point_id = tmp[i].id;
		//	overlap_point.way = 3;
		//	points_with_overlap.push_back(overlap_point);
		//	
		//}
		 coord.push_back(tmp[i]);
	}
	}
	for (int i = 0; i < n_tmp.size(); i++)
		nvtr.push_back(n_tmp[i]);
}

void net_generator::generatTrenchGround(int iters){

	Point A,B,C,TMP;
	T_Point t_pnt;
	int start_id = coord.size(), id = 0;
	int coord_layer;
	real b = 0.6 * path[0].r;
	real a_x_step = Xh*(trench_size[1] - trench_size[0]) / n;
	real a_z_step = Zh*(trench_size[5] - trench_size[4]) / n;
	real b_step =   2* b/ n;
	real dx, dz, s, c;

	// Посчитаем все точки и построем сетку на всех слоях
	for (int i = 0; i < circle_centers.size(); i++){

		//Верхняя сторона
		for (int k = 0; k < n; k++)
		{
			// Точка с внешнего квадрата
			A.x = X0 + trench_size[1] * Xh - k*a_x_step;
			A.z = Z0 + trench_size[5] * Zh;

			// Точка с внутреннего квадрата
			C.x = circle_centers[i].O.x + b - k*b_step;
			C.z = circle_centers[i].O.z + b;
			
			c = (A.x - C.x) / sqrt((A.x - C.x)*(A.x - C.x) + (A.z - C.z)*(A.z - C.z));
			s = (A.z - C.z) / sqrt((A.x - C.x)*(A.x - C.x) + (A.z - C.z)*(A.z - C.z));

			// Точка с радиуса
			B.x = circle_centers[i].O.x + c*circle_centers[i].R;
			B.z = circle_centers[i].O.z + s*circle_centers[i].R;

			dx = (A.x - B.x) / (p + 1);
			dz = (A.z - B.z) / (p + 1);

			for (int j = 0; j <= p; j++){
				if (!j) {
					if (k) {
						TMP.x = A.x - dx;
						t_pnt.way = 1;
					}
					else
					{
						t_pnt.way = 3;
						TMP.x = A.x;
					}
					t_pnt.point_id = start_id + id;
					t_points.push_back(t_pnt);
				}
				else
					TMP.x = A.x - j*dx;
				TMP.y = circle_centers[i].O.y;
				TMP.z = A.z - j*dz;
				TMP.id = start_id + id;
				coord.push_back(TMP);
			
				if (j && (!i && initial_zone || i == circle_centers.size() - 1 && end_zone))
				{
					Face_Point A(TMP.id);
					face_points.push_back(A);
				}
				id++;
			}
		}
		//Левая сторона
		for (int k = 0; k < n; k++)
		{
			// Точка с внешнего квадрата
			A.x = X0 + trench_size[0] * Xh;
			A.z = Z0 + trench_size[5] * Zh - k*a_z_step;
			// Точка с внутреннего квадрата
			C.x = circle_centers[i].O.x - b;
			C.z = circle_centers[i].O.z + b - k*b_step;

			c = (C.x - A.x) / sqrt((A.x - C.x)*(A.x - C.x) + (A.z - C.z)*(A.z - C.z));
			s = (A.z - C.z) / sqrt((A.x - C.x)*(A.x - C.x) + (A.z - C.z)*(A.z - C.z));

			// Точка с радиуса
			B.x = circle_centers[i].O.x - c*circle_centers[i].R;
			B.z = circle_centers[i].O.z + s*circle_centers[i].R;

			dx = (A.x - B.x) / (p + 1);
			dz = (A.z - B.z) / (p + 1);

			for (int j = 0; j <= p; j++){
				if (!j) {
					if (k) {
						TMP.z = A.z - dz;
						t_pnt.way = 2;
					}
					else {
						t_pnt.way = 3;
						TMP.z = A.z;
					}

					t_pnt.point_id = start_id + id;
					t_points.push_back(t_pnt);
				}
				else
					TMP.z = A.z - j*dz;
				TMP.x = A.x - j*dx;
				TMP.y = circle_centers[i].O.y;
				TMP.id = start_id + id;
				coord.push_back(TMP);
				id++;

				if (j && (!i && initial_zone || i == circle_centers.size() - 1 && end_zone))
				{
					Face_Point A(TMP.id);
					face_points.push_back(A);
				}
			}
		}
		// Нижняя сторона
		for (int k = 0; k < n; k++)
		{
			// Точка с внешнего квадрата
			A.x = X0 + trench_size[0] * Xh + k*a_x_step;
			A.z = Z0 + trench_size[4] * Zh;
			// Точка с внутреннего квадрата
			C.x = circle_centers[i].O.x - b + k*b_step;
			C.z = circle_centers[i].O.z - b;

			c = (A.x - C.x) / sqrt((A.x - C.x)*(A.x - C.x) + (A.z - C.z)*(A.z - C.z));
			s = (A.z - C.z) / sqrt((A.x - C.x)*(A.x - C.x) + (A.z - C.z)*(A.z - C.z));

			// Точка с радиуса
			B.x = circle_centers[i].O.x + c*circle_centers[i].R;
			B.z = circle_centers[i].O.z + s*circle_centers[i].R;

			dx = (A.x - B.x) / (p + 1);
			dz = (A.z - B.z) / (p + 1);

			for (int j = 0; j <= p; j++){

				if (!j) {
					if (k) {
						TMP.x = A.x - dx;
						t_pnt.way = 1;
					}
					else {
						t_pnt.way = 3;
						TMP.x = A.x;
					}

					t_pnt.point_id = start_id + id;
					t_points.push_back(t_pnt);

				}
				else
					TMP.x = A.x - j*dx;
				TMP.y = circle_centers[i].O.y;
				TMP.z = A.z - j*dz;
				TMP.id = start_id + id;
				coord.push_back(TMP);
				id++;

				if (j && (!i && initial_zone || i == circle_centers.size() - 1 && end_zone))
				{
					Face_Point A(TMP.id);
					face_points.push_back(A);
				}
			}
		}
		//Правая сторона
		for (int k = 0; k < n; k++)
		{
			// Точка с внешнего квадрата
			A.x = X0 + trench_size[1] * Xh;
			A.z = Z0 + trench_size[4] * Zh + k*a_z_step;
			// Точка с внутреннего квадрата
			C.x = circle_centers[i].O.x + b;
			C.z = circle_centers[i].O.z - b + k*b_step;

			c = (A.x - C.x) / sqrt((A.x - C.x)*(A.x - C.x) + (A.z - C.z)*(A.z - C.z));
			s = (A.z - C.z) / sqrt((A.x - C.x)*(A.x - C.x) + (A.z - C.z)*(A.z - C.z));

			// Точка с радиуса
			B.x = circle_centers[i].O.x + c*circle_centers[i].R;
			B.z = circle_centers[i].O.z + s*circle_centers[i].R;

			dx = (A.x - B.x) / (p + 1);
			dz = (A.z - B.z) / (p + 1);

			for (int j = 0; j <= p; j++){
				if (!j) {
					if (k) {
						TMP.z = A.z - dz;
						t_pnt.way = 2;
					}
					else{ 
						t_pnt.way = 3; 
						TMP.z = A.z;
					}

					t_pnt.point_id = start_id + id;
					t_points.push_back(t_pnt);
				}
				else
					TMP.z = A.z - j*dz;
				TMP.x = A.x - j*dx;
				TMP.y = circle_centers[i].O.y;
				TMP.id = start_id + id;
				coord.push_back(TMP);
				id++;

				if (j && (!i && initial_zone || i == circle_centers.size() - 1 && end_zone))
				{
					Face_Point A(TMP.id);
					face_points.push_back(A);
				}
			}
		}
		if (i == 0) coord_layer = id;
	}

	//Построим сетку

	int material = 0;
	id = 1;
	//Сетка по окружностям
	for (int q = 1; q < iters; q++)
	for (int k = 0; k < 4 * n; k++){
		material = GROUND;
		for (int i = 0; i < p; i++)
		{ 
			//Склейка конца с началом
			if (k == 4 * n - 1){
				NVTR tmp_FE(
					start_id + (p + 1)*k + i + (q - 1)*coord_layer + 1,
					start_id + (p + 1)*k + i + 1 + (q - 1)*coord_layer + 1,
					start_id + i + (q - 1)*coord_layer + 1,
					start_id + i + 1 + (q - 1)*coord_layer + 1,
					start_id + (p + 1)*k + i + q *coord_layer + 1,
					start_id + (p + 1)*k + i + 1 + q *coord_layer + 1,
					start_id + i + q * coord_layer + 1,
					start_id + i + 1 + q *coord_layer + 1,
					nvtr[nvtr.size() - 1].id + id,
					material
					);
				nvtr.push_back(tmp_FE);
				id++;
				// Добавим рёбра
				if (q == 1 && initial_zone || q == iters - 1 && end_zone) {
					int c = 0;
					if (q == iters - 1 && end_zone) c = 4;
					Edge A(tmp_FE.n[0 + c] - 1, tmp_FE.n[1 + c] - 1, 1);
					Edge B(tmp_FE.n[2 + c] - 1, tmp_FE.n[3 + c] - 1, 1);
					Edge C(tmp_FE.n[0 + c] - 1, tmp_FE.n[2 + c] - 1, 2);
					Edge D(tmp_FE.n[1 + c] - 1, tmp_FE.n[3 + c] - 1, 3);
					edges.push_back(A); edges.push_back(B); edges.push_back(C); edges.push_back(D);
				}
			}
			else{
				NVTR tmp_FE(
					start_id + (p + 1)*k + i + (q - 1)*coord_layer + 1,
					start_id + (p + 1)*k + i + 1 + (q - 1)*coord_layer + 1,
					start_id + (p + 1)*(k + 1) + i + (q - 1)*coord_layer + 1,
					start_id + (p + 1)*(k + 1) + i + 1 + (q - 1)*coord_layer + 1,
					start_id + (p + 1)*k + i + q *coord_layer + 1,
					start_id + (p + 1)*k + i + 1 + q *coord_layer + 1,
					start_id + (p + 1)*(k + 1) + i + q * coord_layer + 1,
					start_id + (p + 1)*(k + 1) + i + 1 + q *coord_layer + 1,
					nvtr[nvtr.size() - 1].id + id,
					material
					);
				nvtr.push_back(tmp_FE);
				id++;
				// Добавим рёбра
				if (q == 1 && initial_zone || q == iters - 1  && end_zone) {
					int c = 0;
					if (q == iters - 1 && end_zone) c = 4;
					if (k / n == 0 || k / n == 2) {
						Edge A(tmp_FE.n[0 + c] - 1, tmp_FE.n[1 + c] - 1, 2);
						Edge B(tmp_FE.n[2 + c] - 1, tmp_FE.n[3 + c] - 1, 2);
						Edge C(tmp_FE.n[0 + c] - 1, tmp_FE.n[2 + c] - 1, 1);
						Edge D(tmp_FE.n[1 + c] - 1, tmp_FE.n[3 + c] - 1, 1);
						edges.push_back(A); edges.push_back(B); edges.push_back(C); edges.push_back(D);
					}
					else {
						Edge A(tmp_FE.n[0 + c] - 1, tmp_FE.n[1 + c] - 1, 1);
						Edge B(tmp_FE.n[2 + c] - 1, tmp_FE.n[3 + c] - 1, 1);
						Edge C(tmp_FE.n[0 + c] - 1, tmp_FE.n[2 + c] - 1, 2);
						Edge D(tmp_FE.n[1 + c] - 1, tmp_FE.n[3 + c] - 1, 2);
						edges.push_back(A); edges.push_back(B); edges.push_back(C); edges.push_back(D);
					}
				}
			}
		}
		//Стыковка траншеи с трубой
		//Склейка конца с началом
		if (k == 4 * n - 1){
			NVTR tmp_FE(
				start_id + (p + 1)*k + p + (q - 1)*coord_layer + 1,
				(l + 2)*k + (q - 1)*coor_on_layer + 1,
				start_id + p + (q - 1)*coord_layer + 1,
				(q - 1)*coor_on_layer + 1,
				start_id + (p + 1)*k + p + q *coord_layer + 1,
				(l + 2)*k + q *coor_on_layer + 1,
				start_id + p + q * coord_layer + 1,
				q *coor_on_layer + 1,
				nvtr[nvtr.size() - 1].id + id,
				material
				);
			nvtr.push_back(tmp_FE);
			id++;
			//Добавим рёбра
			if (q == 1 && initial_zone || q == iters - 1 && end_zone) {
				int c = 0;
				if (q == iters - 1 && end_zone) c = 4;
				Edge A(tmp_FE.n[0 + c] - 1, tmp_FE.n[1 + c] - 1, 1);
				Edge B(tmp_FE.n[2 + c] - 1, tmp_FE.n[3 + c] - 1, 1);
				Edge C(tmp_FE.n[0 + c] - 1, tmp_FE.n[2 + c] - 1, 2);
				Edge D(tmp_FE.n[1 + c] - 1, tmp_FE.n[3 + c] - 1, 3);
				edges.push_back(A); edges.push_back(B); edges.push_back(C); edges.push_back(D);
			}
		}
		else{
			NVTR tmp_FE(
				start_id + (p + 1)*k + p + (q - 1)*coord_layer + 1,
				(l + 2)*k + (q - 1)*coor_on_layer + 1,
				start_id + (p + 1)*(k + 1) + p + (q - 1)*coord_layer + 1,
				(l + 2)*(k + 1) + (q - 1)*coor_on_layer + 1,
				start_id + (p + 1)*k + p + q *coord_layer + 1,
				(l + 2)*k + q *coor_on_layer + 1,
				start_id + (p + 1)*(k + 1) + p + q * coord_layer + 1,
				(l + 2)*(k + 1) + q *coor_on_layer + 1,
				nvtr[nvtr.size() - 1].id + id,
				material
				);
			nvtr.push_back(tmp_FE);
			id++;
			//Добавим рёбра
			if (q == 1 && initial_zone || q == iters - 1 && end_zone) {
				int c = 0;
				if (q == iters - 1 && end_zone) c = 4;
				if (k / n == 0 || k / n == 2) {
					Edge A(tmp_FE.n[0 + c] - 1, tmp_FE.n[1 + c] - 1, 2);
					Edge B(tmp_FE.n[2 + c] - 1, tmp_FE.n[3 + c] - 1, 2);
					Edge C(tmp_FE.n[0 + c] - 1, tmp_FE.n[2 + c] - 1, 1);
					Edge D(tmp_FE.n[1 + c] - 1, tmp_FE.n[3 + c] - 1, 1);
					edges.push_back(A); edges.push_back(B); edges.push_back(C); edges.push_back(D);
				}
				else {
					Edge A(tmp_FE.n[0 + c] - 1, tmp_FE.n[1 + c] - 1, 1);
					Edge B(tmp_FE.n[2 + c] - 1, tmp_FE.n[3 + c] - 1, 1);
					Edge C(tmp_FE.n[0 + c] - 1, tmp_FE.n[2 + c] - 1, 2);
					Edge D(tmp_FE.n[1 + c] - 1, tmp_FE.n[3 + c] - 1, 2);
					edges.push_back(A); edges.push_back(B); edges.push_back(C); edges.push_back(D);
				}
			}
		}

	}
	
	// Добавим "большой" элемент в конец траншеи
	if (path[0].O.y >Y0 + Yh * trench_size[2]){
		Point A(X0 + trench_size[0] * Xh, Y0 + trench_size[2] * Yh, Z0 + trench_size[4] * Zh, coord.size()),
			B(X0 + trench_size[1] * Xh, Y0 + trench_size[2] * Yh, Z0 + trench_size[4] * Zh, coord.size() + 1),
			C(X0 + trench_size[0] * Xh, Y0 + trench_size[2] * Yh, Z0 + trench_size[5] * Zh, coord.size() + 2),
			D(X0 + trench_size[1] * Xh, Y0 + trench_size[2] * Yh, Z0 + trench_size[5] * Zh, coord.size() + 3);
		coord.push_back(A);	coord.push_back(B);	coord.push_back(C);	coord.push_back(D);
		int q, w, e, r;
		Point tmp1(X0 + trench_size[0] * Xh, path[0].O.y, Z0 + trench_size[4] * Zh);
		q = findPoint(tmp1);
		Point tmp2(X0 + trench_size[1] * Xh, path[0].O.y, Z0 + trench_size[4] * Zh);
		w = findPoint(tmp2);
		Point tmp3(X0 + trench_size[0] * Xh, path[0].O.y, Z0 + trench_size[5] * Zh);
		e = findPoint(tmp3);
		Point tmp4(X0 + trench_size[1] * Xh, path[0].O.y, Z0 + trench_size[5] * Zh);
		r = findPoint(tmp4);
		NVTR tmp_FE(
			D.id+1,B.id+1,C.id+1,A.id+1,
			r + 1, w + 1, e + 1, q + 1,
			nvtr[nvtr.size() - 1].id + 1,
			GROUND
			);
		nvtr.push_back(tmp_FE);
	}
	
	// Добавим "большой" элемент в начало траншеи
	if (path[n_path - 1].O.y < Y0 + Yh * trench_size[3]){
		Point A(X0 + trench_size[0] * Xh, Y0 + trench_size[3] * Yh, Z0 + trench_size[4] * Zh, coord.size()),
			B(X0 + trench_size[1] * Xh, Y0 + trench_size[3] * Yh, Z0 + trench_size[4] * Zh, coord.size() + 1),
			C(X0 + trench_size[0] * Xh, Y0 + trench_size[3] * Yh, Z0 + trench_size[5] * Zh, coord.size() + 2),
			D(X0 + trench_size[1] * Xh, Y0 + trench_size[3] * Yh, Z0 + trench_size[5] * Zh, coord.size() + 3);
		coord.push_back(A);	coord.push_back(B);	coord.push_back(C);	coord.push_back(D); int q, w, e, r;
		Point tmp1(X0 + trench_size[0] * Xh, path[n_path - 1].O.y, Z0 + trench_size[4] * Zh);
		q = findPoint(tmp1);
		Point tmp2(X0 + trench_size[1] * Xh, path[n_path - 1].O.y, Z0 + trench_size[4] * Zh);
		w = findPoint(tmp2);
		Point tmp3(X0 + trench_size[0] * Xh, path[n_path - 1].O.y, Z0 + trench_size[5] * Zh);
		e = findPoint(tmp3);
		Point tmp4(X0 + trench_size[1] * Xh, path[n_path - 1].O.y, Z0 + trench_size[5] * Zh);
		r = findPoint(tmp4);
		NVTR tmp_FE(
			r + 1, w + 1, e + 1, q + 1,
			D.id + 1, B.id + 1, C.id + 1, A.id + 1,
			nvtr[nvtr.size() - 1].id + 1,
			GROUND
			);
		nvtr.push_back(tmp_FE);
	}
	
	
}

void net_generator::setMaterial(){
	for (int i = 0; i < nvtr.size(); i++)
		if (nvtr[i].material == GROUND && coord[nvtr[i].n[2] - 1].z > 10)
			nvtr[i].material = 4;
	
}

// Задание граничных условий
void net_generator::setBoundaryClause(){
	FILE *bound = fopen("boundary1.txt", "w");
	for (int i = 0; i < coord.size(); i++)
		if (coord[i].x == X0 || coord[i].x == Xn ||
			coord[i].y == Y0 || coord[i].y == Yn ||
			coord[i].z == Z0 || coord[i].z == Zn)
			fprintf(bound, "%d ", coord[i].id);
	fclose(bound);
}

//---------------------------------------------------------------------------------------
// Подпрограмма построения Т-матрицы
void net_generator::buildT_Matrix(){
	
	if (checkTheOverlap())
		fixTheOverlap();
	
	findFacePointNeighbour();
	setPointsInfo();
	renumbePoints();
	setTMatrix();
	t_matrix.writeMatrixInfo();
}
// Подпрограмма перенумерации точек
void net_generator::renumbePoints(){
	int nc = n_uzlov - t_points.size();
	int id_t = 0,id_n_t =0;

	// Сохраним текущие ячейки
	vector <NVTR> tmp_nvtr;
	// Сохраним текущие терминальные узлы
	vector <T_Point> tmp_points = t_points;
	tmp_nvtr = nvtr;
	int j = 0;
	int new_id;
	// Для всех узлов
	for (int i = 0; i < n_uzlov; i++){

		bool flag = true;
		// Проверяем нахождение узла в списке терминальных
		for (j = id_t; flag && j < t_points.size(); j++)
			if (coord[i].id == t_points[j].point_id) flag = false;
		//Даём новый индекс не терминальным узлам
		if (flag){
			new_id = id_n_t;
			id_n_t++;
		}
		//если узел терминальный
		else { 
			new_id = nc + id_t;
			id_t++;
		}
		// Заменяем данную вершину в списке ячеек
		for (int p = 0;  p < tmp_nvtr.size(); p++)
			for (int t = 0; t < 8; t++)
				if (tmp_nvtr[p].n[t] == coord[i].id + 1){
					nvtr[p].n[t] = new_id + 1;
				}
		for (int i = 0; i < edges.size(); i++) {
			if (edges[i].a == coord[i].id) edges[i].a = new_id;
			if (edges[i].b == coord[i].id) edges[i].b = new_id;
		}

		for (int l = 0; l < t_points.size(); l++){
			if (tmp_points[l].point_id == coord[i].id) t_points[l].point_id = new_id;
			if (tmp_points[l].A == coord[i].id) t_points[l].A = new_id;
			if (tmp_points[l].B == coord[i].id) t_points[l].B = new_id;
			if (tmp_points[l].C == coord[i].id) t_points[l].C = new_id;
			if (tmp_points[l].D == coord[i].id) t_points[l].D = new_id;

		}
		coord[i].id = new_id;
	}
	//Сортируем точки
	std::sort(coord.begin(), coord.end(), sort_coord_vect);
}

void net_generator::findFacePointNeighbour()
{
	// Выполним swap, чтобы a - меньший номер 
	for (int i = 0; i < edges.size(); i++)
		if (edges[i].b < edges[i].a) {
			int c = edges[i].b;
			edges[i].b = edges[i].a;
			edges[i].a = c;
		}
	/*
	std::sort(edges.begin(), edges.end(), sort_coord_ed);*/
	//Удалим повторяющиеся рёбра
	vector<Edge> tmp_edge;
	for (int i = 0; i < edges.size(); i++) {
		//Флаг означает, что ребро повторяющееся
		bool flag = false;
		for (int j = 0; j < tmp_edge.size() && !flag; j++) {
			if ((edges[i].a == tmp_edge[j].a && edges[i].b == tmp_edge[j].b ||
				edges[i].b == tmp_edge[j].a && edges[i].a == tmp_edge[j].b )&& tmp_edge[j].axis == edges[i].axis)
				flag = true;
		}
		if (!flag) tmp_edge.push_back(edges[i]);
	}

	edges.clear();
	for (int j = 0; j < tmp_edge.size(); j++) 	edges.push_back(tmp_edge[j]);
	// Теперь повторяющихся рёбер нет
	std::sort(edges.begin(), edges.end(), sort_coord_ed);
	std::sort(face_points.begin(), face_points.end(), sort_faces);

	for (int i = 0; i < face_points.size(); i++) {
		for (int j = 0; j < edges.size(); j++) {
			//Определим какая наша, какая соседняя
			int our_id, neib_id;
			if (face_points[i].point_id == edges[j].a || face_points[i].point_id == edges[j].b) {
				if (face_points[i].point_id == edges[j].a)
				{
					our_id = edges[j].a;
					neib_id = edges[j].b;
				}
				else
				{
					our_id = edges[j].b;
					neib_id = edges[j].a;
				}

				bool flag = false;
				//соседи по горизонтали
				if (edges[j].axis == 1)
					if (coord[our_id].x < coord[neib_id].x) {
						if (face_points[i].right == -1) {
							face_points[i].right = neib_id;
							flag = true;
						}
					}
					else {
						if (face_points[i].left == -1) {
							face_points[i].left = neib_id;
							flag = true;
						}
					}
				//соседи по верткали
				if (edges[j].axis == 2)
					if (coord[our_id].z < coord[neib_id].z) {
						if (face_points[i].up == -1) {
							face_points[i].up = neib_id;
							flag = true;
						}
					}
					else {
						if (face_points[i].down == -1) {
							face_points[i].down = neib_id;
							flag = true;
						}
					}
				if (!flag) {
					if (coord[our_id].x < coord[neib_id].x && face_points[i].right == -1)
					{
						face_points[i].right = neib_id;
						flag = true;
					}
					if (coord[our_id].x > coord[neib_id].x && !flag && face_points[i].left == -1)
					{
						face_points[i].left = neib_id;
						flag = true;
					}
					if (coord[our_id].z > coord[neib_id].z && !flag && face_points[i].down == -1)
					{
						face_points[i].down = neib_id;
						flag = true;
					}
					if (coord[our_id].z < coord[neib_id].z && !flag && face_points[i].up == -1)
					{
						face_points[i].up = neib_id;
						flag = true;
					}
				}
			}
		}
	}
	/*
	//На момент коммента - 0

	// Посчитаем чило вершин с проблемами
	int problems(0); 
	vector<Face_Point> pr;
	for (int i = 0; i < face_points.size(); i++)
		if (face_points[i].left == -1 || face_points[i].right == -1 ||
			face_points[i].up == -1 || face_points[i].down == -1) {
			problems++;
			pr.push_back(face_points[i]);
		}
	problems++;
	*/
}

// Подпрограмма проверки на перехлёст
bool net_generator::checkTheOverlap(){ 
	if (points_with_overlap.size()){
		for (int i = 0; i < points_with_overlap.size(); i++){

			switch (points_with_overlap[i].way)
			{
				// движение по Х
			case 1:
			{
				real a = (trench_size[1] - trench_size[0] )* Xh / n;
				bool flag = true;
				int j = 0;
				for (j = 0; flag && j < n; j++){
					if (coord[points_with_overlap[i].point_id].x >= X0 + trench_size[0] * Xh + j*a &&
						coord[points_with_overlap[i].point_id].x <= X0 + trench_size[0] * Xh + (j + 1)*a )
						flag = false;
				}
				int steps;
				flag = true;
				if (coord[points_with_overlap[i].point_id].z == Z0 + trench_size[4] * Zh) { steps = trench_size[0]; j--; flag = false; }
				else{ steps = trench_size[1];  j = n - j; }
				
				Point A(X0 + steps*Xh, coord[points_with_overlap[i].point_id].y, coord[points_with_overlap[i].point_id].z);
				if (flag){
					points_with_overlap[i].A = findPoint(A) + (p + 1)*j;
					points_with_overlap[i].B = points_with_overlap[i].A + p + 1;
				}
				else{
					points_with_overlap[i].B = findPoint(A) + (p + 1)*j;
					points_with_overlap[i].A = points_with_overlap[i].B + p + 1;
				}
				break;
			}
				// движение по Z
			case 2:
			{
				real  a = (trench_size[5] - trench_size[4])*Zh/n;
				bool flag = true;
				int j = 0;
				for (j = 0; flag && j < n; j++){
					if (coord[points_with_overlap[i].point_id].z >= Z0 + trench_size[4] * Zh + j*a &&
						coord[points_with_overlap[i].point_id].z <= Z0 + trench_size[4] * Zh + (j + 1)*a)
						flag = false;
				}
				flag = true;
				int steps;
				if (coord[points_with_overlap[i].point_id].x == X0 + trench_size[0] * Xh) { steps = trench_size[5]; j = n - j; flag = false; }
				else{ steps = trench_size[4];  j--; 
				}

				Point A(coord[points_with_overlap[i].point_id].x, coord[points_with_overlap[i].point_id].y, Z0 + steps*Zh );
				if (!flag){
					points_with_overlap[i].A = findPoint(A) + (p + 1)*j;
					points_with_overlap[i].B = points_with_overlap[i].A + p + 1;
				}
				else{
					points_with_overlap[i].B = findPoint(A) + (p + 1)*j;
					points_with_overlap[i].A = points_with_overlap[i].B + p + 1;
				}
				break;
			}
				// движение по Y
			/*case 3:
			{
				
				bool flag = true;
				int j = 0;
				for (j = 0; flag && j < circle_centers.size(); j++){
					if (coord[points_with_overlap[i].point_id].y >= circle_centers[j].y &&
						coord[points_with_overlap[i].point_id].y <= circle_centers[j + 1].y)
						flag = false;
				}
				j--;
				Point A(coord[points_with_overlap[i].point_id].x, circle_centers[j + 1].y, coord[points_with_overlap[i].point_id].z);
				points_with_overlap[i].A = findPoint(A);

				Point B(coord[points_with_overlap[i].point_id].x, circle_centers[j].y, coord[points_with_overlap[i].point_id].z);
				points_with_overlap[i].B = findPoint(B);

				break;
			}*/
			default:
				break;
			}
		}
		return true;
	}else
	return false;
}

// Подпрограмма устранеия перехлёста
void net_generator::fixTheOverlap(){
	for (int i = 0; i < points_with_overlap.size(); i++){
		bool flag = true;
		int p;
		int l;

		// Определим элемент, в котором лежит прямая, на которой находится точка перехлёста
		flag = true;
		for (p = start_elem; flag && p < nvtr.size(); p++){
			for (int a = 0; flag && a < 7; a++)
				for (int b = a + 1; flag && b < 8; b++)
					if (nvtr[p].n[a] == points_with_overlap[i].A + 1 && nvtr[p].n[b] == points_with_overlap[i].B + 1 ||
						nvtr[p].n[b] == points_with_overlap[i].A + 1 && nvtr[p].n[a] == points_with_overlap[i].B + 1)
						flag = false;
		} p--;

		switch (points_with_overlap[i].way)
		{
		case 1:
			// Определим вершину, которую нужно изменить
			flag = true;
			for (l = 0; flag && l < 8; l++)
				if (nvtr[p].n[l] != points_with_overlap[i].B + 1 &&
					coord[nvtr[p].n[l] - 1].x == coord[points_with_overlap[i].B].x &&
					coord[nvtr[p].n[l] - 1].y == coord[points_with_overlap[i].B].y
					) flag = false;
			l -- ;
			// Заменим данную вершину
			coord[points_with_overlap[i].B].x = coord[points_with_overlap[i].point_id].x;
			coord[nvtr[p].n[l]- 1 ].x = coord[points_with_overlap[i].point_id].x;

			break;
		case 2:
			// Определим вершину, которую нужно изменить
			flag = true;
			for (l = 0; flag && l < 8; l++)
				if (nvtr[p].n[l] != points_with_overlap[i].B + 1 &&
					coord[nvtr[p].n[l] - 1].z == coord[points_with_overlap[i].B].z &&
					coord[nvtr[p].n[l] - 1].y == coord[points_with_overlap[i].B].y
					) flag = false;
			l--;
			// Заменим данную вершину
			coord[points_with_overlap[i].B].z = coord[points_with_overlap[i].point_id].z;
			coord[nvtr[p].n[l] - 1].z = coord[points_with_overlap[i].point_id].z;
			break;
		//case 3:
		//	// Определим вершину, которую нужно изменить
		//	flag = true;
		//	for (l = 0; flag && l < 8; l++)
		//		if (nvtr[p].n[l] != points_with_overlap[i].B + 1 &&
		//			coord[nvtr[p].n[l] - 1].x == coord[points_with_overlap[i].B].x &&
		//			coord[nvtr[p].n[l] - 1].y == coord[points_with_overlap[i].B].y
		//			) flag = false;
		//	l--;
		//	// Заменим данную вершину
		//	coord[points_with_overlap[i].B].y = coord[points_with_overlap[i].point_id].y;
		//	coord[nvtr[p].n[l] - 1].y = coord[points_with_overlap[i].point_id].y;
		//	break;
		default:
			break;
		}

	}

}

// Подпрограмма задания полной информации для каждого терминального узла
void net_generator::setPointsInfo(){
	for (int i = 0; i < t_points.size(); i++){
		// Определим номера вершин, образующих ребро
		switch ( t_points[i].way )
		{
		// движение по Х
		case 1:
		{
			int  n_x_step = trench_size[1] - trench_size[0];
			bool flag = true; 
			int j = 0;
			for ( j = 0; flag && j < n_x_step; j++){
				if (coord[t_points[i].point_id].x >= X0 + (trench_size[0] + j)*Xh &&
					coord[t_points[i].point_id].x <= X0 + (trench_size[0] + j + 1)*Xh )
					flag = false;
			}
			j--;

			Point A(X0 + (trench_size[0] + j + 1)*Xh, coord[t_points[i].point_id].y, coord[t_points[i].point_id].z);
			t_points[i].A = findPoint(A);

			Point B(X0 + (trench_size[0] + j)*Xh, coord[t_points[i].point_id].y, coord[t_points[i].point_id].z);
			t_points[i].B = findPoint(B);

			if (t_points[i].point_id == t_points[i].A || t_points[i].point_id == t_points[i].B)
			{
				t_points.erase(t_points.begin() + i); i--;
			}
			else{

				if (t_points[i].A == -1){
					A.x = X0 + trench_size[1] * Xh;
					t_points[i].A = findPoint(A);
				}
				if (t_points[i].B == -1){
					B.x = X0 + trench_size[0] * Xh;
					t_points[i].B = findPoint(B);
				}

				real den = coord[t_points[i].A].x - coord[t_points[i].B].x;

				t_points[i].A_value = (coord[t_points[i].point_id].x - coord[t_points[i].B].x) / den;

				t_points[i].B_value = (coord[t_points[i].A].x - coord[t_points[i].point_id].x) / den;
			}
			break;
		}
		// движение по Z
		case 2:
		{
			int  n_z_step = trench_size[5] - trench_size[4];
			bool flag = true;
			int j = 0;
			for (j = 0; flag && j < n_z_step; j++){
				if (coord[t_points[i].point_id].z >= Z0 + (trench_size[4] + j)*Zh &&
					coord[t_points[i].point_id].z <= Z0 + (trench_size[4] + j + 1)*Zh )
					flag = false;
			}
			j--;

			Point A(coord[t_points[i].point_id].x, coord[t_points[i].point_id].y, Z0 + (trench_size[4] + j + 1)*Zh);
			t_points[i].A = findPoint(A);

			Point B(coord[t_points[i].point_id].x, coord[t_points[i].point_id].y, Z0 + (trench_size[4] + j)*Zh);
			t_points[i].B = findPoint(B);

			if (t_points[i].point_id == t_points[i].A || t_points[i].point_id == t_points[i].B)
			{
				t_points.erase(t_points.begin() + i); i--;
			}
			else{
				if (t_points[i].A == -1){
					A.z = Z0 + trench_size[5] * Zh;
					t_points[i].A = findPoint(A);
				}
				if (t_points[i].B == -1){
					B.z = Z0 + trench_size[4] * Zh;
					t_points[i].B = findPoint(B);
				}
				real den = coord[t_points[i].A].z - coord[t_points[i].B].z;

				t_points[i].A_value = (coord[t_points[i].point_id].z - coord[t_points[i].B].z) / den;

				t_points[i].B_value = (coord[t_points[i].A].z - coord[t_points[i].point_id].z) / den;
			}
			break;
		}
		// движение по Y
		case 3:
		{
			int  n_y_step = trench_size[3] - trench_size[2];
			bool flag = true;
			int j = 0;
			for (j = 0; flag && j < n_y_step; j++){
				if (coord[t_points[i].point_id].y >= Y0 + (trench_size[2] + j)*Yh &&
					coord[t_points[i].point_id].y <= Y0 + (trench_size[2] + j + 1)*Yh )
					flag = false;
			}
			j--;
			Point A(coord[t_points[i].point_id].x, Y0 + (trench_size[2] + j + 1)*Yh, coord[t_points[i].point_id].z);
			t_points[i].A = findPoint(A);

			Point B(coord[t_points[i].point_id].x, Y0 + (trench_size[2] + j)*Yh, coord[t_points[i].point_id].z);
			t_points[i].B = findPoint(B);
			if (t_points[i].point_id == t_points[i].A || t_points[i].point_id == t_points[i].B)
			{
				t_points.erase(t_points.begin() + i); i--;
			}
			else{
				real den = coord[t_points[i].A].y - coord[t_points[i].B].y;

				t_points[i].A_value = (coord[t_points[i].point_id].y - coord[t_points[i].B].y) / den;

				t_points[i].B_value = (coord[t_points[i].A].y - coord[t_points[i].point_id].y) / den;
			}
			break;
		}
		default:
			break;
		}
	}
	facetPointsBuilder();
}

//рекурсивное заполнение солонки Т-матрицы
void net_generator::recColomBuild(int id, real val, vector <Matrix_Elem> &column,bool horizontal) {
	bool flag = true; int i;
	for (i = 0; flag && i < t_points.size(); i++)
		if (t_points[i].point_id == id) flag = false;
	i--;

	//Если вершина не терминальнная
	if (flag) {
		Matrix_Elem tmp;
		tmp.id = id;
		tmp.val = val;
		column.push_back(tmp);
	}
	else {
		if ((horizontal && t_points[i].A != -1 && t_points[i].B != -1 )|| (!horizontal && (t_points[i].C == -1 || t_points[i].D == -1))) {
		
			//Обслужим А и В
			recColomBuild(t_points[i].A, val*t_points[i].A_value, column, true);
			recColomBuild(t_points[i].B, val*t_points[i].B_value, column, true);
		}
		else {
			//Обслужим C и D
			recColomBuild(t_points[i].C, val*t_points[i].C_value, column, false);
			recColomBuild(t_points[i].D, val*t_points[i].C_value, column, false);
		}

	}

}

// Подпрограмма заполнения матрицы
void net_generator::setTMatrix(){
	t_matrix.ig.push_back(1);
	vector <Matrix_Elem> column;
	for (int i = 0; i < t_points.size(); i++){
		if (t_points[i].A != -1 ) recColomBuild(t_points[i].A, t_points[i].A_value, column, true);
		if (t_points[i].B != -1) recColomBuild(t_points[i].B, t_points[i].B_value, column, true);
		if (t_points[i].C != -1) recColomBuild(t_points[i].C, t_points[i].C_value, column, false);
		if (t_points[i].D != -1) recColomBuild(t_points[i].D, t_points[i].D_value, column, false);
		/*
		// Проверим точку А
		bool flag = true; int j;
		for (j = 0; flag && j < t_points.size(); j++)
			if (t_points[i].A == t_points[j].point_id) flag = false;
		j--;
		if (!flag){
			Matrix_Elem tmp;
			tmp.id = t_points[j].A;
			tmp.val = t_points[j].A_value * t_points[i].A_value;
			column.push_back(tmp);

			tmp.id = t_points[j].B;
			tmp.val = t_points[j].B_value * t_points[i].A_value;
			column.push_back(tmp);
		}
		else{
			Matrix_Elem tmp;
			tmp.id = t_points[i].A;
			tmp.val = t_points[i].A_value;
			column.push_back(tmp);
		}
		// Проверим точку В
		flag = true;
		// Зачем бегать по циклу, когда можно прверить на больше или...
		for (j = 0; flag && j < t_points.size(); j++)
			if (t_points[i].B == t_points[j].point_id) flag = false;
		j--;
		if (!flag){
			Matrix_Elem tmp;
			tmp.id = t_points[j].A;
			tmp.val = t_points[j].A_value * t_points[i].B_value;
			column.push_back(tmp);

			tmp.id = t_points[j].B;
			tmp.val = t_points[j].B_value * t_points[i].B_value;
			column.push_back(tmp);
		}
		else{
			Matrix_Elem tmp;
			tmp.id = t_points[i].B;
			tmp.val = t_points[i].B_value;
			column.push_back(tmp);
		}
		if(t_points[i].C != -1){

			// Проверим точку С
			flag = true;
			// Зачем бегать по циклу, когда можно прверить на больше или...
			for (j = 0; flag && j < t_points.size(); j++)
				if (t_points[i].C == t_points[j].point_id) flag = false;
			j--;
			if (!flag && t_points[j].C != -1 && t_points[j].D != -1) {
				Matrix_Elem tmp;
				tmp.id = t_points[j].C;
				tmp.val = t_points[j].C_value * t_points[i].C_value;
				column.push_back(tmp);

				tmp.id = t_points[j].D;
				tmp.val = t_points[j].D_value * t_points[i].C_value;
				column.push_back(tmp);
			}
			else {
				Matrix_Elem tmp;
				tmp.id = t_points[i].C;
				tmp.val = t_points[i].C_value;
				column.push_back(tmp);
			}
			
		}
		if (t_points[i].D != -1) {

			// Проверим точку С
			flag = true;
			// Зачем бегать по циклу, когда можно прверить на больше или...
			for (j = 0; flag && j < t_points.size(); j++)
				if (t_points[i].C == t_points[j].point_id) flag = false;
			j--;
			if (!flag && t_points[j].C != -1 && t_points[j].D != -1) {
				Matrix_Elem tmp;
				tmp.id = t_points[j].C;
				tmp.val = t_points[j].C_value * t_points[i].D_value;
				column.push_back(tmp);

				tmp.id = t_points[j].D;
				tmp.val = t_points[j].D_value * t_points[i].D_value;
				column.push_back(tmp);
			}
			else {
				Matrix_Elem tmp;
				tmp.id = t_points[i].D;
				tmp.val = t_points[i].D_value;
				column.push_back(tmp);
			}

		}

		*/
		//Сортируем точки
		std::sort(column.begin(), column.end(), sort_coord_col);
		//Заносим элементы в матрицу
		for (int k = 0; k < column.size(); k++){
			t_matrix.jg.push_back(column[k].id);
			t_matrix.gg.push_back(column[k].val);
		}
		t_matrix.ig.push_back(t_matrix.ig[t_matrix.ig.size() - 1] + column.size());
		column.clear();
	}
}

void net_generator::facetPointsBuilder()
{
	for (int i = 0; i < face_points.size(); i++) {
		T_Point A;
		A.point_id = face_points[i].point_id;
		A.A = face_points[i].left;
		A.B = face_points[i].right;
		A.C = face_points[i].up;
		A.D = face_points[i].down;
		A.way = 0;
		calculateTVertexValue(A);
		t_points.push_back(A);
	}
	std::sort(t_points.begin(), t_points.end(), sort_t_points);
}
void net_generator::calculateTVertexValue(T_Point &P){
	
	double lenX = sqrt((coord[P.A].x - coord[P.B].x)*(coord[P.A].x - coord[P.B].x) +
		(coord[P.A].z - coord[P.B].z)*(coord[P.A].z - coord[P.B].z));
	double lenZ = sqrt((coord[P.C].x - coord[P.D].x)*(coord[P.C].x - coord[P.D].x) +
		(coord[P.C].z - coord[P.D].z)*(coord[P.C].z - coord[P.D].z));
	//Получим значение
	if (P.A_value == -1) {
		double lenXA = sqrt((coord[P.point_id].x - coord[P.B].x)*(coord[P.point_id].x - coord[P.B].x) +
			(coord[P.point_id].z - coord[P.B].z)*(coord[P.point_id].z - coord[P.B].z));
		P.A_value = lenXA / lenX;
		/*
		double r;
		bool flag = false;
		for (int i = 0; !flag && i < face_points.size(); i++)
			if (face_points[i].point_id == P.A)
			{
				flag = true;
				if (t_points[i].A_value == -1)
					calculateTVertexValue(t_points[i]);
				r = t_points[i].A_value;
			}
		if (flag)
			P.A_value = r * T;
		else
			P.A_value = T;*/
	}
	if (P.B_value == -1) {
		double lenXB = sqrt((coord[P.point_id].x - coord[P.A].x)*(coord[P.point_id].x - coord[P.A].x) +
			(coord[P.point_id].z - coord[P.A].z)*(coord[P.point_id].z - coord[P.A].z));
		P.B_value = lenXB / lenX;
	/*
		double r;
		bool flag = false;
		for (int i = 0; !flag && i < t_points.size(); i++)
			if (t_points[i].point_id == P.B)
			{
				flag = true;
				if (t_points[i].B_value == -1)
					calculateTVertexValue(t_points[i]);
				r = t_points[i].B_value;
			}
		if (flag)
			P.B_value = r * T;
		else
			P.B_value = T;*/
	}
	if (P.C_value == -1) {
		double lenXC = sqrt((coord[P.point_id].x - coord[P.D].x)*(coord[P.point_id].x - coord[P.D].x) +
			(coord[P.point_id].z - coord[P.D].z)*(coord[P.point_id].z - coord[P.D].z));
		P.C_value = lenXC / lenZ;
		/*
		double r;
		bool flag = false;
		for (int i = 0; !flag && i < t_points.size(); i++)
			if (t_points[i].point_id == P.C)
			{
				flag = true;
				if (t_points[i].C_value == -1)
					calculateTVertexValue(t_points[i]);
				r = t_points[i].C_value;
			}
		if (flag)
			P.C_value = r * T;
		else
			P.C_value = T;*/
	}
	if (P.D_value == -1) {
		double lenXD = sqrt((coord[P.point_id].x - coord[P.C].x)*(coord[P.point_id].x - coord[P.C].x) +
			(coord[P.point_id].z - coord[P.C].z)*(coord[P.point_id].z - coord[P.C].z));
		P.D_value = lenXD / lenZ;
	}
		/*
		double r;
		bool flag = false;
		for (int i = 0; !flag && i < t_points.size(); i++)
			if (t_points[i].point_id == P.D)
			{
				flag = true;
				if (t_points[i].D_value == -1)
					calculateTVertexValue(t_points[i]);
				r = t_points[i].D_value;
			}
		if (flag)
			P.D_value = r * T;
		else
			P.D_value = T;*/
	
	/*
	for (int i = 0; i < t_points.size(); i++) {
		if (t_points[i].point_id == P.point_id) 
		{
			P.A_value = t_points[i].A_value;
			P.B_value = t_points[i].B_value;
			return true;
		}
	}*/
	/*
	for (int i = 0; i < face_points.size(); i++) {
		if()
	}*/


}


//---------------------------------------------------------------------------------------

int net_generator::findPoint(Point A){
	int i;
	for (i = 0; i < coord.size(); i++)
		if (A.x == coord[i].x && A.y == coord[i].y && A.z == coord[i].z)
			return coord[i].id;
	if (i == coord.size()) return -1;
}