#include "links.h"

class vect{
public:
	real x, y, z;
	//Происходит поворот по, либо против часовой стрелки
	bool clockwiseZ, clockwiseX, clockwiseY;
	// поворот точки на вектор
	Point rotatePoint(Point A, Point v){
		Point res;
		real alfaZ, alfaX, alfaY(1);
		// косинус угола в плоскости XY
		alfaZ = (x*v.x + y*v.y) / (sqrt(v.x*v.x + v.y * v.y) * sqrt(x*x + y*y));
		//косинус угола в плоскости YZ
		alfaX = (z*v.z + y*v.y) / (sqrt(v.z*v.z + v.y * v.y) * sqrt(y*y + z*z));
		res = A;
		/*
		

		if (x!= 0 && y == 0 && z!= 0){
			//косинус угола в плоскости XZ
			
			// поворот исходной точки вокруг оси X
			if (clockwiseZ)
				res = rotateOXB(res, z);
			else
				res = rotateOXF(res, z);

			// поворот исходной точки вокруг оси Z
			if (clockwiseX)
				res = rotateOyF(res, 0.0);
			else
				res = rotateOyB(res, 0.0);
			
		
			res.id = A.id;
			return res;
		}
		*/

		// поворот исходной точки вокруг оси Х
		if (!isnan(alfaX) && alfaX!=-1)
		if (clockwiseZ)
			res = rotateOXB(A, alfaX);
		else 
			res = rotateOXF(A, alfaX);
			
		// поворот полученной точки вокруг оси Z
		if (!isnan(alfaZ) && alfaZ != -1)
		if (clockwiseX)
			res = rotateOzF(res, alfaZ);
		else 
			res = rotateOzB(res, alfaZ);

		res.id = A.id;
		return res;

	}
	
	real length(){
		real l = sqrt(x*x + y*y + z*z);
		return l;
	}
	vect(){};
	// x,y,z, вокруг Z,вокруг X 
	vect(real x1, real x2, real x3,bool A, bool B, bool C){
		x = x1;
		y = x2;
		z = x3;
		clockwiseZ = A;
		clockwiseX = B;
		clockwiseY = C;
	}
private:
	
	// поворото вокруг оси OX против часовой стрелке
	Point rotateOXB(Point A, real cos){
		real s = sqrt(1-cos*cos);
		real c = cos;
		Point R;
		R.x = A.x;
		R.y = A.y*c - A.z*s;
		R.z = A.y*s + A.z*c;
		return R;
	}
	// поворото вокруг оси OY против часовой стрелке
	Point rotateOyB(Point A, real cos){
		real s = sqrt(1 - cos*cos);
		real c = cos;
		Point R;
		R.y = A.y;
		R.x = A.x*c + A.z*s;
		R.z = -A.x*s + A.z*c;
		return R;
	}
	// поворот вокруг оси OZ против часовой стрелке
	Point rotateOzB(Point A, real cos){
		real s = sqrt(1 - cos*cos);
		real c = cos;
		Point R;
		R.x = A.x * c - s * A.y;
		R.y = A.x * s + A.y*c;
		R.z = A.z;
		return R;
	}

	// поворото вокруг оси OX по часовой стрелке
	Point rotateOXF(Point A, real cos){
		real s = sqrt(1 - cos*cos);
		real c = cos;
		Point R;
		R.x = A.x;
		R.y = A.y*c + A.z*s;
		R.z = -A.y*s + A.z*c;
		return R;
	}
	// поворото вокруг оси OY по часовой стрелке
	Point rotateOyF(Point A, real cos){
		real s = sqrt(1 - cos*cos);
		real c = cos;
		Point R;
		R.y = A.y;
		R.x = A.x*c - A.z*s;
		R.z = A.x*s + A.z*c;
		return R;
	}
	// поворот вокруг оси OZ по часовой стрелке
	Point rotateOzF(Point A, real cos){
		real s = sqrt(1 - cos*cos);
		real c = cos;
		Point R;
		R.x = A.x * c + s * A.y;
		R.y =  - A.x * s + A.y*c;
		R.z = A.z;
		return R;
	}
};