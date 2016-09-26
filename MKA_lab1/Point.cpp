class Point {
public:
	double x, y;
	Point() {
		x = 0; y = 0;
	}
	Point(double new_x, double new_y) {
		x = new_x; y = new_y;
	}
	Point& operator &=(Point& point) {
		x = point.x; y = point.y;
		return *this;
	}
};