//Name: Alex Lin
//Date: 9/19/2019
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <math.h>

#define IMAGE_SIZE 800

using namespace std;

struct Point {
	double x;
	double y;
	Point() {
		this->x = (double) rand() / RAND_MAX;
		this->y = (double) rand() / RAND_MAX;
	}
        Point(double x, double y) {
		this->x = x;
		this->y = y;
	}
};

double distance(Point p1, Point p2) {
	return sqrt(pow(p2.x - p1.x, 2.0) + pow(p2.y - p1.y, 2.0));
}

Point calculateIncirclePoint(Point A, Point B, Point C, double a, double b, double c) {
	Point incircle ((a * A.x + b * B.x + c * C.x) / (a + b + c),
				(a * A.y + b * B.y + c * C.y) / (a + b + c));
	return incircle;
}

Point calculateCircumcirclePoint(Point p1, Point p2, Point p3) {
	Point   midpoint1 ((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0),
		midpoint2 ((p1.x + p3.x) / 2.0, (p1.y + p3.y) / 2.0);
	//Finding the Perpendicular Angle Bisector, -bx + ay = c
	double  a1 = -(p1.x - p2.x),
		b1 = p2.y - p1.y,
		c1 = a1 * midpoint1.x + b1 * midpoint1.y,

		a2 = -(p1.x - p3.x),
		b2 = p3.y - p1.y,
		c2 = a2 * midpoint2.x + b2 * midpoint2.y,
		determinant = a1 * b2 - a2 * b1;

	Point circumcircle_center ((b2 * c1 - b1 * c2) / determinant, (a1 * c2 - a2 * c1) / determinant);
	return circumcircle_center;
}

double calculateIncircleRadius(double a, double b, double c) {
	double s = (a + b + c) / 2.0;
	return sqrt(s * (s - a) * (s - b) * (s - c)) / s;
}

double calculateCircumcircleRadius(double a, double b, double c) {
	double s = (a + b + c) / 2.0;
	return (a * b * c) / (4.0 * sqrt(s * (s - a) * (s - b) * (s - c)));
}

//datastructure for storing a pixel's rgb values
struct Color {
	bool r, g, b;
	Color() {
		this->r = false;
		this->g = false;
		this->b = false;
	}
	Color(bool r, bool g, bool b) {
		this->r = r;
		this->g = g;
		this->b = b;
	}
};

void drawLine(Color (&pixMap)[IMAGE_SIZE][IMAGE_SIZE], double x1, double y1, double x2, double y2, Color color) {
	//Bresemham's integer algorithm doesn't draw multiple pixels per column so if its steep it should swap the x and y values to fix this issue
    bool steep = abs(y2 - y1) > abs(x2 - x1);
	
	if (steep) {
		swap(x1, y1);
		swap(x2, y2);
	}
	
	//makes sure that x1 is always to the left of x2
	if (x1 > x2) {
		swap(x1, x2);
		swap(y1, y2);
	}

	//changes the x and y values from 0-1 to 0-IMAGE_SIZE
	x1 *= IMAGE_SIZE;
	y1 *= IMAGE_SIZE;
	x2 *= IMAGE_SIZE;
	y2 *= IMAGE_SIZE;

	int	dx = (int) abs(x2 - x1),
		dy = (int) abs(y2 - y1),
		e = (int) (dy - dx),
		ystep = (y1 < y2) ? 1 : -1,
		y = (int) y1,
		maxX = (int) x2;

	for (int x = (int) x1; x < maxX; x++) {
		if (x >= 0 && x < IMAGE_SIZE && y >= 0 && y < 800)
            //x and y values were swapped if steep
			if (steep)
				pixMap[x][y] = color;
			else
				pixMap[y][x] = color;
		if (e >= 0) {
			y += ystep;
			e -= dx;
		}
		e += dy;
	}		
}

void drawCircle(Color (&pixMap)[IMAGE_SIZE][IMAGE_SIZE], Point center, double radius, Color color) {
	int	x,
		y = (int) (radius * IMAGE_SIZE),
		maxX = (int) ((radius * IMAGE_SIZE) / sqrt(2)),
		y2 = y * y,
		ty = (2 * y) - 1,
		y2_new = y2,
	
		center_x = (int) (center.x * IMAGE_SIZE),
		center_y = (int) (center.y * IMAGE_SIZE);

	for (x = 0; x <= maxX; x++) {
		if ((y2 - y2_new) >= ty) {
			y2 -= ty;
			y -= 1;
			ty -= 2;
		}
	
        //only draws to array if in bounds
		if (x + center_x >= 0 && x + center_x < IMAGE_SIZE && y + center_y >= 0 && y + center_y < IMAGE_SIZE)
			pixMap[y + center_y][x + center_x] = color;
		if (x + center_x >= 0 && x + center_x < IMAGE_SIZE && -y + center_y >= 0 && -y + center_y < IMAGE_SIZE)
			pixMap[-y + center_y][x + center_x] = color;
		if (-x + center_x >= 0 && -x + center_x < IMAGE_SIZE && y + center_y >= 0 && y + center_y < IMAGE_SIZE)
			pixMap[y + center_y][-x + center_x] = color;
		if (-x + center_x >= 0 && -x + center_x < IMAGE_SIZE && -y + center_y >= 0 && -y + center_y < IMAGE_SIZE)
			pixMap[-y + center_y][-x + center_x] = color;


		if (y + center_x >= 0 && y + center_x < IMAGE_SIZE && x + center_y >= 0 && x + center_y < IMAGE_SIZE)
			pixMap[x + center_y][y + center_x] = color;
		if (y + center_x >= 0 && y + center_x < IMAGE_SIZE && -x + center_y >= 0 && -x + center_y < IMAGE_SIZE)
	        	pixMap[-x + center_y][y + center_x] = color;
		if (-y + center_x >= 0 && -y + center_x < IMAGE_SIZE && x + center_y >= 0 && x + center_y < IMAGE_SIZE)
	        	pixMap[x + center_y][-y + center_x] = color;
		if (-y + center_x >= 0 && -y + center_x < IMAGE_SIZE && -x + center_y >= 0 && -x + center_y < IMAGE_SIZE)
			pixMap[-x + center_y][-y + center_x] = color;

		y2_new -= (2 * x) - 3;
	}
}

void generateImage(Color (&pixMap)[IMAGE_SIZE][IMAGE_SIZE]) {
	ofstream imageFile;
	imageFile.open("image.ppm");

	imageFile << "P3 " << IMAGE_SIZE << " " << IMAGE_SIZE << " 1" << endl;

	for (int y = IMAGE_SIZE - 1; y >= 0; y--) {
		for (int x = 0; x < IMAGE_SIZE; x++)
			imageFile << (int) pixMap[y][x].r << " " << (int) pixMap[y][x].g << " " << (int) pixMap[y][x].b << " ";
		imageFile << endl;
	}
	imageFile.close();
}

int main() {
    //random seed for rand
	srand((unsigned)time(NULL));

	Point point_a, point_b, point_c;
	
	Point   midpoint_ab ((point_a.x + point_b.x) / 2, (point_a.y + point_b.y) / 2),
		midpoint_bc ((point_b.x + point_c.x) / 2, (point_b.y + point_c.y) / 2),
		midpoint_ac ((point_a.x + point_c.x) / 2, (point_a.y + point_c.y) / 2);
	
    //sides for triangle and inner-triangle
	double 	side_a = distance(point_b, point_c),
        	side_b = distance(point_a, point_c),
		side_c = distance(point_a, point_b),
		side_ab = distance(midpoint_bc, midpoint_ac),
		side_bc = distance(midpoint_ab, midpoint_ac),
		side_ac = distance(midpoint_ab, midpoint_bc);

	Point	incircle_center = calculateIncirclePoint(point_a, point_b, point_c, side_a, side_b, side_c),
		circumcircle_center = calculateCircumcirclePoint(point_a, point_b, point_c),
		ninepointcircle_center = calculateCircumcirclePoint(midpoint_ab, midpoint_bc, midpoint_ac);

	double	incircle_radius = calculateIncircleRadius(side_a, side_b, side_c),
		circumcircle_radius = calculateCircumcircleRadius(side_a, side_b, side_c),
		ninepointcircle_radius = calculateCircumcircleRadius(side_ab, side_bc, side_ac),

		eulerline_slope = (incircle_center.y - circumcircle_center.y) / (incircle_center.x - circumcircle_center.x),
		eulerline_yint = -eulerline_slope * incircle_center.x + incircle_center.y,
		eulerline_xint = -eulerline_yint / eulerline_slope,
		eulerline_maxVal = (IMAGE_SIZE - 1) / (double) IMAGE_SIZE;

	Color  base_color(true, true, true),	
        triangle_color(false, false, false),
		eulerline_color(true, true, false),
		incircle_color(true, false, false),
		circumcircle_color(false, true, false),
		ninepointcircle_color(false, false, true);
        
	Color pixMap[IMAGE_SIZE][IMAGE_SIZE];	

    for (int i = 0; i < IMAGE_SIZE; i++)
        for (int j = 0; j < IMAGE_SIZE; j++)
            pixMap[i][j] = base_color;
    
    //triangle
	drawLine(pixMap, point_a.x, point_a.y, point_b.x, point_b.y, triangle_color);
	drawLine(pixMap, point_b.x, point_b.y, point_c.x, point_c.y, triangle_color);
	drawLine(pixMap, point_a.x, point_a.y, point_c.x, point_c.y, triangle_color);
	
    //eulerline
	drawLine(pixMap, 0.0, eulerline_yint, eulerline_maxVal, eulerline_slope * eulerline_maxVal + eulerline_yint, eulerline_color);

    //circles
	drawCircle(pixMap, incircle_center, incircle_radius, incircle_color);
	drawCircle(pixMap, circumcircle_center, circumcircle_radius, circumcircle_color);
	drawCircle(pixMap, ninepointcircle_center, ninepointcircle_radius, ninepointcircle_color);

	generateImage(pixMap);

	cout	<< "Points" 
		<< "\n\tA: " << point_a.x << ", " << point_a.y
		<< "\n\tB: " << point_b.x << ", " << point_b.y
		<< "\n\tC: " << point_c.x << ", " << point_c.y
		<< "\n\tMidpoint AB: " << midpoint_ab.x << ", " << midpoint_ab.y
		<< "\n\tMidpoint BC: " << midpoint_bc.x << ", " << midpoint_bc.y
		<< "\n\tMidpoint AC: " << midpoint_ac.x << ", " << midpoint_ac.y
		<< "\nSides"
		<< "\n\ta: " << side_a
		<< "\n\tb: " << side_b
		<< "\n\tc: " << side_c
		<< "\n\tab: " << side_ab
		<< "\n\tbc: " << side_bc
		<< "\n\tac: " << side_ac
		<< "\nIncircle" 
			<< "\n\tPoint: " << incircle_center.x << ", " << incircle_center.y
			<< "\n\tRadius: " << incircle_radius	
		<< "\nCircumcircle" 
			<< "\n\tPoint: " << circumcircle_center.x << ", " << circumcircle_center.y
			<< "\n\tRadius: " << circumcircle_radius	
		<< "\nNine Point Circle" 
			<< "\n\tPoint: " << ninepointcircle_center.x << ", " << ninepointcircle_center.y
			<< "\n\tRadius: " << ninepointcircle_radius
		<<"\nEuler Line"
			<< "\n\tSlope: " << eulerline_slope
			<< "\n\tX-Intercept: " << eulerline_xint
			<< "\n\tY-Intercept: " << eulerline_yint
		<< endl;
	return 0;
}
