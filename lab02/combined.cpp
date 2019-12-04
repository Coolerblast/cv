//Name: Alex Lin
//Date: 10/23/2019

#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>
#include <algorithm>
#include <chrono>

#define IMAGE_SIZE 800
#define NUM_POINTS 1000

using namespace std;

class Point {
    private:
        double x;
        double y;
    public:
        Point() {
            this->x = (double) rand() / RAND_MAX;
            this->y = (double) rand() / RAND_MAX;
        }
        Point(double x, double y) {
            this->x = x;
            this->y = y;
        }
        double getX() {
            return x;
        }
        double getY() {
            return y;
        }
        friend ostream& operator<<(ostream &strm, const Point &a) {
            return strm << "X: " << a.x << " Y: " << a.y;
        }

};

bool compareX(Point a, Point b) {
    return a.getX() < b.getX();
}

bool compareY(Point a, Point b) {
    return a.getY() < b.getY();
}

struct Color {
    bool r, b, g;
    Color () {
        this->r = true;
        this->g = true;
        this->b = true;
    }
    Color(bool r, bool g, bool b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
};

void drawPoint(Color pixMap[], Point points[], Point closestPair[]) {
    for (int i = 0; i < NUM_POINTS; i++) {
        int x = (int) (points[i].getX() * IMAGE_SIZE);
        int y = (int) (points[i].getY() * IMAGE_SIZE);
        pixMap[y * IMAGE_SIZE + x] = Color (false, false, false);
    }
    for (int i = 0; i < 2; i++) {
        int x = (int) (closestPair[i].getX() * IMAGE_SIZE);
        int y = (int) (closestPair[i].getY() * IMAGE_SIZE);
        pixMap[y * IMAGE_SIZE + x] = Color (true, false, false);
    }
}

void generateImage(Color pixMap[]) {
    ofstream imageFile;
    imageFile.open("image.ppm");

    imageFile << "P3 " << IMAGE_SIZE << " " << IMAGE_SIZE << " 1" << endl;

    for (int y = IMAGE_SIZE - 1; y >= 0; y--) {
        for (int x = 0; x < IMAGE_SIZE; x++) 
            imageFile << (int) pixMap[y * IMAGE_SIZE + x].r << " " << (int) pixMap[y * IMAGE_SIZE + x].g << " " << (int) pixMap[y * IMAGE_SIZE + x].b << " ";
        imageFile << endl;
   }
}

double distance(Point p1, Point p2) {
    return sqrt((p2.getX() - p1.getX()) * (p2.getX() - p1.getX()) + (p2.getY() - p1.getY()) * (p2.getY() - p1.getY()));
}

double min(double x, double y) {
    return x < y ? x : y;
}

Point* closestPairGBL;
double* minD;

double bruteForce(Point points[], int n) {
    Point* closestPair = new Point[2];
    double min_distance = 1;
    double temp;
    for (int i = 0; i <= n; i++)
        for (int j = i + 1; j <= n; j++) {
            temp = distance(points[i], points[j]);
            if (temp < min_distance) {
                min_distance = temp;
                closestPair[0] = points[i];
                closestPair[1] = points[j];
            }
        }
    if (min_distance < *minD) {
        closestPairGBL[0] = closestPair[0];
        closestPairGBL[1] = closestPair[1];
        *minD = min_distance;
    }
    delete[] closestPair;
    return min_distance;
}

double stripClosest(Point strip[], int n, double d) {
    Point* closestPair = new Point[2];
    double min_distance = d;
    double temp;

    sort(strip, strip + n, compareY);

    int count;
    for (int i = 0; i < n; i++) {
        count = 0;
        for (int j = i + 1; j < n && (strip[j].getY() - strip[i].getY()) < d; j++) {
            temp = distance(strip[i], strip[j]);
            if (temp < min_distance) {
                min_distance = temp;
                closestPair[0] = strip[i];
                closestPair[1] = strip[j];
            }
            if (++count > 7)
                break;
        }
    }
    if (min_distance < *minD) {
        closestPairGBL[0] = closestPair[0];
        closestPairGBL[1] = closestPair[1];
        *minD = min_distance;
    }
    delete[] closestPair;
    return min_distance;
}

double closestPairStrip(Point points[], int n) {
    if (n <= 3)
        return bruteForce(points, n);
    
    int mid = n / 2;
    Point mid_point = points[mid];

    double min_distance = min(closestPairStrip(points, mid), closestPairStrip(points + mid, n - mid));

    Point* strip = new Point[n];
    int count = 0;
    for (int i = 0; i < n; i++)
        if (abs(points[i].getX() - mid_point.getX()) < min_distance) {
            strip[count] = points[i];
            count++;
        }
    min_distance = min(stripClosest(strip, count, min_distance), min_distance);
    delete[] strip;
    return min_distance;
}

double closestPairBrute(Point points[], int n) {
    if (n <= 3)
        return bruteForce(points, n);

    int mid = n / 2;
    Point mid_point = points[mid];

    double min_distance = min(closestPairBrute(points, mid), closestPairBrute(points + mid, n - mid));

    Point* strip = new Point[n];
    int count = 0;
    for (int i = 0; i < n; i++)
        if (abs(points[i].getX() - mid_point.getX()) < min_distance) {
            strip[count] = points[i];
            count++;
        }
    min_distance = min(bruteForce(strip, count), min_distance);
    delete[] strip;
    return min_distance;
}

void closestStrip(Point points[]) {
    minD = new double(1);
    closestPairGBL = new Point[2];
    sort(points, points + NUM_POINTS, compareX);
    closestPairStrip(points, NUM_POINTS);
}

void closestBrute(Point points[]) {
    minD = new double(1);
    closestPairGBL = new Point[2];
    sort(points, points + NUM_POINTS, compareX);
    closestPairBrute(points, NUM_POINTS);
}

int main() {
    srand((unsigned)time(NULL));

    Point* points = new Point[NUM_POINTS];

    Color* pixMap = new Color[IMAGE_SIZE * IMAGE_SIZE];

    cout << "Number of Points:\t" << NUM_POINTS << endl;

    std::chrono::high_resolution_clock::time_point begin, end;

    //Min Distance is started as 1 because the unit square is from 0 to a max of 1
    minD = new double(1);
    //New array created on the heap
    closestPairGBL = new Point[2];
    begin = std::chrono::high_resolution_clock::now();
    bruteForce(points, NUM_POINTS);
    end = std::chrono::high_resolution_clock::now();
	cout << "Brute Force:\t\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << endl;
    //free memory allocated on the heap
    delete minD;
    delete[] closestPairGBL;

    begin = std::chrono::high_resolution_clock::now();
    closestBrute(points);
    end = std::chrono::high_resolution_clock::now();
	cout << "Recursive (BF):\t\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << endl;
    delete minD;
    delete[] closestPairGBL;

    begin = std::chrono::high_resolution_clock::now();
    closestStrip(points);
    end = std::chrono::high_resolution_clock::now();
	cout << "Recursive (nlogn):\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << endl;

    drawPoint(pixMap, points, closestPairGBL);

    delete minD;
    delete[] closestPairGBL;

    generateImage(pixMap);

    return 0;
}