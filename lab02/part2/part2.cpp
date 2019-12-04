#include <iostream>
#include <chrono>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <vector>
#include <algorithm>

#define IMAGE_SIZE 200
#define NUM_POINTS 200

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

void drawPoint(Color (&pixMap)[IMAGE_SIZE][IMAGE_SIZE], vector<Point> points, Point closestPair[2]) {
    for (int i = 0; i < NUM_POINTS; i++) {
        int x = (int) (points[i].getX() * IMAGE_SIZE);
        int y = (int) (points[i].getY() * IMAGE_SIZE);
        pixMap[y][x] = Color (false, false, false);
    }
    for (int i = 0; i < 2; i++) {
        int x = (int) (closestPair[i].getX() * IMAGE_SIZE);
        int y = (int) (closestPair[i].getY() * IMAGE_SIZE);
        pixMap[y][x] = Color (true, false, false);
    }
}

void generateImage(Color (&pixMap)[IMAGE_SIZE][IMAGE_SIZE]) {
    ofstream imageFile;
    imageFile.open("image.ppm");

    imageFile << "P3 " << IMAGE_SIZE << " " << IMAGE_SIZE << " 1" << endl;

    for (int y = IMAGE_SIZE - 1; y >= 0; y--) {
        for (int x = IMAGE_SIZE - 1; x >= 0; x--) 
            imageFile << (int) pixMap[y][x].r << " " << (int) pixMap[y][x].g << " " << (int) pixMap[y][x].b << " ";
        imageFile << endl;
   }
}

double distance(Point p1, Point p2) {
    return sqrt((p2.getX() - p1.getX()) * (p2.getX() - p1.getX()) + (p2.getY() - p1.getY()) * (p2.getY() - p1.getY()));
}

double min(double x, double y) {
    return x < y ? x : y;
}

Point closestPairGBL[2];
double minD = DBL_MAX;

double bruteForce(vector<Point> points, int left, int right) {
    Point closestPair[2];
    double min_distance = DBL_MAX;
    double temp;
    for (int i = left; i <= right; i++)
        for (int j = i + 1; j <= right; j++) {
            temp = distance(points[i], points[j]);
            if (temp < min_distance) {
                min_distance = temp;
                closestPair[0] = points[i];
                closestPair[1] = points[j];
            }
        }
    if (min_distance < minD) {
        closestPairGBL[0] = closestPair[0];
        closestPairGBL[1] = closestPair[1];
        minD = min_distance;
    }
    return min_distance;
}

double stripClosest(vector<Point> strip, double d) {
    Point closestPair[2];
    double min_distance = d;
    double temp;
    for (int i = 0; i < strip.size(); i++)
        for (int j = i + 1; j < strip.size(); j++) {
            temp = distance(strip[i], strip[j]);
            if (temp < min_distance) {
                min_distance = temp;
                closestPair[0] = strip[i];
                closestPair[1] = strip[j];
            }
        }
    if (min_distance < minD) {
        closestPairGBL[0] = closestPair[0];
        closestPairGBL[1] = closestPair[1];
        minD = min_distance;
    }
    return min_distance;
}

double closestPair(vector<Point> points, int left, int right) {
    if (right - left <= 3)
        return bruteForce(points, left, right);

    int mid = left + (right - left) / 2;
    Point mid_point = points[mid];

    double min_distance = min(closestPair(points, left, mid), closestPair(points, mid + 1, right));

    vector<Point> strip;
    int count = 0;
    for (int i = left; i < right; i++)
        if (abs(points[i].getX() - mid_point.getX() < min_distance)) {
            strip.push_back(points[i]);
            count++;
        }
    stripClosest(strip, count);
    
    // if (min_distance < minD) {
    //     closestPairGBL = closestPair;
    //     minD = min_distance;
    // }
    return min_distance;
}

void closest(vector<Point> points) {
    sort(points.begin(), points.end(), compareX);
    closestPair(points, 0, NUM_POINTS);
}

int main() {
    srand((unsigned)time(NULL));

    vector<Point> points;
    for (int i = 0; i < NUM_POINTS; i++)
        points.push_back(Point());

    Color pixMap[IMAGE_SIZE][IMAGE_SIZE];

    std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();

    closest(points);

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();

    cout << closestPairGBL[0] << " " << closestPairGBL[1] << endl;
	cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << std::endl;
    
    drawPoint(pixMap, points, closestPairGBL);

    generateImage(pixMap);

    // sort(points, points + NUM_POINTS, compareX);

    // for (int i = 0; i < NUM_POINTS; i++)
    //     cout << points[i] << endl;
    // cout << "Median X at Index " << NUM_POINTS / 2 << ": " << points[NUM_POINTS / 2];

    return 0;
}