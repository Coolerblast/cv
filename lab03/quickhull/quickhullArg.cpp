#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <chrono>

#define IMAGE_SIZE 400
int NUM_POINTS;

using namespace std;

class Point {
   public:
    double x;
    double y;
    Point() {
        this->x = rand() / (RAND_MAX + 1.0);
        this->y = rand() / (RAND_MAX + 1.0);
    }
    Point(const double& x, const double& y) {
        this->x = x;
        this->y = y;
    }
    friend ostream& operator<<(ostream& strm, const Point& a) {
        return strm << "X: " << a.x << " Y: " << a.y;
    }
};

struct Color {
    bool r, b, g;
    Color() {
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

void drawLine(vector<Color>& pixMap, const Point& a, const Point& b, const Color& color) {
    double x1 = a.x * IMAGE_SIZE,
           y1 = a.y * IMAGE_SIZE,
           x2 = b.x * IMAGE_SIZE,
           y2 = b.y * IMAGE_SIZE;
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
    int dx = (int)abs(x2 - x1),
        dy = (int)abs(y2 - y1),
        e = (int)(dy - dx),
        ystep = (y1 < y2) ? 1 : -1,
        y = (int)y1,
        maxX = (int)x2;

    for (int x = (int)x1; x < maxX; x++) {
        if (x >= 0 && x < IMAGE_SIZE && y >= 0 && y < 800)
            //x and y values were swapped if steep
            if (steep)
                pixMap[x * IMAGE_SIZE + y] = color;
            else
                pixMap[y * IMAGE_SIZE + x] = color;
        if (e >= 0) {
            y += ystep;
            e -= dx;
        }
        e += dy;
    }
}

void drawPoints(vector<Color>& pixMap, const vector<Point>& points, const vector<const Point*>& hull) {
    for (int i = 0; i < points.size(); i++) {
        int x = (int)(points[i].x * IMAGE_SIZE);
        int y = (int)(points[i].y * IMAGE_SIZE);
        pixMap[y * IMAGE_SIZE + x] = Color(false, false, false);
    }
    for (int i = 0; i < hull.size(); i++) {
        int x = (int)((*hull[i]).x * IMAGE_SIZE);
        int y = (int)((*hull[i]).y * IMAGE_SIZE);
        pixMap[y * IMAGE_SIZE + x] = Color(true, false, false);
    }
}

void generateImage(vector<Color>& pixMap) {
    ofstream imageFile;
    imageFile.open("image.ppm");

    imageFile << "P3 " << IMAGE_SIZE << " " << IMAGE_SIZE << " 1" << endl;

    for (int y = IMAGE_SIZE - 1; y >= 0; y--) {
        for (int x = 0; x < IMAGE_SIZE; x++)
            imageFile << (int)pixMap[y * IMAGE_SIZE + x].r << " " << (int)pixMap[y * IMAGE_SIZE + x].g << " " << (int)pixMap[y * IMAGE_SIZE + x].b << " ";
        imageFile << endl;
    }
}

inline double lineDistance(const Point& a, const Point& b, const Point& c) {
    return (c.y - a.y) * (b.x - a.x) - (b.y - a.y) * (c.x - a.x);
}

inline int findSide(const double& dist) {
    return dist < 0 ? -1 : dist > 0 ? 1 : 0;
}

vector<const Point*> hull;

void findHull(const vector<Point>& points, const Point& left, const Point& right, int side) {
    int max = -1;
    double tempDist = 0;
    double maxDist = 0;
    for (int i = 0; i < points.size(); i++) {
        tempDist = lineDistance(left, right, points[i]);
        if (findSide(tempDist) == side && maxDist < abs(tempDist)) {
            max = i;
            maxDist = abs(tempDist);
        }
    }
    if (max == -1) {
        hull.emplace_back(&left);
        hull.emplace_back(&right);
        return;
    }
    findHull(points, points[max], left, -findSide(lineDistance(points[max], left, right)));
    findHull(points, points[max], right, -findSide(lineDistance(points[max], right, left)));
}

vector<Point*>& quickHull(vector<Point>& points) {
    Point* left = &points[0];
    Point* right = &points[0];
    for (auto il = points.begin() + 1; il != points.end(); il++) {
        if (il->x < left->x)
            left = &*il;
        else if (il->x > right->x)
            right = &*il;
    }
    findHull(points, *left, *right, 1);
    findHull(points, *left, *right, -1);
}

int main(int argc, char** argv) {
    NUM_POINTS = atoi(argv[1]);
    srand((unsigned)time(0));
    vector<Point>* points = new vector<Point>(NUM_POINTS);
    auto begin = chrono::steady_clock::now();
    quickHull(*points);
    auto end = chrono::steady_clock::now();
    cout << "Number of Points: " << NUM_POINTS << "\n"
         << "Elapsed Time: " << chrono::duration_cast<chrono::milliseconds>(end - begin).count()
         << "\n";

    /*
    vector<Color>* pixMap = new vector<Color>(IMAGE_SIZE * IMAGE_SIZE);
    for (auto il = hull.begin(); il != hull.end(); il += 2) {
        drawLine(*pixMap, **il, **(il + 1), Color(true, false, true));
    }
    drawPoints(*pixMap, *points, hull);
    generateImage(*pixMap);
    */
    return 0;
}