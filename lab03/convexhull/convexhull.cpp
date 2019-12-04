#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#define IMAGE_SIZE 400
#define NUM_POINTS 75

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

void drawLine(vector<Color>& pixMap, const Point& p1, const Point& p2, const Color& color) {
    double x1 = p1.x * IMAGE_SIZE, y1 = p1.y * IMAGE_SIZE, x2 = p2.x * IMAGE_SIZE,
           y2 = p2.y * IMAGE_SIZE;
    // Bresemham's integer algorithm doesn't draw multiple pixels per column so if its steep it
    // should swap the x and y values to fix this issue
    bool steep = abs(y2 - y1) > abs(x2 - x1);

    if (steep) {
        swap(x1, y1);
        swap(x2, y2);
    }

    // makes sure that x1 is always to the left of x2
    if (x1 > x2) {
        swap(x1, x2);
        swap(y1, y2);
    }

    // changes the x and y values from 0-1 to 0-IMAGE_SIZE
    int dx = (int)abs(x2 - x1), dy = (int)abs(y2 - y1), e = (int)(dy - dx),
        ystep = (y1 < y2) ? 1 : -1, y = (int)y1, maxX = (int)x2;

    for (int x = (int)x1; x < maxX; x++) {
        if (x >= 0 && x < IMAGE_SIZE && y >= 0 && y < 800)
            // x and y values were swapped if steep
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

void drawPoints(vector<Color>& pixMap, const vector<Point>& points,
                const vector<Point*>& convex_points) {
    for (int i = 0; i < points.size(); i++) {
        int x = (int)(points[i].x * IMAGE_SIZE);
        int y = (int)(points[i].y * IMAGE_SIZE);
        pixMap[y * IMAGE_SIZE + x] = Color(false, false, false);
    }
    for (int i = 0; i < convex_points.size(); i++) {
        int x = (int)((*convex_points[i]).x * IMAGE_SIZE);
        int y = (int)((*convex_points[i]).y * IMAGE_SIZE);
        pixMap[y * IMAGE_SIZE + x] = Color(true, false, false);
    }
}

void generateImage(vector<Color>& pixMap) {
    ofstream imageFile;
    imageFile.open("image.ppm");

    imageFile << "P3 " << IMAGE_SIZE << " " << IMAGE_SIZE << " 1" << endl;

    for (int y = IMAGE_SIZE - 1; y >= 0; y--) {
        for (int x = 0; x < IMAGE_SIZE; x++)
            imageFile << (int)pixMap[y * IMAGE_SIZE + x].r << " "
                      << (int)pixMap[y * IMAGE_SIZE + x].g << " "
                      << (int)pixMap[y * IMAGE_SIZE + x].b << " ";
        imageFile << endl;
    }
    imageFile.close();
}

inline double crossProduct(const Point& a, const Point& b) { return a.x * b.y - a.y * b.x; }

// check if angle is convex
inline bool convex(Point a, Point b, Point c) {
    return crossProduct(Point(b.x - a.x, b.y - a.y), Point(c.x - b.x, c.y - b.y)) > 0;
}

Point* origin;

bool compareAngles(const Point& a, const Point& b) {
    if (convex(*origin, a, b)) return true;
    if (convex(*origin, b, a)) return false;
    //compare the distances if angle is the same
    return abs(a.x - origin->x) + abs(a.y - origin->y) <
           abs(b.x - origin->x) + abs(b.y - origin->y);
}
vector<Point*>& convexHull(vector<Point>& points) {
    int n = points.size();
    if (n < 3) return *new vector<Point*>{&points[0], &points[1]};

    int lowest = 0;
    for (int i = 1; i < n; i++) {
        if (points[i].y > points[lowest].y) {
            lowest = i;
        }
    }
    swap(points[0], points[lowest]);
    origin = &points[0];

    // sort all points but origin
    sort(points.begin() + 1, points.end(), compareAngles);

    // used as stack
    vector<Point*>* st = new vector<Point*>;
    st->push_back(&points[0]);
    st->push_back(&points[1]);

    int i = 2;
    while (i < n) {
        Point* a = (*st)[st->size() - 2];
        Point* b = (*st)[st->size() - 1];
        Point* c = &points[i];
        if (convex(*a, *b, *c)) {
            st->push_back(c);
            i++;
        } else {
            st->pop_back();
            if (st->size() < 2) {
                st->push_back(c);
                i++;
            }
        }
    }

    // remove mid if 3 collinear points are found
    Point* a = (*st)[st->size() - 2];
    Point* b = (*st)[st->size() - 1];
    Point* c = (*st)[0];
    if (!convex(*a, *b, *c)) {
        st->pop_back();
    }
    return *st;
}

int main() {
    srand((unsigned)time(0));
    vector<Point>* points = new vector<Point>(NUM_POINTS);
    auto begin = chrono::steady_clock::now();
    vector<Point*>& convex_points = convexHull(*points);
    auto end = chrono::steady_clock::now();
    cout << "Number of Points: " << NUM_POINTS << "\n"
         << "Elapsed Time: " << chrono::duration_cast<chrono::milliseconds>(end - begin).count()
         << "\n";

    vector<Color>* pixMap = new vector<Color>(IMAGE_SIZE * IMAGE_SIZE);
    for (auto il = convex_points.begin(); il != convex_points.end() - 1; il += 1) {
        drawLine(*pixMap, **il, **(il + 1), Color(true, false, true));
    }
    drawLine(*pixMap, *convex_points[0], *convex_points[convex_points.size() - 1],
             Color(true, false, true));
    drawPoints(*pixMap, *points, convex_points);
    generateImage(*pixMap);

    return 0;
}