//Name: Alex Lin
//Date: 11/10/2019
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_map>

#define IMAGE_SIZE 800
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
    Point(double x, double y) {
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

void drawPoint(Color pixMap[], const Point points[], const Point closestPair[]) {
    for (int i = 0; i < NUM_POINTS; i++) {
        int x = (int)(points[i].x * IMAGE_SIZE);
        int y = (int)(points[i].y * IMAGE_SIZE);
        pixMap[y * IMAGE_SIZE + x] = Color(false, false, false);
    }
    for (int i = 0; i < 2; i++) {
        int x = (int)(closestPair[i].x * IMAGE_SIZE);
        int y = (int)(closestPair[i].y * IMAGE_SIZE);
        pixMap[y * IMAGE_SIZE + x] = Color(true, false, false);
    }
}

void generateImage(const Color pixMap[]) {
    ofstream imageFile;
    imageFile.open("image.ppm");

    imageFile << "P3 " << IMAGE_SIZE << " " << IMAGE_SIZE << " 1" << endl;

    for (int y = IMAGE_SIZE - 1; y >= 0; y--) {
        for (int x = 0; x < IMAGE_SIZE; x++)
            imageFile << (int)pixMap[y * IMAGE_SIZE + x].r << " " << (int)pixMap[y * IMAGE_SIZE + x].g << " " << (int)pixMap[y * IMAGE_SIZE + x].b << " ";
        imageFile << endl;
    }
}

inline double distance(const Point& p1, const Point& p2) {
    return sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
}

template <typename T>
void random_shuffle(T data[], int n) {
    for (int i = 0; i < n - 1; i++)
        swap(data[i], data[rand() % (i + 1)]);
}

void makeDictionary(unordered_map<int, const Point*>& dict, const Point* points, int n, int num_squares, const double& side_length) {
    for (int i = 0; i <= n; i++) {
        dict[(int)(points[i].x / side_length) * num_squares + (int)(points[i].y / side_length)] = &points[i];
    }
}

Point* closestDictionary(const Point* points) {
    double min_distance = distance(points[0], points[1]);
    Point* closestPair = new Point[2]{points[0], points[1]};

    if (min_distance != 0) {
        double side_length = min_distance / 2;
        int num_squares = (int)(1 / side_length);
        unordered_map<int, const Point*>* dict = new unordered_map<int, const Point*>();
        makeDictionary(*dict, points, 1, num_squares, side_length);

        int x, y, pos;
        double temp_distance;

        const Point* temp;

        bool update = false;
        for (int i = 2; i < NUM_POINTS; i++) {
            x = (int)(points[i].x / side_length);
            y = (int)(points[i].y / side_length);
            for (int ty = y - 2; ty <= y + 2; ty++) {
                pos = ty * num_squares + x;
                for (int loc = pos - 2; loc <= pos + 2; loc++) {
                    if (dict->find(loc) != dict->end()) {
                        temp = dict->at(loc);
                        temp_distance = distance(points[i], *temp);
                        if (temp_distance < min_distance) {
                            min_distance = temp_distance;
                            update = true;
                        }
                    }
                }
            }
            if (update) {
                update = false;
                side_length = min_distance / 2;
                num_squares = (int)(1 / side_length);
                closestPair[0] = points[i];
                closestPair[1] = *temp;
                if (min_distance == 0)
                    break;
                dict->clear();
                makeDictionary(*dict, points, i, num_squares, side_length);
            } else {
                (*dict)[y * num_squares + x] = &points[i];
            }
        }
        delete dict;
    }
    return closestPair;
}

int main(int argc, char** argv) {
    NUM_POINTS = atoi(argv[1]);
    srand((unsigned)time(0));

    Point* points = new Point[NUM_POINTS];

    cout << "Number of Points:\t" << NUM_POINTS << endl;

    random_shuffle(points, NUM_POINTS);

    chrono::steady_clock::time_point begin, end;

    begin = chrono::steady_clock::now();

    Point* closestPair = closestDictionary(points);

    end = chrono::steady_clock::now();
    cout << "Closest Pair(s):\t" << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "ms" << endl;

    cout << closestPair[0] << " " << closestPair[1] << endl;

    // Color* pixMap = new Color[IMAGE_SIZE * IMAGE_SIZE];
    // drawPoint(pixMap, points, closestPair);
    // generateImage(pixMap);

    return 0;
}