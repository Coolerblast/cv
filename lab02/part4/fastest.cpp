//Name: Alex Lin Period: 4
//Date: 11/12/2019
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

using namespace std;

class Point {
   public:
    double x;
    double y;
    Point(const double& x = rand() / (RAND_MAX + 1.0), const double& y = rand() / (RAND_MAX + 1.0)) {
        this->x = x;
        this->y = y;
    }
    Point(istream& istrm) {
        istrm >> this->x >> this->y;
    }
    friend ostream& operator<<(ostream& ostrm, const Point& a) {
        return ostrm << "X: " << a.x << " Y: " << a.y;
    }
};

inline double distance(const Point& p1, const Point& p2) {
    return sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
}

template <typename T>
void random_shuffle(vector<T>& data) {
    for (int i = 0; i < data.size() - 1; i++)
        swap(data[i], data[rand() % (i + 1)]);
}

void FindClosestPair(vector<Point> points) {
    random_shuffle(points);
    double min_distance = distance(points[0], points[1]);
    vector<const Point*> closestPair{&(points[0]), &(points[1])};
    if (min_distance != 0) {
        double side_length = min_distance / 2;
        int num_squares = (int)(1 / side_length);
        unordered_map<int, const Point*> dict;
        dict[(int)(points[0].x / side_length) * num_squares + (int)(points[0].y / side_length)] = &points[0];
        dict[(int)(points[1].x / side_length) * num_squares + (int)(points[1].y / side_length)] = &points[1];

        int x, y, pos;
        double temp_distance;

        const Point* temp;

        bool update = false;
        for (int i = 2; i < points.size(); i++) {
            x = (int)(points[i].x / side_length);
            y = (int)(points[i].y / side_length);
            for (int ty = y - 2; ty <= y + 2; ty++) {
                pos = ty * num_squares + x;
                for (int loc = pos - 2; loc <= pos + 2; loc++) {
                    if (dict.find(loc) != dict.end()) {
                        temp = dict.at(loc);
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
                closestPair[0] = &points[i];
                closestPair[1] = temp;
                if (min_distance == 0)
                    break;
                dict.clear();
                for (int j = 0; j <= i; j++) {
                    dict[(int)(points[j].x / side_length) * num_squares + (int)(points[j].y / side_length)] = &points[j];
                }
            } else {
                dict[y * num_squares + x] = &points[i];
            }
        }
    }
    cout << *closestPair[0] << "\n"
         << *closestPair[1] << "\n";
}

int main() {
    srand((unsigned)time(NULL));
    ifstream file("points.txt", ifstream::in);
    vector<Point> points;
    int numberpoints;
    file >> numberpoints;
    for (int i = 0; i < numberpoints; i++)
        points.emplace_back(Point(file));
    FindClosestPair(points);
    return 0;
}