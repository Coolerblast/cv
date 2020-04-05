// Name: Alex Lin
// Date: 3/15/2020

#define _USE_MATH_DEFINES
#include <fstream>
#include <opencv2/opencv.hpp>

using namespace std;

struct point_comp {
    bool operator()(const pair<int, int>& a, const pair<int, int>& b) const {
        return less_comp(minmax(a.first, a.second), minmax(b.first, b.second));
    }
    less<std::pair<int, int>> less_comp;
};

double inline distance(cv::Point3f a, cv::Point3f b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

void gluPerspective(const float& fov, const float& aspectRatio, const float& n, const float& f,
                    float& t, float& b, float& l, float& r) {
    float scale = tan(fov * 0.5 * M_PI / 180) * n;
    r = aspectRatio * scale, l = -r;
    t = scale, b = -t;
}

void glFrustum(const float& t, const float& b, const float& l, const float& r, const float& n,
               const float& f, cv::Mat& M) {
    M = (cv::Mat1f(4, 4) << 2 * n / (r - l), 0, 0, 0,                                     //
         0, 2 * n / (t - b), 0, 0,                                                        //
         (r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -2 * f * n / (f - n),  //
         0, 0, -1, 0);
}

void multPointMatrix(const cv::Point3f& in, cv::Point3f& out, const cv::Mat& M) {
    cv::Mat temp = (cv::Mat1f(4, 1) << in.x, in.y, in.z, 1);
    temp = M * temp;
    out = cv::Point3f(temp.at<float>(0, 0), temp.at<float>(1, 0), temp.at<float>(2, 0));
    if (temp.at<float>(0, 3) != 1) {
        out.x /= temp.at<float>(3, 0);
        out.y /= temp.at<float>(3, 0);
        out.z /= temp.at<float>(3, 0);
    }
    // cout << temp.at<float>(3, 0) << " " << out << endl;
}

int main(int argc, char** argv) {
    cout << "Input name of platonic solid: ";
    string name;
    cin >> name;

    ifstream infile(name + ".txt", ios::in);

    vector<cv::Point3f> points;
    double ix, iy, iz;
    while (infile >> ix >> iy >> iz) points.emplace_back(cv::Point3f(ix, iy, iz));

    double mindist = INFINITY;
    for (int i = 1; i < points.size(); ++i) {
        mindist = min(mindist, distance(points[0], points[i]));
        // cout << mindist << endl;
    }

    set<pair<int, int>, point_comp> edges;
    for (int i = 0; i < points.size(); ++i)
        for (int j = 0; j < points.size(); ++j)
            if (distance(points[i], points[j]) - mindist < 0.0001)
                edges.emplace(pair<int, int>(i, j));
    // for (auto a : edges) {
    //     cout << a.first << " " << a.second << endl;
    //     cout << points[a.first] << " " << points[a.second] << endl;
    // }

    cv::Mat Mproj, worldCamera = (cv::Mat1f(4, 4) << 1, 0, 0, 0,  //
                                  0, 1, 0, 0,                     //
                                  0, 0, -1, -4,                   //
                                  0, 0, 0, 1                      //
                   );

    int width = 512, height = 512;
    cv::Size size(width, height);
    cv::VideoWriter outputVideo("output.avi", cv::VideoWriter::fourcc('M','J','P','G'), 20,
                                size);

    if (!outputVideo.isOpened()) {
        cout << "Could not open the output video for write: " << "output.avi" << endl;
        return -1;
    }

    cv::Scalar color(255, 0, 0);

    float fov = 90;
    float near = 0.1;
    float far = 100;
    float aspectRatio = width / (float)height;
    float b, t, l, r;
    gluPerspective(fov, aspectRatio, near, far, t, b, l, r);
    glFrustum(t, b, l, r, near, far, Mproj);


    cout << "Generating video..." << endl;
    for (float angle = 0; angle < 30; angle+=0.1) {
        cv::Mat frame = cv::Mat::zeros(size, CV_8UC3);
        vector<cv::Point2i> projected(points.size());
        for (int i = 0; i < points.size(); i++) {
            cv::Point3f vertCam, projVert;
            multPointMatrix(points[i], vertCam, worldCamera);
            multPointMatrix(vertCam, projVert, Mproj);
            projected[i] =
                cv::Point2i(min(width - 1, (int)((projVert.x + 1) * 0.5 * width)),
                            min(height - 1, (int)((1 - (projVert.y + 1) * 0.5) * height)));
        }

        for (auto elem : edges) {
            cv::Point3f vertCamera, projVert;
            // cout << projected[elem.first] << projected[elem.second] << endl;
            cv::line(frame, projected[elem.first], projected[elem.second], color, 3);
        }

        float rad = angle * M_PI / 180.0;
        for (auto& vert : points) {
            float old = vert.x;
            vert.x = vert.x * cos(rad) - vert.z * sin(rad);
            vert.z = vert.z * cos(rad) + old * sin(rad);
        }

        // cv::imshow("Hello", frame);
        // cv::waitKey(0);

        outputVideo << frame;
    }

    outputVideo.release();

    cout << "Done!" << endl;

    return 0;
}