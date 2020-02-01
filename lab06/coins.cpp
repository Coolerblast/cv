#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

int main(int argc, char** argv) {
    const char* filename = argc >= 2 ? argv[1] : "coinseasy.jpg";
    Mat src = imread(samples::findFile(filename), IMREAD_COLOR);
    // Check if image is loaded fine
    if (src.empty()) {
        printf(" Error opening image\n");
        printf(" Program Arguments: [image_name -- default %s] \n", filename);
        return -1;
    }
    Mat gray;
    cvtColor(src, gray, COLOR_BGR2GRAY);
    GaussianBlur(gray, gray, Size(5, 5), 5);

    vector<Vec3f> circles;
    HoughCircles(gray, circles, HOUGH_GRADIENT, 1, gray.cols / 50, 100, 30, gray.cols / 80,
                 gray.cols / 16);
    for (const Vec3i& c : circles) {
        Point center = Point(c[0], c[1]);
        circle(src, center, 1, Scalar(0, 100, 100), 3, LINE_AA);
        circle(src, center, c[2], Scalar(255, 0, 255), 3, LINE_AA);
    }
    // Show results
	imshow("Canny Circles", gray);
    imshow("Detected Circles", src);
    // Wait and Exit
    waitKey();
    return 0;
}