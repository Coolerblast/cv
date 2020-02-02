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
    GaussianBlur(gray, gray, Size(5, 5), 20);

    vector<Vec3f> circles;
    HoughCircles(gray, circles, HOUGH_GRADIENT, 1, gray.cols / 80, 100, 30, gray.cols / 80,
                 gray.cols / 16);

    Mat radii(circles.size(), 1, CV_32F), coin_labels, centers;
    for (int i = 0; i < circles.size(); i++) {
        radii.at<float>(i) = circles[i][2];
    }
    kmeans(radii, 3, coin_labels,
           TermCriteria(TermCriteria::MAX_ITER | TermCriteria::EPS, 10, .01), 10,
           KMEANS_PP_CENTERS, centers);

    for (int i = 0; i < coin_labels.rows; i++) {
        cout << coin_labels.at<int>(i) << endl;
    }
    float total = 0.0f;
    for (int i = 0; i < circles.size(); i++) {
        Point center = Point(circles[i][0], circles[i][1]);
        circle(src, center, 1, Scalar(0, 100, 100), 3, LINE_AA);
        Scalar clr;
        switch (coin_labels.at<int>(i)) {
            case 0:
                clr = Scalar(0, 0, 255);
                total += 0.01f;
                break;
            case 1:
                clr = Scalar(0, 255, 255);
                total += 0.05f;
                break;
            case 2:
                clr = Scalar(255, 0, 255);
                total += 0.25f;
                break;
            case 3:
                clr = Scalar(0, 0, 0);
                total += 0.50f;
        }
        printf("Total: %.2f\n", total);
        circle(src, center, circles[i][2], clr, 3, LINE_AA);
    }
    printf("Total: %.2f", total);
    // Show results
    imshow("Canny Circles", gray);
    imshow("Detected Circles", src);
    // Wait and Exit
    waitKey();
    return 0;
}