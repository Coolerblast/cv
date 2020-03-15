// Name: Alex Lin
// Date: 02/23/2020

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

    //imshow("Originial", src);

    // Mat morph = src.clone();
    // for (int r = 1; r <= 2; r++) {
    //     Mat kernel = getStructuringElement(MORPH_ELLIPSE, Size(2 * r + 1, 2 * r + 1));
    //     morphologyEx(morph, morph, MORPH_CLOSE, kernel);
    //     morphologyEx(morph, morph, MORPH_OPEN, kernel);
    // }
    // imshow("Smoothing", morph);

    // // Mat shifted;
    // // pyrMeanShiftFiltering(src, shifted, 21, 51);
    // // imshow("Mean shifting", shifted);

    // Mat gray;
    // cvtColor(morph, gray, COLOR_BGR2GRAY);

    // /* take morphological gradient */
    // Mat mgrad;
    // Mat kernel = getStructuringElement(MORPH_ELLIPSE, Size(3, 3));
    // morphologyEx(morph, mgrad, MORPH_GRADIENT, kernel);
    // imshow("Grad", mgrad);
    // Mat ch[3], merged;
    // /* split the gradient image into channels */
    // split(mgrad, ch);
    // /* apply Otsu threshold to each channel */
    // threshold(ch[0], ch[0], 0, 255, THRESH_BINARY | THRESH_OTSU);
    // threshold(ch[1], ch[1], 0, 255, THRESH_BINARY | THRESH_OTSU);
    // threshold(ch[2], ch[2], 0, 255, THRESH_BINARY | THRESH_OTSU);
    // /* merge the channels */
    // merge(ch, 3, merged);
    // imshow("Merged", merged);
    // Mat gray, binary;
    // cvtColor(merged, gray, COLOR_BGR2GRAY);
    // threshold(gray, binary, 1, 255, THRESH_BINARY);
    // imshow("Binary", binary);

    // morphologyEx(binary, binary, MORPH_CLOSE, getStructuringElement(MORPH_RECT, Size(5,5)));
    // imshow("MorphologyEx CLOSE", binary);

    // waitKey();
    Mat gray;
    cvtColor(src, gray, COLOR_BGR2GRAY);
    GaussianBlur(gray, gray, Size(5, 5), 20, 20);
    // medianBlur(gray, gray, 5);
    //imshow("Gaussian Blur", gray);

    // Mat thresh;
    // adaptiveThreshold(gray, thresh, 255, ADAPTIVE_THRESH_GAUSSIAN_C, THRESH_BINARY_INV, 11, 1);
    // imshow("Adaptive Threshold", thresh);
    // Mat morph;
    // morphologyEx(thresh, morph, MORPH_CLOSE, getStructuringElement(MORPH_ELLIPSE, Size(7, 7)));
    // imshow("Morphology CLOSE", morph);

    // waitKey();

    // vector<vector<Point>> contours;
    // vector<Vec4i> hierarchy;
    // findContours(morph, contours, hierarchy, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);

    // Mat filled = Mat::zeros(src.rows, src.cols, CV_8UC1);
    // drawContours(filled, contours, -1, Scalar(255), -1);
    // imshow("Contour filling", filled);

    // // TO DO LATER
    // Mat dist;
    // distanceTransform(filled, dist, DIST_L2, 3);
    // normalize(dist, dist, 0.0, 1.0, NORM_MINMAX);
    // imshow("Distance", dist);

    vector<Vec3f> circles;
    HoughCircles(gray, circles, HOUGH_GRADIENT, 1, src.cols / 55, 80, 38, src.cols / 80,
                 src.cols / 16);

    float min = circles[0][2], max;
    Mat radii(circles.size(), 1, CV_32F), coin_labels, centers;
    radii.at<float>(0) = circles[0][2];
    for (int i = 1; i < circles.size(); i++) {
        radii.at<float>(i) = circles[i][2];
        if (circles[i][2] < min)
            min = circles[i][2];
        else if (circles[i][2] > max)
            max = circles[i][2];
    }
    int k = max / min > 2.2 ? 4 : 3;
    kmeans(radii, k, coin_labels,
           TermCriteria(TermCriteria::MAX_ITER | TermCriteria::EPS, 10, .01), 10,
           KMEANS_PP_CENTERS, centers);

    int k_rad_index[3];
    map<float, int> k_rad_sorted;
    // cout << coin_labels.total() << endl;
    for (int i = 0; i < k; i++) {
        k_rad_index[i] =
            find(coin_labels.ptr<int>(0), coin_labels.ptr<int>(coin_labels.total()), i) -
            coin_labels.ptr<int>(0);
        k_rad_sorted[radii.at<float>(k_rad_index[i])] = i;
        // printf("Index: %d, radius: %f\n", k_rad_index[i], radii.at<float>(k_rad_index[i]));
    }
    // for (auto t : k_rad_sorted) {
    //     printf("K: %d, radius: %.2f\n", t.second, t.first);
    // }


    float total = 0.0f;
    for (int i = 0; i < circles.size(); i++) {
        Point center = Point(circles[i][0], circles[i][1]);
        circle(src, center, 1, Scalar(0, 100, 100), 3, LINE_AA);
        Scalar clr;
        map<float, int>::iterator iter = k_rad_sorted.begin();
        if (coin_labels.at<int>(i) == iter++->second) {
            clr = Scalar(0, 0, 255);
            total += 0.01f;
        } else if (coin_labels.at<int>(i) == iter++->second) {
            clr = Scalar(0, 255, 255);
            total += 0.05f;
        } else if (coin_labels.at<int>(i) == iter++->second) {
            clr = Scalar(255, 0, 255);
            total += 0.25f;
        } else if (coin_labels.at<int>(i) == iter++->second) {
            clr = Scalar(0, 0, 0);
            total += 0.50f;
        }
        // printf("Total: %.2f\n", total);
        circle(src, center, circles[i][2], clr, 3, LINE_AA);
    }

    printf("Total: %.2f\n", total);
    imwrite("detected_coins.jpg", src);
    // Show results
    imshow("Detected Circles", src);
    // Wait and Exit
    waitKey();
    return 0;
}