// Name: Alex Lin
// Date: 4/2/2020

#define _USE_MATH_DEFINES
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

int main(int argc, char** argv) {
    VideoCapture cap("withChessBoard.mov");
    if (!cap.isOpened()) {
        cout << "Error opening video stream or file" << endl;
        return -1;
    }

    Size videoSize(cap.get(CAP_PROP_FRAME_WIDTH), cap.get(CAP_PROP_FRAME_HEIGHT));
    VideoWriter outputVideo("output.avi", VideoWriter::fourcc('M', 'J', 'P', 'G'),
                            cap.get(CAP_PROP_FPS), videoSize);

    if (!outputVideo.isOpened()) {
        cout << "Could not open the output video for write: "
             << "output.avi" << endl;
        return -1;
    }

    vector<Point3f> object_points{
        Point3f(-1, -1, 0),  // 0
        Point3f(1, -1, 0),   // 1
        Point3f(-1, 1, 0),   // 2
        Point3f(1, 1, 0),    // 3
        Point3f(-1, -1, 2),  // 4
        Point3f(1, -1, 2),   // 5
        Point3f(-1, 1, 2),   // 6
        Point3f(1, 1, 2),    // 7
    };

    vector<pair<int, int>> object_lines{
        pair<int, int>(0, 1),  //
        pair<int, int>(1, 3),  //
        pair<int, int>(3, 2),  //
        pair<int, int>(0, 2),  //
        pair<int, int>(4, 5),  //
        pair<int, int>(5, 7),  //
        pair<int, int>(7, 6),  //
        pair<int, int>(4, 6),  //
        pair<int, int>(0, 4),  //
        pair<int, int>(1, 5),  //
        pair<int, int>(2, 6),  //
        pair<int, int>(3, 7),
    };

    Size chessboardSize(7, 7);
    vector<Point3f> checker_object_pts;
    for (int i = -chessboardSize.width / 2; i <= chessboardSize.width / 2; i++)
        for (int j = -chessboardSize.height / 2; j <= chessboardSize.height / 2; j++)
            checker_object_pts.push_back(Point3f(i, j, 0));

    int frames = cap.get(CAP_PROP_FRAME_COUNT);
    for (int i = 1; i <= frames / 3; i++) {
        printf_s("\rGenerating frame %d of %d...", i, frames);
        Mat frame, gray;
        cap >> frame;
        if (frame.empty()) break;
        flip(frame, frame, 0);
        cvtColor(frame, gray, COLOR_BGR2GRAY);

        vector<Point2f> corners;
        bool found = findChessboardCorners(frame, chessboardSize, corners,
                                           CALIB_CB_ADAPTIVE_THRESH | CALIB_CB_NORMALIZE_IMAGE);
        if (found) {
            Mat cameraMatrix, distCoeffs, R, T;

            cv::TermCriteria criteria(TermCriteria::EPS | TermCriteria::MAX_ITER, 30, 0.001);
            cv::cornerSubPix(gray, corners, cv::Size(11, 11), cv::Size(-1, -1), criteria);
            drawChessboardCorners(frame, chessboardSize, corners, found);
            calibrateCamera(vector<vector<Point3f>>{checker_object_pts},
                            vector<vector<Point2f>>{corners}, cv::Size(gray.rows, gray.cols),
                            cameraMatrix, distCoeffs, R, T);
            Mat rvec;
            Mat tvec;

            solvePnP(checker_object_pts, corners, cameraMatrix, distCoeffs, rvec, tvec);

            vector<Point2f> output_points;
            projectPoints(checker_object_pts, rvec, tvec, cameraMatrix, distCoeffs, output_points);
            // for (const pair<int, int>& lines : object_lines) {
            //     line(frame, output_points[lines.first], output_points[lines.second],
            //          Scalar(0, 255, 0), 10);
            // }
            for (Point2f pt : output_points)
                circle(frame, pt, 5, Scalar(0, 255, 0), 2);
        }

        imshow("Frame", frame);
        char c = (char)waitKey(25);
        if (c == 27) break;
        outputVideo << frame;
    }

    cap.release();
    outputVideo.release();
    cv::destroyAllWindows();

    return 0;
}