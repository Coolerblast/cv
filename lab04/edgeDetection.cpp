#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

using namespace std;

int WIDTH, HEIGHT;

void generateImage(vector<unsigned char>& pixMap, string filename) {
    ofstream imageFile;
    imageFile.open(filename);

    imageFile << "P3 " << WIDTH << " " << HEIGHT << " 255" << endl;

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++)
            imageFile << (int)pixMap[y * WIDTH + x] << " " << (int)pixMap[y * WIDTH + x] << " "
                      << (int)pixMap[y * WIDTH + x] << " ";
        imageFile << endl;
    }
    imageFile.close();
}

// Usinzg unsigned char because it stores values from 0 to 255 and RGB values can only go to 255
vector<unsigned char> loadGrayScaleImage(const char* fileDir) {
    ifstream imageFile(fileDir, ios::in | ios::binary);
    string imageType;
    short int colorSize;
    imageFile >> imageType >> WIDTH >> HEIGHT >> colorSize;
    vector<unsigned char> image(WIDTH * HEIGHT);
    unsigned char r, g, b;
    int i = 0;
    if (imageType.compare("P3") == 0)
        while (i < image.size())
            if (imageFile >> r >> g >> b)
                image[i++] = (unsigned char)((r + g + b) / 3);
            else
                cerr << "ERROR: PPM image could not be read correctly.\n";
    else if (imageType.compare("P6") == 0) {
        unsigned char buffer[3];
        while (i < image.size())
            if (imageFile.read((char*)buffer, 3))
                image[i++] = (unsigned char)((buffer[0] + buffer[1] + buffer[2]) / 3);
            else
                cerr << "ERROR: PPM image could not be read correctly.\n";
    }
    colorSize += 1;
    // This normalizes all the input PPMs, scaling it up to 8 bit color values
    if (colorSize != 256)
        for (int i = 0; i < image.size(); i++) image[i] = (image[i] + 1) * (256 / colorSize) - 1;
    imageFile.close();
    return image;
}

vector<double> generateGaussianFilter(const int& n, const double& sigma) {
    vector<double> filter(n * n);
    double d = 2 * sigma * sigma;
    int mid = n / 2;
    int min = mid - n + 1;
    double sum = 0;
    for (int y = min; y <= mid; y++)
        for (int x = min; x <= mid; x++) {
            filter[(y - min) * n + x - min] = exp(-(x * x + y * y) / d) / (M_PI * d);
            sum += filter[(y - min) * n + x - min];
        }
    for (int i = 0; i < filter.size(); i++) filter[i] /= sum;

    // for (int y = min; y <= mid; y++) {
    //     for (int x = min; x <= mid; x++) cout << filter[(y - min) * n + x - min] << " ";
    //     cout << endl;
    // }
    return filter;
}

void gaussianBlur(vector<unsigned char>& image, const int& n, const double& sigma) {
    vector<double> filter = generateGaussianFilter(n, sigma);
    vector<unsigned char> blur(image.size());

    int mid = n / 2;
    int min = mid - n + 1;
    double sum, color;

    for (int y = 0; y < HEIGHT; y++)
        for (int x = 0; x < WIDTH; x++) {
            sum = 0;
            color = 0;
            for (int fy = min; fy <= mid; fy++)
                for (int fx = min; fx <= mid; fx++)
                    if (y + fy >= 0 && y + fy < HEIGHT && x + fx >= 0 && fx < WIDTH) {
                        sum += filter[(fy - min) * n + fx - min];
                        color +=
                            filter[(fy - min) * n + fx - min] * image[(y + fy) * WIDTH + x + fx];
                    }
            color /= sum;
            blur[y * WIDTH + x] = (unsigned char)color;
        }
    image = blur;
}

int sobelKernalX[] = {
    -1, 0, 1,  //
    -2, 0, 2,  //
    -1, 0, 1   //
};

int sobelKernalY[] = {
    1,  2,  1,  //
    0,  0,  0,  //
    -1, -2, -1  //
};

vector<double> applySobelOperator(vector<unsigned char>& image, vector<double>& angles) {
    vector<double> g(image.size());
    short int gX, gY;
    double gMax = 0;
    for (int y = 1; y < HEIGHT - 1; y++)
        for (int x = 1; x < WIDTH - 1; x++) {
            gX = 0;
            gY = 0;
            for (int fy = -1; fy <= 1; fy++)
                for (int fx = -1; fx <= 1; fx++) {
                    gX += sobelKernalX[(fy + 1) * 3 + (fx + 1)] * image[(y + fy) * WIDTH + x + fx];
                    gY += sobelKernalY[(fy + 1) * 3 + (fx + 1)] * image[(y + fy) * WIDTH + x + fx];
                }
            g[y * WIDTH + x] = sqrt(gX * gX + gY * gY);
            if (g[y * WIDTH + x] > gMax) gMax = g[y * WIDTH + x];
            angles[y * WIDTH + x] = atan2(gY, gX);
        }
    for (int i = 0; i < image.size(); i++) image[i] = (unsigned char)((g[i] / gMax) * 255);
    return g;
}

void applyNonMaxSuppression(vector<unsigned char>& image, const vector<double>& angles) {
    vector<unsigned char> edges(image.size());
    vector<float> dir(angles.size());
    for (int i = 0; i < angles.size(); i++) {
        dir[i] = fmod(angles[i] + M_PI, M_PI) / M_PI * 8;
    }

    for (int y = 1; y < HEIGHT - 1; y++)
        for (int x = 1; x < WIDTH - 1; x++) {
            for (int fy = -1; fy <= 1; fy++)
                for (int fx = -1; fx <= 1; fx++) {
                    int i = y * WIDTH + x;
                    int tc = i - WIDTH;  // top-center
                    int tl = tc - 1;     // top-left
                    int tr = tc + 1;     // top-right
                    int cl = i - 1;      // center-left
                    int cr = i + 1;      // center-right
                    int bc = i + WIDTH;  // bottom-center
                    int bl = bc - 1;     // bottom-left
                    int br = bc + 1;     // bottom-right

                    if (((dir[i] <= 1 || dir[i] > 7) && image[i] > image[cr] &&
                         image[i] > image[cl]) ||  // 0 deg
                        ((dir[i] <= 3 && dir[i] > 1) && image[i] > image[tr] &&
                         image[i] > image[bl]) ||  // 45 deg
                        ((dir[i] <= 5 && dir[i] > 3) && image[i] > image[tc] &&
                         image[i] > image[bc]) ||  // 90 deg
                        ((dir[i] <= 7 && dir[i] > 5) && image[i] > image[tl] &&
                         image[i] > image[br]))  // 135 deg
                        edges[i] = image[i];
                    else
                        edges[i] = 0;
                }
        }
    image = edges;
}

void calculateThresholdValues(vector<double>& g, double& tMinRatio, double& tMaxRatio) {
    double sum = 0;
    for (double d : g) sum += d;
    double avg = sum / g.size();
    double var = 0;
    for (double d : g) var += (d - avg) * (d - avg);
    var /= g.size();
    double sd = sqrt(var);
    tMinRatio = (sd) / 255;
    // tMaxRatio = (255 - sd * 2) / 255;
    tMaxRatio = (avg + sd / 2) / 255;
    cout << "Average: " << avg << endl
         << "Variance: " << var << endl
         << "Standard Deviation: " << sd << endl
         << "tMin: " << tMinRatio << endl
         << "tMax: " << tMaxRatio << endl;
}

vector<int> threshold(vector<unsigned char>& image, const double& tMinRatio,
                      const double& tMaxRatio) {
    int maxValue = 255;
    int med = maxValue / 2;

    int tMax = maxValue * tMaxRatio;
    int tMin = tMax * tMinRatio;

    vector<int> strong;

    for (int i = 0; i < image.size(); i++)
        if (image[i] > tMax) {
            image[i] = maxValue;
            strong.emplace_back(i);
        } else if (image[i] > tMin)
            image[i] = med;
        else
            image[i] = 0;
    return strong;
}

void defineEdges(vector<unsigned char>& image, vector<int> strong, const int& max) {
    vector<unsigned char> edges(image.size());
    int i;
    vector<int> nedges;
    for (int c : strong) {
        if (edges[c] == 0) {
            edges[c] = max;
            nedges.emplace_back(c);
            do {
                i = nedges.back();
                nedges.pop_back();
                int neighbors[8] = {
                    i - WIDTH - 1,  // top-left
                    i - WIDTH,      // top-center
                    i - WIDTH + 1,  // top-right
                    i - 1,          // center-left
                    i + 1,          // center-right
                    i + WIDTH - 1,  // bottom-left
                    i + WIDTH,      // bottom-center
                    i + WIDTH + 1   // bottom-right
                };
                for (int n : neighbors) {
                    if (image[n] == 127 && edges[n] == 0) {
                        edges[n] = max;
                        nedges.emplace_back(n);
                    }
                }
            } while (!nedges.empty());
        }
    }
    image = edges;
}

class InputParser {
   public:
    InputParser(const int& argc, char** argv) {
        for (int i = 1; i < argc; ++i) this->tokens.push_back(string(argv[i]));
        this->re = regex("-\\w+");
    }

    bool getCmdOption(const string& option, initializer_list<string*> out) const {
        vector<string>::const_iterator itr =
            find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end())
            if (++itr != this->tokens.end())
                for (auto o = out.begin(); o != out.end(); ++o) {
                    if (itr != this->tokens.end() && !regex_match(*itr, this->re))
                        **o = *(itr++);
                    else {
                        cerr << "ERROR: Too few arguments for the option \"" << option << '\"'
                             << endl;
                        return false;
                    }
                }
            else {
                cerr << "ERROR: Too few arguments for the option \"" << option << '\"' << endl;
                return false;
            }
        return true;
    }

    bool cmdOptionExists(const string& option) const {
        return find(this->tokens.begin(), this->tokens.end(), option) != this->tokens.end();
    }

   private:
    vector<string> tokens;
    regex re;
};

bool detectEdges(const int& argc, char** argv, vector<unsigned char>& image,
                 string& outputFilename) {
    InputParser input(argc, argv);
    string stMinRatio, stMaxRatio, sgSize, sgSigma, stage;

    if (!(input.getCmdOption("-g", {&sgSize, &sgSigma}) &&
          input.getCmdOption("-t", {&stMinRatio, &stMaxRatio}) &&
          input.getCmdOption("-o", {&outputFilename})))
        return false;

    if (argc > 1 && regex_match(argv[1], regex("-\\w+")))
        image = loadGrayScaleImage(argv[1]);
    else
        image = loadGrayScaleImage("keyboard.ppm");

    int gSize = 5;
    double gSigma = 1.0;
    if (!sgSize.empty() && !sgSigma.empty()) {
        gSize = stoi(sgSize);
        gSigma = stod(sgSigma);
    }

    gaussianBlur(image, gSize, gSigma);

    vector<double> angles(image.size());
    vector<double> g = applySobelOperator(image, angles);

    if (input.cmdOptionExists("-sobel")) return true;

    applyNonMaxSuppression(image, angles);

    if (input.cmdOptionExists("-nonmax")) return true;

    double tMinRatio, tMaxRatio;
    if (!stMinRatio.empty() && !stMaxRatio.empty()) {
        tMinRatio = stod(stMinRatio);
        tMaxRatio = stod(stMaxRatio);
    } else
        calculateThresholdValues(g, tMinRatio, tMaxRatio);

    vector<int> strong = threshold(image, tMinRatio, tMaxRatio);

    if (input.cmdOptionExists("-threshold")) return true;

    defineEdges(image, strong, 255);

    return true;
}

int main(int argc, char** argv) {
    string outputFilename = "output.ppm";
    vector<unsigned char> image;
    if (!detectEdges(argc, argv, image, outputFilename)) return -1;
    generateImage(image, outputFilename);
    return 0;
}