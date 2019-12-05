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

string helpMessage =
    "Syntax:\tprogram.exe [input file] -arguments...\n"
    "Options and arguments:\n"
    "-o dir\t\t: specify file output location\n"
    "-g int dbl\t: specify gradient size and sigma value [gradient sigma]\n"
    "-t dbl dbl\t: specify threshold min and max percentages [min max]\n"
    "-h\t\t: displays this help message\n"
    "-grayscale\t: output grayscale image\n"
    "-blur\t\t: output image after gaussian blur has been applied\n"
    "-sobel\t\t: output image after sobel operator has been applied\n"
    "-nonmax\t\t: output image after non max suppression has been applied\n"
    "-threshold\t: output image after values have been thresholded\n";

struct Color {
    unsigned char r, g, b;
    Color(unsigned char r = 0, unsigned char g = 0, unsigned char b = 0) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
};

void generateImage(vector<Color>& pixMap, string filename) {
    ofstream imageFile;
    imageFile.open(filename);

    imageFile << "P3 " << WIDTH << " " << HEIGHT << " 255" << endl;

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++)
            imageFile << (int)pixMap[y * WIDTH + x].r << " " << (int)pixMap[y * WIDTH + x].g << " "
                      << (int)pixMap[y * WIDTH + x].b << " ";
        imageFile << endl;
    }
    imageFile.close();
}

// Usinzg unsigned char because it stores values from 0 to 255 and RGB values can only go to 255
bool loadImage(vector<Color>& image, const string& fileDir) {
    if (fileDir.substr(fileDir.find_last_of(".") + 1) != "ppm") {
        cerr << "ERROR: " << fileDir << " is not a supported image type. Supported types: ppm"
             << endl;
        return false;
    }
    ifstream imageFile(fileDir, ios::in | ios::binary);
    string imageType;
    short int colorSize;
    imageFile >> imageType >> WIDTH >> HEIGHT >> colorSize;
    image.resize(WIDTH * HEIGHT);
    unsigned char r, g, b;
    int i = 0;
    if (imageType.compare("P3") == 0)
        while (i < image.size())
            if (imageFile >> r >> g >> b)
                image[i++] = Color(r, g, b);
            else
                cerr << "ERROR: PPM image could not be read correctly.\n";
    else if (imageType.compare("P6") == 0) {
        unsigned char buffer[3];
        while (i < image.size())
            if (imageFile.read((char*)buffer, 3))
                image[i++] = Color(buffer[0], buffer[1], buffer[2]);
            else
                cerr << "ERROR: PPM image could not be read correctly.\n";
    }
    colorSize += 1;
    // This normalizes all the input PPMs, scaling it up to 8 bit color values
    if (colorSize != 256)
        for (int i = 0; i < image.size(); i++) {
            image[i].r = (image[i].r + 1) * (256 / colorSize) - 1;
            image[i].g = (image[i].g + 1) * (256 / colorSize) - 1;
            image[i].b = (image[i].b + 1) * (256 / colorSize) - 1;
        }
    imageFile.close();
    return true;
}

void grayScaleImage(vector<Color>& image) {
    for (int i = 0; i < image.size(); i++) {
        unsigned char gray = (unsigned char)((image[i].r + image[i].g + image[i].b) / 3);
        image[i].r = gray;
        image[i].g = gray;
        image[i].b = gray;
    }
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

void gaussianBlur(vector<Color>& image, const int& n, const double& sigma) {
    vector<double> filter = generateGaussianFilter(n, sigma);
    vector<Color> blur(image.size());

    int mid = n / 2;
    int min = mid - n + 1;
    double sum, rVal, gVal, bVal, filterVal;
    Color color;

    for (int y = 0; y < HEIGHT; y++)
        for (int x = 0; x < WIDTH; x++) {
            sum = rVal = gVal = bVal = 0;
            for (int fy = min; fy <= mid; fy++)
                for (int fx = min; fx <= mid; fx++)
                    if (y + fy >= 0 && y + fy < HEIGHT && x + fx >= 0 && fx < WIDTH) {
                        filterVal = filter[(fy - min) * n + fx - min];
                        color = image[(y + fy) * WIDTH + x + fx];
                        sum += filterVal;
                        rVal += filterVal * color.r;
                        gVal += filterVal * color.g;
                        bVal += filterVal * color.b;
                    }
            rVal /= sum;
            gVal /= sum;
            bVal /= sum;
            blur[y * WIDTH + x] =
                Color((unsigned char)rVal, (unsigned char)gVal, (unsigned char)bVal);
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

vector<double> applySobelOperator(vector<Color>& image, vector<unsigned char>& out,
                                  vector<double>& g, vector<double>& angles) {
    short int sX, sY, gX, gY, gXr, gYr, gXg, gYg, gXb, gYb;
    Color color;
    double gMax = 0;
    for (int y = 1; y < HEIGHT - 1; y++)
        for (int x = 1; x < WIDTH - 1; x++) {
            gX = gY = gXr = gYr = gXg = gYg = gXb = gYb = 0;
            for (int fy = -1; fy <= 1; fy++)
                for (int fx = -1; fx <= 1; fx++) {
                    color = image[(y + fy) * WIDTH + x + fx];
                    sX = sobelKernalX[(fy + 1) * 3 + (fx + 1)];
                    sY = sobelKernalY[(fy + 1) * 3 + (fx + 1)];
                    gXr += sX * color.r;
                    gYr += sY * color.r;
                    gXg += sX * color.g;
                    gYg += sY * color.g;
                    gXb += sX * color.b;
                    gYb += sY * color.b;
                }
            gX = abs(gXr) > abs(gXg) ? gXr : abs(gXg) > abs(gXb) ? gXg : gXb;
            gY = abs(gYr) > abs(gYg) ? gYr : abs(gYg) > abs(gYb) ? gYg : gYb;
            g[y * WIDTH + x] = sqrt(gX * gX + gY * gY);
            if (g[y * WIDTH + x] > gMax) gMax = g[y * WIDTH + x];
            angles[y * WIDTH + x] = atan2(gY, gX);
        }
    for (int i = 0; i < out.size(); i++) out[i] = (unsigned char)((g[i] / gMax) * 255);
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
    tMaxRatio = (avg) / 255;
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

    bool getInputFile(string& inputFileDir) const {
        if (!regex_match(tokens[0], re)) {
            inputFileDir = tokens[0];
            return true;
        }
        return false;
    }

   private:
    vector<string> tokens;
    regex re;
};

bool detectEdges(const InputParser& input, vector<Color>& image, string& outputFilename) {
    string stMinRatio, stMaxRatio, sgSize, sgSigma, stage;

    if (!(input.getCmdOption("-g", {&sgSize, &sgSigma}) &&
          input.getCmdOption("-t", {&stMinRatio, &stMaxRatio}) &&
          input.getCmdOption("-o", {&outputFilename})))
        return false;

    //! change this later to null by default!
    string inputFileDir = "keyboard.ppm";
    input.getInputFile(inputFileDir);
    if (!loadImage(image, inputFileDir)) return false;

    if (input.cmdOptionExists("-grayscale")) {
        grayScaleImage(image);
        return true;
    }

    int gSize = 5;
    double gSigma = 1.0;
    if (!sgSize.empty() && !sgSigma.empty()) {
        gSize = stoi(sgSize);
        gSigma = stod(sgSigma);
    }

    gaussianBlur(image, gSize, gSigma);
    if (input.cmdOptionExists("-blur")) return true;

    vector<unsigned char> imageGS(image.size());
    vector<double> g(image.size());
    vector<double> angles(image.size());
    applySobelOperator(image, imageGS, g, angles);
    if (!input.cmdOptionExists("-sobel")) {
        applyNonMaxSuppression(imageGS, angles);
        if (!input.cmdOptionExists("-nonmax")) {
            double tMinRatio, tMaxRatio;
            if (!stMinRatio.empty() && !stMaxRatio.empty()) {
                tMinRatio = stod(stMinRatio);
                tMaxRatio = stod(stMaxRatio);
            } else
                calculateThresholdValues(g, tMinRatio, tMaxRatio);

            vector<int> strong = threshold(imageGS, tMinRatio, tMaxRatio);
            if (!input.cmdOptionExists("-threshold")) defineEdges(imageGS, strong, 255);
        }
    }
    for (int i = 0; i < image.size(); i++) image[i] = Color(imageGS[i], imageGS[i], imageGS[i]);

    return true;
}

int main(int argc, char** argv) {
    InputParser input(argc, argv);
    if (input.cmdOptionExists("-h")) {
        cout << helpMessage << endl;
        return 0;
    }
    string outputFilename = "output.ppm";
    vector<Color> image;

    if (!detectEdges(input, image, outputFilename)) return -1;
    generateImage(image, outputFilename);
    return 0;
}