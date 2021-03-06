#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

using namespace std;

int WIDTH, HEIGHT, COLOR_SIZE = 255;

string helpMessage =
    "Syntax:\tprogram.exe [input file] -arguments...\n"
    "Options and arguments:\n"
    "-o dir\t\t: specify file output location\n"
    "-b int dbl\t: specify gaussian blur filter size and sigma value [size sigma]\n"
    "-t dbl dbl\t: specify threshold min and max percentages [min max]\n"
    "-c r|g|b\t: specify which color channel(s) to use algorithm on | default is grayscale\n"
    "-h\t\t: displays this help message\n"
    "-compressed\t: output a P6 PPM instead of a P3 PPM\n"
    "-grayscale\t: output a grayscale version of the image\n"
    "-blur\t\t: output image after gaussian blur has been applied\n"
    "-sobel\t\t: output image after sobel operator has been applied\n"
    "-gradient mdl\t: output image of the color gradient using a specified color model [hsl|hsv] "
    "| default is hsv\n"
    "-nonmax\t\t: output image after non max suppression has been applied\n"
    "-threshold\t: output image after values have been thresholded\n";

void HSVtoRGB(const double& h, const double& s, const double& v, unsigned char& r,
              unsigned char& g, unsigned char& b) {
    double angle = fmod(h * (180 / M_PI) + 360, 360);

    auto f = [](int n, const double& h, const double& s,
                const double& v) {  // lambda expression for readablility
        double k = fmod(n + h / 60.0, 6);
        return v - v * s * max(min(k, min(4 - k, 1.0)), 0.0);
    };

    r = f(5, angle, s, v) * 255;
    g = f(3, angle, s, v) * 255;
    b = f(1, angle, s, v) * 255;
}

void HSLtoRGB(const double& h, const double& s, const double& l, unsigned char& r,
              unsigned char& g, unsigned char& b) {
    double angle = fmod(h * (180 / M_PI) + 360, 360);

    auto f = [](int n, const double& h, const double& s,
                const double& l) {  // lambda expression for readablility
        double k = fmod(n + h / 30.0, 12);
        double a = s * min(l, 1 - l);
        return l - a * max(min(k - 3, min(9 - k, 1.0)), -1.0);
    };

    r = f(0, angle, s, l) * 255;
    g = f(8, angle, s, l) * 255;
    b = f(4, angle, s, l) * 255;
}

bool loadImage(vector<unsigned char>& r, vector<unsigned char>& g, vector<unsigned char>& b,
               const string& fileDir) {
    if (fileDir.substr(fileDir.find_last_of(".") + 1) != "ppm") {
        cerr << "ERROR: " << fileDir << " is not a supported image type. Supported types: ppm\n";
        return false;
    }
    ifstream imageFile(fileDir, ios::in | ios::binary);
    string imageType;
    short int colorSize;
    imageFile >> imageType >> WIDTH >> HEIGHT >> colorSize;
    r.resize(WIDTH * HEIGHT);
    g.resize(WIDTH * HEIGHT);
    b.resize(WIDTH * HEIGHT);
    int i = 0;
    if (imageType.compare("P3") == 0) {
        int tr, tb, tg;
        while (i < r.size())
            if (imageFile >> tr >> tb >> tg) {  // checks if all 3 RGB values could be read
                r[i] = tr;
                b[i] = tb;
                g[i++] = tg;
            } else
                cerr << "ERROR: PPM image could not be read correctly.\n";
    } else if (imageType.compare("P6") == 0) {
        imageFile.get();  // skips the whitespace after the PPM header
        unsigned char buffer[3];
        while (i < r.size())
            if (imageFile.read((char*)buffer, 3)) {  // checks if all 3 RGB values could be read
                r[i] = buffer[0];
                g[i] = buffer[1];
                b[i++] = buffer[2];
            } else
                cerr << "ERROR: PPM image could not be read correctly.\n";
    }
    // this normalizes all the input PPMs, scaling it up to 8 bit color values
    if (++colorSize != 256)
        for (int i = 0; i < r.size(); i++) {
            r[i] = (r[i] + 1) * (256 / colorSize) - 1;
            g[i] = (g[i] + 1) * (256 / colorSize) - 1;
            b[i] = (b[i] + 1) * (256 / colorSize) - 1;
        }
    imageFile.close();
    return true;
}

bool generateImage(const vector<vector<unsigned char>*>& pixMap, string filename, bool compression,
                   string options = "") {
    ofstream imageFile(filename, ios::out | ios::binary);
    if (!imageFile) return false;

    imageFile << (compression ? "P6 " : "P3 ") << WIDTH << " " << HEIGHT << " 255\n";

    bool grayscale = options.empty();
    bool r = options.find("r") != string::npos;
    bool g = options.find("g") != string::npos;
    bool b = options.find("b") != string::npos;

    if (grayscale) {
        if (compression) {
            int j = 0;
            while (j < pixMap[0]->size()) {
                unsigned char buffer[3] = {(*pixMap[0])[j], (*pixMap[0])[j], (*pixMap[0])[j++]};
                imageFile.write((char*)buffer, 3);
            }
        } else {
            short int gray;
            for (int y = 0; y < HEIGHT; y++) {
                for (int x = 0; x < WIDTH; x++) {
                    gray = (short int)(*pixMap[0])[y * WIDTH + x];
                    imageFile << gray << " " << gray << " " << gray << " ";
                }
                imageFile << '\n';  // for formating purposes
            }
        }
    } else {
        short int i;
        unsigned char buffer[3];
        if (compression) {
            int j = -1;  // starts from -1 because loop iterates j in the conditional statement
            while (++j < pixMap[0]->size()) {
                i = 0;
                buffer[0] = r ? (*pixMap[i++])[j] : (unsigned char)0;
                buffer[1] = g ? (*pixMap[i++])[j] : (unsigned char)0;
                buffer[2] = b ? (*pixMap[i])[j] : (unsigned char)0;
                imageFile.write((char*)buffer, 3);
            }
        } else {
            for (int y = 0; y < HEIGHT; y++) {
                for (int x = 0; x < WIDTH; x++) {
                    i = 0;
                    buffer[0] = r ? (*pixMap[i++])[y * WIDTH + x] : (unsigned char)0;
                    buffer[1] = g ? (*pixMap[i++])[y * WIDTH + x] : (unsigned char)0;
                    buffer[2] = b ? (*pixMap[i])[y * WIDTH + x] : (unsigned char)0;
                    imageFile << (int)buffer[0] << " " << (int)buffer[1] << " " << (int)buffer[2]
                              << " ";
                }
                imageFile << '\n';  // for formating purposes
            }
        }
    }
    imageFile.close();
    return true;
}

void grayScaleImage(vector<unsigned char>& gray, const vector<unsigned char>& r,
                    const vector<unsigned char>& g, const vector<unsigned char>& b) {
    for (int i = 0; i < r.size(); i++) {
        unsigned char a = (unsigned char)((r[i] + g[i] + b[i]) / 3);  // averages the 3 rgb values
        gray[i] = a;
    }
}

vector<double> generateGaussianFilter(const int& n, const double& sigma) {
    vector<double> filter(n * n);  // filter is a square
    double d = 2 * sigma * sigma;
    int mid = n / 2;
    int min = mid - n + 1;
    double sum = 0;
    for (int y = min; y <= mid; y++)
        for (int x = min; x <= mid; x++) {
            filter[(y - min) * n + x - min] = exp(-(x * x + y * y) / d) / (M_PI * d);
            sum += filter[(y - min) * n + x - min];
        }
    for (double& d : filter) d /= sum;  // turns the values into percentages
    return filter;
}

void gaussianBlur(vector<vector<unsigned char>*>& channels, const int& n, const double& sigma) {
    vector<double> filter = generateGaussianFilter(n, sigma);
    vector<vector<unsigned char>> blur(channels.size());
    int mid = n / 2;
    int min = mid - n + 1;
    double sum, filterVal;
    vector<double> val(channels.size());
    int index;

    for (int y = 0; y < HEIGHT; y++)
        for (int x = 0; x < WIDTH; x++) {
            sum = 0;
            for (double& d : val) d = 0;         // reset val vector to 0
            for (int fy = min; fy <= mid; fy++)  // go through the filter
                for (int fx = min; fx <= mid; fx++)
                    // check if the index is inbounds
                    if (y + fy >= 0 && y + fy < HEIGHT && x + fx >= 0 && fx < WIDTH) {
                        filterVal = filter[(fy - min) * n + fx - min];
                        index = (y + fy) * WIDTH + x + fx;
                        sum += filterVal;
                        for (int i = 0; i < val.size(); i++)
                            val[i] += filterVal * (*channels[i])[index];
                    }
            for (int i = 0; i < val.size(); i++) {
                val[i] /= sum;  // scales all the values to 0-255
                blur[i].emplace_back(val[i]);
            }
        }
    for (int i = 0; i < channels[0]->size(); i++)  // copy blur to channel
        for (int j = 0; j < channels.size(); j++) (*channels[j])[i] = blur[j][i];
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

vector<double> applySobelOperator(const vector<vector<unsigned char>*>& channels,
                                  vector<unsigned char>& out, vector<double>& grad,
                                  vector<double>& angles) {
    short int sX, sY, gX, gY;
    vector<short int> gXc(channels.size()), gYc(channels.size());
    int index;
    double gMax = 0;
    for (int y = 1; y < HEIGHT - 1; y++)
        for (int x = 1; x < WIDTH - 1; x++) {
            gX = gY = 0;
            for (int i = 0; i < channels.size(); i++)
                gXc[i] = gYc[i] = 0;  // reset all the values to 0
            for (int fy = -1; fy <= 1; fy++)
                for (int fx = -1; fx <= 1; fx++) {
                    index = (y + fy) * WIDTH + x + fx;
                    sX = sobelKernalX[(fy + 1) * 3 + (fx + 1)];
                    sY = sobelKernalY[(fy + 1) * 3 + (fx + 1)];
                    for (short int i = 0; i < channels.size(); i++) {
                        gXc[i] += sX * (*channels[i])[index];
                        gYc[i] += sY * (*channels[i])[index];
                    }
                }
            for (short int i = 0; i < channels.size(); i++) {
                // find the value with highest intensity
                if (abs(gXc[i]) > abs(gX)) gX = gXc[i];
                if (abs(gYc[i]) > abs(gY)) gY = gYc[i];
            }
            index = y * WIDTH + x;
            grad[index] = sqrt(gX * gX + gY * gY);
            if (grad[index] > gMax) gMax = grad[index];
            angles[index] = atan2(gY, gX);
        }
    for (int i = 0; i < out.size(); i++)
        // copy the values to the out vector with a range of 0-255
        out[i] = (unsigned char)((grad[i] / gMax) * 255);
}

void applyNonMaxSuppression(vector<unsigned char>& image, const vector<double>& angles) {
    vector<unsigned char> edges(image.size());
    vector<float> dir(angles.size());
    // convert radian values to a float between 0-8
    for (int i = 0; i < angles.size(); i++) dir[i] = fmod(angles[i] + M_PI, M_PI) / M_PI * 8;

    int i, tc, tl, tr, cl, cr, bc, bl, br;
    for (int y = 1; y < HEIGHT - 1; y++)
        for (int x = 1; x < WIDTH - 1; x++) {
            i = y * WIDTH + x;
            tc = i - WIDTH;  // top-center
            tl = tc - 1;     // top-left
            tr = tc + 1;     // top-right
            cl = i - 1;      // center-left
            cr = i + 1;      // center-right
            bc = i + WIDTH;  // bottom-center
            bl = bc - 1;     // bottom-left
            br = bc + 1;     // bottom-right

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
    image = edges;
}

void calculateThresholdValues(vector<double>& grad, double& tMinRatio, double& tMaxRatio) {
    double sum = 0;
    for (const double& d : grad) sum += d;
    double avg = sum / grad.size();
    double var = 0;
    for (const double& d : grad) var += (d - avg) * (d - avg);
    var /= grad.size();
    double sd = sqrt(var);
    tMinRatio = (sd) / 255;
    tMaxRatio = (avg + sd / 8) / 255;
    cout << "Average: " << avg << "\nVariance: " << var << "\nStandard Deviation: " << sd
         << "\ntMin: " << tMinRatio << "\ntMax: " << tMaxRatio;
}

vector<int> threshold(vector<unsigned char>& image, const double& tMinRatio,
                      const double& tMaxRatio) {
    int med = COLOR_SIZE / 2;
    int tMax = COLOR_SIZE * tMaxRatio;
    int tMin = tMax * tMinRatio;

    vector<int> strong;

    for (int i = 0; i < image.size(); i++)
        // if value it larger than max threshold, its an edge
        if (image[i] > tMax) {
            image[i] = COLOR_SIZE;
            strong.emplace_back(i);
            // if it's inbetween it could be an edge
        } else if (image[i] > tMin)
            image[i] = med;
        else
            // if value it larger than min threshold, its not an edge
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
                // look through neighbors and see if its touching any pixels that could be an edge
                for (int n : neighbors) {
                    if (image[n] == 127 && edges[n] == 0) {
                        edges[n] = max;
                        // recurse with added edge
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
        if (!(tokens.empty() || regex_match(tokens[0], re))) {
            inputFileDir = tokens[0];
            return true;
        }
        return false;
    }

   private:
    vector<string> tokens;
    regex re;
};

bool detectEdges(const InputParser& input) {
    string inputFileDir = "";
    if (!input.getInputFile(inputFileDir)) {
        cerr << "ERROR: Input file not specified.\n";
        return false;
    }

    string outputFileDir = "output.ppm";
    string stMinRatio, stMaxRatio, sgSize, sgSigma, stage, channelOptions;
    if (!(input.getCmdOption("-b", {&sgSize, &sgSigma}) &&
          input.getCmdOption("-t", {&stMinRatio, &stMaxRatio}) &&
          input.getCmdOption("-o", {&outputFileDir}) &&
          input.getCmdOption("-c", {&channelOptions})))
        return false;
    bool compressed = input.cmdOptionExists("-compressed");

    vector<unsigned char> r, g, b;
    if (!loadImage(r, g, b, inputFileDir)) return false;
    vector<unsigned char> gray(r.size());

    vector<vector<unsigned char>*> channels;
    if (channelOptions.empty() || input.cmdOptionExists("-grayscale")) {
        grayScaleImage(gray, r, g, b);
        channels.emplace_back(&gray);
        if (input.cmdOptionExists("-grayscale"))
            return generateImage(channels, outputFileDir, compressed, channelOptions);
    } else {
        if (channelOptions.find("r") != string::npos) channels.emplace_back(&r);
        if (channelOptions.find("g") != string::npos) channels.emplace_back(&g);
        if (channelOptions.find("b") != string::npos) channels.emplace_back(&b);
    }

    int gSize = 5;
    double gSigma = 1.0;
    if (!sgSize.empty() && !sgSigma.empty()) {
        gSize = stoi(sgSize);
        gSigma = stod(sgSigma);
    }

    gaussianBlur(channels, gSize, gSigma);
    if (input.cmdOptionExists("-blur"))
        return generateImage(channels, outputFileDir, compressed, channelOptions);

    vector<double> grad(r.size()), angles(r.size());
    applySobelOperator(channels, gray, grad, angles);

    string gradientType;
    if (input.cmdOptionExists("-gradient") && input.getCmdOption("-gradient", {&gradientType})) {
        vector<vector<unsigned char>*> RGB{&r, &g, &b};
        if (gradientType.compare("hsv") == 0 || gradientType.empty())
            for (int i = 0; i < r.size(); i++)
                HSVtoRGB(angles[i], 1.0, gray[i] / 255.0, r[i], g[i], b[i]);
        else if (gradientType.compare("hsl") == 0)
            for (int i = 0; i < r.size(); i++)
                HSLtoRGB(angles[i], 1.0, gray[i] / 255.0, r[i], g[i], b[i]);
        return generateImage(RGB, outputFileDir, compressed, "rgb");
    }
    if (input.cmdOptionExists("-sobel"))
        return generateImage(vector<vector<unsigned char>*>{&gray}, outputFileDir, compressed);

    applyNonMaxSuppression(gray, angles);
    if (input.cmdOptionExists("-nonmax"))
        return generateImage(vector<vector<unsigned char>*>{&gray}, outputFileDir, compressed);

    double tMinRatio, tMaxRatio;
    if (!stMinRatio.empty() && !stMaxRatio.empty()) {
        tMinRatio = stod(stMinRatio);
        tMaxRatio = stod(stMaxRatio);
    } else
        calculateThresholdValues(grad, tMinRatio, tMaxRatio);

    vector<int> strong = threshold(gray, tMinRatio, tMaxRatio);
    if (!input.cmdOptionExists("-threshold")) defineEdges(gray, strong, 255);
    return generateImage(vector<vector<unsigned char>*>{&gray}, outputFileDir, compressed);
}

int main(int argc, char** argv) {
    InputParser input(argc, argv);
    if (input.cmdOptionExists("-h")) {
        cout << helpMessage << endl;
        return 0;
    }
    return detectEdges(input);
}