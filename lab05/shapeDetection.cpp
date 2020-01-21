#define _USE_MATH_DEFINES

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <vector>

using namespace std;

int WIDTH, HEIGHT, COLOR_SIZE = 255;

string helpMessage = "this message is helpful!\n";

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
}

vector<int> applyThreshold(vector<unsigned char>& image, const double& tMinRatio,
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

struct Edge {
    int x, y;
    Edge(const int& x, const int& y) {
        this->x = x;
        this->y = y;
    }
};

void findAllEdges(const vector<unsigned char>& image, vector<Edge>& edges) {
    for (int y = 0; y < HEIGHT; y++)
        for (int x = 0; x < WIDTH; x++)
            if (image[y * WIDTH + x] == 255) edges.emplace_back(Edge(x, y));
}

inline double lineDistance(const Edge& a, const Edge& b, const Edge& c) {
    return ((b.y - a.y) * c.x - (b.x - a.x) * c.y + b.x * a.y - b.y * a.x) /
           sqrt((b.y - a.y) * (b.y - a.y) + (b.x - a.x) * (b.x - a.x));
}

void drawPixel(vector<unsigned char>& r, vector<unsigned char>& g, vector<unsigned char>& b, int x,
               int y, unsigned char rgb[3]) {
    if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) return;
    r[y * WIDTH + x] = rgb[0];
    g[y * WIDTH + x] = rgb[1];
    b[y * WIDTH + x] = rgb[2];
}

void drawCircle(vector<unsigned char>& r, vector<unsigned char>& g, vector<unsigned char>& b,
                int cx, int cy, int rad, unsigned char rgb[3]) {
    int x = rad;
    int y = 0;
    int err = 0;
    while (x >= y) {
        drawPixel(r, g, b, cx + x, cy + y, rgb);
        drawPixel(r, g, b, cx + y, cy + x, rgb);
        drawPixel(r, g, b, cx - y, cy + x, rgb);
        drawPixel(r, g, b, cx - x, cy + y, rgb);
        drawPixel(r, g, b, cx - x, cy - y, rgb);
        drawPixel(r, g, b, cx - y, cy - x, rgb);
        drawPixel(r, g, b, cx + y, cy - x, rgb);
        drawPixel(r, g, b, cx + x, cy - y, rgb);
        if (err <= 0) {
            y++;
            err += 2 * y + 1;
        } else {
            x--;
            err -= 2 * x + 1;
        }
    }
}

void accumPixel(vector<int>& accum, int x, int y) {
    if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) return;
    accum[y * WIDTH + x]++;
}

void accumCircle(vector<int>& accum, int cx, int cy, int rad) {
    int x = rad;
    int y = 0;
    int err = 0;
    while (x >= y) {
        accumPixel(accum, cx + x, cy + y);
        accumPixel(accum, cx + y, cy + x);
        accumPixel(accum, cx - y, cy + x);
        accumPixel(accum, cx - x, cy + y);
        accumPixel(accum, cx - x, cy - y);
        accumPixel(accum, cx - y, cy - x);
        accumPixel(accum, cx + y, cy - x);
        accumPixel(accum, cx + x, cy - y);
        if (err <= 0) {
            y++;
            err += 2 * y + 1;
        } else {
            x--;
            err -= 2 * x + 1;
        }
    }
}

bool isInCicles(const map<int, vector<pair<int, int>>>& allradius, int x, int y, int r) {
    for (const pair<int, vector<pair<int, int>>>& cc : allradius) {
        for (const pair<int, int>& p : cc.second) {
            int distSq = (p.first - x) * (p.first - x) + (p.second - y) * (p.second - y);
            if (distSq <= cc.first * cc.first / 2) {
                return true;
            }
        }
    }
    return false;
}

void houghTransform(vector<unsigned char>& r, vector<unsigned char>& g, vector<unsigned char>& b,
                    const vector<Edge>& edges) {
    int tHeight = sqrt(WIDTH * WIDTH + HEIGHT * HEIGHT) / 2;
    int tWidth = 360;
    int rMin = min(WIDTH, HEIGHT) / 40;
    int rMax = min(WIDTH, HEIGHT) / 16;
    cout << rMax << endl;
    int rDepth = rMax - rMin;

    vector<unsigned char> out(r.size(), 0);
    map<int, vector<pair<int, int>>> circleCollector;

    for (int rad = rMax; rad >= rMin; --rad) {
        vector<int> accum(WIDTH * HEIGHT, 0);
        // cout << rad << '\n';
        for (Edge edge : edges) {
            accumCircle(accum, edge.x, edge.y, rad);
        }
        // unsigned int i, t, l, r, b;
        unsigned int threshold = 10.0 / 3.0 * rad;
        for (int y = 1; y < HEIGHT - 1; y++)
            for (int x = 1; x < WIDTH - 1; x++) {
                if (accum[y * WIDTH + x] > threshold) {
                    int index = y * WIDTH + x;
                    // Check if it is a local max
                    if (accum[index] > accum[index + 1] &&      //
                        accum[index] > accum[index - 1] &&      //
                        accum[index] > accum[index + WIDTH] &&  //
                        accum[index] > accum[index - WIDTH])
                        if (!isInCicles(circleCollector, x, y, rad)) {
                            circleCollector[rad].emplace_back(make_pair(x, y));
                        }
                }
            }
    }

    int sum = 0, numCircles = 0;
    for (const pair<int, vector<pair<int, int>>>& cc : circleCollector) {
        sum += cc.first * cc.second.size();
        numCircles += cc.second.size();
    }
    double avg = (double)sum / numCircles;
    double var = 0;
    for (const pair<int, vector<pair<int, int>>>& cc : circleCollector)
        var += (cc.first - avg) * (cc.first - avg) * cc.second.size();
    var /= numCircles;
    double sd = sqrt(var);
    double lower = avg - sd * 2;
    double upper = avg + sd * 2;

    cout << "sum:\t\t" << sum << "\nnum:\t\t" << numCircles << "\naverage:\t" << avg << "\nsd:\t\t"
         << sd << "\nlower:\t\t" << lower << "\nupper:\t\t" << upper << endl;

    while (circleCollector.begin()->first < lower) {
        circleCollector.erase(circleCollector.begin());
    }

    for (const pair<int, vector<pair<int, int>>>& cc : circleCollector) {
        cout << cc.first << ":";
        if (cc.first < lower || cc.first > upper) continue;
        for (const pair<int, int> p : cc.second) {
            cout << " <" << p.first << ", " << p.second << ">";
        }
        cout << endl;
    }

    // diameter of the coins in inches, smallest to largest
    static float sizes[] = {.705f, .750f, .835f, .955f, 1.205f};
    // corresponding values
    static float values[] = {.10f, .01f, .05f, .25f, .50f};

    int baseSize = circleCollector.begin()->first;
    int startingCenterIndex = circleCollector.begin()->second[0].second * WIDTH +
                              circleCollector.begin()->second[0].first;  // Y * WIDTH + X
    int startingCoin = 0;

    #define RED_RATIO 1.1
    // if red is the most prominent color of the starting coin
    if ((float) r[startingCenterIndex] / g[startingCenterIndex] > RED_RATIO && (float) r[startingCenterIndex] / b[startingCenterIndex] > RED_RATIO )
        startingCoin++;
    int coin = 1;  // change later to starting coin
    vector<int> sizeValue;

    sizeValue.push_back(baseSize * sizes[coin + 1] / sizes[coin++]);
    while (coin < 4) {
        sizeValue.push_back(sizeValue.back() * sizes[coin + 1] / sizes[coin]);
        coin++;
    }

    float total = 0;

    int count[5] = {0, 0, 0, 0, 0};
    int index = 0;
    int indexBeforeNickelCutoff = 2 - startingCoin;
    unsigned char colors[5][3] = {
        {255, 255, 255}, {255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}};
    for (auto const& cc : circleCollector) {
        cout << sizeValue[index] << " " << values[index] << endl;
        while (cc.first > sizeValue[index]) index++;
        if (index <= indexBeforeNickelCutoff) {
            for (auto const& p : cc.second) {
                int position = p.second * WIDTH + p.first;
                cout << (int)r[position] << " " << (int)g[position] << " " << (int)b[position]
                     << endl;

                if ((float) r[position] / g[position] > RED_RATIO && (float) r[position] / b[position] > RED_RATIO ) {
                    cout << "PENNY!\n";
                    // its a penny is r is most prominent color
                    total += 0.01;
                    drawCircle(r, g, b, p.first, p.second, cc.first, colors[1]);
                    count[1]++;
                } else {
                    total += values[index + startingCoin];
                    drawCircle(r, g, b, p.first, p.second, cc.first, colors[index + startingCoin]);
                    count[index + startingCoin]++;
                }
            }

        } else {
            count[index + startingCoin] += cc.second.size();
            total += cc.second.size() * values[index + startingCoin];
            for (auto const& p : cc.second)
                drawCircle(r, g, b, p.first, p.second, cc.first, colors[index + startingCoin]);
        }
    }

    cout << "PENNIES:\t" << count[1] << endl;
    cout << "DIMES:\t\t" << count[0] << endl;
    cout << "NICKELS:\t" << count[2] << endl;
    cout << "QUARTERS:\t" << count[3] << endl;
    cout << "HALF DOLLARS:\t" << count[4] << endl;
    cout << "Total: " << total << endl;
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

bool detectCoins(const InputParser& input) {
    string inputFileDir = "";
    if (!input.getInputFile(inputFileDir)) {
        cerr << "ERROR: Input file not specified.\n";
        return false;
    }

    string outputFileDir = "output.ppm";
    string temp;
    if (!(input.getCmdOption("-o", {&outputFileDir}))) return false;
    bool compressed = input.cmdOptionExists("-compressed");

    /**********************************************************
     * * * * * * * * * * * CANNY EDGE * * * * * * * * * * * * *
     **********************************************************/

    vector<unsigned char> r, g, b;
    if (!loadImage(r, g, b, inputFileDir)) return false;
    vector<unsigned char> gray(r.size());

    vector<vector<unsigned char>*> channels;
    grayScaleImage(gray, r, g, b);
    channels.emplace_back(&gray);

    gaussianBlur(channels, 5, 1.0);

    vector<double> grad(r.size()), angles(r.size());
    applySobelOperator(channels, gray, grad, angles);

    applyNonMaxSuppression(gray, angles);

    double tMinRatio, tMaxRatio;
    calculateThresholdValues(grad, tMinRatio, tMaxRatio);

    vector<int> strong = applyThreshold(gray, tMinRatio, tMaxRatio);
    defineEdges(gray, strong, 255);

    /**********************************************************/

    vector<Edge> edges;

    findAllEdges(gray, edges);
    cout << gray.size() << endl;
    cout << edges.size() << endl;
    houghTransform(r, g, b, edges);

    return generateImage(vector<vector<unsigned char>*>{&r, &g, &b}, outputFileDir, compressed,
                         "rgb");
}

int main(int argc, char** argv) {
    InputParser input(argc, argv);
    if (input.cmdOptionExists("-h")) {
        cout << helpMessage << endl;
        return 0;
    }
    return detectCoins(input);
}