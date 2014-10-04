#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <vector>
#include <tuple>
#include <queue>
#include <cmath>
#include <algorithm>

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::tuple;
using std::tie;
using std::make_tuple;
using std::queue;

#include "io.h"
#include "matrix.h"

typedef tuple<uint, uint, uint, uint> Rect;

class MedianFilter
{
public:
    tuple<uint, uint, uint> operator () (const Image &m) const
    {
        uint size = 2 * radius + 1;
        vector <uint> r(size * size), g(size * size), b(size * size);
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r[size * i + j], g[size * i + j], b[size * i + j]) = m(i, j);
            }
        }
        sort(r.begin(), r.end());
        sort(g.begin(), g.end());
        sort(b.begin(), b.end());

        return make_tuple(r[size*size / 2], g[size*size / 2], b[size*size / 2]);
    }
    // Radius of neighbourhoud, which is passed to that operator
    static const int radius = 2;
};


void autolevels(Image &img)
{
    uint r, g, b;
    uint max_r(0), max_g(0), max_b(0);
    uint min_r(255), min_g(255), min_b(255);

    for (uint i = 0; i < img.n_rows; i++)
        for (uint j = 0; j < img.n_cols; j++) {
            tie(r, g, b) = img(i, j);
            if (r > max_r)
                max_r = r;
            if (g > max_g)
                max_g = g;
            if (b > max_b)
                max_b = b;
            if (r < min_r)
                min_r = r;
            if (g < min_g)
                min_g = g;
            if (b < min_b)
                min_b = b;
        }

    for (uint i = 0; i < img.n_rows; i++)
        for (uint j = 0; j < img.n_cols; j++) {
            tie(r, g, b) = img(i, j);
            r = 255 * (r - min_r) / (max_r - min_r);
            g = 255 * (g - min_g) / (max_g - min_g);
            b = 255 * (b - min_b) / (max_b - min_b);
            img(i, j) = make_tuple(r, g, b);
        }
}


// Finds connected components
uint bfs(vector <vector <uint>> &used, const Image &binimg)
{
    uint m, n;
    uint k = 0;
    queue <tuple <uint, uint>> que;

    for (uint i = 0; i < used.size(); i++)
        for (uint j = 0; j < used[i].size(); j++) {
            if (!used[i][j]) {
                que.push(make_tuple(i, j));
                used[i][j] = ++k;
                while (!que.empty()) {
                    tie(m, n) = que.front();
                    que.pop();
                    if (m >= 1 && !used[m-1][n] && binimg(m,n) == binimg(m-1, n)) {
                        used[m-1][n] = k;
                        que.push(make_tuple(m-1, n));
                    }
                    if (m+1 < used.size() && !used[m+1][n] && binimg(m,n) == binimg(m+1, n)) {
                        used[m+1][n] = k;
                        que.push(make_tuple(m+1, n));
                    }
                    if (n >= 1 && !used[m][n-1] && binimg(m,n) == binimg(m, n-1)) {
                        used[m][n-1] = k;
                        que.push(make_tuple(m, n-1));
                    }
                    if (n+1 < used[i].size() && !used[m][n+1] && binimg(m,n) == binimg(m, n+1)) {
                        used[m][n+1] = k;
                        que.push(make_tuple(m, n+1));
                    }
                }
            }
        }
    return k;
}


// Returns threshold for grayscale img
uint otsuThreshold(vector <uint> &image)
{
    uint min = image[0], max = image[0];
    uint temp;

    /**** Histogram construction ****/
    /* Find max and min grayscale */
    for(uint i = 0; i < image.size(); i++) {
        temp = image[i];
        if (temp < min)
            min = temp;
        if (temp > max)
            max = temp;
    }

    vector <uint> hist(max - min + 1);
    for (uint i = 0; i < hist.size(); i++)
        hist[i] = 0;
    
    for(uint i = 0; i < image.size(); i++)
        hist[image[i] - min]++;

    // for(uint i = 0; i < image.size(); i++)
        // hist[image[i] - min] /= image.size();
    /**** Histogram has constructed ****/

    uint temp1 = temp = 0;
    uint alpha(0), beta(0), threshold(0);
    double sigma, maxSigma = -1;
    
    /* Counting expectation */
    for(uint i = 0; i <= (max - min); i++) {
        temp += i * hist[i];
        temp1 += hist[i];
    }

    /* Main cycle for finding threshold
    find grayscale with minimal dispersion */
    for(uint i = 0; i < (max - min); i++) {
        alpha += i * hist[i];
        beta += hist[i];
        double w1 = static_cast<double>(beta) / static_cast<double>(temp1);
        double a = static_cast<double>(alpha) / static_cast<double>(beta) - 
                   static_cast<double>((temp - alpha)) / static_cast<double>((temp1 - beta));

        sigma = w1 * (1 - w1) * a * a;

        if (sigma > maxSigma) {
            maxSigma = sigma;
            threshold = i;
        }
    }
    return threshold + min;
}


uint BHThreshold(vector <uint> &image)
{
    uint min = image[0], max = image[0];
    uint temp;

    /**** Histogram construction ****/
    /* Find max and min grayscale */
    for(uint i = 0; i < image.size(); i++) {
        temp = image[i];
        if (temp < min)
            min = temp;
        if (temp > max)
            max = temp;
    }

    vector <uint> hist(max - min + 1);
    for (uint i = 0; i < hist.size(); i++)
        hist[i] = 0;
    
    for(uint i = 0; i < image.size(); i++)
        hist[image[i] - min]++;

    /**** Histogram has constructed ****/

    uint i_s(0), i_e(hist.size()-1);
    uint i_m = static_cast<uint>((i_s + i_e) / 2.0f);

    uint w_l(0), w_r(0);
    for (uint i = i_s; i < i_m + 1; i++)
        w_l += hist[i];

    for (uint i = i_m + 1; i < i_e + 1; i++)
        w_r += hist[i];

    while (i_s <= i_e) {
        if (w_r > w_l) { // правая часть тяжелее
            w_r -= hist[i_e--];
            if (((i_s + i_e) / 2) < i_m) {
                w_r += hist[i_m];
                w_l -= hist[i_m--];
            }
        } else { // левая часть тяжелее
            w_l -= hist[i_s++]; 
            if (((i_s + i_e) / 2) > i_m) {
                w_l += hist[i_m + 1];
                w_r -= hist[i_m + 1];
                i_m++;
            }
        }
    }
    return i_m;
}


// Makes img binary by luminance
void make_binarization(Image &img)
{
    uint r, g, b;
    vector <uint> gs_img;
    for (uint i = 0; i < img.n_rows; i++)
        for (uint j = 0; j < img.n_cols; j++) {
            tie(r, g, b) = img(i ,j);
            gs_img.push_back(0.299 * r + 0.587 * g + 0.114 * b);
        }

    double threshold = 33;
    // cout << BHThreshold(gs_img);
    // cout << otsuThreshold(gs_img);
    for (uint i = 0; i < img.n_rows; i++)
        for (uint j = 0; j < img.n_cols; j++) {
            tie(r, g, b) = img(i, j);
            uint lumin = 255;
            // RGB Luminance value = 0.299 R + 0.587 G + 0.114 B
            if (0.299 * r + 0.587 * g + 0.114 * b < threshold)
                lumin = 0;
            img(i, j) = make_tuple(lumin, lumin, lumin);
        }
}


vector <double> moment(const vector <vector <uint>> &used, const int m, const int n,
            const vector <uint> &avg_x, const vector <uint> &avg_y)
{
    vector <double> v(avg_x.size());

    for (long long i = 0; i < static_cast<long long>(used.size()); i++)
        for (long long j = 0; j < static_cast<long long>(used[0].size()); j++) {
            v[used[i][j]-1] += pow(j-avg_x[used[i][j]-1], m) * pow(i-avg_y[used[i][j]-1], n);
        }
    return v;
}


void count_geometrical_characteristics(const vector <vector <uint>> &used, vector <vector <tuple<uint, uint>>> &border,
                                       vector <uint> &area, vector <uint> &avg_x,
                                       vector <uint> &avg_y, vector <uint> &perim,
                                       vector <double> &elongation, vector <double> &theta)
{
    // Counting area, perimeter and mass center
    for (uint i = 0; i < used.size(); i++)
        for (uint j = 0; j < used[0].size(); j++) {
            area[used[i][j]-1]++;
            avg_x[used[i][j]-1] += j;
            avg_y[used[i][j]-1] += i;
            if (used[i][j] > 1 &&
                (used[i >= 1 ? i-1 : i][j]-1) * (used[i+1 < used.size() ? i+1 : i][j]-1) *
                (used[i][j >= 1 ? j-1 : j]-1) * (used[i][j+1 < used[0].size() ? j+1 : j]-1) == 0) {
                perim[used[i][j]-1]++;
                border[used[i][j]-1].push_back(make_tuple(i, j));
            }
        }

    for (uint i = 0; i < avg_x.size(); i++) {
        avg_x[i] /= area[i];
        avg_y[i] /= area[i];
    }

    // counting moments and elongations
    vector <double> moment20 = moment(used, 2, 0, avg_x, avg_y);
    vector <double> moment02 = moment(used, 0, 2, avg_x, avg_y);
    vector <double> moment11 = moment(used, 1, 1, avg_x, avg_y);

    for (uint i = 0; i < avg_x.size(); i++) {
        elongation[i] = (moment20[i] + moment02[i] + sqrt(pow((moment20[i] - moment02[i]), 2) + 4*pow(moment11[i], 2))) /
                        (moment20[i] + moment02[i] - sqrt(pow((moment20[i] - moment02[i]), 2) + 4*pow(moment11[i], 2)));
    }

    // counting angles
    for (uint i = 0; i < avg_x.size(); i++)
        theta[i] = atan(2 * moment11[i] / (moment20[i] - moment02[i])) / 2;
}


// finds tip of ptr and Rect in which elements are inscribed
void find_dirs(const Image &in, const Image &img, const vector <vector <tuple<uint, uint>>> &border,
               vector <uint> &avg_x, vector <uint> &avg_y, vector <double> &theta,
               vector <uint> &dir_x, vector <uint> &dir_y, vector <Rect> &rects)
{
    vector <tuple <uint, uint>> ok;
    long long a, b;    
    double EPS = 0.025;
    for (uint i = 1; i < border.size(); i++) {
        uint max_xe = 0;
        uint min_xe = img.n_cols;
        uint max_ye = 0;
        uint min_ye = img.n_rows;
        // find direction dots
        for (uint j = 0; j < border[i].size(); j++) {
            tie(a, b) = border[i][j];

            if (a > max_ye)
                max_ye = a;
            if (a < min_ye)
                min_ye = a;
            if (b > max_xe)
                max_xe = b;
            if (b < min_xe)
                min_xe = b;

            double at;
            if (b == avg_x[i])
                at = M_PI/2;
            else
                at = atan(fabs((a - static_cast<long long>(avg_y[i]))) / fabs((b - static_cast<long long>(avg_x[i]))));
            
            if (at > M_PI/4)
                at -= M_PI/2;
            if (at < -M_PI/4)
                at += M_PI/2;

            if (fabs(at - theta[i]) < EPS || fabs(at + theta[i]) < EPS) {
                ok.push_back(make_tuple(a, b));
            }
        }

        uint r, g, bl;
        long long gr_x(-1), gr_y(-1);
        rects.push_back(make_tuple(min_xe, min_ye, max_xe, max_ye));
        for (uint ii = min_ye; ii <= max_ye; ii++)
            for (uint j = min_xe; j <= max_xe; j++) {
                tie(r, g, bl) = in(ii, j);
                if (r < 150 && bl < 150 && g > 180) {
                    gr_x = j;
                    gr_y = ii;
                }
            }
        
        // ...cont. of finding direction
        tie(a, b) = ok[0];
        long long max_x = b, max_y = a;
        for (uint k = 0; k < ok.size(); k++) {
            tie(a, b) = ok[k];
            long long dgr_x = gr_x - avg_x[i];
            long long dgr_y = gr_y - avg_y[i];
            double dx = b - avg_x[i];
            double dy = a - avg_y[i];
            double m_dx = max_x - avg_x[i];
            double m_dy = max_y - avg_y[i];
            if (dx * dx + dy * dy > m_dy * m_dy + m_dx * m_dx) {
                if (gr_x > 0) {
                    if ((dx > 0 && dgr_x > 0) || (dy > 0 && dgr_y > 0) ||
                    (dx < 0 && dgr_x < 0) || (dy < 0 && dgr_y < 0)) {
                        max_x = b;
                        max_y = a;
                    }
                } else {
                    max_x = b;
                    max_y = a;
                }
            }
        }
        dir_x[i] = max_x;
        dir_y[i] = max_y;
        ok.clear();
    }
}


// draws lines of path and rects of elements on path
void find_way(Image &img, const vector <vector <uint>> &used,
              vector <uint> &avg_x, vector <uint> &avg_y, uint i,
              vector <uint> &dir_x, vector <uint> &dir_y, vector <uint> &perim,
              vector <uint> &area, vector <double> &elongation, vector <Rect> &rects, vector <Rect> &path)
{
    if (perim[i]*perim[i] / area[i] > 17 || perim[i]*perim[i] / area[i] < 13 ||
        elongation[i] < 3.1 || elongation[i] > 4.25) {
        path.push_back(rects[i-1]);
        return;
    }
    path.push_back(rects[i-1]);
    long long l, dx, dy;
    long long xr = abs(dir_x[i] - avg_x[i]);
    long long yr = abs(dir_y[i] - avg_y[i]);
    if (xr > yr)
        l = xr;
    else
        l = yr;
    long long px = (avg_x[i] << 12) + (1 << 11);
    long long py = (avg_y[i] << 12) + (1 << 11);
    long long ex = (dir_x[i] << 12) + (1 << 11);
    long long ey = (dir_y[i] << 12) + (1 << 11);
    if (l != 0) {
        dx = (ex - px) / l;
        dy = (ey - py) / l;
    } else {
        dx = 0;
        dy = 0;
    }

    while ((used[py >> 12][px >> 12] == i+1 || used[py >> 12][px >> 12] == 1)) {
        img(py >> 12, px >> 12) = make_tuple(255, 0, 255);
        px += dx;
        py += dy;
        if ((py >> 12) < 0 || (px >> 12) < 0 ||
           (py >> 12) >= static_cast<long long>(img.n_rows) || (px >> 12) >= static_cast<long long>(img.n_cols)) {
            cout << endl << i << endl;
            return;            
        }
    }

    find_way(img, used, avg_x, avg_y, used[py >> 12][px >> 12] - 1, dir_x, dir_y, perim, area, elongation, rects, path);
}


// finds and retruns red ptr label
uint find_red_ptr(const Image &img, const vector <vector <uint>> &used)
{
    uint r, g, b, pix = 0;
    for (uint i = 1; i < img.n_rows-1; i++)
        for (uint j = 1; j < img.n_cols-1; j++) {
            pix = 0;
            for (uint k = i-1; k <= i+1; k++)
                for (uint l = j-1; l <= j+1; l++) {
                    tie(r, g, b) = img(k, l);
                    if (r > 160 && g < 80 && b < 80)
                        pix++;
                }
            if (pix >= 7)
                return used[i][j];
        }
    throw "No red pointer";
    return 0;
}


void draw_rects(Image &img, vector <Rect> &path)
{
    uint x, y, w, h;
    for (uint i = 0; i < path.size(); i++) {
        tie(x, y, w, h) = path[i];
        uint dx = x;
        uint dy = y;
        while (dx <= w) {
            img(y, dx) = make_tuple(255, 0, 0);
            img(h, dx) = make_tuple(255, 0, 0);
            dx++;
        }
        while (dy <= h) {
            img(dy, x) = make_tuple(255, 0, 0);
            img(dy, w) = make_tuple(255, 0, 0);
            dy++;
        }
    }
}


tuple<vector<Rect>, Image>
find_treasure(const Image& in)
{
    auto path = vector<Rect>();
        
    Image filtered_img = in.deep_copy();

    /* filteering step */
    autolevels(filtered_img);
    filtered_img = filtered_img.unary_map(MedianFilter());
    // return make_tuple(path, filtered_img);

    /* binarization step */
    Image binimg = filtered_img.deep_copy();
    make_binarization(binimg);

    /* binary img filtering step */
    binimg = binimg.unary_map(MedianFilter());
    // return make_tuple(path, binimg);

    vector <vector <uint>> used(in.n_rows);
    for (auto &rows : used)
        rows.resize(in.n_cols);

    uint k = bfs(used, binimg);

    vector <uint> area(k), avg_x(k), avg_y(k), perim(k);
    vector <double> elongation(k);
    vector <double> theta(k);
    vector <vector <tuple <uint, uint>>> border(k);
    count_geometrical_characteristics(used, border, area, avg_x, avg_y, perim, elongation, theta);
    
    /*
    for (uint i = 0; i < used.size(); i++)
        for (uint j = 0; j < used[0].size(); j++)
            if (area[used[i][j]-1] < 350 || area[used[i][j]-1] > 4800)
                used[i][j] = 1;
    */

    /*
    for (uint i = 0; i < k; i++) {
        cout << i+1 << ": " << " \tPerim: " << perim[i] << " \tArea: " << area[i] <<
        " \tCompact: " << perim[i]*perim[i] / area[i] << " \tElongation: " << elongation[i] <<
        " \tAngles: " << theta[i] << " \tAvg point (" << avg_x[i] << ", " << avg_y[i] << ")" << endl;
    }
    */

    Image img = in.deep_copy();
    /*
    for (uint i = 0; i < in.n_rows; i++)
        for (uint j = 0; j < in.n_cols; j++) {
            img(i, j) = make_tuple(255-5*used[i][j], 255-15*used[i][j], 255-30*used[i][j]);
        }
    */
    vector <uint> dir_x(k), dir_y(k);
    auto rects = vector<Rect>();

    find_dirs(filtered_img, img, border, avg_x, avg_y, theta, dir_x, dir_y, rects);
    find_way(img, used, avg_x, avg_y, find_red_ptr(filtered_img, used) - 1, dir_x, dir_y, perim, area, elongation, rects, path);

    vector <Rect> treasure;
    treasure.push_back(path[path.size() - 1]);
    draw_rects(img, treasure);

    return make_tuple(path, img);
}


int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cout << "Usage: " << endl << argv[0]
             << " <in_image.bmp> <out_image.bmp> <out_path.txt>" << endl;
        return 0;
    }

    try {
        Image src_image = load_image(argv[1]);
        ofstream fout(argv[3]);

        vector<Rect> path;
        Image dst_image;
        tie(path, dst_image) = find_treasure(src_image);
        save_image(dst_image, argv[2]);

        uint x, y, width, height;
        for (const auto &obj : path)
        {
            tie(x, y, width, height) = obj;
            fout << x << " " << y << " " << width << " " << height << endl;
        }

    } catch (const string &s) {
        cerr << "Error: " << s << endl;
        return 1;
    }
}
