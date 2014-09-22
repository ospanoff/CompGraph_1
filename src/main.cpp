#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <vector>
#include <tuple>
#include <queue>

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

tuple<vector<Rect>, Image>
find_treasure(const Image& in)
{
    // Base: return Rect of treasure only
    // Bonus: return Rects of arrows and then treasure Rect

    auto path = vector<Rect>();

    // RGB Luminance value = 0.3 R + 0.59 G + 0.11 B
    Image binimg = in.deep_copy();

    uint r, g, b;
    for (uint i = 0; i < in.n_rows; i++)
        for (uint j = 0; j < in.n_cols; j++) {
            tie(r, g, b) = binimg(i, j);
            uint lumin = 0;
            if (1 - (0.299 * r + 0.587 * g + 0.114 * b) / 255 < 0.83)
                lumin = 255;
            binimg(i, j) = make_tuple(lumin, lumin, lumin);
        }

    queue <tuple <uint, uint>> que;
    vector <vector <uint>> used(in.n_rows);
    for (auto &rows : used)
        rows.resize(in.n_cols);

    uint m, n;
    uint k = 0;
    vector<uint> numofel;
    numofel.push_back(0);
    for (uint i = 0; i < used.size(); i++)
        for (uint j = 0; j < used[i].size(); j++) {
            if (!used[i][j]) {
                // cout << k++ << endl;
                que.push(make_tuple(i, j));
                used[i][j] = ++k;
                numofel.push_back(0);
                while (!que.empty()) {
                    tie(m, n) = que.front();
                    que.pop();
                    if (m >= 1 && !used[m-1][n] && binimg(m,n) == binimg(m-1, n)) {
                        used[m-1][n] = k;
                        que.push(make_tuple(m-1, n));
                        numofel[k]++;
                    }
                    if (m+1 < used.size() && !used[m+1][n] && binimg(m,n) == binimg(m+1, n)) {
                        used[m+1][n] = k;
                        que.push(make_tuple(m+1, n));
                        numofel[k]++;
                    }
                    if (n >= 1 && !used[m][n-1] && binimg(m,n) == binimg(m, n-1)) {
                        used[m][n-1] = k;
                        que.push(make_tuple(m, n-1));
                        numofel[k]++;
                    }
                    if (n+1 < used[i].size() && !used[m][n+1] && binimg(m,n) == binimg(m, n+1)) {
                        used[m][n+1] = k;
                        que.push(make_tuple(m, n+1));
                        numofel[k]++;
                    }
                }
            }
        }

    for (auto &a : numofel)
        cout << a << endl;

    for (uint i = 0; i < in.n_rows; i++)
        for (uint j = 0; j < in.n_cols; j++) {
            if (numofel[used[i][j]] < 100)
                used[i][j] = 1;
        }

    return make_tuple(path, binimg/*in.deep_copy()*/);
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
