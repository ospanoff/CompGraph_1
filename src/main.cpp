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
    
    ofstream fout("./out/2.txt");
    auto path = vector<Rect>();

    // RGB Luminance value = 0.3 R + 0.59 G + 0.11 B

    queue <tuple <uint, uint>> que;
    vector <vector <uint>> used(in.n_rows);
    for (auto &rows : used)
        rows.resize(in.n_cols);

    uint m, n;
    uint k = 0;
    for (uint i = 0; i < used.size(); i++)
        for (uint j = 0; j < used[i].size(); j++) {
            if (!used[i][j]) {
                cout << k++;
                que.push(make_tuple(i, j));
                used[i][j] = k;
                while (!que.empty()) {
                    tie(m, n) = que.front();
                    que.pop();
                    if (m >= 1 && !used[m-1][n] && in(m,n) == in(m-1, n)) {
                        used[m-1][n] = k;
                        que.push(make_tuple(m-1, n));
                    }
                    if (m+1 < used.size() && !used[m+1][n] && in(m,n) == in(m+1, n)) {
                        used[m+1][n] = k;
                        que.push(make_tuple(m+1, n));
                    }
                    if (n >= 1 && !used[m][n-1] && in(m,n) == in(m, n-1)) {
                        used[m][n-1] = k;
                        que.push(make_tuple(m, n-1));
                    }
                    if (n+1 < used[i].size() && !used[m][n+1] && in(m,n) == in(m, n+1)) {
                        used[m][n+1] = k;
                        que.push(make_tuple(m, n+1));
                    }
                }
            }
        }

    /*for (auto &rows : used) {
        for (auto &x : rows)
            fout << x;
        fout << endl;
    }*/

    Image img = in.deep_copy();
    for (uint i = 0; i < img.n_rows; i++)
        for (uint j = 0; j < img.n_cols; j++) {
            int r = used[i][j];
            img(i, j) = make_tuple(r, r, r);
        }

    return make_tuple(path, img/*in.deep_copy()*/);
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
