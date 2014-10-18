#include <iostream>
#include <vector>
#include <algorithm>
#include "matrix.h"
#include "io.h"

using std::cout;
using std::endl;

using std::tuple;
using std::vector;
using std::sort;
using std::get;
using std::tie;
using std::make_tuple;

// Matrix usage example
// Also see: matrix.h, matrix.hpp for comments on how filtering works

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

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        cout << "Usage: " << argv[0] << " <input_image> <output_image>" << endl;
        return 0;
    }
    // Image = Matrix<tuple<uint, uint, uint>>
    // tuple is from c++ 11 standard
    Image img = load_image(argv[1]);
    Image img2 = img.unary_map(MedianFilter());

    save_image(img2, argv[2]);
}
