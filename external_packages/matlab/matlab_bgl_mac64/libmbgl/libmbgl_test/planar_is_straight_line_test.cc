#include <include/matlab_bgl.h>
#include <algorithm>
#include <iostream>

int main(int argc, char **argv) {
    mbglIndex n = 4;
    mbglIndex rp[] = {0,3,6,9,12};
    mbglIndex ci[] = {1,2,3,0,2,3,0,1,3,0,1,2};
	double X[] = { 1,2,1,0,2,1,0,1 };
	int is_sldraw = 0;
    is_straight_line_drawing(n, ci, rp, X, &is_sldraw);
    std::cout << "is_sldraw = " << is_sldraw << std::endl;
}
