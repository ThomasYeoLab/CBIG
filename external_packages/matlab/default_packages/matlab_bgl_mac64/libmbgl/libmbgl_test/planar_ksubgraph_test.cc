#include <include/matlab_bgl.h>
#include <algorithm>
#include <iostream>

int main(int argc, char **argv) {
    mbglIndex n = 5;
    mbglIndex rp[] = {0,4,8,12,16,20};
    mbglIndex ci[] = {1,2,3,4,0,2,3,4,0,1,3,4,0,1,2,4,0,1,2,3};
	mbglIndex ki[20];
	mbglIndex kj[20];
	int is_planar = 0;
	mbglIndex nedges = 0;
    boyer_myrvold_planarity_test(n, ci, rp, &is_planar, ki, kj, &nedges,
	  NULL, NULL);
    std::cout << "is_planar = " << is_planar << std::endl;
    std::cout << "nedges = " << nedges << std::endl;
    if (nedges != 10) {
      return -1;
    }
    if (is_planar != 0) {
      return -1;
    }
    return 0;
}
