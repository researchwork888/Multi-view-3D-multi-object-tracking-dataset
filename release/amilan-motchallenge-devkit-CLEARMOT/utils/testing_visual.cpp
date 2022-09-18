#include <matrix.h>
#include "mex.h"
#include <cmath>
#include <omp.h>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <unordered_map>
using namespace std;

double x1;
double y1;
double z1;
double xaxis1;
double yaxis1;
double zaxis1 ;
double x2; 
double y2;
double z2;
double xaxis2;
double yaxis2;
double zaxis2;

    double topright_corner_1x , topright_corner_1y;
    topright_corner_1x = x1 + xaxis1;
    topright_corner_1y = y1 + yaxis1;
    
    double bottomleft_corner_1x , bottomleft_corner_1y;
    bottomleft_corner_1x = x1 - xaxis1;
    bottomleft_corner_1y = y1 - yaxis1;
    
    double z1_top, z1_bottom;
    z1_top = z1 + zaxis1;
    z1_bottom = z1 - zaxis1;
    
    double rec_area1, vol_1;
    rec_area1 = (topright_corner_1x - bottomleft_corner_1x) * (topright_corner_1y - bottomleft_corner_1y);
    vol_1 = rec_area1*(z1_top - z1_bottom);
    
    double topright_corner_2x , topright_corner_2y;
    topright_corner_2x = x2 + xaxis2;
    topright_corner_2y = y2 + yaxis2;
    
    double bottomleft_corner_2x , bottomleft_corner_2y;
    bottomleft_corner_2x = x2 - xaxis2;
    bottomleft_corner_2y = y2 - yaxis2;
    
    double z2_top, z2_bottom;
    z2_top = z2 + zaxis2;
    z2_bottom = z2 - zaxis2;
    
    double rec_area2 ,vol_2;
    rec_area2 = (topright_corner_2x - bottomleft_corner_2x) * (topright_corner_2y - bottomleft_corner_2y);
    vol_2 = rec_area2*(z2_top - z2_bottom);
    
    min(topright_corner_1x,topright_corner_2x)

    
   
