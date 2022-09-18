mex utils/MinCostMatching.cpp -outdir utils CXXFLAGS="$CXXFLAGS --std=c++11"
mex utils/clearMOTMex.cpp -outdir utils CXXFLAGS="$CXXFLAGS --std=c++11"
mex utils/costBlockMex.cpp -outdir utils COMPFLAGS="/openmp $COMPFLAGS" CXXFLAGS="$CXXFLAGS --std=c++11"


mex '-I/usr/local/Cellar/libomp/8.0.0/include/omp.h/' '/Users/jonahong/Google Drive/PHD 2018/Codes/Visual MOT Challenge DevKit/amilan-motchallenge-devkit-aa05c63476c9/utils/clearMOTMex.cpp' 

mex '/Users/jonahong/Google Drive/PHD 2018/Codes/Visual MOT Challenge DevKit/amilan-motchallenge-devkit-aa05c63476c9/utils/costBlockMex.cpp'
 
mex '-I/usr/local/Cellar/libomp/8.0.0/include/omp.h/' '/Users/jonahong/Google Drive/PHD 2018/Codes/Visual MOT Challenge DevKit/amilan-motchallenge-devkit-aa05c63476c9/utils/MinCostMatching.cpp'


matlab

lapacklib = fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft','libmwlapack.lib');
