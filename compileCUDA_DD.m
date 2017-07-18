%------------------------------------------------------------------------- 
% GE Confidential. General Electric Proprietary Data (c) 2015 General Electric Company
% Date: Nov 20, 2015
% Routine: compileCUDA_DD.m
% Author
%	Rui Liu, Lin Fu
% Organization: 
%  GE Global Research and GE Healthcare, General Electric Company
%-------------------------------------------------------------------------
% % 
% nvmex('DD3_GPU_proj_ker.cu');
% matlabroot '/extern/include 
system( '/usr/local/cuda/bin/nvcc -Xcompiler -fopenmp -O3 --use_fast_math --compile   -o DD3_GPU_proj_ker.o  --compiler-options -fPIC  -I"/software/matlab_r2014a/extern/include " -I/usr/local/cuda/include -I/usr/local/cuda/samples/common/inc "DD3_GPU_proj_ker.cu" ' );

% no need to run mex here
% mex DD3_GPU_proj_ker.o -L/usr/local/cuda/lib64 -lcudart -lgomp

% nvmex('DD3_GPU_back_ker.cu');
system( '/usr/local/cuda/bin/nvcc -Xcompiler -fopenmp -O3 --use_fast_math --compile   -o DD3_GPU_back_ker.o  --compiler-options -fPIC  -I"/software/matlab_r2014a/extern/include " -I/usr/local/cuda/include -I/usr/local/cuda/samples/common/inc "DD3_GPU_back_ker.cu" ' );
mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -fopenmp" -L/usr/local/cuda/lib64 -lcudart -lgomp DD3_GPU_mex.cpp DD3_GPU_proj_ker.o DD3_GPU_back_ker.o

% 
% % Compile qGGMRF
% system( '/usr/local/cuda/bin/nvcc -Xcompiler  -O3 --use_fast_math --compile   -o qGGMRF.o  --compiler-options -fPIC  -I"/usr/local/MATLAB/R2015b/extern/include " -I/usr/local/cuda/include -I/usr/local/cuda/samples/common/inc "qGGMRF.cu" ' );
% mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -std=c++11" -L/usr/local/cuda/lib64 -lcudart qGGMRF_MEX.cpp qGGMRF.o
