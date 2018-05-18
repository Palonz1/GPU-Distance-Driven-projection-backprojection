%------------------------------------------------------------------------- 
% Date: MAY 18, 2018
% Routine: compileCUDA_DD.m
% Author
%	Rui Liu
%-------------------------------------------------------------------------
if isunix
    system( 'nvcc -Xcompiler -fopenmp -O3 --use_fast_math --compile -o DD3_GPU_proj_ker.o  --compiler-options -fPIC  -I/usr/local/cuda/include -I/usr/local/cuda/samples/common/inc "DD3_GPU_proj_ker.cu" ' );
    system( 'nvcc -Xcompiler -fopenmp -O3 --use_fast_math --compile -o DD3_GPU_back_ker.o  --compiler-options -fPIC  -I/usr/local/cuda/include -I/usr/local/cuda/samples/common/inc "DD3_GPU_back_ker.cu" ' );
    mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -fopenmp" -L"\usr\local\CUDA\lib64\" -lcudart  DD3_GPU_mex.cpp DD3_GPU_proj_ker.o DD3_GPU_back_ker.o
elseif ispc
    system( 'nvcc.exe -v -Xcompiler -openmp -O3 --compile --use_fast_math -o DD3_GPU_proj_ker.lib "DD3_GPU_proj_ker.cu" ' );
    system( 'nvcc.exe -v -Xcompiler -openmp -O3 --compile --use_fast_math -o DD3_GPU_back_ker.lib "DD3_GPU_back_ker.cu" ' );
    mex -v -largeArrayDims COMPFLAGS="$COMPFLAGS /MT" -lcudart  DD3_GPU_mex.cpp DD3_GPU_proj_ker.lib DD3_GPU_back_ker.lib    
elseif ismac
    % We do not support mac yet.
end


