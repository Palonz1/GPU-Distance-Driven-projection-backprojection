/*
 * GE Confidential. General Electric Proprietary Data (c) 2015 General Electric Company
 * Organization: 
 *  GE Global Research and GE Healthcare, General Electric Company
 * DD3_GPU_mex.cpp
 * Matlab mex gateway routine for the GPU based distance-driven projector
 *
 * author: Rui Liu (GE Global Research)
 * date: 2015.09.01
 * version: 1.0
 */

#include "mex.h"
#include "matrix.h"
#include <cstring>
#include <iostream>
#include "omp.h"

#include "DD3_GPU_proj.h"
#include "DD3_GPU_back.h"
typedef unsigned char byte;

// C interface for GPU projection
extern "C"
void DD3Proj_gpu(
			  float x0, float y0, float z0, // src position in gantry coordinates
			  int DNU, int DNV,  // detector channel#, detector row#
			  float* xds, float* yds, float* zds, // detector center position in x, y, z dimension
			  float imgXCenter, float imgYCenter, float imgZCenter, // img center in world coordinate
			  float* hangs, float* hzPos, int PN, // view angles, src position in Z in world coordinate, view#
			  int XN, int YN, int ZN, // dimension of image volume
			  float* hvol, float* hprj, // img data, projection data
			  float dx, float dz, // img size in xy, img size in z
			  byte* mask, int gpunum, int prjMode);
// C interface for GPU backprojection
extern "C"
void DD3Back_gpu(
			  float x0, float y0, float z0, // src position in gantry coordinates
			  int DNU, int DNV,  // detector channel#, detector row#
			  float* xds, float* yds, float* zds, // detector center position in x, y, z dimension
			  float imgXCenter, float imgYCenter, float imgZCenter, // img center in world coordinate
			  float* hangs, float* hzPos, int PN, // view angles, src position in Z in world coordinate, view#
			  int XN, int YN, int ZN, // dimension of image volume
			  float* hvol, float* hprj, // img data, projection data
			  float dx, float dz, // img size in xy, img size in z
			  byte* mask, int gpunum, int squared, int prjMode);


namespace CTMBIR
{
  void DD3_GPU_help_mex()
  {
    std::cout<<"Usage DD3_GPU_Proj:\n"
             <<"y = function('Proj', x0, y0, z0, nrdetcols, nrdetrows, *xds, *yds, *zds,*xdsl,*ydsl,*xdsr,*ydsr, imgXoffset, imgYoffset, imgZoffset, *viewangles, *zshifts, nrviews, nrcols, nrrows, nrplanes, *pOrig, vox_xy_size, vox_zsize, xy_mask, gpunum)\n"
	     <<"float x0, source x coordinate (before rotating)\n"
	     <<"float y0, source y coordinate (before rotating)\n"
	     <<"float z0, source z coordinate (before rotating)\n"
	     <<"int nrdetcols, number of detector columns (in-plane)\n"
	     <<"int nrdetrows, number of detector rows (in-z)\n"
	     <<"float* xds, nrdetcols detector x coordinates (before rotating)\n"
	     <<"float* yds, nrdetcols detector y coordinates (before rotating)\n"
	     <<"float* zds, nrdetrows detector z coordinates (before rotating)\n"
	     <<"float imgXoffsets, the x-offset of the image center relative to COR\n"
	     <<"float imgYoffsets, the y-offset of the image center relative to COR\n"
	     <<"float imgZoffsets, the z-offset of the image center relative to COR\n"
	     <<"float* viewangles, nrviews rotation angles\n"
	     <<"float* zshifts, nrviews z-increments of the patient/grantry(for helical)\n"
	     <<"int nrviews, number of angles\n"
	     <<"int nrcols, number of columns in image\n"
	     <<"int nrrows, number of rows in image\n"
	     <<"int nrplanes, number of planes in image\n"
	     <<"float* pOrig, least significant index is plane, then row, then col\n" // What is that?
	     <<"float vox_xy_size, voxel size in x and y\n"
	     <<"float vox_z_size, voxel size in z\n"
	     <<"byte* xy_mask, xy plane [mask, mask_trans]\n"
	     <<"int gpunum, multi GPU in views\n"
             <<"int prjMode, Different projection modes"
	     <<"out: [det-col, det-row, view] sinogram\n"
	     <<"no need to input zero sinogram\n"
	     <<"y = function('Back', x0, y0, z0, nrdetcols, nrdetrows, *xds, *yds, *zds, imgXoffset, imgYoffset, imgZoffset, *viewangles,*zshifts,nrviews,*sinogram,nrcols,nrrows,nrplanes,vox_xy_size,vox_z_size,xy_mask,gpunum,projector_type)\n"
	     <<"int projector_type: 0 corresponds to a standard backprojector\n"
	     <<"                    1 corresponds to a squared backprojector\n"
	     <<"float* sinogram, least significant index is det-row, then det-col then view\n"
	     <<"out: [plane row col] image\n"
	     <<"no need to input zero image\n";
  }
}


void DD3_GPU_Proj_mex(
   mxArray* plhs[], // [view det-row det-col]
   const mxArray* mx_x0, //source x coordinate (before rotating)
   const mxArray* mx_y0, //source y coordinate (before rotating)
   const mxArray* mx_z0, //source z coordinate (before rotating)
   const mxArray* mx_nrdetcols, //number of detector columns (in-plane)
   const mxArray* mx_nrdetrows, //number of detector rows (in-Z)
   const mxArray* mx_xds, //nrdetcols detector x coordinates (before rotating)
   const mxArray* mx_yds, //nrdetcols detector y coordinates (before rotating)
   const mxArray* mx_zds, //nrdetrows detector z coordinates (before rotating)
   const mxArray* mx_imgXoffset, //the x-offset of the image center relative to COR
   const mxArray* mx_imgYoffset, //the y-offset of the image center relative to COR
   const mxArray* mx_imgZoffset, //the z-offset of the image center relative to COR
   const mxArray* mx_viewangles, //nrviews rotation angles
   const mxArray* mx_zshifts, //nrviews z-increments of the patient gantry (for helical)
   const mxArray* mx_nrviews, //number of angles
   const mxArray* mx_nrcols, // number of columns in image
   const mxArray* mx_nrrows, // number of rows in image
   const mxArray* mx_nrplanes, //number of planes in image
   const mxArray* mx_pOrig, //least significant index is plane, then row, then col
   const mxArray* mx_vox_xy_size, //voxel size in x and y
   const mxArray* mx_vox_z_size, //voxel size in Z
   const mxArray* mx_xy_mask, //xy plane [mask, mask_trans]
   const mxArray* mx_gpunum,//multi gpu in views
   const mxArray* mx_prjMode)  //projection modes
{
  int nrdetcols = *((int*)mxGetData(mx_nrdetcols));
  int nrdetrows = *((int*)mxGetData(mx_nrdetrows));
  int nrviews = *((int*)mxGetData(mx_nrviews));
  //Create output array of class mxREAL and mxSINGLE_CLASS
  const mwSize dims[] = {nrdetrows, nrdetcols, nrviews};
  plhs[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL);

  float* xds = (float*)mxGetPr(mx_xds);
  float* yds = (float*)mxGetPr(mx_yds);
  float* zds = (float*)mxGetPr(mx_zds);

  DD3Proj_gpu(
    *((float*)mxGetData(mx_x0)),
    *((float*)mxGetData(mx_y0)),
    *((float*)mxGetData(mx_z0)),
    nrdetcols, nrdetrows,
    xds,
    yds,
    zds ,
    *((float*)mxGetData(mx_imgXoffset)),
    *((float*)mxGetData(mx_imgYoffset)),
    *((float*)mxGetData(mx_imgZoffset)),
    (float*) mxGetPr(mx_viewangles),
    (float*) mxGetPr(mx_zshifts),
    *((int*) mxGetData(mx_nrviews)),
    *((int*) mxGetData(mx_nrcols)),
    *((int*) mxGetData(mx_nrrows)),
    *((int*) mxGetData(mx_nrplanes)),
    (float*) mxGetPr(mx_pOrig),
    (float*) mxGetPr(plhs[0]),
    *((float*)mxGetData(mx_vox_xy_size)),
    *((float*)mxGetData(mx_vox_z_size)),
    (byte*)mxGetPr(mx_xy_mask),
    *((int*)mxGetData(mx_gpunum)),
    *((int*)mxGetData(mx_prjMode)));
}



/**
 * DD3_GPU_Back_mex()
 */
void DD3_GPU_Back_mex(
      mxArray* plhs[], //[plane, XN(row), YN(col)]
      const mxArray* mx_x0, //source x coordinate (before rotating)
      const mxArray* mx_y0, //source y coordinate (before rotating)
      const mxArray* mx_z0, //source z coordinate (before translating)
      const mxArray* mx_nrdetcols, //number of detector columns (in-plane)
      const mxArray* mx_nrdetrows, //number of detector rows(in-z)
      const mxArray* mx_xds, //nrdetcols detector x coordinate (before rotating)
      const mxArray* mx_yds, //nrdetcols detector y coordinate (before rotating)
      const mxArray* mx_zds, //nrdetrows detector z coordinate (before translating)
      const mxArray* mx_imgXoffset, // the x-offset of the image center relative to COR
      const mxArray* mx_imgYoffset, // the y-offset of the image center relative to COR
      const mxArray* mx_imgZoffset, // the z-offset of the image center relative to COR
      const mxArray* mx_viewangles, //nrviews rotation angles
      const mxArray* mx_zshifts, //nrviews z-increments ofthe patient/gantry (for helical)
      const mxArray* mx_nrviews, //number of angles,
      const mxArray* mx_nrcols, //# columns in image
      const mxArray* mx_nrrows, //# rows in image
      const mxArray* mx_nrplanes, //# planes in image
      const mxArray* mx_sinogram, //# least significant index is view, then det-row, then det-col
      const mxArray* mx_vox_xy_size, //voxel size in x and y
      const mxArray* mx_vox_z_size, //voxel size in z
      const mxArray* mx_xy_mask, // xy plane [mask/ no trans!]
      const mxArray* mx_gpunum, // how many gpu would you like to use,the maximum number is 2 and the minimum is 1
      const mxArray* mx_projector_type, //0 : standard backprojector; 1: squared backprojector
      const mxArray* mx_prjMode) // we would like to use string/char*
{
  int nrcols = *((int*)mxGetData(mx_nrcols));
  int nrrows = *((int*)mxGetData(mx_nrrows));
  int nrplanes = *((int*)mxGetData(mx_nrplanes));

  int nrdetcols = *((int*)mxGetData(mx_nrdetcols));
  int nrdetrows = *((int*)mxGetData(mx_nrdetrows));
  int nrviews = *((int*)mxGetData(mx_nrviews));

  const mwSize dims[]={nrplanes, nrrows, nrcols};
  plhs[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS, mxREAL);
  
  float* xds = (float*)mxGetPr(mx_xds);
  float* yds = (float*)mxGetPr(mx_yds);
  float* zds = (float*)mxGetPr(mx_zds);

  DD3Back_gpu(
		  *((float*)mxGetData(mx_x0)),
		  *((float*)mxGetData(mx_y0)),
		  *((float*)mxGetData(mx_z0)),
		  nrdetcols, nrdetrows,
		  xds,yds,zds,
		  *((float*)mxGetData(mx_imgXoffset)),
		  *((float*)mxGetData(mx_imgYoffset)),
		  *((float*)mxGetData(mx_imgZoffset)),
		  (float*)mxGetPr(mx_viewangles),
		  (float*)mxGetPr(mx_zshifts),
		  *((int*)mxGetData(mx_nrviews)),
		  *((int*)mxGetData(mx_nrcols)),
		  *((int*)mxGetData(mx_nrrows)),
		  *((int*)mxGetData(mx_nrplanes)),
		  (float*)(mxGetPr(plhs[0])),
		  (float*)mxGetPr(mx_sinogram),
		  *((float*)mxGetData(mx_vox_xy_size)),
		  *((float*)mxGetData(mx_vox_z_size)),
		  (byte*)mxGetPr(mx_xy_mask),
		  *((int*)mxGetData(mx_gpunum)),
	      *((int*)mxGetData(mx_projector_type)),
	      *((int*)mxGetData(mx_prjMode)));
}



/*
 * DD3_GPU_mex()
 */
void DD3_GPU_mex(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  if(nrhs < 1)
  {
    CTMBIR::DD3_GPU_help_mex();
    std::cerr<<"Error: usage\n";
  }

  if(!mxIsChar(prhs[0]))
  {
    std::cerr<<"First argument must be a string\n";
    CTMBIR::DD3_GPU_help_mex();
    std::cerr<<"Error: usage\n";
  }
  
  //Read and copy first in-out string to arg[0]
  int arg0len = mxGetM(prhs[0]) * mxGetN(prhs[0]) + 1;
  char* arg0 = (char*)mxCalloc(arg0len, sizeof(char));
  if(mxGetString(prhs[0],arg0, arg0len))
  {
    std::cerr<<"but with mxGetString\n";
  }

  //Forward projection
  if(!std::strcmp(arg0, "Proj"))
  {
    if(nrhs != 24 || nlhs != 1)
    {
      CTMBIR::DD3_GPU_help_mex();
      std::cerr<<"Error: usage\n";
    }
    
    DD3_GPU_Proj_mex(plhs,prhs[1],prhs[2],prhs[3],prhs[4],prhs[5],prhs[6],prhs[7],prhs[8],prhs[9],prhs[10],prhs[11],prhs[12],prhs[13],prhs[14],prhs[15],prhs[16],prhs[17],prhs[18],prhs[19],prhs[20],prhs[21],prhs[22],prhs[23]);
  }
  else if(!std::strcmp(arg0,"Back"))
  {
    if(nrhs != 25 || nlhs != 1)
    {
      CTMBIR::DD3_GPU_help_mex();
      std::cerr<<"Error: usage\n";
    }
    
    DD3_GPU_Back_mex(plhs,prhs[1],prhs[2],prhs[3],prhs[4],prhs[5],prhs[6],prhs[7],prhs[8],prhs[9],prhs[10],prhs[11],prhs[12],prhs[13],prhs[14],prhs[15],prhs[16],prhs[17],prhs[18],prhs[19],prhs[20],prhs[21],prhs[22],prhs[23],prhs[24]);
  }
  else
  {
    std::cerr<<"Error: unsupported action "<<(const char*)mxGetData(prhs[0])<<"\n";
  }

}





void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  if(!nlhs && !nrhs)
  {
    CTMBIR::DD3_GPU_help_mex();
    return;
  }
  DD3_GPU_mex(nlhs, plhs, nrhs, prhs);
}
