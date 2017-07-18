//------------------------------------------------------------------------- 
// GE Confidential. General Electric Proprietary Data (c) 2015 General Electric Company
// Date: Sep, 2015
// Organization: 
//  GE Global Research and GE Healthcare, General Electric Company
// -------------------------------------------------------------------------


typedef unsigned char byte;

// \brief: GPU backprojection interface
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
