# GPU-Distance-Driven-projection-backprojection

If you are using this library, please cite
"GPU-based Branchless Distance-Driven Projection and Backprojection", IEEE Transactions on Computational Imaging. DOI: 10.1109/TCI.2017.2675705
Author : Rui Liu, Lin Fu, Bruno De Man, Hengyong Yu

Email: ctlab.ruiliu@gmail.com

NOTE: currently, we only support linux platform. If you want to use it in Windows, please contact ctlab.ruiliu@gmail.com or liurui1217@gmail.com.
Thanks.

The interface of projection and backprojection
[out, cfg] = dd3(flag_fw, cfg, in, view_ind, mask, bp_type);

flag_fw: Indicate which mode you want to use. The data type is a string. You can have the following options:

        1. 'fp'                         multi-CPUs based distance-driven forward projection (try to use all of your cores)
        2. 'fp1'                        one core CPU based DD forward projection
        3. 'fp8'                        8 core CPU based dd forward projection
        4. 'bp'                         multi-CPUs based distance-driven backprojection (try to use all of your cores)
        5. 'bp1'                        one core CPU based DD backprojection
        6. 'bp8'                        8 core CPU based DD backprojection
        7. 'fp_gpu_branchless'          GPU based branchless DD projection 
        8. 'fp_gpu_pseudo_dd'           GPU based pseudo DD projection
        9. 'fp_gpu_soft_interp'         GPU based software-interpolation based DD projection 
        10.'fp_gpu_branches'            GPU based Distance Driven projection
        11.'bp_gpu_branchless'          GPU based branchless DD backprojection
        12.'bp_gpu_pseudo_dd'           GPU based pseudo DD backprojection.
        13.'bp_gpu_soft_interp'         GPU baesd software-interpolation based DD backprojection
        14.'bp_gpu_zline'               GPU based z-line DD backprojection

cfg: The geometry configuration, you have to fill in the following fields manually:

        cfg.acq.sid = 541;                                                  %    source to object distance
        cfg.acq.sdd = 949;                                                  %    source to detector distance
        cfg.acq.views_per_rotation = rotN;                                  %    how many views are sampled per rotation
        cfg.acq.num_data_views = rotN;                                      %    Totally how many views are sampled (if helical scan, it will be different from the previous variable)
        cfg.acq.row_count = DNV;                                            %    number of detector cell along transverse plane
        cfg.acq.total_col_count = DNU;                                      %    number of detector cell along bench moving direction
        cfg.acq.col_width = 1.0239 / 888 * DNU;                             %    Detector cell width size along transverse plane
        cfg.acq.col_height = 1.096439 / 64 * DNV;                           %    Detector cell height size along bench moving direction
        cfg.acq.rotation_direction = 1;                                     %    NOTE: this parameter is currently ignored
        cfg.acq.first_view_angle = 0;                                       %    The first view (unit: degree)
        cfg.acq.first_view_zposition = 0;                                   %    First source position in z direction
        cfg.acq.helical_pitch = 0;                                          %    Helical pitch
        cfg.acq.col_height_at_iso = 0.625;                                  %    optional parameters (automatically provided if obtained from real scan files) 

        cfg.recon.col_center = (DNU+1) / 2;                                 %    Center index of the detector cell along transverse direction (1 based)
        cfg.recon.row_center = (DNV+1) / 2;                                 %    Center index of the detector cell along bench moving direction (1 based, NOTE: it is ignored in current version)
        cfg.recon.recon_pixels_z = ZN;                                      %    number of voxels of the image volume to be reconstructed along bench moving direction
        cfg.recon.recon_pixels_y = XN;                                      %    number of voxels of the image volume to be reconstructed along x axis.
        cfg.recon.recon_pixels_x = YN;                                      %    number of voxels of the image volume to be reconstructed along y axis.
        cfg.recon.dfov_mm = 250 ;                                           %    Field of view (in : mm)
        cfg.recon.recon_slice_spacing = 0.625 / 64 * ZN;                    %    image slice thickness along bench moving direction
        cfg.recon.recon_center_z = 0;                                       %    center Z position of the image volume in the world coordinate.

in:  input data

       if forward projection, it is the image volume stored in [Z, X, Y] order.
       if backprojection, it is the projection data stored in [ZN, DN, PN] order : ZN is the number of detector cell along bench moving direction, DN is the number of detector cell along transverse plane and PN is the number of views.

view_ind (optional):

       The indices of the views that will be used for projection/backprojection. It will be helpful in OS technique. If not provided, all views will be applied.

mask (optional):

       The mask used along transverse plane to indicate the ROI. If not provided, the ones(X,Y) mask will be used.

bp_type (optional):

        It can be used for PCG algorithm. However, DO NOT USE IT!!
