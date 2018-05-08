
% load system cfg structure
	% all in unit of mm, except otherwise noted
	cfg.acq.sid = 541;
	cfg.acq.sdd = 949;
	cfg.acq.views_per_rotation = 984;
	cfg.acq.num_data_views = 984;
	cfg.acq.row_count = 64;
	cfg.acq.total_col_count = 888;
	cfg.acq.col_width = 1.0239;
	cfg.acq.col_height = 1.096439;
	cfg.acq.rotation_direction = 1; % NOTE: this parameter is currently ignored in dd3
	cfg.acq.first_view_angle = 0; % in the unit of degree
	cfg.acq.first_view_zposition = 0;
	cfg.acq.helical_pitch = 0;
	% optional parameters (automatically provided if obtained from real scan files)
	cfg.acq.col_height_at_iso = 0.625;

	cfg.recon.col_center = 444.75; % one-based detector index
	cfg.recon.row_center = 8.5; % one-based detector index. NOTE: this parameter is currently ignored in dd3
	cfg.recon.recon_pixels_z = 64;
	cfg.recon.recon_pixels_y = 512;
	cfg.recon.recon_pixels_x = 512;
	cfg.recon.dfov_mm = 250;
	cfg.recon.recon_slice_spacing = 0.625;
	cfg.recon.recon_center_z = 0;

	% create a 3D volume with a line object in z
	x = ones( cfg.recon.recon_pixels_z, cfg.recon.recon_pixels_x, cfg.recon.recon_pixels_y, 'single' );
	x( :, ceil(size(x,2)/4), ceil(size(x,3)/2) ) = 1;
	figure; imagesc( squeeze(x(8,:,:)) );

	% test forward and back projection
    tic;
	y  = dd3( 'fp_gpu_branchless_sat2d', cfg, x );
    toc;
	figure; imagesc( squeeze(y(32,:,:)) );
	tic;
	xx = dd3( 'bp_gpu_branchless_sat2d', cfg, y );
    toc;
	figure; imagesc( squeeze(xx(32,:,:)) );
