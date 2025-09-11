steps = [29]; %change this if you want to run different steps
    % Step 1: creates .exp files for ice shelf extent
    % Step 2: calculates the max buttressing, WITH gap fill interpolation
    % Step 3: calculates the max buttresing with NO gap fill interpolation!
    % Step 4: plots the thickness of the ice shelf, with NO gap fill
    % Step 5: plots the change in thickness compared to 2006
    % Step 6: plotting ALL years of max buttressing - NO MASK THO

    % Step 7: plotting ALL years of thickness change minus 2006
    % Step 8: calculates thickness MEAN
    % Step 9: THICKNESS, NO gap fill, WTH MASK (similar to step 4)
    % Step 10: THICKNESS ANOMALY (steps 8 + 9 MUST be run BEFORE)
    % Step 11: plotting ALL years of thickness anomalies

    % Step 12: KNMAX with MASK applied
        % AND also saving as a georeferenced geotiff
    % Step 13: plotting ALL years of KNMAX with MASK together
    % Step 14: plots of KNMAX change ( MINUS 2006 )
    % Step 15: plotting all KNMAX minus 2006 change together
    % Step 16: calculating KNMAX MEAN for all years
        % and saves as a geotiff
    % Step 17: KNMAX ANOMALY (steps 12 and 16 MUST be run BEFORE this)
    % Step 18: knmax anomalies all together

    % Step 19: plotting velocity for each year
    % Step 20: calculating change in in velocity ( minus 2006)
    % Step 21: plotting all changes in velocity together
    % Step 22: calculating velocity mean
        % and saving as a geotiff
    % Step 23: plotting velocity anomalies (requires steps 19 + 22)
    % Step 24: plotting all velocity anomalies together 

    % Step 25: longitudinal strain rate - thomas params
    % Step 26: plotting all longitudinal strain years together

    % Step 27: GL flux working 1
        % THIS CODE GIVES 1.3 Gt/yr FOR EVERY YEAR
    % Step 28: GL flux working 2
    % Step 29: GL flux working 3
        % USE THIS CODE!!!
    % Step 30: visualisaing the FLUX GATE
        % before running this step, run step for chosen md model
    % Step 31: creating .exp file for GL ONLY (FLUX GATE)

    % Step 32: CHANGE in area

    % Step 33: calculating EXY
    % Step 34: plotting ALL EXY together
    % Step 35: calculating EYY
    % Step 36: plotting ALL EYY together



clustername = 'oshostname';
%% Cluster settings {{{
if strcmpi(clustername,'oshotname'),
	cluster=generic('name',oshostname(),'np',1);
	lock = 1;
	loadonly = 0;
end

% Set paths 
org = organizer('repository', ['Models'], 'prefix', ['wilkes-'], 'steps', steps, 'color', '34;47;2'); clear steps;
datadir='Data/'; %% Jeremy need to change to relevant data directory
%}}}



%	steps 1-XX: calculate fields for each year of observational data
if perform(org, 'Calc_exp'),% {{{

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	for yy = 1:length(years)
		% CONVERT SHP TO EXP
		indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory
		%shp2exp([indir 'ice_extent_' num2str(years(yy)) '.shp'],[indir 'ice_extent_' num2str(years(yy)) '_v2.exp']);
		shpa = shpread([indir 'ice_extent_' num2str(years(yy)) '.shp']);
		if max(shpa.x)<1e6
			[x,y] = ll2xy(shpa.y,shpa.x,-1);
			expa.x = x;
			expa.y = y;
		else
			expa.x = shpa.x;
			expa.y = shpa.y;
		end
		expwrite(expa, [indir 'ice_extent_' num2str(years(yy)) '.exp']);
	end


end%}}}


%%%%%%%%%%%%%%%%% STEP TWO:

if perform(org, 'Buttressing_with_gapfill'),% {{{

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	for yy = 1:length(years)
		indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory

        % indir2 = ['z20250728_buttressing_feedback/ice_extents/exps/']; % this links to different ice extent shapefiles

		in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1);
        % in=ContourToNodes(md.mesh.x,md.mesh.y,[indir2 'ice_extent_' num2str(years(yy)) '.exp'],1);
		md.mask.ocean_levelset(find(in)) = -1;
		md.mask.ice_levelset(find(in)) = -1;

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		vx_in = double(ncread(fname,'VX'));
		vy_in = double(ncread(fname,'VY'));
		thickness_in = double(ncread(fname,'Thickness'));

		% Get the vx, vy onto the ISSM grid
		vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);

		vxnew = md.inversion.vx_obs;
		vynew = md.inversion.vy_obs;
		vxnew(find(in)) = vx(find(in));
		vynew(find(in)) = vy(find(in));

		flag = in & isnan(vxnew);
		pos = find(flag); pos2 = find(~flag);
		vxnew(pos) = griddata(md.mesh.x(pos2),md.mesh.y(pos2), vxnew(pos2), md.mesh.x(pos), md.mesh.y(pos));
		md.inversion.vx_obs(find(in)) = vxnew(find(in));

		flag = in & isnan(vynew);
		pos = find(flag); pos2 = find(~flag);
		vynew(pos) = griddata(md.mesh.x(pos2),md.mesh.y(pos2), vynew(pos2), md.mesh.x(pos), md.mesh.y(pos));
		md.inversion.vy_obs(find(in)) = vynew(find(in));

		md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);

		% Get the thickness onto the ISSM grid
		thickness = InterpFromGrid(x_in,y_in,thickness_in',md.mesh.x,md.mesh.y);
		Hnew = md.geometry.thickness;
		Hnew(find(in)) = thickness(find(in));

		flag = in & isnan(Hnew);
		pos = find(flag); pos2 = find(~flag);
		Hnew(pos) = griddata(md.mesh.x(pos2),md.mesh.y(pos2), Hnew(pos2), md.mesh.x(pos), md.mesh.y(pos));
		md.geometry.thickness(find(in)) = Hnew(find(in));

		% calculate the buttressing numbers
		[Knflow{yy},Nflow{yy},Knmax{yy},Nmax{yy},N0{yy}]=buttressingnumber(md);

		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;

		plotmodel(md,'data',Knmax{yy},'figure', yy,'xlim#all',[0.9e6, 1.2e6],'ylim#all',[-2.2e6, -2e6],'colormap#1','redblue','caxis',[-5 5],'mask#all',mask,'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp'])
		exportgraphics(gcf,['Plots/buttressing_calc_' num2str(years(yy)) '_gapfill.png'])
        % exportgraphics(gcf,['Plots/buttressing_calc_' num2str(years(yy)) '_gapfill_indir2.png'])
	end

	%if want to save all output
	save([indir 'buttressing_alloutfiles_gapfill.mat'], 'Knflow', 'Nflow', 'Knmax', 'Nmax','masks');


end%}}}



%%%%%%%%%%%%%%%%% STEP THREE:

% KNMAX = MAXIMUM BUTTRESSING

if perform(org, 'Buttressing_without_gapfill'),% {{{

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];

	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	for yy = 1:length(years)
		indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory
        % indir2 = ['z20250728_buttressing_feedback/ice_extents/exps/']; % this links to different ice extent shapefiles


		in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1);
        % in=ContourToNodes(md.mesh.x,md.mesh.y,[indir2 'ice_extent_' num2str(years(yy)) '.exp'],1);
		md.mask.ocean_levelset(find(in)) = -1;
		md.mask.ice_levelset(find(in)) = -1;

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		vx_in = double(ncread(fname,'VX'));
		vy_in = double(ncread(fname,'VY'));
		thickness_in = double(ncread(fname,'Thickness'));

		% Get the vx, vy onto the ISSM grid
		vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		md.inversion.vx_obs(find(in)) = vx(find(in));
		md.inversion.vy_obs(find(in)) = vy(find(in));
		md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);

		% Get the thickness onto the ISSM grid
		thickness = InterpFromGrid(x_in,y_in,thickness_in',md.mesh.x,md.mesh.y);
		md.geometry.thickness(find(in)) = thickness(find(in));

		% calculate the buttressing numbers
		[Knflow{yy},Nflow{yy},Knmax{yy},Nmax{yy},N0{yy}]=buttressingnumber(md);

		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;

		%if want to save output per iteration of the for loop
		% plotmodel(md,'data',Knmax{yy},'figure', yy,'xlim#all',[0.9e6, 1.2e6],'ylim#all',[-2.2e6, -2e6],'colormap#1','redblue','caxis',[-5 5],'mask#all',mask,'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp'])
		
        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        plotmodel(md,...
            'data',Knmax{yy},...
            'figure', yy,...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'colormap#1','redblue',...
            'caxis',[-5 5],...
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        cb = colorbar;
        ylabel(cb, 'Maximum Buttressing Number', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['Knmax ' num2str(years(yy))], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        % exportgraphics(gcf,['Plots/buttressing_calc_' num2str(years(yy)) '_nogapfill.png'])
        % exportgraphics(gcf,['Plots/buttressing_calc_' num2str(years(yy)) '_nogapfill_indir2.png'])
        exportgraphics(gcf,['Plots/buttressing_calc_' num2str(years(yy)) '_nogapfill_3.png'])
	end

	%if want to save all output
	save([indir 'buttressing_alloutfiles_no_gapfill.mat'], 'Knflow', 'Nflow', 'Knmax', 'Nmax','masks');

end%}}}



%%%%%%%%%%%%%%%%% STEP FOUR: THICKNESS

if perform(org, 'Thickness_without_gapfill'),% {{{

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	for yy = 1:length(years)
		indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory

		in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1);
		md.mask.ocean_levelset(find(in)) = -1;
		md.mask.ice_levelset(find(in)) = -1;

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		thickness_in = double(ncread(fname,'Thickness'));

		% Get the thickness onto the ISSM grid
		thickness = InterpFromGrid(x_in,y_in,thickness_in',md.mesh.x,md.mesh.y);
		md.geometry.thickness(find(in)) = thickness(find(in));
		H = md.geometry.thickness;

		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;

		%if want to save output per iteration of the for loop
		% plotmodel(md,'data',H,'figure', yy,'xlim#all',[0.9e6, 1.2e6],'ylim#all',[-2.2e6, -2e6],'caxis',[300 1000],'mask#all',mask,'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp'])
		
        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        plotmodel(md,...
            'data',H,...
            'figure', yy,...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'colormap#1','parula',...
            'caxis',[300 1000],...
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        cb = colorbar;
        ylabel(cb, 'Ice Thickness (m)', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['Thickness ' num2str(years(yy))], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        
        exportgraphics(gcf,['Plots/thickness_' num2str(years(yy)) '_nogapfill.png'])
	end

	%if want to save all output
	save([indir 'thickness_alloutfiles.mat'], 'H','masks');

end%}}}




%%%%%%%%%%%%%%%%% STEP FIVE: THICKNESS CHANGE COMPARED TO 2006

if perform(org, 'Thickness_minus_2006'),% {{{

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];

    % Load model and 2006 reference thickness
    md = loadmodel(org, 'Inversion_BHO');
    
    % Read 2006 thickness data
    refyear = 2006;
    indir_ref = ['Data/' num2str(refyear) '/'];
    in_ref = logical(ContourToNodes(md.mesh.x, md.mesh.y, [indir_ref 'ice_extent_' num2str(refyear) '.exp'], 1));
    fname_ref = [indir_ref num2str(refyear) '.nc'];
    x_in = double(ncread(fname_ref, 'x'));
    y_in = double(ncread(fname_ref, 'y'));
    thickness_in = double(ncread(fname_ref, 'Thickness'));
    thickness_ref = InterpFromGrid(x_in, y_in, thickness_in', md.mesh.x, md.mesh.y);
	    velo = md.inversion.vel_obs;

    delta_H_all = {};
    masks = {};

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	
    for yy = 1:length(years)
        thisyear = years(yy);
        indir = ['Data/' num2str(thisyear) '/'];
    
        md.mask.ice_levelset = ones(md.mesh.numberofvertices, 1);
        md.mask.ocean_levelset = ones(md.mesh.numberofvertices, 1);

        in = logical(ContourToNodes(md.mesh.x, md.mesh.y, [indir 'ice_extent_' num2str(thisyear) '.exp'], 1));
        md.mask.ocean_levelset(in) = -1;
        md.mask.ice_levelset(in) = -1;
    
        fname = [indir num2str(thisyear) '.nc'];
        x_in = double(ncread(fname, 'x'));
        y_in = double(ncread(fname, 'y'));
        thickness_in = double(ncread(fname, 'Thickness'));
        thickness = InterpFromGrid(x_in, y_in, thickness_in', md.mesh.x, md.mesh.y);
    
        % Difference from 2006
        delta_H = thickness - thickness_ref;
    
        % Mask
        mask = ones(md.mesh.numberofvertices,1);
        mask(md.mask.ice_levelset > 0) = 0;
        delta_H(mask == 0) = NaN;
    
        delta_H_all{yy} = delta_H;
        masks{yy} = mask;

        % Plot thickness change
        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
    
        plotmodel(md,...
            'data', delta_H,...
            'figure', yy,...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'caxis', [-10 10],...  % adjust ?
            'mask#all', mask,...
            'expdisp', [indir 'ice_extent_' num2str(thisyear) '.exp']);
    
        colormap(redblue(64));
        cb = colorbar;
        ylabel(cb, 'ΔThickness (m)', 'FontName', 'Times New Roman', 'FontSize', 12);
        xlabel('Eastings (m)', 'FontName', 'Times New Roman', 'FontSize', 12);
        ylabel('Northings (m)', 'FontName', 'Times New Roman', 'FontSize', 12);
        title(['Thickness Change: ' num2str(thisyear) ' – 2006'], 'FontName', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 12);
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

    
        exportgraphics(gcf, ['Plots/thickness_change_' num2str(thisyear) '_minus_2006_redblue.png']);
    end


	%if want to save all output
	save([indir 'thickness_change_minus_2006_alloutfiles.mat'], 'delta_H','masks');

end%}}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % After the for loop finishes, save the last year's delta_H as GeoTIFF
% 
% % Choose the last year's delta_H and corresponding mask
% delta_H_last = delta_H_all{end};
% mask_last = masks{end};
% 
% % Coordinates of mesh points
% x_coords = md.mesh.x;
% y_coords = md.mesh.y;
% 
% % Define grid spacing (adjust dx, dy as needed)
% dx = 500; %or 100?
% dy = 500; 
% 
% % Define grid limits based on min/max coordinates (or your xlim/ylim)
% xlim_grid = [min(x_coords) max(x_coords)];
% ylim_grid = [min(y_coords) max(y_coords)];
% 
% % Create regular grid vectors
% xrange = xlim_grid(1):dx:xlim_grid(2);
% yrange = ylim_grid(1):dy:ylim_grid(2);
% 
% % Create meshgrid
% [Xgrid, Ygrid] = meshgrid(xrange, yrange);
% 
% % Interpolate scattered data (delta_H) onto the regular grid
% F = scatteredInterpolant(x_coords, y_coords, delta_H_last, 'linear', 'none');
% Zgrid = F(Xgrid, Ygrid);
% 
% % Define referencing object for GeoTIFF
% R = maprefcells([min(xrange) max(xrange)], [min(yrange) max(yrange)], size(Zgrid));
% 
% 
% 
% % Save GeoTIFF
% geotiffwrite('Data/thickness_2006_2018_delta_500.tif', Zgrid, R, 'CoordRefSysCode', 'EPSG:3031');
% 
% fprintf('GeoTIFF saved: Data/thickness_2006_2018_delta_500.tif\n');


%%%%%%%%%%%%%%%%% STEP SIX: PLOTTING ALL MAX BUTTRESSING TOGETHER

% KNMAX = MAXIMUM BUTTRESSING

if perform(org, 'Buttressing_ALL_years'),% {{{

	files = {};
for yy = 1:length(years)
    files{end+1} = ['Plots/buttressing_calc_' num2str(years(yy)) '_nogapfill_3.png'];
end

% Create montage
figure;
montage(files, 'Size', [3 4], 'BackgroundColor', 'w'); % 3 rows, 4 cols — adjust as needed
title('Knmax for All Years');
exportgraphics(gcf, 'Plots/Knmax_all_years.png', 'Resolution', 300);

end%}}}


%%%%%%%%%%%%%%%%% STEP SEVEN: PLOTTING ALL THICKNESS CHANGE TOGETHER

% KNMAX = MAXIMUM BUTTRESSING

if perform(org, 'Thickness_Change_ALL_years'),% {{{

	files = {};
for yy = 1:length(years)
    files{end+1} = ['Plots/thickness_change_' num2str(years(yy)) '_minus_2006_redblue.png'];
end

% Create montage
figure;
montage(files, 'Size', [3 4], 'BackgroundColor', 'w'); % 3 rows, 4 cols — adjust as needed
title('Thickness Change for All Years');
exportgraphics(gcf, 'Plots/change_thickness_all_years.png', 'Resolution', 300);

end%}}}




%%%%%%%%%%%%%%%%% STEP EIGHT: THICKNESS MEAN


if perform(org, 'Thickness_MEAN'),% {{{


    years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];

    md = loadmodel(org, 'Inversion_BHO');
    velo = md.inversion.vel_obs;

    Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	

    n_nodes = md.mesh.numberofvertices;
    n_years = length(years);
    thickness_all_years = NaN(n_nodes, n_years);

    % Calculating the mean across all years

    for yy = 1:length(years)
		    
            md.mask.ice_levelset(:) = 1;   
            md.mask.ocean_levelset(:) = 1;


            indir = ['Data/' num2str(years(yy)) '/']; 
    
		    in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1); % makes values within the .exp file = 1
		    md.mask.ocean_levelset(find(in)) = -1; % sets ocean_levelset to -1 for nodes inside the ice extent
		    md.mask.ice_levelset(find(in)) = -1; % sets ice_levelset to -1 for nodes inside the ice extent
    
		    % READ THE ITS_LIVE DATA ONTO ISSM MESH
		    fname = [indir num2str(years(yy)) '.nc'];
		    x_in = double(ncread(fname,'x'));
		    y_in = double(ncread(fname,'y'));
		    thickness_in = double(ncread(fname,'Thickness'));
    
		    % Get the thickness onto the ISSM grid
		    thickness = InterpFromGrid(x_in,y_in,thickness_in',md.mesh.x,md.mesh.y);
		    md.geometry.thickness(find(in)) = thickness(find(in));
		    H = md.geometry.thickness;
    
		    % mask= ones(md.mesh.numberofvertices,1);
		    % mask(find(md.mask.ice_levelset>0)) = 0;
		    % masks{yy} = mask;

            mask = ones(md.mesh.numberofvertices,1);
            mask(md.mask.ice_levelset > 0) = 0;
            masks{yy} = mask;
            H(mask == 0) = NaN; % removes thickness values outside the ice extent

            thickness_all_years(:,yy) = H;
    end


    thickness_mean = mean(thickness_all_years, 2, 'omitnan');


    plotmodel(md,...
        'data', thickness_mean,...
        'figure', yy,...
        'xlim#all',[998648.2458, 1171935.9538],...
        'ylim#all', [-2148431.9935, -2033743.9935],...
        'caxis', [0 1000],...  % adjust ?
        'mask#all', mask,...
        'expdisp', [indir 'ice_extent_' num2str(thisyear) '.exp']);
 
    title('Mean Thickness');
    exportgraphics(gcf, 'Plots/MEAN_thickness.png', 'Resolution', 300);

	save([indir 'thickness_MEAN.mat'], 'thickness_mean','masks');


 end%}}}





%%%%%%%%%%%%%%%%% STEP NINE: THICKNESS WITH MASK

if perform(org, 'Thickness_without_gapfill_MASKED'),% {{{

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	for yy = 1:length(years)
		
        md.mask.ice_levelset(:) = 1;   
        md.mask.ocean_levelset(:) = 1; 
        
        indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory

        in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1); % makes values within the .exp file = 1        
        md.mask.ocean_levelset(find(in)) = -1; % sets ocean_levelset to -1 for nodes inside the ice extent
		md.mask.ice_levelset(find(in)) = -1; % sets ice_levelset to -1 for nodes inside the ice extent

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		thickness_in = double(ncread(fname,'Thickness'));

		% Get the thickness onto the ISSM grid
		thickness = InterpFromGrid(x_in,y_in,thickness_in',md.mesh.x,md.mesh.y);
		md.geometry.thickness(find(in)) = thickness(find(in));
		H = md.geometry.thickness;

        mask = ones(md.mesh.numberofvertices,1);
        mask(md.mask.ice_levelset > 0) = 0;
        masks{yy} = mask;
        H(mask == 0) = NaN; % removes thickness values outside the ice extent

		%if want to save output per iteration of the for loop
		% plotmodel(md,'data',H,'figure', yy,'xlim#all',[0.9e6, 1.2e6],'ylim#all',[-2.2e6, -2e6],'caxis',[300 1000],'mask#all',mask,'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp'])
		
        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        plotmodel(md,...
            'data',H,...
            'figure', yy,...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'colormap#1','parula',...
            'caxis',[300 1000],...
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        cb = colorbar;
        ylabel(cb, 'Ice Thickness (m)', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['Thickness ' num2str(years(yy))], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        
        exportgraphics(gcf,['Plots/thickness_' num2str(years(yy)) '_nogapfill_MASKED.png'])
	end

	%if want to save all output
	save([indir 'thickness_alloutfiles_MASKED.mat'], 'H','masks');

end%}}}



%%%%%%%%%%%%%%%%% STEP TEN: THICKNESS ANOMALY 
        %% STEPS 8 AND 9 MUST BE COMPLETED FIRST

if perform(org, 'Thickness_Anomaly'),% {{{

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

    min_vals = NaN(length(years),1); % made to store min and max anom
    max_vals = NaN(length(years),1);

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	for yy = 1:length(years)
		
        md.mask.ice_levelset(:) = 1;   
        md.mask.ocean_levelset(:) = 1; 
        
        indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory

        in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1); % makes values within the .exp file = 1        
        md.mask.ocean_levelset(find(in)) = -1; % sets ocean_levelset to -1 for nodes inside the ice extent
		md.mask.ice_levelset(find(in)) = -1; % sets ice_levelset to -1 for nodes inside the ice extent

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		thickness_in = double(ncread(fname,'Thickness'));

		% Get the thickness onto the ISSM grid
		thickness = InterpFromGrid(x_in,y_in,thickness_in',md.mesh.x,md.mesh.y);
		md.geometry.thickness(find(in)) = thickness(find(in));
		H = md.geometry.thickness;

        mask = ones(md.mesh.numberofvertices,1);
        mask(md.mask.ice_levelset > 0) = 0;
        masks{yy} = mask;
        H(mask == 0) = NaN; % removes thickness values outside the ice extent

		
        % Anomaly Calculation
        anomaly_H = H - thickness_mean;

        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');

        plotmodel(md,...
            'data',anomaly_H,...
            'figure', yy,...
            'colormap#1','parula',...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'caxis',[-5 5],...       %'caxis',[-30 15],... 
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        colormap(redblue(64));
        cb = colorbar;
        ylabel(cb, 'Ice Thickness Anomaly (m)', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['Thickness Anomaly ' num2str(years(yy))], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        
        exportgraphics(gcf,['Plots/thickness_' num2str(years(yy)) '_nogapfill_MASKED_ANOMALY_redblue.png'])

        min_vals(yy) = min(anomaly_H, [], 'omitnan');
        max_vals(yy) = max(anomaly_H, [], 'omitnan');

	end

	%if want to save all output
	save([indir 'thickness_alloutfiles_MASKED_ANOMALY.mat'], 'anomaly_H','H','masks');

end%}}}



%%%%%%%%%%%%%%%%% STEP ELEVEN: PLOTTING ALL THICKNESS ANOMALY TOGETHER

% KNMAX = MAXIMUM BUTTRESSING

if perform(org, 'Thickness_ANOMALY_Change_ALL_years'),% {{{

	files = {};
for yy = 1:length(years)
    files{end+1} = ['Plots/thickness_' num2str(years(yy)) '_nogapfill_MASKED_ANOMALY_redblue.png'];
end

% Create montage
figure;
montage(files, 'Size', [3 4], 'BackgroundColor', 'w'); % 3 rows, 4 cols — adjust as needed
title('Thickness Anomaly Change for All Years');
exportgraphics(gcf, 'Plots/anomaly_change_thickness_all_years_redblue.png', 'Resolution', 300);

end%}}}


%%%%%%%%%%%%%%%%% STEP TWELVE: KNMAX MAX BUTTRESSING with MASK APPLIED


if perform(org, 'KNMAX_no_gap_fill_MASKED'),% {{{

    addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	for yy = 1:length(years)
		
        md.mask.ice_levelset(:) = 1;   
        md.mask.ocean_levelset(:) = 1;

        indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory

        in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1); % makes values within the .exp file = 1        
        md.mask.ocean_levelset(find(in)) = -1; % sets ocean_levelset to -1 for nodes inside the ice extent
		md.mask.ice_levelset(find(in)) = -1; % sets ice_levelset to -1 for nodes inside the ice extent		

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		vx_in = double(ncread(fname,'VX'));
		vy_in = double(ncread(fname,'VY'));
		thickness_in = double(ncread(fname,'Thickness'));

		% Get the vx, vy onto the ISSM grid
		vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		md.inversion.vx_obs(find(in)) = vx(find(in));
		md.inversion.vy_obs(find(in)) = vy(find(in));
		md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);

		% Get the thickness onto the ISSM grid
		thickness = InterpFromGrid(x_in,y_in,thickness_in',md.mesh.x,md.mesh.y);
		md.geometry.thickness(find(in)) = thickness(find(in));

		% calculate the buttressing numbers
		[Knflow{yy},Nflow{yy},Knmax{yy},Nmax{yy},N0{yy}]=buttressingnumber(md);

		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;


        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        plotmodel(md,...
            'data',Knmax{yy},...
            'figure', yy,...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'colormap#1','redblue',...
            'caxis',[-5 5],...
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        colormap(redblue(64));
        cb = colorbar;
        ylabel(cb, 'Maximum Buttressing Number', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['Knmax ' num2str(years(yy))], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        exportgraphics(gcf,['Plots/KNMAX_no_gap_fill_' num2str(years(yy)) '_masked.png'])


        %%%%%%%%%%%%%%%%%%%%%%%% EXPORTING AS A GEOTIFF:

        % grid resol in metres
        dx = 100;
        dy = 100;
        
        xrange = x(1):dx:x(2);
        yrange = y(1):dy:y(2);
        
        [Xgrid, Ygrid] = meshgrid(xrange, yrange);
        
        % Interpolate Knmax onto the regular grid
        F = scatteredInterpolant(md.mesh.x, md.mesh.y, Knmax{yy}, 'linear', 'none');
        Zgrid = F(Xgrid, Ygrid);
        
        % Define referencing object using limits, not full vectors
        R = maprefcells([min(xrange) max(xrange)], [min(yrange) max(yrange)], size(Zgrid));
        
        % Write GeoTIFF
        geotiffwrite(['Data/KNMAX_' num2str(years(yy)) '_unmasked.tif'], Zgrid, R, ...
            'CoordRefSysCode', 'EPSG:3031');



	end

	%if want to save all output
	save([indir 'KNMAX_no_gap_fill_MASKED_all_years.mat'], 'Knflow', 'Nflow', 'Knmax', 'Nmax','masks');

end%}}}




%%%%%%%%%%%%%%%%% STEP THIRTEEN: PLOTTING ALL masked KNMAX TOGETHER

% KNMAX = MAXIMUM BUTTRESSING

if perform(org, 'KNMAX_all_years_masked'),% {{{

	files = {};
for yy = 1:length(years)
    files{end+1} = ['Plots/KNMAX_no_gap_fill_' num2str(years(yy)) '_masked.png'];
end

% Create montage
figure;
montage(files, 'Size', [3 4], 'BackgroundColor', 'w'); % 3 rows, 4 cols — adjust as needed
title('KNMAX all years');
exportgraphics(gcf, 'Plots/KNMAX_no_gap_fill_masked.png', 'Resolution', 300);

end%}}}

%%%%%%%%%%%%%%%%% STEP FOURTEEN: KNMAX change compared to 2006


if perform(org, 'KNMAX_minus_2006'),% {{{

    % Define years to process
    years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];

    % Load model
    md = loadmodel(org, 'Inversion_BHO');

    % Initialize containers
    Knflow = {};
    Nflow = {};
    Knmax = {};
    Nmax = {};
    N0 = {};
    delta_KNMAX = {};
    masks = {};

    % ---- Load and compute reference KNMAX for 2006 ----
    refyear = 2006;
    indir_ref = ['Data/' num2str(refyear) '/'];
    in_ref = ContourToNodes(md.mesh.x, md.mesh.y, [indir_ref 'ice_extent_' num2str(refyear) '.exp'], 1);

    % Reset and apply 2006 mask
    md.mask.ocean_levelset(:) = 1;
    md.mask.ice_levelset(:) = 1;
    md.mask.ocean_levelset(find(in_ref)) = -1;
    md.mask.ice_levelset(find(in_ref)) = -1;

    % Interpolate thickness and velocity fields for 2006
    fname_ref = [indir_ref num2str(refyear) '.nc'];
    x_in = double(ncread(fname_ref, 'x'));
    y_in = double(ncread(fname_ref, 'y'));
    vx_in = double(ncread(fname_ref, 'VX'));
    vy_in = double(ncread(fname_ref, 'VY'));
    thickness_in = double(ncread(fname_ref, 'Thickness'));

    vx_ref = InterpFromGrid(x_in, y_in, vx_in', md.mesh.x, md.mesh.y);
    vy_ref = InterpFromGrid(x_in, y_in, vy_in', md.mesh.x, md.mesh.y);
    thickness_ref = InterpFromGrid(x_in, y_in, thickness_in', md.mesh.x, md.mesh.y);

    md.inversion.vx_obs(find(in_ref)) = vx_ref(find(in_ref));
    md.inversion.vy_obs(find(in_ref)) = vy_ref(find(in_ref));
    md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
    md.geometry.thickness(find(in_ref)) = thickness_ref(find(in_ref));

    % Compute KNMAX for 2006
    [~, ~, Knmax_ref, ~, ~] = buttressingnumber(md);

    % ---- Loop through all other years ----
    for yy = 1:length(years)
        thisyear = years(yy);
        indir = ['Data/' num2str(thisyear) '/'];

        % Reset mask
        md.mask.ocean_levelset(:) = 1;
        md.mask.ice_levelset(:) = 1;

        % Apply ice extent mask
        in = ContourToNodes(md.mesh.x, md.mesh.y, [indir 'ice_extent_' num2str(thisyear) '.exp'], 1);
        md.mask.ocean_levelset(find(in)) = -1;
        md.mask.ice_levelset(find(in)) = -1;

        % Interpolate velocity and thickness
        fname = [indir num2str(thisyear) '.nc'];
        x_in = double(ncread(fname, 'x'));
        y_in = double(ncread(fname, 'y'));
        vx_in = double(ncread(fname, 'VX'));
        vy_in = double(ncread(fname, 'VY'));
        thickness_in = double(ncread(fname, 'Thickness'));

        vx = InterpFromGrid(x_in, y_in, vx_in', md.mesh.x, md.mesh.y);
        vy = InterpFromGrid(x_in, y_in, vy_in', md.mesh.x, md.mesh.y);
        thickness = InterpFromGrid(x_in, y_in, thickness_in', md.mesh.x, md.mesh.y);

        md.inversion.vx_obs(find(in)) = vx(find(in));
        md.inversion.vy_obs(find(in)) = vy(find(in));
        md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
        md.geometry.thickness(find(in)) = thickness(find(in));

        % Compute buttressing numbers
        [Knflow{yy}, Nflow{yy}, Knmax{yy}, Nmax{yy}, N0{yy}] = buttressingnumber(md);

        % Difference from 2006
        delta_KNMAX{yy} = Knmax{yy} - Knmax_ref;

        % Define plot mask
        mask = ones(md.mesh.numberofvertices, 1);
        mask(md.mask.ice_levelset > 0) = 0;
        masks{yy} = mask;

        % Plot
        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);

        plotmodel(md, ...
            'data', delta_KNMAX{yy}, ...
            'figure', yy, ...
            'colormap#1','redblue',...
            'xlim#all', [998648.2458, 1171935.9538], ...
            'ylim#all', [-2148431.9935, -2033743.9935], ...
            'caxis', [-3 3], ...
            'mask#all', mask, ...
            'expdisp', [indir 'ice_extent_' num2str(thisyear) '.exp']);

        colormap(redblue(64));
        cb = colorbar;
        ylabel(cb, 'ΔKN_{max}', 'FontName', 'Times New Roman', 'FontSize', 12);
        xlabel('Eastings (m)', 'FontName', 'Times New Roman', 'FontSize', 12);
        ylabel('Northings (m)', 'FontName', 'Times New Roman', 'FontSize', 12);
        title(['KN_{max} Change: ' num2str(thisyear) ' – 2006'], 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontSize', 12);
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

        exportgraphics(gcf, ['Plots/KNMAX_change_' num2str(thisyear) '_minus_2006.png'])
    end

    % Optional: save results
    save('Data/KNMAX_change_minus_2006_alloutfiles.mat', 'Knmax', 'delta_KNMAX', 'masks');

end % }}}

%%%%%%%%%%%%%%%%% STEP FIFTEEN: PLOTTING ALL masked KNMAX change minus 2006 TOGETHER

% KNMAX = MAXIMUM BUTTRESSING

if perform(org, 'KNMAX_all_years_masked_minus_2006'),% {{{

	files = {};
for yy = 1:length(years)
    files{end+1} = ['Plots/KNMAX_change_' num2str(years(yy)) '_minus_2006.png'];
end

% Create montage
figure;
montage(files, 'Size', [3 4], 'BackgroundColor', 'w'); % 3 rows, 4 cols — adjust as needed
title('KNMAX minus 2006');
exportgraphics(gcf, 'Plots/KNMAX_no_gap_fill_masked_MINUS_2006.png', 'Resolution', 300);

end%}}}



%%%%%%%%%%%%%%%%% STEP SIXTEEN: KNMAX MEAN


if perform(org, 'KNMAX_MEAN'),% {{{


    years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];

    md = loadmodel(org, 'Inversion_BHO');
    velo = md.inversion.vel_obs;

    Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	

    n_nodes = md.mesh.numberofvertices;
    n_years = length(years);
    knmax_all_years = NaN(n_nodes, n_years);

    % Calculating the mean across all years

    for yy = 1:length(years)
		    
            md.mask.ice_levelset(:) = 1;   
            md.mask.ocean_levelset(:) = 1;


            indir = ['Data/' num2str(years(yy)) '/']; 
    
		    in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1); % makes values within the .exp file = 1
		    md.mask.ocean_levelset(find(in)) = -1; % sets ocean_levelset to -1 for nodes inside the ice extent
		    md.mask.ice_levelset(find(in)) = -1; % sets ice_levelset to -1 for nodes inside the ice extent
    
		    % READ THE ITS_LIVE DATA ONTO ISSM MESH
		    fname = [indir num2str(years(yy)) '.nc'];
		    x_in = double(ncread(fname,'x'));
		    y_in = double(ncread(fname,'y'));
		    thickness_in = double(ncread(fname,'Thickness'));
    
		    % Get the thickness onto the ISSM grid
		    thickness = InterpFromGrid(x_in,y_in,thickness_in',md.mesh.x,md.mesh.y);
		    md.geometry.thickness(find(in)) = thickness(find(in));

            [~, ~, Knmax, ~, ~] = buttressingnumber(md);

		    % mask= ones(md.mesh.numberofvertices,1);
		    % mask(find(md.mask.ice_levelset>0)) = 0;
		    % masks{yy} = mask;

            mask = ones(md.mesh.numberofvertices,1);
            mask(md.mask.ice_levelset > 0) = 0;
            masks{yy} = mask;
            
            Knmax(mask == 0) = NaN; % removes knmax values outside the ice extent
            knmax_all_years(:, yy) = Knmax;

    end

    knmax_mean = mean(knmax_all_years, 2, 'omitnan');

    plotmodel(md,...
        'data', knmax_mean,...
        'figure', yy,...
        'colormap#1','redblue',...
        'xlim#all',[998648.2458, 1171935.9538],...
        'ylim#all', [-2148431.9935, -2033743.9935],...
        'caxis', [-5 5],...  % adjust ?
        'mask#all', mask,...
        'expdisp', [indir 'ice_extent_' num2str(years(yy)) '.exp']);

    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'KNMAX', 'FontName', 'Times New Roman', 'FontSize', 12);

    figure(yy);
    xlabel('Eastings (m)');
    ylabel('Northings (m)');
    title(['KNMAX Mean'], 'FontWeight', 'normal');

    box off
    x = xlim;
    y = ylim;
    
    hold on
    plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
    plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
    plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
    plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
 
    title('Mean KNMAX');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    exportgraphics(gcf, 'Plots/MEAN_KNMAX.png', 'Resolution', 300);



    %%%%%%%%%%%%%%%%%%%%%%%% EXPORTING AS A GEOTIFF:

    % grid resol in metres
    dx = 100;
    dy = 100;
    
    xrange = x(1):dx:x(2);
    yrange = y(1):dy:y(2);
    
    [Xgrid, Ygrid] = meshgrid(xrange, yrange);
    
    % Interpolate Knmax mean onto the regular grid
    F = scatteredInterpolant(md.mesh.x, md.mesh.y, knmax_mean, 'linear', 'none');
    Zgrid = F(Xgrid, Ygrid);
    
    % Define referencing object using limits, not full vectors
    R = maprefcells([min(xrange) max(xrange)], [min(yrange) max(yrange)], size(Zgrid));
    
    % Write GeoTIFF
    geotiffwrite(['Data/KNMAX_mean_masked.tif'], Zgrid, R, ...
        'CoordRefSysCode', 'EPSG:3031');




	save([indir 'KNMAX_MEAN.mat'], 'knmax_mean','masks');


 end%}}}



%%%%%%%%%%%%%%%%% STEP SEVENTEEN: KNMAX ANOMALY 
        %% STEPS 12 AND 16 MUST BE COMPLETED FIRST

if perform(org, 'KNMAX_Anomaly'),% {{{

    load('Data/KNMAX_MEAN.mat', 'knmax_mean', 'masks'); % from step 16

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

    min_vals_knmax = NaN(length(years),1); % made to store min and max anom
    max_vals_knmax = NaN(length(years),1);

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	for yy = 1:length(years)
		
        md.mask.ice_levelset(:) = 1;   
        md.mask.ocean_levelset(:) = 1; 
        
        indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory

        in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1); % makes values within the .exp file = 1        
        md.mask.ocean_levelset(find(in)) = -1; % sets ocean_levelset to -1 for nodes inside the ice extent
		md.mask.ice_levelset(find(in)) = -1; % sets ice_levelset to -1 for nodes inside the ice extent

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		thickness_in = double(ncread(fname,'Thickness'));

		% Get the thickness onto the ISSM grid
		thickness = InterpFromGrid(x_in,y_in,thickness_in',md.mesh.x,md.mesh.y);
		md.geometry.thickness(find(in)) = thickness(find(in));

        [~, ~, Knmax, ~, ~] = buttressingnumber(md);

        mask = ones(md.mesh.numberofvertices,1);
        mask(md.mask.ice_levelset > 0) = 0;
        masks{yy} = mask;
        Knmax(mask == 0) = NaN; % removes thickness values outside the ice extent

		
        % Anomaly Calculation
        anomaly_Knmax = Knmax - knmax_mean;

        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');

        plotmodel(md,...
            'data',anomaly_Knmax,...
            'figure', yy,...
            'colormap#1','parula',...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'caxis',[-0.01 0.01],...       %'caxis',[-30 15],... 
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        colormap(redblue(64));
        cb = colorbar;
        ylabel(cb, 'KNMAX Anomaly', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['KNMAX Anomaly ' num2str(years(yy))], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        
        exportgraphics(gcf,['Plots/KNMAX_' num2str(years(yy)) '_ANOMALY_redblue.png'])

        min_vals_knmax(yy) = min(anomaly_Knmax, [], 'omitnan');
        max_vals_knmax(yy) = max(anomaly_Knmax, [], 'omitnan');

	end

	%if want to save all output
	save([indir 'KNMAX_ANOMALIES.mat'], 'anomaly_Knmax','Knmax','masks');

end%}}}



%%%%%%%%%%%%%%%%% STEP EIGEHTEEN: PLOTTING all knmax anomalies together

% KNMAX = MAXIMUM BUTTRESSING

if perform(org, 'KNMAX_anomalies_all_years'),% {{{

	files = {};
for yy = 1:length(years)
    files{end+1} = ['Plots/KNMAX_' num2str(years(yy)) '_ANOMALY_redblue.png'];
end

% Create montage
figure;
montage(files, 'Size', [3 4], 'BackgroundColor', 'w'); % 3 rows, 4 cols — adjust as needed
title('KNMAX Anomalies');
exportgraphics(gcf, 'Plots/KNMAX_anomalies_all_years.png', 'Resolution', 300);

end%}}}



%%%%%%%%%%%%%%%%% STEP NINETEEN: Annual Velocities


if perform(org, 'Velocity_annual'),% {{{

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;
    
    annual_velocity = {};

	masks = {};

	for yy = 1:length(years)
		
        md.mask.ice_levelset(:) = 1;   
        md.mask.ocean_levelset(:) = 1;

        indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory

        in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1); % makes values within the .exp file = 1        
        md.mask.ocean_levelset(find(in)) = -1; % sets ocean_levelset to -1 for nodes inside the ice extent
		md.mask.ice_levelset(find(in)) = -1; % sets ice_levelset to -1 for nodes inside the ice extent		

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		vx_in = double(ncread(fname,'VX'));
		vy_in = double(ncread(fname,'VY'));

		% Get the vx, vy onto the ISSM grid
		vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		md.inversion.vx_obs(find(in)) = vx(find(in));
		md.inversion.vy_obs(find(in)) = vy(find(in));
		md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
        vel = md.inversion.vel_obs


		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
        vel(mask == 0) = NaN;
        masks{yy} = mask;
        annual_velocity{yy} = vel;


        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        plotmodel(md,...
            'data',vel,...
            'figure', yy,...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'colormap#1','redblue',...
            'caxis',[0 1000],...
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        colormap(redblue(64));
        cb = colorbar;
        ylabel(cb, 'Velocity (m/yr)', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['Velocity ' num2str(years(yy))], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        exportgraphics(gcf,['Plots/Velocity_' num2str(years(yy)) '_masked.png'])
	end

	%if want to save all output
    save('Data/surface_velocity_all_years.mat', 'annual_velocity', 'masks');

end%}}}



%%%%%%%%%%%%%%%%% STEP TWENTY: Velocity CHANGE ( MINUS 2006 )


if perform(org, 'Velocity_delta'),% {{{

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
	refyear = 2006;

    md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;
    
    annual_velocity = {};
    delta_velocity = {};

	masks = {};



    indir_ref = ['Data/' num2str(refyear) '/'];
    in_ref = ContourToNodes(md.mesh.x, md.mesh.y, [indir_ref 'ice_extent_' num2str(refyear) '.exp'], 1);

    md.mask.ice_levelset(:) = 1;
    md.mask.ocean_levelset(:) = 1;
    md.mask.ice_levelset(find(in_ref)) = -1;
    md.mask.ocean_levelset(find(in_ref)) = -1;

    fname_ref = [indir_ref num2str(refyear) '.nc'];
    x_in = double(ncread(fname_ref, 'x'));
    y_in = double(ncread(fname_ref, 'y'));
    vx_in = double(ncread(fname_ref, 'VX'));
    vy_in = double(ncread(fname_ref, 'VY'));

    vx_ref = InterpFromGrid(x_in, y_in, vx_in', md.mesh.x, md.mesh.y);
    vy_ref = InterpFromGrid(x_in, y_in, vy_in', md.mesh.x, md.mesh.y);

    % md.inversion.vx_obs(in_ref) = vx_ref(in_ref);
    % md.inversion.vy_obs(in_ref) = vy_ref(in_ref);
    % md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
    % vel_ref = md.inversion.vel_obs;


    md.inversion.vx_obs(find(in_ref)) = vx_ref(find(in_ref));
    md.inversion.vy_obs(find(in_ref)) = vy_ref(find(in_ref));
    md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
    vel_ref = md.inversion.vel_obs;

	for yy = 1:length(years)
		
        md.mask.ice_levelset(:) = 1;   
        md.mask.ocean_levelset(:) = 1;

        indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory

        in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1); % makes values within the .exp file = 1        
        md.mask.ocean_levelset(find(in)) = -1; % sets ocean_levelset to -1 for nodes inside the ice extent
		md.mask.ice_levelset(find(in)) = -1; % sets ice_levelset to -1 for nodes inside the ice extent		

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		vx_in = double(ncread(fname,'VX'));
		vy_in = double(ncread(fname,'VY'));

		% Get the vx, vy onto the ISSM grid
		vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		md.inversion.vx_obs(find(in)) = vx(find(in));
		md.inversion.vy_obs(find(in)) = vy(find(in));
		md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
        vel = md.inversion.vel_obs


        delta_vel = vel - vel_ref


		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
        delta_vel(mask == 0) = NaN;
        masks{yy} = mask;
        annual_velocity{yy} = vel;
        delta_velocity{yy} = delta_vel;


        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        plotmodel(md,...
            'data',delta_vel,...
            'figure', yy,...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'colormap#1','redblue',...
            'caxis',[-150 150],...
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        colormap(redblue(64));
        cb = colorbar;
        ylabel(cb, 'Delta Velocity (m/yr)', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['Delta Velocity ' num2str(years(yy)) ' - 2006'], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        exportgraphics(gcf,['Plots/Velocity_delta_' num2str(years(yy)) '_minus_2006.png'])
	end

	%if want to save all output
    save('Data/delta_velocity.mat', 'delta_velocity', 'masks');

end%}}}




%%%%%%%%%%%%%%%%% STEP TWENTY ONE: PLOTTING all velocity changes together

% minus 2006

if perform(org, 'velocity_changes_all_years'),% {{{

	files = {};
for yy = 1:length(years)
    files{end+1} = ['Plots/Velocity_delta_' num2str(years(yy)) '_minus_2006.png'];
end

% Create montage
figure;
montage(files, 'Size', [3 4], 'BackgroundColor', 'w'); % 3 rows, 4 cols — adjust as needed
title('Velocity change since 2006');
exportgraphics(gcf, 'Plots/velocity_delta_all_years.png', 'Resolution', 300);

end%}}}




%%%%%%%%%%%%%%%%% STEP TWENTY TWO: Velocity MEAN


if perform(org, 'Velocity_MEAN'),% {{{


    years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];

    md = loadmodel(org, 'Inversion_BHO');

    Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	

    n_nodes = md.mesh.numberofvertices;
    n_years = length(years);
    vel_all_years = NaN(n_nodes, n_years);
    
    % Calculating the mean across all years

    for yy = 1:length(years)
		    
            md.mask.ice_levelset(:) = 1;   
            md.mask.ocean_levelset(:) = 1;


            indir = ['Data/' num2str(years(yy)) '/']; 
    
		    in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1); % makes values within the .exp file = 1
		    md.mask.ocean_levelset(find(in)) = -1; % sets ocean_levelset to -1 for nodes inside the ice extent
		    md.mask.ice_levelset(find(in)) = -1; % sets ice_levelset to -1 for nodes inside the ice extent
    
		    % READ THE ITS_LIVE DATA ONTO ISSM MESH
		    fname = [indir num2str(years(yy)) '.nc'];
		    x_in = double(ncread(fname,'x'));
		    y_in = double(ncread(fname,'y'));
		    thickness_in = double(ncread(fname,'Thickness'));
    

            fname = [indir num2str(years(yy)) '.nc'];
		    x_in = double(ncread(fname,'x'));
		    y_in = double(ncread(fname,'y'));
		    vx_in = double(ncread(fname,'VX'));
		    vy_in = double(ncread(fname,'VY'));
    
		    vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		    vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		    md.inversion.vx_obs(find(in)) = vx(find(in));
		    md.inversion.vy_obs(find(in)) = vy(find(in));
		    md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
            vel = md.inversion.vel_obs


            mask = ones(md.mesh.numberofvertices,1);
            mask(md.mask.ice_levelset > 0) = 0;
            masks{yy} = mask;

            vel(mask == 0) = NaN;
            vel_all_years(:, yy) = vel;


    end

    vel_mean = mean(vel_all_years, 2, 'omitnan');

    plotmodel(md,...
        'data', vel_mean,...
        'figure', yy,...
        'colormap#1','redblue',...
        'xlim#all',[998648.2458, 1171935.9538],...
        'ylim#all', [-2148431.9935, -2033743.9935],...
        'caxis', [0 1000],...  % adjust ?
        'mask#all', mask,...
        'expdisp', [indir 'ice_extent_' num2str(years(yy)) '.exp']);

    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Velocity', 'FontName', 'Times New Roman', 'FontSize', 12);

    figure(yy);
    xlabel('Eastings (m)');
    ylabel('Northings (m)');
    title(['Velocity Mean'], 'FontWeight', 'normal');

    box off
    x = xlim;
    y = ylim;
    
    hold on
    plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
    plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
    plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
    plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
 
    title('Mean Velocity');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    exportgraphics(gcf, 'Plots/MEAN_velocity.png', 'Resolution', 300);



    %%%%%%%%%%%%%%%%%%%%%%%% EXPORTING AS A GEOTIFF:

    % grid resol in metres
    dx = 100;
    dy = 100;
    
    xrange = x(1):dx:x(2);
    yrange = y(1):dy:y(2);
    
    [Xgrid, Ygrid] = meshgrid(xrange, yrange);
    
    % Interpolate velocity mean onto the regular grid
    F = scatteredInterpolant(md.mesh.x, md.mesh.y, vel_mean, 'linear', 'none');
    Zgrid = F(Xgrid, Ygrid);
    
    % Define referencing object using limits, not full vectors
    R = maprefcells([min(xrange) max(xrange)], [min(yrange) max(yrange)], size(Zgrid));
    
    % Write GeoTIFF
    geotiffwrite(['Data/VEL_mean_masked.tif'], Zgrid, R, ...
        'CoordRefSysCode', 'EPSG:3031');






	save([indir 'Velocity_MEAN.mat'], 'vel_mean','masks');


 end%}}}




%%%%%%%%%%%%%%%%% STEP TWENTY THREE: VELOCITY ANOMALY 
        %% STEPS 19 AND 22 MUST BE COMPLETED FIRST

if perform(org, 'Vel_Anomaly'),% {{{

    load('Data/Velocity_MEAN.mat', 'vel_mean', 'masks');

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

    min_vals_vel = NaN(length(years),1); % made to store min and max anom
    max_vals_vel = NaN(length(years),1);

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};
	for yy = 1:length(years)
		
        md.mask.ice_levelset(:) = 1;   
        md.mask.ocean_levelset(:) = 1; 
        
        indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory

        in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1); % makes values within the .exp file = 1        
        md.mask.ocean_levelset(find(in)) = -1; % sets ocean_levelset to -1 for nodes inside the ice extent
		md.mask.ice_levelset(find(in)) = -1; % sets ice_levelset to -1 for nodes inside the ice extent

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		vx_in = double(ncread(fname,'VX'));
		vy_in = double(ncread(fname,'VY'));

		% Get the vx, vy onto the ISSM grid
		vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		md.inversion.vx_obs(find(in)) = vx(find(in));
		md.inversion.vy_obs(find(in)) = vy(find(in));
		md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
        vel = md.inversion.vel_obs

        mask = ones(md.mesh.numberofvertices,1);
        mask(md.mask.ice_levelset > 0) = 0;
        masks{yy} = mask;
        vel(mask == 0) = NaN; 

		
        % Anomaly Calculation
        anomaly_vel = vel - vel_mean;

        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');

        plotmodel(md,...
            'data',anomaly_vel,...
            'figure', yy,...
            'colormap#1','parula',...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'caxis',[-30 30],...       %'caxis',[-30 15],... 
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        colormap(redblue(64));
        cb = colorbar;
        ylabel(cb, 'Velocity Anomaly (m/yr)', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['Velocity Anomaly ' num2str(years(yy))], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        
        exportgraphics(gcf,['Plots/velocity_' num2str(years(yy)) '_ANOMALY.png'])

        min_vals_vel(yy) = min(anomaly_vel, [], 'omitnan');
        max_vals_vel(yy) = max(anomaly_vel, [], 'omitnan');

	end

	%if want to save all output
	save([indir 'velocity_anomalies_all.mat'], 'anomaly_vel','vel','vel_mean','masks','min_vals_vel', 'max_vals_vel');

end%}}}




%%%%%%%%%%%%%%%%% STEP TWENTY FOUR: PLOTTING all velocity anomalies together


if perform(org, 'velocity_anomalies_all_years'),% {{{

	files = {};
for yy = 1:length(years)
    files{end+1} = ['Plots/velocity_' num2str(years(yy)) '_ANOMALY.png'];
end

% Create montage
figure;
montage(files, 'Size', [3 4], 'BackgroundColor', 'w'); % 3 rows, 4 cols — adjust as needed
title('Velocity Anomalies');
exportgraphics(gcf, 'Plots/velocity_anomalies_all_years.png', 'Resolution', 300);

end%}}}


%%%%%%%%%%%%%%%%% STEP TWENTY FIVE: thomas params longitudinal strain rate


if perform(org, 'longitudinal_thomasparams'),% {{{

    years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
    md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};

	for yy = 1:length(years)
		
        indir = ['Data/' num2str(years(yy)) '/']; 
        in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1);
	    md.mask.ocean_levelset(find(in)) = -1;
	    md.mask.ice_levelset(find(in)) = -1;
    
	    % READ THE ITS_LIVE DATA ONTO ISSM MESH
	    fname = [indir num2str(years(yy)) '.nc'];
	    x_in = double(ncread(fname,'x'));
	    y_in = double(ncread(fname,'y'));
	    vx_in = double(ncread(fname,'VX'));
	    vy_in = double(ncread(fname,'VY'));
	    thickness_in = double(ncread(fname,'Thickness'));
    
	    % Get the vx, vy onto the ISSM grid
	    vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
	    vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
	    md.inversion.vx_obs(find(in)) = vx(find(in));
	    md.inversion.vy_obs(find(in)) = vy(find(in));
	    md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
    
	    % Get the thickness onto the ISSM grid
	    thickness = InterpFromGrid(x_in,y_in,thickness_in',md.mesh.x,md.mesh.y);
	    md.geometry.thickness(find(in)) = thickness(find(in));
    

       	mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        md = mechanicalproperties(md, vx, vy);
        [~,~,~,exx,~] = thomasparams(md,'eq','Thomas','smoothing',0,'coordsys','longitudinal');
        
        sec_per_year = 365 * 24 * 60 * 60;
        exx_per_year = exx * sec_per_year
        exx_per_year(mask == 0) = NaN;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');

        plotmodel(md,...
            'data',exx_per_year,...
            'figure', yy,...
            'colormap#1','parula',...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'caxis',[-0.02 0.02],...    
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        colormap(redblue(64));
        cb = colorbar;
        ylabel(cb, 'exx (longitudinal) (a-1)', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['exx ' num2str(years(yy))], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        
        exportgraphics(gcf,['Plots/exx_' num2str(years(yy)) '_longitudinal_thomasparams.png'])

    end



end%}}}



%%%%%%%%%%%%%%%%% STEP TWENTY SIX: PLOTTING all exx together


if perform(org, 'exx_together'),% {{{

	files = {};
for yy = 1:length(years)
    files{end+1} = ['Plots/exx_' num2str(years(yy)) '_longitudinal_thomasparams.png'];
end

% Create montage
figure;
montage(files, 'Size', [3 4], 'BackgroundColor', 'w'); % 3 rows, 4 cols — adjust as needed
title('exx longitudinal');
exportgraphics(gcf, 'Plots/exx_longitudinal_all_years.png', 'Resolution', 300);

end%}}}




%%%%%%%%%%%%%%%%% STEP TWENTY SEVEN: GL FLUX test 1
    % THIS CODE GIVES 1.3 Gt/yr FOR EVERY YEAR


if perform(org, 'GL_Flux_test_1'),% {{{


% Add required path for outflux()
addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\ISSM-Windows-MATLAB\ISSM-Windows-MATLAB\src\m\contrib\morlighem\massbalance')

% Define the years you're working with
years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];

% Load base model
md = loadmodel(org, 'Inversion_BHO');

% Initialize flux storage
fluxes = zeros(length(years), 1);

for yy = 1:length(years)

    % Reset mask to floating ice
    md.mask.ice_levelset(:) = 1;   
    md.mask.ocean_levelset(:) = 1;

    % Set input directory for the current year
    indir = ['Data/' num2str(years(yy)) '/'];

    % Read and apply the .exp ice extent mask for current year
    in = ContourToNodes(md.mesh.x, md.mesh.y, [indir 'ice_extent_' num2str(years(yy)) '.exp'], 1);
    md.mask.ocean_levelset(find(in)) = -1; 
	md.mask.ice_levelset(find(in)) = -1;

    % Read ITS_LIVE NetCDF velocity and thickness data
    fname = [indir num2str(years(yy)) '.nc'];
    x_in = double(ncread(fname, 'x'));
    y_in = double(ncread(fname, 'y'));
    vx_in = double(ncread(fname, 'VX'));
    vy_in = double(ncread(fname, 'VY'));
    thickness_in = double(ncread(fname, 'Thickness'));

    % Interpolate velocity onto ISSM mesh
    vx = InterpFromGrid(x_in, y_in, vx_in', md.mesh.x, md.mesh.y);
    vy = InterpFromGrid(x_in, y_in, vy_in', md.mesh.x, md.mesh.y);
    md.inversion.vx_obs(find(in)) = vx(find(in));
	md.inversion.vy_obs(find(in)) = vy(find(in));
    md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);

    % Interpolate thickness onto ISSM mesh
    thickness = InterpFromGrid(x_in, y_in, thickness_in', md.mesh.x, md.mesh.y);
    md.geometry.thickness(find(in)) = thickness(find(in));

    % Compute grounding line flux (Gt/yr)
    fluxes(yy) = outflux(md);
    fprintf('Year %d: Flux = %.10f Gt/yr\n', years(yy), fluxes(yy));

    fprintf('Year %d: vx range = [%.2f, %.2f], thickness range = [%.2f, %.2f]\n', ...
    years(yy), min(vx(find(in))), max(vx(find(in))), min(thickness(find(in))), max(thickness(find(in))));


end

% Save results to .mat file
save('Data/GL_flux_timeseries.mat', 'years', 'fluxes');


% addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\ISSM-Windows-MATLAB\ISSM-Windows-MATLAB\src\m\contrib\morlighem\massbalance')
% flux=outflux(md);


end%}}}



%%%%%%%%%%%%%%%%% STEP TWENTY EIGHT: GL FLUX test 2


if perform(org, 'GL_Flux_test_2'),% {{{

    % Add required path for outflux()
addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\ISSM-Windows-MATLAB\ISSM-Windows-MATLAB\src\m\contrib\morlighem\massbalance')

% Define the years you're working with
years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];

% Load base model
md = loadmodel(org, 'Inversion_BHO');

% Initialize flux storage
fluxes = zeros(length(years), 1);

for yy = 1:length(years)

    % Reset mask to floating ice
    md.mask.ice_levelset(:) = 1;   
    md.mask.ocean_levelset(:) = 1;

    % Set input directory for the current year
    indir = ['Data/' num2str(years(yy)) '/'];

    % Read and apply the .exp ice extent mask for current year
    in = ContourToNodes(md.mesh.x, md.mesh.y, [indir 'ice_extent_' num2str(years(yy)) '.exp'], 1);
    md.mask.ocean_levelset(find(in)) = -1; 
    md.mask.ice_levelset(find(in)) = -1;

    % Read ITS_LIVE NetCDF velocity and thickness data
    fname = [indir num2str(years(yy)) '.nc'];
    x_in = double(ncread(fname, 'x'));
    y_in = double(ncread(fname, 'y'));
    vx_in = double(ncread(fname, 'VX'));
    vy_in = double(ncread(fname, 'VY'));
    thickness_in = double(ncread(fname, 'Thickness'));

    % Interpolate full domain
    vx = InterpFromGrid(x_in, y_in, vx_in', md.mesh.x, md.mesh.y);
    vy = InterpFromGrid(x_in, y_in, vy_in', md.mesh.x, md.mesh.y);
    thickness = InterpFromGrid(x_in, y_in, thickness_in', md.mesh.x, md.mesh.y);

    % Replace NaNs
    vx(isnan(vx)) = 0;
    vy(isnan(vy)) = 0;
    thickness(isnan(thickness)) = 1;  % or use min(valid thickness) if needed
    

    % Assign to model
    md.inversion.vx_obs = vx;
    md.inversion.vy_obs = vy;
    md.inversion.vel_obs = sqrt(vx.^2 + vy.^2);
    md.geometry.thickness = thickness;

    % Ensure ice is unmasked at the gate (optional but recommended)
    gate = logical(ContourToNodes(md.mesh.x, md.mesh.y, 'Data/GL_Shapefile/GL_2003_UNION.exp', 1));
    md.mask.ice_levelset(gate) = -1;
    md.mask.ocean_levelset(gate) = 1;

    % Compute flux across the gate
    fluxes(yy) = outfluxgate(md, 'Data/GL_Shapefile/GL_2003_UNION.exp');
    fprintf('Year %d: Flux across gate = %.2f Gt/yr\n', years(yy), fluxes(yy));

    % Debugging: check ranges
    fprintf('Year %d: vx range = [%.2f, %.2f], thickness range = [%.2f, %.2f]\n', ...
        years(yy), min(vx(find(in))), max(vx(find(in))), ...
        min(thickness(find(in))), max(thickness(find(in))));
end


% Save results to .mat file
save('Data/GL_flux_timeseries.mat', 'years', 'fluxes');


% addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\ISSM-Windows-MATLAB\ISSM-Windows-MATLAB\src\m\contrib\morlighem\massbalance')
% flux=outflux(md);


end%}}}


%%%%%%%%%%%%%%%%% STEP TWENTY NINE: GL FLUX test 3
    % THIS WORKS!!!


if perform(org, 'GL_Flux_test_3'),% {{{

    % Add required path for outflux()
addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\ISSM-Windows-MATLAB\ISSM-Windows-MATLAB\src\m\contrib\morlighem\massbalance')

% Define the years you're working with
years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];

% Load base model
md = loadmodel(org, 'Inversion_BHO');

% Initialize flux storage
fluxes = zeros(length(years), 1);

for yy = 1:length(years)

    % Reset mask to floating ice
    md.mask.ice_levelset(:) = 1;   
    md.mask.ocean_levelset(:) = 1;

    % Set input directory for the current year
    indir = ['Data/' num2str(years(yy)) '/'];

    % Read and apply the .exp ice extent mask for current year
    in = ContourToNodes(md.mesh.x, md.mesh.y, [indir 'ice_extent_' num2str(years(yy)) '.exp'], 1);
    md.mask.ocean_levelset(find(in)) = -1; 
	md.mask.ice_levelset(find(in)) = -1;

    % Read ITS_LIVE NetCDF velocity and thickness data
    fname = [indir num2str(years(yy)) '.nc'];
    x_in = double(ncread(fname, 'x'));
    y_in = double(ncread(fname, 'y'));
    vx_in = double(ncread(fname, 'VX'));
    vy_in = double(ncread(fname, 'VY'));
    thickness_in = double(ncread(fname, 'Thickness'));

    % Interpolate full domain
    vx = InterpFromGrid(x_in, y_in, vx_in', md.mesh.x, md.mesh.y);
    vy = InterpFromGrid(x_in, y_in, vy_in', md.mesh.x, md.mesh.y);
    thickness = InterpFromGrid(x_in, y_in, thickness_in', md.mesh.x, md.mesh.y);
    
    % Replace NaNs — use zeros or another default value
    vx(isnan(vx)) = 0;
    vy(isnan(vy)) = 0;
    thickness(isnan(thickness)) = 600; % this is ~~ average GL border thickness

    vx(find(abs(vx)>5e3)) = 0;
    vy(find(abs(vy)>5e3)) = 0;
    
    % Assign to model
    md.inversion.vx_obs = vx;
    md.inversion.vy_obs = vy;
    % md.inversion.vel_obs = sqrt(vx.^2 + vy.^2);
    md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
    md.geometry.thickness = thickness;


    % Compute grounding line flux (Gt/yr)
    
    fluxes(yy) = outflux(md);

    % return

    fprintf('Year %d: Flux = %.2f Gt/yr\n', years(yy), fluxes(yy));

    fprintf('Year %d: vx range = [%.2f, %.2f], thickness range = [%.2f, %.2f]\n', ...
    years(yy), min(vx(find(in))), max(vx(find(in))), min(thickness(find(in))), max(thickness(find(in))));


end

% Save results to .mat file
save('Data/GL_flux_timeseries.mat', 'years', 'fluxes');


% addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\ISSM-Windows-MATLAB\ISSM-Windows-MATLAB\src\m\contrib\morlighem\massbalance')
% flux=outflux(md);


end%}}}




%%%%%%%%%%%%%%%%% STEP THIRTY: VISUALISING the flux gate
    % To make sure the grounding line is being used correctly as the gate
    % gate is in yellow

    %% RUN WHATEVER MD I WANT TO APPLY THIS TO FIRST!!


if perform(org, 'flux_gate'),% {{{

% Assumes 'md' is already loaded (e.g. from loadmodel)
figure;
hold on; axis equal;
title('Grounding Line (Flux Gate) Visualization');

% Plot mesh
patch('Faces', md.mesh.elements, 'Vertices', [md.mesh.x md.mesh.y], ...
    'FaceColor', 'none', 'EdgeColor', [0.7 0.7 0.7]);

% Identify grounded and floating elements
grounded = md.mask.ice_levelset > 0;
floating = md.mask.ice_levelset < 0;

% Plot grounded nodes
scatter(md.mesh.x(grounded), md.mesh.y(grounded), 10, 'red', 'filled');
% Plot floating nodes
scatter(md.mesh.x(floating), md.mesh.y(floating), 10, 'blue', 'filled');

% Highlight grounding line elements (where levelset sign changes in an element)
gl_elems = find(any(md.mask.ice_levelset(md.mesh.elements) > 0, 2) & ...
                any(md.mask.ice_levelset(md.mesh.elements) < 0, 2));

% Extract the element nodes and draw their edges
for i = 1:length(gl_elems)
    elem = md.mesh.elements(gl_elems(i), :);
    x = md.mesh.x(elem);
    y = md.mesh.y(elem);
    patch(x, y, 'k', 'EdgeColor', 'yellow', 'FaceColor', 'none', 'LineWidth', 1.5);
end

legend('Mesh edges', 'Grounded nodes', 'Floating nodes', 'Grounding line');
xlabel('X (m)');
ylabel('Y (m)');


end%}}}



%%%%%%%%%%%%%%%%% STEP THIRTY ONE: extracting .exp from GL shapefile

if perform(org, 'GL_exp'),% {{{

	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};

    indir = 'Data/GL_Shapefile/'; % can change this to wherever my .shp is

    shpa = shpread([indir 'GL_2003_UNION.shp']);
	
    if max(shpa.x)<1e6
		[x,y] = ll2xy(shpa.y,shpa.x,-1);
		expa.x = x;
		expa.y = y;
	else
		expa.x = shpa.x;
		expa.y = shpa.y;
	end
	expwrite(expa, [indir 'GL_2003_UNION.exp']);

end%}}}



%%%%%%%%%%%%%%%%% STEP THIRTY TWO: CHANGE in area


if perform(org, 'change_in_area'),% {{{

    addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');

	years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;


    areas_per_year = zeros(length(years),1);


	masks = {};

	for yy = 1:length(years)
		
        md.mask.ice_levelset(:) = 1;   
        md.mask.ocean_levelset(:) = 1;

        indir = ['Data/' num2str(years(yy)) '/']; 

        in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1); % makes values within the .exp file = 1        
        md.mask.ocean_levelset(find(in)) = -1; % sets ocean_levelset to -1 for nodes inside the ice extent
		md.mask.ice_levelset(find(in)) = -1; % sets ice_levelset to -1 for nodes inside the ice extent		

		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;


        % calculating area:
        areas = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y);
        A = averaging(md, areas, 0);
        A_masked = A(find(in)); % in metres squared?
        A_masked_kilometres_squared = A_masked / 1e6;

        areas_per_year(yy) = sum(A_masked_kilometres_squared);


	end

	%if want to save all output
	save([indir 'areas_per_year.mat'], 'areas_per_year');


    % to visualise the area calculated, to make sure the domain is right,
    % copy and paste the following into the command window:

        % A_plot = A;          % copy the full-length vertex-based area vector
        % A_plot(~in) = NaN;   % set vertices outside polygon to NaN
        % 
        % figure;
        % plotmodel(md, 'data', A_plot / 1e6, ... % km²
        %           'title', ['Ice-covered area for ' num2str(years(yy))], ...
        %           'edgecolor', 'k');
        % colorbar;
        % ylabel(colorbar, 'Area (km²)');
        % 
        % total_area_km2 = nansum(A_plot) / 1e6; 

end%}}}



%%%%%%%%%%%%%%%%% STEP THIRTY THREE: thomas params EXY


if perform(org, 'EXY_thomasparams'),% {{{

    years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
    md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};

	for yy = 1:length(years)
		
        indir = ['Data/' num2str(years(yy)) '/']; 
        in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1);
	    md.mask.ocean_levelset(find(in)) = -1;
	    md.mask.ice_levelset(find(in)) = -1;
    
	    % READ THE ITS_LIVE DATA ONTO ISSM MESH
	    fname = [indir num2str(years(yy)) '.nc'];
	    x_in = double(ncread(fname,'x'));
	    y_in = double(ncread(fname,'y'));
	    vx_in = double(ncread(fname,'VX'));
	    vy_in = double(ncread(fname,'VY'));
    
	    % Get the vx, vy onto the ISSM grid
	    vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
	    vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
	    md.inversion.vx_obs(find(in)) = vx(find(in));
	    md.inversion.vy_obs(find(in)) = vy(find(in));
	    md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
    

       	mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        md = mechanicalproperties(md, vx, vy);
        [~,~,~,exy,~] = thomasparams(md,'eq','Thomas','smoothing',0,'coordsys','longitudinal');
        
        sec_per_year = 365 * 24 * 60 * 60;
        exy_per_year = exy * sec_per_year
        exy_per_year(mask == 0) = NaN;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');

        plotmodel(md,...
            'data',exy_per_year,...
            'figure', yy,...
            'colormap#1','parula',...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'caxis',[-0.02 0.02],...    
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        colormap(redblue(64));
        cb = colorbar;
        ylabel(cb, 'exy (longitudinal) (a-1)', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['exy ' num2str(years(yy))], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        
        exportgraphics(gcf,['Plots/exy_' num2str(years(yy)) '_thomasparams.png'])

    end



end%}}}



%%%%%%%%%%%%%%%%% STEP THIRTY FOUR: PLOTTING ALL EXY TOGETHER

% KNMAX = MAXIMUM BUTTRESSING

if perform(org, 'EXY_all_years'),% {{{

	files = {};
for yy = 1:length(years)
    files{end+1} = ['Plots/exy_' num2str(years(yy)) '_thomasparams.png'];
end

% Create montage
figure;
montage(files, 'Size', [3 4], 'BackgroundColor', 'w'); % 3 rows, 4 cols — adjust as needed
title('EXY');
exportgraphics(gcf, 'Plots/EXY_all_years.png', 'Resolution', 300);

end%}}}




%%%%%%%%%%%%%%%%% STEP THIRTY FIVE: thomas params EYY


if perform(org, 'EYY_thomasparams'),% {{{

    years = [2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
    md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};

	for yy = 1:length(years)
		
        indir = ['Data/' num2str(years(yy)) '/']; 
        in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1);
	    md.mask.ocean_levelset(find(in)) = -1;
	    md.mask.ice_levelset(find(in)) = -1;
    
	    % READ THE ITS_LIVE DATA ONTO ISSM MESH
	    fname = [indir num2str(years(yy)) '.nc'];
	    x_in = double(ncread(fname,'x'));
	    y_in = double(ncread(fname,'y'));
	    vx_in = double(ncread(fname,'VX'));
	    vy_in = double(ncread(fname,'VY'));
    
	    % Get the vx, vy onto the ISSM grid
	    vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
	    vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
	    md.inversion.vx_obs(find(in)) = vx(find(in));
	    md.inversion.vy_obs(find(in)) = vy(find(in));
	    md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);
    

       	mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        md = mechanicalproperties(md, vx, vy);
        [~,~,~,eyy,~] = thomasparams(md,'eq','Thomas','smoothing',0,'coordsys','longitudinal');
        
        sec_per_year = 365 * 24 * 60 * 60;
        eyy_per_year = eyy * sec_per_year
        eyy_per_year(mask == 0) = NaN;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        width_in = 160 / 25.4;
        height_in = 160 / 25.4;
        
        fig = figure(yy);
        set(fig, 'Units', 'inches', 'Position', [1 1 width_in height_in]);
        
        addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');

        plotmodel(md,...
            'data',eyy_per_year,...
            'figure', yy,...
            'colormap#1','parula',...
            'xlim#all',[998648.2458, 1171935.9538],...
            'ylim#all', [-2148431.9935, -2033743.9935],...
            'caxis',[-0.02 0.02],...    
            'mask#all',mask,...
            'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp']);

        colormap(redblue(64));
        cb = colorbar;
        ylabel(cb, 'eyy (longitudinal) (a-1)', 'FontName', 'Times New Roman', 'FontSize', 12);

        figure(yy);
        xlabel('Eastings (m)');
        ylabel('Northings (m)');
        title(['eyy ' num2str(years(yy))], 'FontWeight', 'normal');

        box off
        x = xlim;
        y = ylim;
        
        hold on
        plot([x(1) x(2)], [y(1) y(1)], 'k'); % bottom
        plot([x(1) x(2)], [y(2) y(2)], 'k'); % top
        plot([x(1) x(1)], [y(1) y(2)], 'k'); % left
        plot([x(2) x(2)], [y(1) y(2)], 'k'); % right
        
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        
        
        exportgraphics(gcf,['Plots/eyy_' num2str(years(yy)) '_thomasparams.png'])

    end



end%}}}



%%%%%%%%%%%%%%%%% STEP THIRTY SIX: PLOTTING ALL EYY TOGETHER


if perform(org, 'EYY_all_years'),% {{{

	files = {};
for yy = 1:length(years)
    files{end+1} = ['Plots/eyy_' num2str(years(yy)) '_thomasparams.png'];
end

% Create montage
figure;
montage(files, 'Size', [3 4], 'BackgroundColor', 'w'); % 3 rows, 4 cols — adjust as needed
title('EYY');
exportgraphics(gcf, 'Plots/EYY_all_years.png', 'Resolution', 300);

end%}}}








% TO QUICKLY PLOT:
% 
% plotmodel(md, 'data', md.geometry.thickness, ...
%               'title', 'Ice Thickness (m)', ...
%               'caxis', [0 max(md.geometry.thickness)], ...
%               'colorbar', 'on');
% 
% vel = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);  % or initialization.*
% plotmodel(md, 'data', vel, ...
%               'title', 'Observed Velocity Magnitude (m/yr)', ...
%               'caxis', [0 max(vel)], ...
%               'colorbar', 'on');
