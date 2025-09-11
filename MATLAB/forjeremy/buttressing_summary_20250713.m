steps = [2]; %change this if you want to run different steps
    %step 1 = % and absolute change in buttressing number
    %step 2 = absolute and percentage change in strain rate
    %step 3 = percentage change in velocity
    %step 4 = components of buttressing (% and absolute): thickness, xx,
                 %xy, yy, velocity
    %steo 5 + 6 = calculating the 'dynamic' and 'geometric' term for the
        %changes in buttressing - see google docs for info on how equation was
        %derived from the original BN calc.

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

%	steps 1
if perform(org, 'Reparam'),% {{{

	years = [2008, 2018];    % CHANGE THE YEARS AS I LIKE
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
		shp2exp([indir 'ice_extent_' num2str(years(yy)) '.shp'],[indir 'ice_extent_' num2str(years(yy)) '.exp']);

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

		vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		vel = sqrt(vx.^2+vy.^2);
		md.inversion.vx_obs(find(in)) = vx(find(in));
		md.inversion.vy_obs(find(in)) = vy(find(in));
		md.inversion.vel_obs(find(in)) = vel(find(in));

		thickness = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		md.geometry.thickness(find(in)) = thickness(find(in));

		% calculate the buttressing numbers
		[Knflow{yy},Nflow{yy},Knmax{yy},Nmax{yy},N0{yy}]=buttressingnumber(md);

		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;

		%% if want to save output per iteration of the for loop
		% Knflow_temp = Knflow{yy};
		% save([indir 'buttressing_outfiles.mat'], 'varnames');

        %% these two lines of code save each individual year's data
        % outfilename = ['buttressing_outfiles_' num2str(years(yy)) '.mat'];
        % save(outfilename, 'Knflow', 'Nflow', 'Knmax', 'Nmax','masks', 'x_in', 'y_in');
	end

	%if want to save all output
	% save([indir 'buttressing_alloutfiles.mat'], 'Knflow', 'Nflow', 'Knmax', 'Nmax','masks', 'x_in', 'y_in'); % HOW DO I ADD X AND Y CORDS TO THIS?

	%% TODO: FIX THE BELOW BY ADDING THE NUMBRER OF KNFLOW DATASETS YOU WANT TO PLOT
	% plotmodel(md,'data',md.inversion.vel_obs, 'data',md.inversion.vel_obs-velo,'caxis#2',[-50 50],'colormap#2','redblue',...
	% 	'data',Knflow{1},'mask#3',masks{1},'caxis#3',[0 2],...
	% 	'data',Knflow{2},'mask#4',masks{2},'caxis#4',[0 2],...
	% 	'figure',1,'title#1','new velocity', 'title#2','new velocity - old velocity',...
	% 	'title#3','buttressing number for year XX'...
	% 	);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SAVES EACH INDIVIDUAL YEAR AS A JPG
    % f = figure;
    % plotmodel(md, 'data', Knflow{yy}, 'mask#all', mask, 'caxis', [0 2], 'title', [num2str(years(yy))]);
    % xlim([998648.2458 1171935.9538]); 
    % ylim([-2148431.9935 -2033743.9935]); 
    % cb = colorbar;
    % ylabel(cb, 'Buttressing Number');
    % 
    % % f = gcf; % defines f to save plot
    % exportgraphics(f, 'Plots/2020.jpg', 'Resolution', 300);
    % % exportgraphics(f,['Plots/2001.pdf']);   % actually saves the plot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% VISUALISING THE PERCENTAGE CHANGE BETWEN TWO YERAS

percent_change = 100 * (Knflow{2} - Knflow{1}) ./ Knflow{1};

% to 'guard' against division by 0 which will give NaN or Inf values:
zero_nans_BN1 = Knflow{1};
zero_nans_BN1(zero_nans_BN1 == 0) = NaN;
percent_change = 100 * (Knflow{2} - zero_nans_BN1) ./ zero_nans_BN1;

mask = masks{1};  % using ice mask from 2008
percent_change(mask == 0) = NaN;

f = figure;
plotmodel(md, ...
    'data', percent_change, ...
    'caxis', [-100 100], ...  
    'title', 'Percent Change in Buttressing Number (2008 to 2018)');

addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');
colormap(redblue(64)); % THIS COLORMAP is downloaded from online. white = 0, negative = red, positive = blue. opposite to matlab inbuilt 'redblue'.

    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]); 
    cb = colorbar;
    ylabel(cb, 'Change (%)');

% exportgraphics(f, 'Plots/2008_2018_percent_change.jpg', 'Resolution', 300);

% save([indir '2008_2018_percent_change.mat'], 'percent_change', 'x_in', 'y_in');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Absolute change in buttressing number between two years

absolute_change_BN = Knflow{2} - Knflow{1};
mask = masks{2}; % uses the 2018 mask, which includes CWIS
absolute_change_BN(mask == 0) = NaN;

f = figure;
plotmodel(md, ...
    'data', absolute_change_BN, ...
    'caxis', [-5 5], ...  
    'title', 'Absolute Change in Buttressing Number (2008 to 2018)');

addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');
colormap(redblue(64)); % THIS COLORMAP is downloaded from online. white = 0, negative = red, positive = blue. opposite to matlab inbuilt 'redblue'.

    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]); 
    cb = colorbar;
    ylabel(cb, 'Change in Buttressing (BN)');

end%}}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%	2: calculating PERCENT CHANGE in strain rate WITH MASK !!!!!!!!!!!!!
if perform(org, 'TESTINGGGGG'),% {{{

	years = [2008, 2018];    % CHANGE THE YEARS AS I LIKE
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
		shp2exp([indir 'ice_extent_' num2str(years(yy)) '.shp'],[indir 'ice_extent_' num2str(years(yy)) '.exp']);

		in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1);

        md.mask.ice_levelset  = +1 * ones(md.mesh.numberofvertices,1);   %  +1 = grounded by default
        md.mask.ocean_levelset = +1 * ones(md.mesh.numberofvertices,1);

		md.mask.ocean_levelset(find(in)) = -1;
		md.mask.ice_levelset(find(in)) = -1;

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		vx_in = double(ncread(fname,'VX'));
		vy_in = double(ncread(fname,'VY'));
		thickness_in = double(ncread(fname,'Thickness'));

		vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		vel = sqrt(vx.^2+vy.^2);
		md.inversion.vx_obs(find(in)) = vx(find(in));
		md.inversion.vy_obs(find(in)) = vy(find(in));
		md.inversion.vel_obs(find(in)) = vel(find(in));

		thickness = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		md.geometry.thickness(find(in)) = thickness(find(in));

        % md=mechanicalproperties(md,md.inversion.vx_obs,md.inversion.vy_obs);
        % [alpha,beta,theta,exx,sigxx]=thomasparams(md,'eq','Thomas','smoothing',2,'coordsys','longitudinal');

		% calculate the buttressing numbers
		[Knflow{yy},Nflow{yy},Knmax{yy},Nmax{yy},N0{yy}]=buttressingnumber(md);

		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;


        md = mechanicalproperties(md, vx, vy);
        [~,~,~,exx,~] = thomasparams(md,'eq','Thomas','smoothing',0,'coordsys','longitudinal');
        exx_all{yy} = exx;

		%% if want to save output per iteration of the for loop
		% Knflow_temp = Knflow{yy};
		% save([indir 'buttressing_outfiles.mat'], 'varnames');

        %% these two lines of code save each individual year's data
        % outfilename = ['buttressing_outfiles_' num2str(years(yy)) '.mat'];
        % save(outfilename, 'Knflow', 'Nflow', 'Knmax', 'Nmax','masks', 'x_in', 'y_in');
	end

    % PERCENTAGE CHANGE IN STRAIN RATE

    sec_per_year = 365.25 * 24 * 60 * 60;
    exx_2008 = exx_all{1} * sec_per_year;
    exx_2018 = exx_all{2} * sec_per_year;


    exx_ref = exx_2008;
    absolute_change_strain = exx_2018 - exx_2008;
    percent_change_strain = 100 * (exx_2018 - exx_ref) ./ exx_ref;


    f = figure;
	plotmodel(md, 'data', percent_change_strain, 'title', 'Percent Change in Strain Rate (2008 to 2018)', 'mask#all', mask, 'caxis', [-100 100]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Change in Strain Rate (%)');
    
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ABSOLUTE CHANGE IN STRAIN RATE:

    absolute_change_strain = exx_2018 - exx_2008;

    f = figure;
	plotmodel(md, 'data', absolute_change_strain, 'title', 'Change in Strain Rate (2008 to 2018)', 'mask#all', mask, 'caxis', [-0.01 0.01]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Change in Strain Rate (a-1)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');
% colormap(redblue(64)); % THIS COLORMAP is downloaded from online. white = 0, negative = red, positive = blue. opposite to matlab inbuilt 'redblue'.

%exportgraphics(f, 'Plots/2008_2018_percent_change_STRAIN.jpg', 'Resolution', 300);

% save([indir '2008_2018_percent_change_STRAIN.mat'], 'percent_change_strain', 'x_in', 'y_in');


end%}}}


%% AND THEN COPY AND PASTE THE FOLLOWING CODE INTO THE COMMAND WINDOW IF PLOT SHOWS BLANK:

    % sec_per_year = 365.25 * 24 * 60 * 60;
    % exx_2008 = exx_all{1} * sec_per_year;
    % exx_2018 = exx_all{2} * sec_per_year;
    % 
    % 
    % exx_ref = exx_2008;
    % percent_change_strain = 100 * (exx_2018 - exx_ref) ./ exx_ref;
    % 
    % 
    % f = figure;
	% plotmodel(md, 'data', percent_change_strain, 'title', 'Percent Change in Strain Rate (2008 to 2018)', 'mask#all', mask, 'caxis', [-100 100]);
    % colormap(redblue(64));
    % cb = colorbar;
    % ylabel(cb, 'Change in Strain Rate (%)');
    % 
    % xlim([998648.2458 1171935.9538]); 
    % ylim([-2148431.9935 -2033743.9935]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%	3: CHANGE IN VELOCITY PERCENTAGE
if perform(org, 'velocity_change'),% {{{

	years = [2008, 2018];    % CHANGE THE YEARS AS I LIKE
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
		shp2exp([indir 'ice_extent_' num2str(years(yy)) '.shp'],[indir 'ice_extent_' num2str(years(yy)) '.exp']);

		in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1);

        md.mask.ice_levelset  = +1 * ones(md.mesh.numberofvertices,1);   %  +1 = grounded by default
        md.mask.ocean_levelset = +1 * ones(md.mesh.numberofvertices,1);

		md.mask.ocean_levelset(find(in)) = -1;
		md.mask.ice_levelset(find(in)) = -1;

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		vx_in = double(ncread(fname,'VX'));
		vy_in = double(ncread(fname,'VY'));
		thickness_in = double(ncread(fname,'Thickness'));

		vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		vel = sqrt(vx.^2+vy.^2);
		md.inversion.vx_obs(find(in)) = vx(find(in));
		md.inversion.vy_obs(find(in)) = vy(find(in));
		md.inversion.vel_obs(find(in)) = vel(find(in));

		thickness = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		md.geometry.thickness(find(in)) = thickness(find(in));



		% calculate the buttressing numbers
		[Knflow{yy},Nflow{yy},Knmax{yy},Nmax{yy},N0{yy}]=buttressingnumber(md);

		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;

        vel_all{yy} = vel;


	end

    vel_change = vel_all{2} - vel_all{1};
    vel_change_percent = 100 * vel_change ./ vel_all{1};


    f = figure;
	plotmodel(md, 'data', vel_change_percent, 'title', 'Percent Change in Velocity (2008 to 2018)', 'mask#all', mask, 'caxis', [-100 100]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Change in Strain Rate (%)');
    
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);

% addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');
% colormap(redblue(64)); % THIS COLORMAP is downloaded from online. white = 0, negative = red, positive = blue. opposite to matlab inbuilt 'redblue'.

%exportgraphics(f, 'Plots/2008_2018_percent_change_STRAIN.jpg', 'Resolution', 300);

% save([indir '2008_2018_percent_change_STRAIN.mat'], 'percent_change_strain', 'x_in', 'y_in');



%%IN ORDER TO FIND THE MAGNITUDE OF CHANGE AT A GIVEN POINT:
% Specify your coordinate of interest (in the same projection as md.mesh.x/y)
x_target = 1051730;     % example x coordinate
y_target = -2105720;     % example y coordinate

% Find the nearest mesh node
dist2 = (md.mesh.x - x_target).^2 + (md.mesh.y - y_target).^2;
[~, node_id] = min(dist2);

% Extract the velocity change values at that node
vchange_value = vel_change(node_id);
vchange_percent = vel_change_percent(node_id);

% Print or display
fprintf('Velocity change at (%.1f, %.1f): %.2f m/yr\n', x_target, y_target, vchange_value);
fprintf('Percent change at (%.1f, %.1f): %.2f %%\n', x_target, y_target, vchange_percent);




end%}}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	4: trying to find components of buttressing
if perform(org, 'buttressing_comp'),% {{{

	years = [2008, 2018];    % CHANGE THE YEARS AS I LIKE
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
		shp2exp([indir 'ice_extent_' num2str(years(yy)) '.shp'],[indir 'ice_extent_' num2str(years(yy)) '.exp']);

		in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1);

        md.mask.ice_levelset  = +1 * ones(md.mesh.numberofvertices,1);   %  +1 = grounded by default
        md.mask.ocean_levelset = +1 * ones(md.mesh.numberofvertices,1);

		md.mask.ocean_levelset(find(in)) = -1;
		md.mask.ice_levelset(find(in)) = -1;

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		vx_in = double(ncread(fname,'VX'));
		vy_in = double(ncread(fname,'VY'));
		thickness_in = double(ncread(fname,'Thickness'));

		vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		vel = sqrt(vx.^2+vy.^2);
		md.inversion.vx_obs(find(in)) = vx(find(in));
		md.inversion.vy_obs(find(in)) = vy(find(in));
		md.inversion.vel_obs(find(in)) = vel(find(in));

		thickness = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		md.geometry.thickness(find(in)) = thickness(find(in));


		% calculate the buttressing numbers
		[Knflow{yy},Nflow{yy},Knmax{yy},Nmax{yy},N0{yy},xxa{yy},xya{yy}, yya{yy},nxf{yy},nyf{yy}]=buttressingnumber(md);

		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;

        vel_all{yy} = vel;

	end


% COPY AND PASTE EACH DESIRED SECTION INTO COMMAND WINDOW:
% COPY AND PASTE EACH DESIRED SECTION INTO COMMAND WINDOW:
% COPY AND PASTE EACH DESIRED SECTION INTO COMMAND WINDOW:
% COPY AND PASTE EACH DESIRED SECTION INTO COMMAND WINDOW:
% COPY AND PASTE EACH DESIRED SECTION INTO COMMAND WINDOW:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes in N0 (thickness)
    delta_N0 = N0{2} - N0{1};

    f = figure;
	plotmodel(md, 'data', delta_N0, 'title',  'Change in N0 (2008 to 2018)', 'mask#all', mask);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Delta N0 (Pa)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERCENTchanges in N0 (thickness)
    delta_N0 = N0{2} - N0{1};
    delta_N0_percent = 100 * delta_N0 ./ N0{1};

    f = figure;
	plotmodel(md, 'data', delta_N0_percent, 'title',  'Percent Change in N0 (2008 to 2018)', 'mask#all', mask, 'caxis', [-100 100]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'N0 Change (%)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes in average xx stress component:

    delta_xxa = xxa{2} - xxa{1};

    f = figure;
	plotmodel(md, 'data', delta_xxa, 'title',  'Change in xxa (2008 to 2018)', 'mask#all', mask, 'caxis', [-1e6 1e6]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Delta xxa (Pa)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes in PERCENT average xx stress component:

    delta_xxa = xxa{2} - xxa{1};
    delta_xxa_percent = 100 * delta_xxa ./ xxa{1};

    f = figure;
	plotmodel(md, 'data', delta_xxa_percent, 'title',  'Percent Change in xxa (2008 to 2018)', 'mask#all', mask,'caxis', [-100 100]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'xxa Change (%)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes in PERCENT average xy stress component:

    delta_xya = xya{2} - xya{1};
    delta_xya_percent = 100 * delta_xya ./ xya{1};

    f = figure;
	plotmodel(md, 'data', delta_xya_percent, 'title',  'Percent Change in xya (2008 to 2018)', 'mask#all', mask,'caxis', [-100 100]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'xya Change (%)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes in PERCENT average yy stress component:

    delta_yya = yya{2} - yya{1};
    delta_yya_percent = 100 * delta_yya ./ yya{1};

    f = figure;
	plotmodel(md, 'data', delta_yya_percent, 'title',  'Percent Change in yya (2008 to 2018)', 'mask#all', mask,'caxis', [-100 100]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'yya Change (%)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes in PERCENT average nxf flow component: SHOULD REALLY PLOT
% ABSOLUTE CHANGE OR ANGUALR CHANGE INSTEAD

    delta_nxf = nxf{2} - nxf{1};
    delta_nxf_percent = 100 * delta_nxf ./ nxf{1};

    f = figure;
	plotmodel(md, 'data', delta_nxf_percent, 'title',  'Percent Change in nxf (2008 to 2018)', 'mask#all', mask,'caxis', [-100 100]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'nxf Change (%)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes in PERCENT VELOCITY:

    delta_vel = vel_all{2} - vel_all{1};
    delta_vel_percent = 100 * delta_vel ./ vel_all{1};

    f = figure;
	plotmodel(md, 'data', delta_vel_percent, 'title',  'Percent Change in Velocity (2008 to 2018)', 'mask#all', mask,'caxis', [-100 100]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Velocity Change (%)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes in absolute velocity:

f = figure;
	plotmodel(md, 'data', delta_vel, 'title',  'Change in Velocity (2008 to 2018)', 'mask#all', mask, 'caxis', [-150 150]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Velocity Change (m/yr)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');
% colormap(redblue(64)); % THIS COLORMAP is downloaded from online. white = 0, negative = red, positive = blue. opposite to matlab inbuilt 'redblue'.

%exportgraphics(f, 'Plots/2008_2018_percent_change_STRAIN.jpg', 'Resolution', 300);

% save([indir '2008_2018_percent_change_STRAIN.mat'], 'percent_change_strain', 'x_in', 'y_in');


end%}}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	5: trying to quantify the reasons for change in buttressing 
if perform(org, 'buttressing_comp_quantifying'),% {{{

	years = [2008, 2018];    % CHANGE THE YEARS AS I LIKE
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

    out = repmat(struct('Knflow', [], 'N0', [], 'xxa', [], 'xya', [], 'yya', [], 'vel', []), length(years), 1);

	Knflow = {};
	Nflow = {};
	Knmax = {};
	Nmax = {};
	N0 = {};
	masks = {};

    exx = cell(1, length(years));

	for yy = 1:length(years)
		% CONVERT SHP TO EXP
		indir = ['Data/' num2str(years(yy)) '/']; %% Jeremy: change this to relevant data directory
		shp2exp([indir 'ice_extent_' num2str(years(yy)) '.shp'],[indir 'ice_extent_' num2str(years(yy)) '.exp']);

		in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_' num2str(years(yy)) '.exp'],1);

        md.mask.ice_levelset  = +1 * ones(md.mesh.numberofvertices,1);   %  +1 = grounded by default
        md.mask.ocean_levelset = +1 * ones(md.mesh.numberofvertices,1);

		md.mask.ocean_levelset(find(in)) = -1;
		md.mask.ice_levelset(find(in)) = -1;

		% READ THE ITS_LIVE DATA ONTO ISSM MESH
		fname = [indir num2str(years(yy)) '.nc'];
		x_in = double(ncread(fname,'x'));
		y_in = double(ncread(fname,'y'));
		vx_in = double(ncread(fname,'VX'));
		vy_in = double(ncread(fname,'VY'));
		thickness_in = double(ncread(fname,'Thickness'));

		vx = InterpFromGrid(x_in,y_in,vx_in',md.mesh.x,md.mesh.y);
		vy = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		vel = sqrt(vx.^2+vy.^2);
		md.inversion.vx_obs(find(in)) = vx(find(in));
		md.inversion.vy_obs(find(in)) = vy(find(in));
		md.inversion.vel_obs(find(in)) = vel(find(in));

		thickness = InterpFromGrid(x_in,y_in,vy_in',md.mesh.x,md.mesh.y);
		md.geometry.thickness(find(in)) = thickness(find(in));


		% calculate the buttressing numbers
		[Knflow{yy},Nflow{yy},Knmax{yy},Nmax{yy},N0{yy},xxa{yy},xya{yy}, yya{yy},nxf{yy},nyf{yy}]=buttressingnumber(md);

		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;

        out(yy).mask = mask;

        vel_all{yy} = vel;

        md = mechanicalproperties(md, vx, vy);
        [~,~,~,exx,~] = thomasparams(md,'eq','Thomas','smoothing',0,'coordsys','longitudinal');
        exx_all{yy} = exx;

        out(yy).Knflow = Knflow{yy};   % buttressing number field
        out(yy).N0     = N0{yy};       % hydrostatic baseline stress
        out(yy).xxa    = xxa{yy};      % ⟨σ'xx⟩
        out(yy).xya    = xya{yy};      % ⟨σ'xy⟩
        out(yy).yya    = yya{yy};      % ⟨σ'yy⟩
        out(yy).vel    = vel_all{yy};          % speed field
        out(yy).exx    = exx_all{yy};
        out(yy).Nflow = Nflow{yy};
        % out(yy).thickness = thickness{yy};

    end

% CALCULATING CHANGE FIELDS (2018 minus 2008)
    dKn   = out(2).Knflow - out(1).Knflow;
    dN0   = out(2).N0     - out(1).N0;
    dxxa  = out(2).xxa    - out(1).xxa;
    dxya  = out(2).xya    - out(1).xya;
    dyya  = out(2).yya    - out(1).yya;
    dvel  = out(2).vel    - out(1).vel;
    dexx  = out(2).exx    - out(1).exx;   % if available
    % dH = out(2).thickness - out(1).thickness;

    mask_common = out(1).mask & out(2).mask;   % 1 where floating in both years
    flds = {'dKn','dN0','dxxa','dxya','dyya','dvel','dexx'};
    for f = 1:numel(flds)
        tmp = eval(flds{f});
        tmp(~mask_common) = NaN;
        assignin('caller', flds{f}, tmp);      % overwrite
    end


end%}}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%STEP 6:

if perform(org, 'buttressing_comp_quantifying_PART_TWO'),% {{{

% qantify thinning (geom) vs stress (dyn)
    % see google docs for how these equations were derived from the
         % original buttressing equation
    dNflow = out(2).Nflow - out(1).Nflow;          % compute the same way as dxxa …
    term_dyn   = - dNflow ./ out(1).N0;            % dynamic contribution
    term_geom  =  (out(2).Nflow .* dN0) ./ (out(1).N0.^2);   % geometric contribution

% Map of ΔKnflow
    figure;
    plotmodel(md, 'data', dKn, 'title','ΔKn_{flow}', 'caxis',[-0.5 0.5], ...
                   'mask#all', mask_common); 
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Knflow Change in Buttressing Number (dimensionless)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);

% Dynamic term
    figure;
    plotmodel(md, 'data', term_dyn, 'title','Dynamic Term (Strain Influence)', ...
                   'caxis',[-0.5 0.5], 'mask#all', mask_common);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Dynamic Change in Buttressing Number (dimensionless)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);

% Geometric term
    figure;
    plotmodel(md, 'data', term_geom, 'title','Geometric Term (Thinning Influence)', ...
                   'caxis',[-0.5 0.5], 'mask#all', mask_common);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Geometric Change in Buttressing Number (dimensionless)');
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);


end%}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






