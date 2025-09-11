steps = [2]; %change this if you want to run different steps
    % step 1 = calculating strain rate for a single year
    % step 2 = crevasse formation for 0.001 - 0.163 a-1
    % step 3 = crevasse formation for 0.004 - 0.163 a-1
    % step 4 = crevasse formation for => 0.01 a-1

    % step 5 = calculating deviatoric stress for a single year
    % step 6 = crevasse formation for 90 to 320 kPa ** do not use


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	1: calculating strain rate:
if perform(org, 'strain_reparam'),% {{{

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

        % md=mechanicalproperties(md,md.inversion.vx_obs,md.inversion.vy_obs);
        % [alpha,beta,theta,exx,sigxx]=thomasparams(md,'eq','Thomas','smoothing',2,'coordsys','longitudinal');

		% calculate the buttressing numbers
		[Knflow{yy},Nflow{yy},Knmax{yy},Nmax{yy},N0{yy}]=buttressingnumber(md);

		mask= ones(md.mesh.numberofvertices,1);
		mask(find(md.mask.ice_levelset>0)) = 0;
		masks{yy} = mask;

        %% these two lines of code save each individual year's data
        % outfilename = ['buttressing_outfiles_' num2str(years(yy)) '.mat'];
        % save(outfilename, 'Knflow', 'Nflow', 'Knmax', 'Nmax','masks', 'x_in', 'y_in');
	end

    %% FROM FELICITY DISCORD FOR STRAIN RATE?

    md=mechanicalproperties(md,md.inversion.vx_obs,md.inversion.vy_obs);
    [alpha,beta,theta,exx,sigxx]=thomasparams(md,'eq','Thomas','smoothing',0,'coordsys','longitudinal');

    % — 2. Convert strain‑rate units:  s‑1  →  a‑1  (years‑1) —
    sec_per_year = 365.25 * 24 * 60 * 60;   % = 31 557 600 s
    exx_yr = exx * sec_per_year;            % now in yr‑1 (a‑1)

    % mask = masks{2};
    mask = double(in);
    exx_yr(mask == 0) = NaN;


    addpath('C:\Users\jjls0\OneDrive\Documents\MATLAB\redblue_v101');
    colormap(redblue(64)); % THIS COLORMAP is downloaded from online. white = 0, negative = red, positive = blue. opposite to matlab inbuilt 'redblue'.


    f = figure;
	plotmodel(md, 'data', exx_yr, 'title', 'exx 2018', 'mask#all', mask, 'caxis', [-0.02 0.02]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Strain Rate (a-1)');
    
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);


%exportgraphics(f, 'Plots/2018_EXX.jpg', 'Resolution', 300);

% save([indir '2008_2018_percent_change.mat'], 'percent_change', 'x_in', 'y_in');


end%}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	2: strain rate between 0.001 and 0.163 for crevasse formation
            % as suggested by campbell et al., 2017
if perform(org, 'crevasse_classification_1'),% {{{

color_data = nan(size(exx_yr));   

color_data(exx_yr < 0.001) = 0;           % below range → blue
color_data(exx_yr > 0.163) = 0;           % above range → blue
color_data(exx_yr >= 0.001 & exx_yr <= 0.163) = 1;  % within range → red

cmap = [0 0 1;   % blue
        1 0 0];  % red

f = figure;
plotmodel(md, ...
    'data', color_data, ...
    'mask#all', mask, ...
    'title', 'exx 2018: crevasse classification (0.001 to 0.163 a-1)', ...
    'caxis', [0 1]);

colormap(cmap);
cb = colorbar;
cb.Ticks = [0, 1];
cb.TickLabels = {'Outside range', '0.001–0.163'};
ylabel(cb, 'Strain Rate (a-1)');

xlim([998648.2458 1171935.9538]); 
ylim([-2148431.9935 -2033743.9935]);


end%}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	3: strain rate between 0.004 and 0.163 for crevasse formation
            % as suggested byHambrey and Muller 1978
if perform(org, 'crevasse_classification_2'),% {{{

color_data = nan(size(exx_yr));   

color_data(exx_yr < 0.004) = 0;           % below range → blue
color_data(exx_yr > 0.163) = 0;           % above range → blue
color_data(exx_yr >= 0.004 & exx_yr <= 0.163) = 1;  % within range → red

cmap = [0 0 1;   % blue
        1 0 0];  % red

f = figure;
plotmodel(md, ...
    'data', color_data, ...
    'mask#all', mask, ...
    'title', 'exx 2018: crevasse classification (0.004 to 0.163 a-1)', ...
    'caxis', [0 1]);

colormap(cmap);
cb = colorbar;
cb.Ticks = [0, 1];
cb.TickLabels = {'Outside range', '0.004–0.163'};
ylabel(cb, 'Strain Rate (a-1)');

xlim([998648.2458 1171935.9538]); 
ylim([-2148431.9935 -2033743.9935]);


end%}}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	4: strain rate greater than 0.01 for crevasse formation
            % as suggested by Holdsworth, 1969
if perform(org, 'crevasse_classification_3'),% {{{

color_data = nan(size(exx_yr));   

color_data(exx_yr < 0.01) = 0;       % below threshold → blue
color_data(exx_yr >= 0.01) = 1;      % ≥ 0.01 → red

cmap = [0 0 1;   % blue
        1 0 0];  % red

f = figure;
plotmodel(md, ...
    'data', color_data, ...
    'mask#all', mask, ...
    'title', 'exx 2018: crevasse classification (=> 0.163 a-1)', ...
    'caxis', [0 1]);

colormap(cmap);
cb = colorbar;
cb.Ticks = [0, 1];
cb.TickLabels = {'< 0.01', '≥ 0.01'};
ylabel(cb, 'Strain Rate Classification');

xlim([998648.2458 1171935.9538]); 
ylim([-2148431.9935 -2033743.9935]);

end %}}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	5: calculating stress for the given year in step 1 (kPa)
if perform(org, 'crevasse_classification_4'),% {{{

    md=mechanicalproperties(md,md.inversion.vx_obs,md.inversion.vy_obs);
    [alpha,beta,theta,exx,sigxx]=thomasparams(md,'eq','Thomas','smoothing',0,'coordsys','longitudinal');

    clean_sigxx = sigxx / 1000

    % mask = masks{2};
    mask = double(in);
    clean_sigxx(mask == 0) = NaN;

    f = figure;
	plotmodel(md, 'data', clean_sigxx, 'title', 'sigxx 2018', 'mask#all', mask, 'caxis', [-300 300]);
    colormap(redblue(64));
    cb = colorbar;
    ylabel(cb, 'Stress (kPa)');
    
    xlim([998648.2458 1171935.9538]); 
    ylim([-2148431.9935 -2033743.9935]);


end %}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	6: stress within 90 to 320 kPa for crevasse formation
            % as suggested by Vaughan, 1993
        %% THIS DOES NOT WORK, BECAUSE ONLY DEVIATORIC STRESS is calculated, but we need the TOTAL STRESS TENSOR instead.
if perform(org, 'crevasse_classification_5'),% {{{

color_data = nan(size(clean_sigxx));   

color_data(clean_sigxx < 90) = 0;           % below range → blue
color_data(clean_sigxx > 320) = 0;           % above range → blue
color_data(clean_sigxx >= 90 & clean_sigxx <= 320) = 1;  % within range → red

cmap = [0 0 1;   % blue
        1 0 0];  % red

f = figure;
plotmodel(md, ...
    'data', color_data, ...
    'mask#all', mask, ...
    'title', 'sigxx 2018: crevasse classification (90 to 320 kPa)', ...
    'caxis', [0 1]);

colormap(cmap);
cb = colorbar;
cb.Ticks = [0, 1];
cb.TickLabels = {'Outside range', '90 to 320'};
ylabel(cb, 'Deviatoric Stress (kPa)');

xlim([998648.2458 1171935.9538]); 
ylim([-2148431.9935 -2033743.9935]);


end %}}}





