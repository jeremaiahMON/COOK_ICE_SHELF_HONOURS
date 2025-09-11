steps = [1]; %change this if you want to run different steps

% Step 1: just 2015 (might not work if i deleted the 2015 nc file)
% Step 2: plotting velocities, maximum buttressing, and buttressing number
%   for all the years listed
% Step 3: plotting the buttressing number for all years in one figure



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
if perform(org, 'Reparam_2015'),% {{{

   md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

	% CONVERT SHP TO EXP
	indir = 'Data/2015/'; %% Jeremy: change this to relevant data directory
	shp2exp([indir 'ice_extent_2015.shp'],[indir 'ice_extent_2015.exp']);

	in=ContourToNodes(md.mesh.x,md.mesh.y,[indir 'ice_extent_2015.exp'],1);
	md.mask.ocean_levelset(find(in)) = -1;
	md.mask.ice_levelset(find(in)) = -1;

	% READ THE ITS_LIVE DATA ONTO ISSM MESH
	fname = [indir '2015.nc'];
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
	[Knflow2015,Nflow2015,Knmax2015,Nmax2015,N02015]=buttressingnumber(md);

	mask2015 = ones(md.mesh.numberofvertices,1);
	mask2015(find(md.mask.ice_levelset>0)) = 0;
	plotmodel(md,'data',md.inversion.vel_obs, 'data',md.inversion.vel_obs-velo,'caxis#2',[-50 50],'colormap#2','redblue',...
		'data',Knmax2015, 'data',Knflow2015, 'figure',5,'caxis#3,4',[0 2],'mask#all',mask2015,...
		'title#1','new velocity', 'title#2','new velocity - old velocity',...
		'title#3','maximum buttressing', 'title#4','buttressing number, flow calc')


end%}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP TWO:

if perform(org, 'Reparam'),% {{{

	years = [2001, 2006]; % JEZ CHANGE THIS
	%years = [2015, 2016, 2020, 2022];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

    

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

		plotmodel(md,'data',md.inversion.vel_obs, 'data',md.inversion.vel_obs-velo,'caxis#2',[-50 50],'colormap#2','redblue',...
			'data',Knmax{yy}, 'data',Knflow{yy}, 'figure',yy,'caxis#3,4',[0 2],'mask#all',mask,...
			'title#1','new velocity', 'title#2','new velocity - old velocity',...
			'title#3','maximum buttressing', 'title#4','buttressing number, flow calc')

		%ax = gcf;
		%exportgraphics(ax,['Plots/buttressing_calc_' num2str(years(yy)) '.pdf'])
	end

    f = gcf; % defines f to save plot
    exportgraphics(f,['Plots/2006_test.pdf']);   % actually saves the plot

end%}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ORIGINAL:
		% plotmodel(md,'data',md.inversion.vel_obs, 'data',md.inversion.vel_obs-velo,'caxis#2',[-50 50],'colormap#2','redblue',...
		% 	'data',Knmax{yy}, 'data',Knflow{yy}, 'figure',yy,'caxis#3,4',[0 2],'mask#all',mask,...
		% 	'title#1','new velocity', 'title#2','new velocity - old velocity',...
		% 	'title#3','maximum buttressing', 'title#4','buttressing number, flow calc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % PLOTTING A SINGLE YEAR ?:
        % plotmodel(md, ...
        %     'data', Knflow{yy}, 'figure',yy,'caxis',[0 2],'mask#all',mask, 'title', ['buttressing number ' num2str(years(yy))])
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNCOMMENT AND PASTE INTO TERMINAL FOR ALL PLOTS:

% plot_opts = cell(1,0);
% for yy = 1:length(years)
%     plot_opts = [plot_opts, { ...
%         'data',   Knflow{yy}, ...
%         'title',  ['Buttressing Number ' num2str(years(yy))]}];
% end
% 
% plotmodel(md, plot_opts{:}, ...
%           'nlines', 4, 'ncols', 4, ...
%           'mask#all', mask);
% 
% % Get all axes and set uniform caxis, xlim, ylim
% fig = gcf;
% allAxes = findall(fig, 'type', 'axes');
% for ax = allAxes'
%     caxis(ax, [0 2]);
%     xlim(ax, [998648.2458 1171935.9538]); % <-- Set your desired x-limits
%     ylim(ax, [-2148431.9935 -2033743.9935]); % <-- Set your desired y-limits
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP THREE

if perform(org, 'new_param_all_years'),% {{{

	years = [2001, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020]; % JEZ CHANGE THIS
	%years = [2015, 2016, 2020, 2022];
	md = loadmodel(org, 'Inversion_BHO');
	velo = md.inversion.vel_obs;

    

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

		%ax = gcf;
		%exportgraphics(ax,['Plots/buttressing_calc_' num2str(years(yy)) '.pdf'])
	end

% PLOTTING ALL THE YEARS IN ONE FIGURE:

plot_opts = cell(1,0);
for yy = 1:length(years)
    plot_opts = [plot_opts, { ...
        'data',   Knflow{yy}, ...
        'title',  ['Buttressing Number ' num2str(years(yy))]}];
end

plotmodel(md, plot_opts{:}, ...
          'nlines', 4, 'ncols', 4, ...
          'mask#all', mask);

% Get all axes and set uniform caxis, xlim, ylim
fig = gcf;
allAxes = findall(fig, 'type', 'axes');
for ax = allAxes'
    caxis(ax, [0 2]);
    xlim(ax, [998648.2458 1171935.9538]); % <-- Set your desired x-limits
    ylim(ax, [-2148431.9935 -2033743.9935]); % <-- Set your desired y-limits
end


end%}}}










