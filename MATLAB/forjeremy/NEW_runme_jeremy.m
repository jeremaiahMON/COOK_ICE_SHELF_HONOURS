steps = [1]; %change this if you want to run different steps

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
if perform(org, 'Reparam'),% {{{

	years = [2015, 2016, 2020];   % CHANGE THE YEARS AS I LIKE
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
		exportgraphics(f,['Plots/buttressing_calc_' num2str(years(yy)) '.pdf'])

		%if want to save output per iteration of the for loop
		Knflow_temp = Knflow{yy};
		save([indir 'buttressing_outfiles.mat'], 'varnames');
	end

	%if want to save all output
	save([indir 'buttressing_alloutfiles.mat'], 'Knflow', 'Nflow', 'Knmax', 'Nmax');


end%}}}

