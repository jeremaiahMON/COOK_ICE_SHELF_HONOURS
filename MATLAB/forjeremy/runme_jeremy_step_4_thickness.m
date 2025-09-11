steps = [4]; %change this if you want to run different steps

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

	years = [2001, 2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
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
if perform(org, 'Buttressing_with_gapfill'),% {{{

	years = [2006:1:2020];
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
	end

	%if want to save all output
	save([indir 'buttressing_alloutfiles.mat'], 'Knflow', 'Nflow', 'Knmax', 'Nmax','masks');


end%}}}
if perform(org, 'Buttressing_without_gapfill'),% {{{

	years = [2006:1:2020];
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
		plotmodel(md,'data',Knmax{yy},'figure', yy,'xlim#all',[0.9e6, 1.2e6],'ylim#all',[-2.2e6, -2e6],'colormap#1','redblue','caxis',[-5 5],'mask#all',mask,'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp'])
		exportgraphics(gcf,['Plots/buttressing_calc_' num2str(years(yy)) '_nogapfill.png'])
	end

	%if want to save all output
	save([indir 'buttressing_alloutfiles.mat'], 'Knflow', 'Nflow', 'Knmax', 'Nmax','masks');

end%}}}
if perform(org, 'Thickness_without_gapfill'),% {{{

	years = [2001, 2006, 2008, 2009, 2010, 2011, 2014, 2015, 2016, 2017, 2018];
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
		plotmodel(md,'data',H,'figure', yy,'xlim#all',[0.9e6, 1.2e6],'ylim#all',[-2.2e6, -2e6],'caxis',[300 1000],'mask#all',mask,'expdisp',[indir 'ice_extent_' num2str(years(yy)) '.exp'])
		exportgraphics(gcf,['Plots/thickness_' num2str(years(yy)) '_nogapfill.png'])
	end

	%if want to save all output
	% save([indir 'thickness_alloutfiles.mat'], 'H','masks');

end%}}}

