%% To be used in conjunction with 'WAOM_3031.ipynb'
% Calculating / Prepping Neutral Densities

% ENSURE that you have downloaded these codes: https://www.teos-10.org/preteos10_software/neutral_density.html 

% Path to the files from the link above ----------^^^
addpath('C:\Users\jjls0\OneDrive\Documents\HONOURS\PROGRAMS\zOTHER DATASETS\Calculation - modes of melt\eos80_legacy_gamma_n');
addpath('C:\Users\jjls0\OneDrive\Documents\HONOURS\PROGRAMS\zOTHER DATASETS\Calculation - modes of melt\eos80_legacy_gamma_n\library');

% CHANGE the step number depending on what cross section you want to
% calculate the neutral densities for.
steps = [3]; 


org = organizer('steps', steps);

%	1: Cross Section A to A'
if perform(org, 'Section_1_Neutral_Density'),% {{{

    load('gamma_n_inputs_section_1.mat')

    % Checking that all the sizes are the same 3D structure
    disp(size(SP));
    disp(size(t));
    disp(size(p));
    
    long = 152.5;
    lat = -68.5;
    
    disp(size(long));
    disp(size(lat));

    [gamma_n] = eos80_legacy_gamma_n(SP, t, p, long, lat);

    save('gamma_n_outputs_section_1.mat', 'gamma_n')

end%}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	2: Cross Section B to B'
if perform(org, 'Section_2_Neutral_Density'),% {{{

    load('gamma_n_inputs_section_2.mat')

    % Checking that all the sizes are the same 3D structure
    disp(size(SP));
    disp(size(t));
    disp(size(p));
    
    long = 152.5;
    lat = -68.5;
    
    disp(size(long));
    disp(size(lat));

    [gamma_n] = eos80_legacy_gamma_n(SP, t, p, long, lat);

    save('gamma_n_outputs_section_2.mat', 'gamma_n')

end%}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	3: Cross Section C to C'
if perform(org, 'Section_3_Neutral_Density'),% {{{

    load('gamma_n_inputs_section_3.mat')

    % Checking that all the sizes are the same 3D structure
    disp(size(SP));
    disp(size(t));
    disp(size(p));
    
    long = 152.5;
    lat = -68.5;
    
    disp(size(long));
    disp(size(lat));

    [gamma_n] = eos80_legacy_gamma_n(SP, t, p, long, lat);

    save('gamma_n_outputs_section_3.mat', 'gamma_n')

end%}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






