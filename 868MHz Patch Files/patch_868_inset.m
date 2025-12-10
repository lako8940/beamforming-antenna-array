%% Inset-Fed 868 MHz Patch Antenna
%
% Based on: Simple Patch Antenna Tutorial from openEMS
% Modified for:
%   - f0 = 868 MHz
%   - Patch width  = 63 mm (x-direction)
%   - Patch length = 80  mm (y-direction)
%   - Inset feed depth = 25 mm
%   - Microstrip feed line with 3 mm width, 1mm gap in notch
%
% Tested style:
%  - Matlab / Octave
%  - openEMS >= v0.0.35
%
% (C) adapted for your use

close all
clear
clc

%% Setup the Simulation
physical_constants;
unit = 1e-3; % all length in mm

%% Patch geometry
patch.width   = 63;    % x-direction (reduced from 63mm to lower input impedance)
patch.length  = 80;  % y-direction (tuned for 868 MHz with εᵣ=4.3)
patch.inset   = 25;    % inset depth along +y from bottom edge (decreased +2mm to lower input impedance)

%% Substrate (FR-4-ish)
% FR-4 εᵣ varies with frequency:
%   1 MHz:   ~4.7
%   1 GHz:   ~4.3-4.35
%   At 868 MHz: typically 4.3-4.35 (manufacturer dependent)
% Conservative choice: 4.3 (closer to measured values at UHF)
substrate.epsR      = 4.3;  % Changed from 4.4 to more accurate value at 868 MHz
tandelta            = 0.02;
f0                  = 868e6;             % center frequency
substrate.width     = patch.width  + 40; % some margin
substrate.length    = patch.length + 40;
substrate.thickness = 1.6;               % mm
substrate.cells     = 4;                 % mesh cells along thickness (reduced from 6)
substrate.kappa     = tandelta * 2*pi*f0 * EPS0 * substrate.epsR; % conductive loss

%% Feed (microstrip inset)
feed.R       = 50;  % lumped port impedance
feed.width   = 3; % microstrip width [mm] (50Ω line on 1.6mm FR4 with εr=4.3)
feed.gap     = 1.0; % gap to patch in notch [mm] (increased from 0.5 for better meshing)
feed.extra   = 20;  % line length extending outside patch [mm]

% Notch width (slot cut into patch)
notch.width = feed.width + 2*feed.gap;   % [mm]

%% Size of the simulation box (air box)
% Reduced from 400x400x300 to speed up simulation
% Still gives ~100mm margin around patch (which is 102x80mm + 40mm substrate margin)
SimBox = [300 300 200];   % [mm] x, y, z extent

%% Setup FDTD Parameter & Excitation Function
fc = 400e6;  % 20 dB corner frequency (bandwidth)
% Increased timesteps significantly - 868 MHz needs longer simulation
% Relaxed end criteria from -50dB to -30dB for practical convergence
FDTD = InitFDTD( 'NrTs', 500000, 'EndCriteria', 1e-3 );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% Setup CSXCAD Geometry & Mesh
CSX = InitCSX();

% initialize the mesh with the "air-box" dimensions
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)/3 SimBox(3)*2/3];

%% Create Substrate
CSX = AddMaterial( CSX, 'substrate' );
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa );

start = [-substrate.width/2  -substrate.length/2  0];
stop  = [ substrate.width/2   substrate.length/2  substrate.thickness];
CSX = AddBox( CSX, 'substrate', 0, start, stop );

% add extra cells to discretize the substrate thickness
mesh.z = [linspace(0,substrate.thickness,substrate.cells+1) mesh.z];

%% Create Ground (same size as substrate)
CSX = AddMetal( CSX, 'gnd' ); % perfect electric conductor (PEC)
start(3) = 0;
stop(3)  = 0;
CSX = AddBox( CSX, 'gnd', 10, start, stop );

%% Create Patch + Inset Notch + Microstrip Feed
CSX = AddMetal( CSX, 'patch' ); % PEC for patch + feed

z_patch = substrate.thickness;

% Patch coordinates (centered at origin)
x_left   = -patch.width/2;
x_right  =  patch.width/2;
y_bottom = -patch.length/2;
y_top    =  patch.length/2;

% Notch coordinates (cut from bottom edge into patch by patch.inset)
x_notch_left  = -notch.width/2;
x_notch_right =  notch.width/2;
y_notch_start = y_bottom;
y_notch_end   = y_bottom + patch.inset;

% 1) Top patch section above the notch
start = [x_left  y_notch_end  z_patch];
stop  = [x_right y_top        z_patch];
CSX = AddBox( CSX, 'patch', 10, start, stop );

% 2) Left bottom wing
start = [x_left        y_bottom  z_patch];
stop  = [x_notch_left  y_notch_end z_patch];
CSX = AddBox( CSX, 'patch', 10, start, stop );

% 3) Right bottom wing
start = [x_notch_right y_bottom  z_patch];
stop  = [x_right       y_notch_end z_patch];
CSX = AddBox( CSX, 'patch', 10, start, stop );

% 4) Microstrip feed line (same metal 'patch', so it merges electrically)
y_feed_end   = y_notch_end;               % up to the bottom of top patch
y_feed_start = y_bottom - feed.extra;     % extend beyond substrate edge

x_feed_left  = -feed.width/2;
x_feed_right =  feed.width/2;

start = [x_feed_left  y_feed_start z_patch];
stop  = [x_feed_right y_feed_end   z_patch];
CSX = AddBox( CSX, 'patch', 10, start, stop );

%% Apply the Excitation & Port as a Lumped Source
% Place the lumped port along the feed line,
% vertically from ground (z=0) to patch/feed metal (z=substrate.thickness).

feed.pos_x = 0;                                % centered with feed line
feed.pos_y = y_bottom - feed.extra/2;          % middle of extended feed line

% Port spans from ground to top of substrate
port_start = [feed.pos_x feed.pos_y 0];
port_stop  = [feed.pos_x feed.pos_y substrate.thickness];

% Debug: Display port location
disp(['Port location: x=' num2str(feed.pos_x) ', y=' num2str(feed.pos_y) ', z=0 to ' num2str(substrate.thickness)]);

[CSX port] = AddLumpedPort(CSX, 5, 1, feed.R, port_start, port_stop, [0 0 1], true);

%% Finalize the Mesh
% detect all edges except of the patch
mesh = DetectEdges(CSX, mesh, 'ExcludeProperty', 'patch');

% detect and set a special 2D metal edge mesh for the patch (including feed)
% Relaxed from /50 to /30 to avoid extremely fine mesh at metal edges
mesh = DetectEdges(CSX, mesh, 'SetProperty', 'patch', ...
                   '2D_Metal_Edge_Res', c0/(f0+fc)/unit/30);

% generate a smooth mesh with max. cell size: lambda_min / 15 (coarser for faster sim)
% At 868 MHz, lambda ~= 345 mm, so max cell ~= 23 mm
max_res = c0/(f0+fc)/unit/15;  % Relaxed from /20 to /15 for coarser mesh
mesh = SmoothMesh(mesh, max_res);

CSX = DefineRectGrid(CSX, unit, mesh);

%% Field dump (optional example)
CSX = AddDump(CSX, 'Hf', 'DumpType', 11, 'Frequency', [f0]);
CSX = AddBox(CSX, 'Hf', 10, ...
             [-substrate.width -substrate.length -10*substrate.thickness], ...
             [ substrate.width  substrate.length  10*substrate.thickness]);

%% Time-domain field dumps for animation
% DumpType 0 = E-field, DumpType 1 = H-field, DumpType 2 = current density
% SubSampling: [10,10,10] reduces file size by taking every 10th timestep
% OptResolution: mesh resolution for the dump (use coarser for smaller files)

% XY-plane cut through the patch (at substrate top surface)
CSX = AddDump(CSX, 'Et_xy', 'DumpType', 0, 'SubSampling', '10,10,10');
CSX = AddBox(CSX, 'Et_xy', 10, ...
             [-substrate.width/2, -substrate.length/2, substrate.thickness], ...
             [ substrate.width/2,  substrate.length/2, substrate.thickness]);

% XZ-plane cut through the feed line (y=0)
CSX = AddDump(CSX, 'Et_xz', 'DumpType', 0, 'SubSampling', '10,10,10');
CSX = AddBox(CSX, 'Et_xz', 10, ...
             [-substrate.width/2, 0, 0], ...
             [ substrate.width/2, 0, SimBox(3)*2/3]);

% YZ-plane cut through the center (x=0)
CSX = AddDump(CSX, 'Et_yz', 'DumpType', 0, 'SubSampling', '10,10,10');
CSX = AddBox(CSX, 'Et_yz', 10, ...
             [0, -substrate.length/2, 0], ...
             [0,  substrate.length/2, SimBox(3)*2/3]);

%% NF2FF box; size is 3 cells away from MUR boundary
start = [mesh.x(4)     mesh.y(4)     mesh.z(4)];
stop  = [mesh.x(end-3) mesh.y(end-3) mesh.z(end-3)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

%% Prepare and Run Simulation
Sim_Path = 'tmp_Patch_868_Inset';
Sim_CSX  = 'patch_868_inset.xml';

% create an empty working directory
[status, message, messageid] = rmdir( Sim_Path, 's' ); %#ok<ASGLU>
[status, message, messageid] = mkdir( Sim_Path );       %#ok<ASGLU>

% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX );

%% Postprocessing & Plots
freq = linspace( max([f0-fc, 400e6]), f0+fc, 501 );

% Debug: Check port structure before calcPort
disp('Port structure before calcPort:');
disp(port);

% Port files are HDF5 format (no extension) - just call calcPort directly
port = calcPort(port, Sim_Path, freq);

% Debug: Check port structure after calcPort
disp('Port structure after calcPort:');
disp(['  port is array: ' num2str(length(port)) ' elements']);
disp(['  port class: ' class(port)]);
if isfield(port, 'uf')
    disp('  port.uf exists');
    if isfield(port.uf, 'ref')
        disp(['    port.uf.ref size: ' num2str(size(port.uf.ref))]);
    end
else
    disp('  ERROR: port.uf does not exist!');
end

if isempty(port) || ~isfield(port, 'uf')
    error('calcPort failed. Port structure is invalid.');
end

% Check if port needs to be accessed as port{1} (cell array) instead of port
disp('Attempting to plot reflection coefficient...');

%% Smith chart port reflection
try
    plotRefl(port, 'threshold', -10)
catch err
    disp('plotRefl failed with error:');
    disp(err.message);
    disp('Trying alternative approach without plotRefl...');
    % Skip Smith chart plot for now
end
title( 'reflection coefficient' );

% plot feed point impedance
Zin = port.uf.tot ./ port.if.tot;
figure
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
s11 = port.uf.ref ./ port.uf.inc;
figure
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}| [dB]' );

drawnow

%% NFFF Plots
% find resonance frequency from s11
[~, f_res_ind] = min(abs(s11));   % robust index
f_res = freq(f_res_ind);

disp(['Resonant frequency: ' num2str(f_res/1e6) ' MHz']);
disp(['S11 at resonance: ' num2str(20*log10(abs(s11(f_res_ind)))) ' dB']);

% calculate the far field at phi=0 degrees and at phi=90 degrees
disp( 'calculating far field at phi=[0 90] deg...' );

try
    nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, ...
                      (-180:2:180)*pi/180, [0 90]*pi/180);
catch err
    disp('Far-field calculation failed with error:');
    disp(err.message);
    disp('Skipping far-field plots. S11 and impedance data are still valid.');
    return;  % Exit script but keep plots
end

% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt'] );
disp( ['directivity: Dmax = ' num2str(nf2ff.Dmax) ...
       ' (' num2str(10*log10(nf2ff.Dmax)) ' dBi)'] );
disp( ['efficiency: nu_rad = ' ...
        num2str(100*nf2ff.Prad./port.P_inc(f_res_ind)) ' %'] );

% normalized directivity as polar plot
figure
polarFF(nf2ff,'xaxis','theta','param',[1 2],'normalize',1)

% log-scale directivity plot
figure
plotFFdB(nf2ff,'xaxis','theta','param',[1 2])

drawnow

% Show 3D pattern
disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
thetaRange = (0:2:180);
phiRange   = (0:2:360) - 180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, ...
                  thetaRange*pi/180, phiRange*pi/180, ...
                  'Verbose', 1, 'Outfile', '3D_Pattern.h5');

figure
plotFF3D(nf2ff,'logscale',-20);

E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'], E_far_normalized, ...
           thetaRange, phiRange, 'scale', 1e-3);
