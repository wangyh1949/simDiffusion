%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 9/24/2024
    Last update date: 9/29/2024

Description: This code simulation the diffusion following fractional
Brownian Motion (fBM), this allows subdiffusion simulation

fBM core: 10.1103/PhysRevE.110.014105 (Eq.1)

Blurring effects (locErr & motion blur) are not considered in this version

---------------------------------------------------------------------------
%}

clear, clc, close all

simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)
strain = 'sim';

% imaging parameters
nTracks = 50000; 
nFrames = 100;      % frane number: 100 frames
frameT = 10e-3;    % frame interval: 100 ms
dt = frameT;         % simulation time step: 10 ms

% diffusion parameters
alpha = 0.4;        % subdiffusion exponent
D = 0.005;          % um^2/s
locErr = 0e-3;   % 20 nm, unit: um

% secondary parameters
nAvg = frameT/dt;       % number of averaging frames
nSteps = nFrames* nAvg; % simulation steps

fprintf( '~~~~ Simulation Starts ~~~~\n')
fprintf( '   D = %.1g,  dt = %d ms,  nFrame = %d,  locErr = %.1g um\n\n', D, dt*1e3, nFrames, locErr)


%% Subdiffusion simulation (fBM)

t = dt*( 1: nSteps);    covMat = zeros( nSteps);

tStart = tic;

% creating the covariance matrix
for i = 1: nSteps
    for j = 1: i
        % symmetric matrix, only calculate lower half
        covMat( i, j) = t(i)^alpha + t(j)^alpha - abs( t(i)-t(j))^alpha;        
    end
end

% Cholesky decomposition
L = chol( D* covMat, 'lower');


% tf( nTracks, 1).traj = []; % create tf structure
% traj = wfbm( alpha/2, 100); % matlab function to generate fbm tracks

EnsMSD = nan( nTracks, nSteps-1); 

for i = 1: nTracks

    traj = L* randn( nSteps, 3); % generate subdiffusive tracks
    
    % calculate EA-MSD
    dr = traj( 2:end, :) - traj( 1, :); 
    EnsMSD( i, :) = sum( dr.^2, 2); % squared displacement (t-1 timelags)
end

fprintf( '~~~~  Time: %.2f s (%d tracks, %d steps)  ~~~~\n\n', toc( tStart), nTracks, nSteps)


%% plot MSD & Fit

figure( 'Position', [1000 400 420 400])

fitR = 1: 3; % fitting region of MSD
dim = 3; % dimension of the system

    eaMSD = mean( EnsMSD, 1);
    time = dt*( 1:length( eaMSD));

    % linear fit in log-log scale 
    f = polyfit( log( time( fitR)), log( eaMSD( fitR)), 1);
    alphaFit = f(1);    DFit = exp( f(2))/ (2*dim);
    
    scatter( time, eaMSD, 40, 'LineWidth', 2, ...
        'MarkerEdgeAlpha', 0.5, 'HandleVisibility','off'), hold on

    tFit = time( 1: round( end/2));  MSDFit =  2* dim* DFit* tFit.^ alphaFit;
    plot( tFit, MSDFit, 'LineWidth', 2, 'DisplayName', sprintf( 'Fit: \\alpha=%.2g, D=%.3g', alphaFit, DFit))


% plot setting
figure( gcf)
set( gca, 'FontSize', 12)
xlabel( 'Time (s)', 'FontSize', 14)
ylabel( 'EA-MSD (µm^2)', 'FontSize', 14)
legend( 'Location', 'northwest', 'FontSize', 12)
title( '3D Subdiffusion Simulation', 'FontSize', 13)
subtitle( sprintf( 'Input: \\alpha=%g, D=%.0g, \\sigma=%.0fnm', alpha, D, locErr*1e3), 'FontSize', 13)
set( gca, 'Xscale', 'log', 'YScale', 'log')

simTxt = sprintf( '%d tracks\n  dt = %dms', nTracks, dt*1e3);
text( 0.55, 0.2, simTxt, 'FontSize', 13, 'Units', 'normalized')



