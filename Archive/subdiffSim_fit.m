%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 9/25/2024
    Last update date: 9/26/2024

Description: This code simulation the diffusion following fractional
Brownian Motion (fBM)

this version: to test theoretical fitting using comprehensive equation
compare between linear vs non-linear fitting

---------------------------------------------------------------------------
%}

clear, clc, close all

simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)
strain = 'sim';

% imaging parameters
nTracks = 10000; 
nFrames = 100;      % frane number: 100 frames
frameT = 100e-3;    % frame interval: 100 ms
dt = 100e-3;         % simulation time: 10 ms

% diffusion parameters
alpha = 0.4;        % subdiffusion exponent
D = 0.005;          % um^2/s
locErr = 50e-3;   % 20 nm, unit: um

% secondary parameters
nAvg = frameT/dt;       % number of averaging frames
nSteps = nFrames* nAvg; % simulation steps

fprintf( '~~~~ Simulation Starts ~~~~\n')
fprintf( '   D = %.1e,  dt = %d ms,  step = %.2f um,  locErr = %.2f um\n', D, dt*1e3, sqrt( 2*D*dt), locErr)
fprintf( '   nFrame = %d, total step = %.2f um\n', nFrames, sqrt( 2*D*dt*nFrames))

%% Subdiffusion simulation (fBM)
tic

t = dt*( 1: nSteps);    covMat = zeros( nSteps-1);

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
EnsMSD = nan( nTracks, nSteps-1); 

for i = 1: nTracks

    traj = L* randn( nSteps, 3); % subdiffusive tracks

    locE = sqrt( 2*locErr^2)* randn( nSteps, 3); % unit: um    
    traj = traj + locE;
            
    dr = traj( 2:end, :) - traj( 1, :);
    
    MSD = sum( dr.^2, 2); % squared displacement (t-1 timelags)
    EnsMSD( i, :) = MSD; % store it in ensemble MSD matrix (fixed length) for EA-MSD plotting
end

    eaMSD = mean( EnsMSD, 1);
    time = dt*( 1:length( eaMSD));

toc

%% plotting & fitting

close all

figure( 'Position', [1000 400 420 400])

fitR = 1: 3;
fitR2 = 1: 20;
dim = 3; % dimension of the system

    % plot MSD
    scatter( time, eaMSD, 30, 'LineWidth', 2, ...
        'MarkerEdgeAlpha', 0.5, 'HandleVisibility','off'), hold on


    % linear fit in log-log scale 
    f = polyfit( log( time( fitR)), log( eaMSD( fitR)), 1);
    alphaFit = f(1);    Dapp = exp( f(2))/ (2*dim);
    

    % plot fitting 
    tFit = time( 1: round( end/4));  MSDFit =  2* dim* Dapp* tFit.^ alphaFit;
    plot( tFit, MSDFit, 'LineWidth', 2, 'DisplayName',...
        sprintf( '\\alpha=%.2f, D=%.1e', alphaFit, Dapp))
    

    % fit using a comprehensive form
    fun = fittype( 'log( 6*a*(x)^b+6*c^2)');
    x0 = [ Dapp, alphaFit, 0];
    xmin = [ 0, 0, 0];
    xmax = [ inf, 2, inf];

    f = fit( time(fitR2)', log( eaMSD(fitR2))', fun, 'StartPoint', x0, 'Lower', xmin, 'Upper', xmax);
    Dalpha = f.a;     alphaFit2 = f.b;     locErrFit = f.c; % unit: um
            
    t =  [ time( 1: round( end/2))];    MSDFit = 2*dim*( Dalpha* t.^ alphaFit2 + locErrFit^2);
    plot( t, MSDFit, 'LineWidth', 2, 'DisplayName', ...
        sprintf( '\\alpha=%.2f, D=%.1e, \\sigma=%.0fnm', alphaFit2, Dalpha, locErrFit*1e3))
        
    % MSDtheo = 2*dim*( D* t.^alpha + locErr^2);
    % plot( t, MSDtheo, 'LineWidth', 2, 'DisplayName', 'theory')


% plot setting
figure( gcf)
set( gca, 'FontSize', 12)
xlabel( 'Time (s)', 'FontSize', 14)
ylabel( 'EA-MSD (µm^2)', 'FontSize', 14)
legend( 'Location', 'northwest', 'FontSize', 11)
title( '3D subdiffusion Simulation', 'FontSize', 12)
subtitle( sprintf( 'Input: \\alpha=%g, D=%.0e, \\sigma=%.0fnm', alpha, D, locErr*1e3))
set( gca, 'Xscale', 'log', 'YScale', 'log')

simTxt = sprintf( '%d tracks\n  dt = %dms', nTracks, dt*1e3);
text( 0.6, 0.2, simTxt, 'FontSize', 13, 'Units', 'normalized')



