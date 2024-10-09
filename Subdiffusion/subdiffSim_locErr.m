%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 9/24/2024
    Last update date: 10/9/2024

Description: This code simulation the diffusion following fractional
Brownian Motion (fBM), this allows subdiffusion simulation

In this version, I considered locErr as blurring, but I didn't consider
effect of motion blurring (averaging). There're two fitting methods to
choose from: 1) linear fitting in log scale. 2) non-linear fitting in log
scale

fBM core: 10.1103/PhysRevE.110.014105 (Eq.1)

---------------------------------------------------------------------------
%}

clear, clc, close all

simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)
strain = 'sim';

% imaging parameters
nTracks = 50000;
nFrames = 100;      % frane number: 100 frames
frameT = 100e-3;    % frame interval: 100 ms
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
fprintf( '~~~~  Time: %.2f s (Covariance Matrix Created)  ~~~~\n\n', toc( tStart))


%% Create Tracks & MSD Plotting

close all

figure( 'Position', [1000 400 420 400])
colorList = get( gca,'colororder');
cList = repmat( colorList, [2, 1]);

fitR = 1: 3; % fitting region of MSD
dim = 3; % dimension of the system

theoFlag = true; % flag for plotting the theoretical MSD
linearFitFlag = false; % 1: linear fit, 0: non-linear fit

locErrList = [ 0 20 30 40]* 1e-3;
% locErrList = 20e-3;

fitResult = nan( length( locErrList), 2);

for c = 1: length( locErrList)

    locErr = locErrList( c);
    EnsMSD = nan( nTracks, nSteps-1);

    tStart = tic;
    for i = 1: nTracks

        traj = L* randn( nSteps, 3); % generate subdiffusive tracks
        locE = locErr* randn( nSteps, 3); % unit: um
        
        trajLoc = traj + locE;

        % calculate EA-MSD
        dr = trajLoc( 2:end, :) - trajLoc( 1, :);
        EnsMSD( i, :) = sum( dr.^2, 2); % squared displacement (t-1 timelags)
    end

    eaMSD = mean( EnsMSD, 1);
    time = dt*( 1: nSteps-1);

    note = sprintf( '\\sigma=%d', locErr*1e3);
    
    % plot MSD
    scatter( time, eaMSD, 40, 'LineWidth', 2, 'MarkerEdgeColor', colorList( c,:),...
        'MarkerEdgeAlpha', 0.5, 'HandleVisibility','off'), hold on

    % linear fit in log-log scale
    f = polyfit( log( time( fitR)), log( eaMSD( fitR)), 1);
    alphaFit = f(1);    DFit = exp( f(2))/ (2*dim);
    
    if logical( theoFlag)
        % plot theoretical MSD
        tFit = linspace( time(1), time( round( nFrames/4)), 1000);
        MSDtheo = 2*dim* D* tFit.^ alpha + 2*dim* locErr^2;
        plot( tFit, MSDtheo, 'LineWidth', 2, 'color', colorList( c,:), 'DisplayName', sprintf( 'theory [%s]', note))

    elseif logical( linearFitFlag)
        % plot linear fit
        tFit = time( 1: round( end/6));  MSDFit =  2* dim* DFit* tFit.^ alphaFit;
        plot( tFit, MSDFit, 'LineWidth', 2, 'color', colorList( c,:), ...
            'DisplayName', sprintf( '\\alpha=%.2f, D=%.1e [%s]', alphaFit, DFit, note))
    else
        % non-linear fit using a comprehensive form: MSD = 6Dt^a + 6*locE^2 (problematic form for locE, overestimate)
        fun = fittype( 'log( 6*a*(x)^b+6*c)');
        x0 = [ DFit, alphaFit, 0];
        xmin = [ 0, 0, 0];
        xmax = [ inf, 2, inf];
        fitR2 = 1: 20;

        f = fit( time(fitR2)', log( eaMSD(fitR2))', fun, 'StartPoint', x0, 'Lower', xmin, 'Upper', xmax);
        DFit = f.a;  alphaFit = f.b;  locErrFit = sqrt( f.c); % unit: um
        
        tFit = linspace( time(1), time( round( nFrames/2)), 1000);    
        MSDFit = 2*dim*( DFit* tFit.^ alphaFit + f.c);
        plot( tFit, MSDFit, 'LineWidth', 2, 'color', colorList( c,:), 'DisplayName', ...
            sprintf( '\\alpha=%.2f, D=%.1e, \\sigma=%.0f', alphaFit, DFit, locErrFit*1e3))
    end

    fitResult( c, :) = [alphaFit, DFit];
    fprintf( '~~~~  Time: %.2f s (%d tracks, %d steps)  ~~~~\n', toc( tStart), nTracks, nSteps)

end

% plot setting
figure( gcf)
set( gca, 'FontSize', 12)
xlabel( 'Time (s)', 'FontSize', 14)
ylabel( 'EA-MSD (µm^2)', 'FontSize', 14)
legend( 'Location', 'northwest', 'FontSize', 10)
legend( 'Location', 'southeast', 'FontSize', 10)
title( '3D Subdiffusion Simulation', 'FontSize', 13)
subtitle( sprintf( 'Input: \\alpha=%g, D=%.0g, \\sigma varies', alpha, D), 'FontSize', 13)
set( gca, 'Xscale', 'log', 'YScale', 'log')

simTxt = sprintf( '%d tracks\n  dt = %dms', nTracks, dt*1e3);
% text( 0.5, 0.2, simTxt, 'FontSize', 13, 'Units', 'normalized')
txt = text( 0.1, 0.85, simTxt, 'FontSize', 13, 'Units', 'normalized');


