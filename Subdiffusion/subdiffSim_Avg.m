%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 9/24/2024
    Last update date: 10/9/2024

Description: This code simulation the diffusion following fractional
Brownian Motion (fBM), this allows subdiffusion simulation

This version examine the motion blur (averaging) effect solely. 

To speed up the calculation, I first generated tracks in a structure, and
then performed calculation using cellfun without a for loop

There're two fitting methods to choose from: 1) linear fitting in log
scale. 2) non-linear fitting in log scale

fBM core: 10.1103/PhysRevE.110.014105 (Eq.1)
---------------------------------------------------------------------------
%}

clear, clc, close all

simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)
strain = 'sim';

% imaging parameters
nTracks = 20000; 
nFrames = 100;      % frane number: 100 frames
frameT = 100e-3;    % frame interval: 100 ms
% expT = 100e-3;      % exposure time: 100 ms
dt = 5e-3;         % simulation step time: 20 ms

% secondary parameters
nInterval = frameT/ dt;         % number of simulation steps per frame
nSteps = nFrames* nInterval;    % simulation steps, simulation track length
time = dt* nInterval*( 1:nFrames-1); % readout time points

% diffusion parameters
alpha = 0.4;    % subdiffusion exponent
D = 0.005;      % um^2/s
locErr = 0e-3;  % unit: um


fprintf( '~~~~ Simulation Starts ~~~~\n')
fprintf( '   D = %.1g,  dt = %d ms,  frameT = %d ms,  nFrame = %d,  locErr = %.1g nm\n\n', D, dt*1e3, frameT*1e3, nFrames, locErr*1e3)


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


% create tracks as tf structure
tf( nTracks, 1).traj = []; % initialization

for i = 1: nTracks
    tf(i).traj = L* randn( nSteps, 3); % subdiffusive tracks
    % tf(i).locE = locErr* randn( nFrames, 3); % unit: um
end
traj = { tf.traj}';
% locE = { tf.locE}';

fprintf( '~~~~  Time: %.2f s (%d tracks generated)  ~~~~\n\n', toc( tStart), nTracks)

%% Compute MSD & Plotting

close all

figure( 'Position', [1000 400 420 400])
colorList = get( gca,'colororder');
cList = repmat( colorList, [2, 1]);

fitR = 1: 3; % fitting region of MSD
dim = 3; % dimension of the system

theoFlag = true; % flag for plotting the theoretical MSD
linearFitFlag = false; % 1: linear fit, 0: non-linear fit

expTList = [0 20 50 100]* 1e-3; % exposure time: 100 ms

avgFlag = true( 1, length( expTList));
avgFlag( 1) = false;

fitResult = nan( length( expTList), 2);

for c = 1: length( expTList)

    expT = expTList( c);
    nAvg = expT/ dt;        % number of averaging steps per exposure
    
    fprintf( '    expT = %d ms, locErr^2 = %d nm^2, motion blur = %.2f nm^2\n',...
        expT*1e3, locErr^2* 1e6, 2*D*expT^alpha/ (1+alpha)/ (2+alpha)* 1e6)

    if logical( avgFlag( c))

        % divide simulation steps by each frame (nInterval steps) into the 3rd dimension
        % dimension: [x, y, z] - [3, nInterval, nFrames]
        tmp = cellfun( @(x) reshape( x', 3, nInterval, []), traj, 'UniformOutput', false);

        % take average of steps during the exposure time window 
        % (average of 1:nAvg in the nInterval steps), dimension: [x, y] - [nFrames, 3]
        trajAvg = cellfun( @(x) squeeze( mean( x(:,1:nAvg,:),2))', tmp, 'UniformOutput', false);

    else
        % no blur (no averaging & locErr), dimension: [x, y] - [nFrames, 3]
        trajAvg = cellfun( @(x) x( 1: nInterval: end, :), traj, 'UniformOutput', false);
    end

    % trajLoc = cellfun( @plus, trajAvg, locE, 'UniformOutput', false);
    
    EnsMSD = cell2mat( cellfun( @(x) sum( (x(2:end,:)- x(1,:)).^2, 2), trajAvg, 'UniformOutput', false)')';
    
    eaMSD = mean( EnsMSD, 1);
    time = dt* nInterval*( 1:nFrames-1); % readout time points

    if logical( avgFlag( c))
        note = sprintf( 'expT=%dms', expT*1e3);
    else
        note = 'no blur';
    end
    
    % plot MSD
    scatter( time, eaMSD, 40, 'LineWidth', 2, 'MarkerEdgeColor', colorList( c,:),...
        'MarkerEdgeAlpha', 0.5, 'HandleVisibility','off'), hold on


    % linear fit in log-log scale 
    f = polyfit( log( time( fitR)), log( eaMSD( fitR)), 1);
    alphaFit = f(1);    DFit = exp( f(2))/ (2*dim);
    
    if logical( theoFlag)
        % plot theoretical MSD
        tFit = linspace( time(1), time( round( nFrames/4)), 1000);
        MSDtheo = 2*dim* D* tFit.^ alpha - 4*dim* D* expT^alpha/ (( 1+alpha)*( 2+alpha));
        plot( tFit, MSDtheo, 'LineWidth', 2, 'color', colorList( c,:), 'DisplayName', sprintf( 'theory [%s]', note))

    elseif logical( linearFitFlag)
        % plot linear fit
        tFit = time( 1: round( end/10));  MSDFit =  2* dim* DFit* tFit.^ alphaFit;
        plot( tFit, MSDFit, 'LineWidth', 2, 'color', colorList( c,:), ...
            'DisplayName', sprintf( '\\alpha=%.2f, D=%.1e [%s]', alphaFit, DFit, note))
    else
        % non-linear fit using a comprehensive form: MSD = 6Dt^a + 6*locE^2 (problematic form for locE, overestimate)
        fun = fittype( 'log( 6*a*(x)^b+6*c)');
        x0 = [ DFit, alphaFit, 0];
        xmin = [ 0, 0, -inf];
        xmax = [ inf, 2, inf];
        fitR2 = 1: 20;

        f = fit( time(fitR2)', log( eaMSD(fitR2))', fun, 'StartPoint', x0, 'Lower', xmin, 'Upper', xmax);
        DFit = f.a;  alphaFit = f.b;  locErrFit = sign( f.c)* sqrt( abs( f.c)); % unit: um

        tFit = linspace( time(1), time( round( nFrames/2)), 1000);
        MSDFit = 2*dim*( DFit* tFit.^ alphaFit + f.c);
        plot( tFit, MSDFit, 'LineWidth', 2, 'color', colorList( c,:), 'DisplayName', ...
            sprintf( '\\alpha=%.2f, D=%.1e, \\sigma=%.0fnm', alphaFit, DFit, locErrFit*1e3))
    end

    fitResult( c, :) = [alphaFit, DFit];
end

% plot setting
figure( gcf)
set( gca, 'FontSize', 12)
xlabel( 'Time (s)', 'FontSize', 14)
ylabel( 'EA-MSD (µm^2)', 'FontSize', 14)
legend( 'Location', 'northwest', 'FontSize', 10)

title( '3D Subdiffusion Simulation', 'FontSize', 13)
subtitle( sprintf( 'Input: \\alpha=%g, D=%.0g, \\sigma=0', alpha, D), 'FontSize', 13)
set( gca, 'Xscale', 'log', 'YScale', 'log')

simTxt = sprintf( '%d tracks\nframeT = %dms\ndt = %dms', nTracks, frameT*1e3, dt*1e3);
    
txt = text( 0.5, 0.2, simTxt, 'FontSize', 13, 'Units', 'normalized');
% txt = text( 0.1, 0.85, simTxt, 'FontSize', 13, 'Units', 'normalized');


