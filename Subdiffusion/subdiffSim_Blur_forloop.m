%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 9/24/2024
    Last update date: 

Description: This code simulation the diffusion following fractional
Brownian Motion (fBM), this allows subdiffusion simulation

I added the effect of locErr & motion blur compared to subdiffSim.m

Old: I used for loop to genearte tracks, instead of stroing in a structure

There're two fitting methods to choose from: 1) linear fitting in log
scale. 2) non-linear fitting in log scale
fBM core: 10.1103/PhysRevE.110.014105 (Eq.1)
---------------------------------------------------------------------------
%}

clear, clc, close all

simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)
strain = 'sim';

% imaging parameters
nTracks = 10000; 
nFrames = 100;      % frane number: 100 frames
frameT = 100e-3;    % frame interval: 500 ms
expT = 100e-3;      % exposure time: 100 ms
dt = expT/ 5;       % simulation step time: 20 ms

% secondary parameters
nInterval = frameT/ dt;         % number of simulation steps per frame
nAvg = expT/ dt;                % number of averaging steps per exposure
nSteps = nFrames* nInterval;    % simulation steps, simulation track length
time = dt* nInterval*( 1:nFrames-1); % readout time points

% diffusion parameters
alpha = 0.4;        % subdiffusion exponent
D = 0.005;          % um^2/s
locErr = 20e-3;     % 20 nm, unit: um


fprintf( '~~~~ Simulation Starts ~~~~\n')
fprintf( '   D = %.1g,  frameT = %d ms,  nFrame = %d,  locErr = %.1g nm\n', D, frameT*1e3, nFrames, locErr*1e3)
fprintf( '   dt = %d ms,  nAvg = %d steps,  overall movement: %.2f um\n', dt*1e3, nAvg, sqrt( 2*D*frameT*nFrames))
fprintf( '   locErr^2 = %d nm^2, motion blur = %.2f nm^2\n\n', locErr^2* 1e6, 2*D*expT^alpha/ (1+alpha)/ (2+alpha)* 1e6)

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
colorList = get( gca,'colororder');     colorList = repmat( colorList, [2, 1]);

fitR = 1: 3; % fitting region of MSD
dim = 3; % dimension of the system

linearFitFlag = 1; % 1: linear fit, 2: non-linear fit

locErrList = [ 0 20 40 70 0]* 1e-3;
% locErrList = 20e-3;

avgFlag = true( 1, length( locErrList));
avgFlag( end) = false;

fitResult = nan( length( locErrList), 2);

for c = 1: length( locErrList)

    locErr = locErrList( c);
    EnsMSD = nan( nTracks, nFrames-1);
    
    tStart = tic;
    for i = 1: nTracks
        
        traj = L* randn( nSteps, 3); % subdiffusive tracks
        locE = sqrt( 2*locErr^2)* randn( nSteps, 3); % unit: um
        
        trajLoc = traj + locE;
        
        if logical( avgFlag( c))
            % motion blur: average position of every nAvg frames 
            tmp = reshape( trajLoc', 3, nInterval, []); % divide by each frame
            trajAvg = squeeze( mean( tmp( :, 1:nAvg, :), 2))'; % average over the exposure
        else
            trajAvg = trajLoc( 1: nInterval: end, :);
        end
            
        dr = trajAvg( 2:end, :) - trajAvg( 1, :);

        MSD = sum( dr.^2, 2); % squared displacement (t-1 timelags)
        EnsMSD( i, :) = MSD; % store it in ensemble MSD matrix (fixed length) for EA-MSD plotting
    end

    eaMSD = mean( EnsMSD, 1);
    time = (1: nFrames-1)* dt*nInterval; % readout time points

    if logical( avgFlag( c))
        note = sprintf( '\\sigma=%dnm', locErr*1e3);
    else
        note = 'no avg';
    end
    
    % plot MSD
    scatter( time, eaMSD, 40, 'LineWidth', 2, 'MarkerEdgeColor', colorList( c,:),...
        'MarkerEdgeAlpha', 0.5, 'HandleVisibility','off'), hold on


    % linear fit in log-log scale 
    f = polyfit( log( time( fitR)), log( eaMSD( fitR)), 1);
    alphaFit = f(1);    DFit = exp( f(2))/ (2*dim);
    
    if logical( linearFitFlag)
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

    % MSDtheo = 2* dim* Dapp* tFit.^ alphaFit;
    fitResult( c, :) = [alphaFit, DFit];
    fprintf( '~~~~  Time: %.2f s (%d tracks, %d steps)  ~~~~\n', toc( tStart), nTracks, nSteps)
end

% plot setting
figure( gcf)
set( gca, 'FontSize', 12)
xlabel( 'Time (s)', 'FontSize', 14)
ylabel( 'EA-MSD (µm^2)', 'FontSize', 14)
legend( 'Location', 'northwest', 'FontSize', 10)

title( '3D Subdiffusion Simulation', 'FontSize', 13)
subtitle( sprintf( 'Input: \\alpha=%g, D=%.0g, \\sigma varies', alpha, D), 'FontSize', 13)
set( gca, 'Xscale', 'log', 'YScale', 'log')

simTxt = sprintf( '%d tracks\nframeT = %dms\ndt = %dms', nTracks, frameT*1e3, dt*1e3);
    
txt = text( 0.5, 0.2, simTxt, 'FontSize', 13, 'Units', 'normalized');
% txt = text( 0.1, 0.85, simTxt, 'FontSize', 13, 'Units', 'normalized');


