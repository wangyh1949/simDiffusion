%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 9/24/2024
    Last update date: 12/18/2025

Description: This code simulation the subdiffusion following fractional
Brownian Motion (fBM), with the effect of locErr & motion averaging.

This code plots mutliple MSD curves with different locErr. Motion blur can
be turned of via 'blurFlag' (optional). 

~~~~~~~~ compared to subdiffSim.m ~~~~~~~~


There're two fitting methods to choose from (see 'linearFitFlag')
    1) linear fitting in log scale
    2) non-linear fitting in log scale
There's also option to plot the theoretical MSD based on the diffusion
parameters via 'theoFlag'

fBM core: 10.1103/PhysRevE.110.014105 (Eq.1)
---------------------------------------------------------------------------
%}

clear, clc, close all

simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)
strain = 'sim';

% imaging parameters
nTracks = 20000; 
nFrames = 100;      % frane number: 100 frames
frameT = 100e-3;    % frame interval: 500 ms
expT = frameT;      % exposure time: 100 ms
dt = expT/ 10;      % simulation step time: 10 ms

% secondary parameters
nInterval = frameT/ dt;         % number of simulation steps per frame
nAvg = expT/ dt;                % number of averaging steps per exposure
nSteps = nFrames* nInterval;    % simulation steps, simulation track length
time = frameT*( 1:nFrames-1);   % readout time points

% diffusion parameters
alpha = 0.4;        % subdiffusion exponent
D = 0.005;          % um^2/s
locErr = 20e-3;     % 20 nm, unit: um

dim = 3; % dimension of the system


fprintf( '~~~~ Simulation Starts ~~~~\n')
fprintf( '   D = %.1g,  frameT = %d ms,  nFrame = %d,  locErr = %.1g nm\n', D, frameT*1e3, nFrames, locErr*1e3)
fprintf( '   dt = %d ms,  nAvg = %d steps,  overall movement: %.2f um\n', dt*1e3, nAvg, sqrt( 2*D*frameT*nFrames))
fprintf( '   locErr^2 = %d nm^2, motion blur = %.2f nm^2\n\n', locErr^2* 1e6, 2*D*expT^alpha/ (1+alpha)/ (2+alpha)* 1e6)

%% Subdiffusion simulation (fBM) & generate tracks

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


% generate tracks as a tf structure
tf( nTracks, 1).traj = []; % initialization

for i = 1: nTracks
    tf(i).traj = L* randn( nSteps, dim); % subdiffusive tracks
    % tf(i).locE = locErr* randn( nSteps, dim); % unit: um
end
traj = { tf.traj}';
% locE = { tf.locE}';

fprintf( '~~~~  Time: %.2f s (%d tracks generated)  ~~~~\n\n', toc( tStart), nTracks)

%% Add locErr & MSD Plotting

close all

figure( 'Position', [1000 400 420 400])
colorList = get( gca,'colororder');     colorList = repmat( colorList, [2, 1]);

theoFlag = false; % flag for plotting the theoretical MSD
linearFitFlag = false; % 1: linear fit, 0: non-linear fit

fitR = 1: 3;    % linear fitting range
fitR2 = 1: 40;  % non-linear fitting range

% list for locErr & flag for motion blurring
locErrList = [ 0 20 30 40 50]* 1e-3; % unit: um
blurFlag = true( 1, length( locErrList));

% locErrList = [ 0 20 40 50 0]* 1e-3; % unit: um
% blurFlag = [ true( 1, length( locErrList)) false];

fitResult = nan( length( locErrList), 2); % [alphaFit, DFit] for each condition

for c = 1: length( locErrList)
    
    % generate locErr for each frame
    locErr = locErrList( c);        
    for i = 1: nTracks        
        tf(i).locE = locErr* randn( nFrames, dim);
    end    
    locE = { tf.locE}'; % unit: um

    % apply motion blurring
    if logical( blurFlag( c))

        % divide simulation steps by each frame (nInterval steps) into the 3rd dimension
        % dimension: [x, y, z] - [dim, nInterval, nFrames]
        tmp = cellfun( @(x) reshape( x', dim, nInterval, []), traj, 'UniformOutput', false);

        % take average of steps during the exposure time window 
        % (average of 1:nAvg in the nInterval steps), dimension: [x, y] - [nFrames, 3]
        trajAvg = cellfun( @(x) squeeze( mean( x(:,1:nAvg,:),2))', tmp, 'UniformOutput', false);
        trajLoc = cellfun( @plus, trajAvg, locE, 'UniformOutput', false);

        note = sprintf( '\\sigma=%dnm', locErr*1e3);
    else
        % no blur (no averaging & locErr)
        trajLoc = cellfun( @(x) x( 1: nInterval: end, :), traj, 'UniformOutput', false);
        note = 'no blur';
    end
    
    EnsMSD = cell2mat( cellfun( @(x) sum( (x(2:end,:)- x(1,:)).^2, 2), trajLoc, 'UniformOutput', false)')';
    eaMSD = mean( EnsMSD, 1);
    
    % plot MSD
    scatter( time, eaMSD, 40, 'LineWidth', 2, 'MarkerEdgeColor', colorList( c,:),...
        'MarkerEdgeAlpha', 0.5, 'HandleVisibility','off'), hold on


    % linear fit in log-log scale 
    f = polyfit( log( time( fitR)), log( eaMSD( fitR)), 1);
    alphaFit = f(1);    DFit = exp( f(2))/ (2*dim);
    
    if logical( theoFlag)
        % plot theoretical MSD
        tFit = linspace( time(1), time( round( nFrames/4)), 1000);
        MSDtheo = 2*dim* D* tFit.^ alpha + 2*dim* locErr^2 - 4*dim* D* expT^alpha/ (( 1+alpha)* (2+alpha));

        if ~logical( blurFlag( c))
            MSDtheo = 2*dim* D* tFit.^alpha;
        end

        plot( tFit, MSDtheo, 'LineWidth', 2, 'color', colorList( c,:), 'DisplayName', ...
            sprintf( 'theory, %s', note))
    elseif logical( linearFitFlag)
        % plot linear fit
        tFit = time( 1: round( end/10));  MSDFit =  2* dim* DFit* tFit.^ alphaFit;
        plot( tFit, MSDFit, 'LineWidth', 2, 'color', colorList( c,:), ...
            'DisplayName', sprintf( '\\alpha=%.2f, D=%.1e [%s]', alphaFit, DFit, note))
    else
        % non-linear fit using a comprehensive form: MSD = 6Dt^a + 6*locE^2 - 12D*dt^a/(1+a)(2+a) 
        fun = fittype( 'log( a*(x)^b+c)');
        x0 = [ DFit, alphaFit, 0];
        xmin = [ 0, 0, -inf];
        xmax = [ inf, 2, inf];

        f = fit( time(fitR2)', log( eaMSD(fitR2))', fun, 'StartPoint', x0, 'Lower', xmin, 'Upper', xmax);
        DFit = f.a/ (2*dim);  alphaFit = f.b;  %locErrFit = sign( f.c)* sqrt( abs( f.c)); % unit: um

        motionBlur = 4*dim* DFit* expT^alphaFit/ (( 1+alphaFit)* (2+alphaFit));
        loc2 = ( f.c + motionBlur)/ (2*dim);
        locErrFit = sign( loc2)* sqrt( abs( loc2)); % unit: um

        % tFit = linspace( time(1), time( round( nFrames/2)), 1000);
        tFit = linspace( time(1), time( fitR2(end)), 1000);
        MSDFit = 2*dim*( DFit* tFit.^ alphaFit) + f.c;
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
legend( 'Location', 'southeast', 'FontSize', 10)

title( '3D Subdiffusion Simulation', 'FontSize', 13)
subtitle( sprintf( 'Input: \\alpha=%g, D=%.0g, \\sigma varies', alpha, D), 'FontSize', 13)
set( gca, 'Xscale', 'log', 'YScale', 'log')

% simTxt = sprintf( '%d tracks\nframeT = %dms\nexpT = %dms\ndt = %dms', nTracks, frameT*1e3, expT*1e3, dt*1e3);
simTxt = sprintf( '%d tracks\nexpT = %dms\ndt = %dms', nTracks, expT*1e3, dt*1e3);
    
% txt = text( 0.5, 0.2, simTxt, 'FontSize', 13, 'Units', 'normalized');
txt = text( 0.1, 0.8, simTxt, 'FontSize', 13, 'Units', 'normalized');


