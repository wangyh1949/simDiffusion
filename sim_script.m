%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 4/15/2024
    last updated date: 7/7/2025

Description: this script simulates the diffusion in cell with different
input parameters, and then compare the MSD by plotting them together.

This script call the simCell function

-------------------------------------------------------------
%}

clear, clc, close all

% track parameters
nTracks = 500; 
nFrames = 100;
frameT = 20e-3; % frame interval
% frameT = 100e-3; % frame interval
expT = frameT;  % continuous exposure
dt = expT/20;   % 10 ms, unit: s (simulation timestep)

% Averaging parameters (motion blur)
nInterval = frameT/dt;          % number of simulation steps per frame
nSteps = nFrames* nInterval;    % total simulation steps
nAvg = expT/dt;                 % number of averaging steps (exposure)
time = (1: nFrames-1)* frameT;  % readout time points

if nTracks <= 100, plotTrackFlag = true; 
else, plotTrackFlag = false; end

% Diffusion parameters
D = 0.05;       % um^2/s
locErr = 0e-2;  % 40 nm, unit: um

% Cell geometry
% cellWid = 0.6;    cellLength = 2; % cell long axis: y & short axis: x-z
cellWid = 1;    cellLength = 3; % cell long axis: y & short axis: x-z, unit: um

fitFlag = 2; % 0: no fit, 1: linear fit, 2: non-linear fit


%% Run Simulation & Plot MSD
close all

fitR = 1:3; fitTxt = '1:3 fit';
dim = 3; % dimension of the system

Dlist = [ 0.001, 0.01, 0.1, 1];
% Dlist = [ 0.0003, 0.001, 0.003, 0.01, 0.03];

figure( 'Position', [1000 400 420 400])
% colorList = get( gca,'colororder');     colorList = repmat( colorList, [2, 1]);

colorList = flip( winter( length( Dlist)+1)); 

tStart = tic;
for c = 1: length( Dlist)

    D = Dlist( c);
    EnsMSD = nan( nTracks, nFrames-1);
    EnsTAMSD = nan( nTracks, nFrames-1);
    
    % 3D simulation in cell
    tf = simCell( nTracks, nFrames, frameT, expT, dt, D, locErr, cellWid, cellLength);
    
    note = sprintf( 'D=%g', D);

    for i = 1: nTracks
        traj = tf(i).traj;

        % calculate MSD
        dr = traj( 2:end, :) - traj( 1, :);
        MSD = sum( dr.^2, 2); % squared displacement (t-1 timelags)
        EnsMSD( i, :) = MSD; % store it in ensemble MSD matrix (fixed length) for EA-MSD plotting
        
        % calculate TA-MSD
        taMSD = nan( 1, nFrames-1);
        for tau = 1: nFrames-1
            dr = traj( tau+1: end, :) - traj( 1: end-tau, :); % displacment of pairs with this timelag
            dr2 = sum( dr.^2, 2);
            taMSD( tau) = mean( dr2, 'omitnan'); % sum over all pairs with this timelag tau
        end
        EnsTAMSD( i, :) = taMSD;
    end
    eaMSD = mean( EnsMSD, 1);
    eataMSD = mean( EnsTAMSD, 1);
    
    if fitFlag == 0
        % plot EATA-MSD
        scatter( time, eataMSD, 30, 'LineWidth', 2, 'MarkerEdgeColor', colorList( c,:),...
            'MarkerEdgeAlpha', 0.5, 'DisplayName', sprintf( '%s', note)), hold on
        fitTxt = '';
    else

        % plot EATA-MSD
        scatter( time, eataMSD, 30, 'LineWidth', 2, 'MarkerEdgeColor', colorList( c,:),...
            'MarkerEdgeAlpha', 0.5, 'HandleVisibility','off'), hold on
    
        % linear fit in log-log scale 
        f = polyfit( log( time( fitR)), log( eataMSD( fitR)), 1);
        alphaFit = f(1);    DFit = exp( f(2))/ (2*dim);    
        
        if fitFlag == 1
            % plot linear fit
            tFit = time( 1: round( end/10));  
            tFit = time( 1:10);
            MSDFit =  2* dim* DFit* tFit.^ alphaFit;
            plot( tFit, MSDFit, 'LineWidth', 2, 'color', colorList( c,:), ...
                'DisplayName', sprintf( '\\alpha=%.2f, D=%.2g [%s]', alphaFit, DFit, note))
        else
            % non-linear fit using a comprehensive form: MSD = 6Dt^a + 6*locE^2 - 12D*dt^a/(1+a)(2+a) 
            fun = fittype( 'log( 6*a*(x)^b+c)');
            x0 = [ DFit, alphaFit, 0];
            xmin = [ 0, 0, -inf];
            xmax = [ inf, 2, inf];

            % xmin = [ 0, 0, -1];   xmax = [ 10, 2, 1];
            fitR2 = 1: 20; fitTxt = '1:20 fit';
            % fitR2 = 1: floor( nFrames/4); fitTxt = sprintf( '1:%d fit', max( fitR2));
    
            f = fit( time(fitR2)', log( eataMSD(fitR2))', fun, 'StartPoint', x0, 'Lower', xmin, 'Upper', xmax);
            DFit = f.a;  alphaFit = f.b;  %locErrFit = sign( f.c)* sqrt( abs( f.c)); % unit: um
    
            motionBlur = 4*dim* DFit* expT^alphaFit/ (( 1+alphaFit)* (2+alphaFit));
            % locErrFit = sqrt( ( f.c + motionBlur)/ (2*dim)); % unit: um
            loc2 = ( f.c + motionBlur)/ (2*dim);
            locErrFit = sign( loc2)* sqrt( abs( loc2)); % unit: um
    
            % tFit = linspace( time(1), time( round( nFrames/2)), 1000);
            tFit = linspace( time(1), time( fitR2( end)), 1000);
            MSDFit = 2*dim*( DFit* tFit.^ alphaFit) + f.c;
            plot( tFit, MSDFit, 'LineWidth', 2, 'color', colorList( c,:), 'DisplayName', ...
                sprintf( '\\alpha=%.2f, D=%.1e, \\sigma=%.0fnm', alphaFit, DFit, locErrFit*1e3))
        end

        % MSDtheo = 2* dim* Dapp* tFit.^ alphaFit;
        fitResult( c, :) = [alphaFit, DFit];
    end

    fprintf( '~~~~  Time: %.2f s (%d tracks, %d steps)  ~~~~\n\n', toc( tStart), nTracks, nSteps)
end


%% plot setting
figure( gcf)
set( gca, 'FontSize', 14)
xlabel( 'Time (s)', 'FontSize', 14)
ylabel( 'EATA-MSD (µm^2)', 'FontSize', 14)
legend( 'Location', 'northwest', 'FontSize', 10)

title( '3D Brownian Diffusion in Cell', 'FontSize', 14)
% subtitle( sprintf( 'Input: D=%.0g, \\sigma varies', D), 'FontSize', 13)
subtitle( sprintf( 'Input: \\sigma=%.0g, D varies  %s', locErr, fitTxt), 'FontSize', 13)
set( gca, 'Xscale', 'log', 'YScale', 'log')

simTxt = sprintf( '%d tracks\nframeT = %dms\nsim dt = %dms', nTracks, frameT*1e3, dt*1e3);
% txt = text( 0.55, 0.15, simTxt, 'FontSize', 13, 'Units', 'normalized');
legend( 'Location', 'southeast', 'FontSize', 10)
txt = text( 0.05, 0.85, simTxt, 'FontSize', 12, 'Units', 'normalized');

%%
figure( gcf)
ylim( [1e-5 20]), xlim( [1e-2 3])
% ylim( [1e-4 10]), xlim( [8e-2 100])
% xlim( [0.05 20]), ylim( [1e-4, 20])
% ylim( [1e-3 2]), xlim( [1 200])

% ylim( [1e-5 200]), xlim( [1e-2 30]) % for 1000 frame, 20ms
% xlim auto, ylim auto



