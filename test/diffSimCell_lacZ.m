%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 12/5/2023
    last updated date: 9/4/2025

    ------- this script is adapted from diffSimCell_Blur.m -------

Description: this script simulates the diffusion of a free particle in a
confined cell shape in 3D with locError & motion blur

======Input=======
nTracks, nFrames, frameT, dt, D, locErr

======Output=======
tracksFinal stucture with fields: traj & MSD
    EnsMSD, EnsTAMSD

MSD fitting parameter: 
    diff, locErr, dalpha, alpha (unit adjusted for plotting)

it simulate the diffusion of lacZ protein in a cell, to test whether
confinment could change normal diffusion into subdiffusion with alpha=0.6

for LacZ (experiments)
    D = 0.17,  alpha=0.61,  locErr = 50 nm

simulation data that could reproduce the data
    D = 0.85,  alpha = 1,   locErr = 65 nm, cellWid = 0.6,  cellLen = 2.2
-------------------------------------------------------------
%}

%% 0. Set up basic diffusion info

clear, clc, close all

% simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)
strain = 'sim';

% Tracking parameters
nTracks = 10; 
nFrames = 100;
frameT = 20e-3; % frame interval, continuous exposure
dt = frameT/5;     % 10 ms, unit: s (simulation time)

if nTracks <= 100, plotTrackFlag = true; 
else, plotTrackFlag = false; end

% Diffusion parameters
D = 0.15;       % um^2/s
locErr = 65e-3;  % 40 nm, unit: um

% Averaging parameters (motion blur)
nAvg = frameT/ dt; % number of averaging frames (exposure)
nSteps = nFrames* nAvg; % simulation steps


% Cell geometry
cellWid = 0.6;    cellLength = 2.2; % cell long axis: y & short axis: x-z, unit: um
% cellWid = 1;    cellLength = 1;
r = cellWid/2;    l = ( cellLength - cellWid)/ 2;

fprintf( '~~~~ Simulation Starts ~~~~\n')
fprintf( '   D = %.2f, dt = %d ms, step = %.2f um, locErr = %.2f um\n', D, dt*1e3, sqrt( 4*D*dt), locErr)
fprintf( '   nFrame = %d, total step = %.2f um\n', nFrames, sqrt( 4*D*dt*nFrames))
fprintf( '   Average frames = %d, total simulation steps = %d\n', nAvg, nSteps)


%% 1. generate randomly distributed origin points inside the cell
    
rng(0) % to make the result repeatable

	% first generate 2*nTracks and then discard ones that are outside 
    ori = (rand( nTracks*4, 3)- 0.5).* [ cellWid, cellLength, cellWid];
    [ x, y, z] = deal( ori(:,1), ori(:,2), ori(:,3));

    inCyl = abs( y) < l & ( x.^2 + z.^2) < r^2; % in the cylinder region
    inCap = abs( y) >= l & ( x.^2 + z.^2 + ( y- l*sign(y)).^2) < r^2; % in the cap region

    ori = ori( inCyl | inCap, :);
    if length( ori) < nTracks
        warning( ' ~~~~~ Not Enough Origin Points !!! ~~~~~')
    end
    
    
%% 2. Diffusion Simulation

clearvars tracksFinal
tracksFinal( nTracks, 1).traj = []; % create tf structure

tic
for i = 1: nTracks

    pos = ori( i, :); % initial position, should be homogeneously generated in cell
    traj = nan( nSteps, 3);    traj( 1, :) = pos;

    % generate random steps, std = sqrt(2Ddt)
    % jumps = normrnd( 0, sqrt( 2*D*dt), [nSteps-1, 3]); % unit: um
    jumps = sqrt( 2*D*dt)* randn( nSteps-1, 3); % unit: um

    % after each step, check whether it goes outside the cell
    for j = 1: nSteps-1
        pos = cellShape( pos+ jumps(j, :), r, l); % corrects pos
        traj( j+1, :) = pos;
    end
    
    % time-averaging to get output traj by frame interval
    trajAvg = nan( nFrames, 3);
    for k = 1: 3
        trajAvg( :,k) = mean( reshape( traj(:,k), nAvg, []), 1);
    end
    
    locE = locErr* randn( nFrames, 3); % unit: um
    trajLoc = trajAvg + locE;
    
    % save track X & Y coordinates into tracksFinal structure
    tracksFinal(i).traj = trajLoc; % x & y coordiate, unit: um

    if plotTrackFlag % plot tracks
        plot3( trajLoc(:,1), trajLoc(:,2), trajLoc(:,3), 'LineWidth', 1.5), hold on        
    end
        
    % ~~~~ Calculate Single Steps ~~~~
    steps3D = sqrt( sum(( trajLoc( 2:end, :)- trajLoc( 1:end-1, :)).^2, 2)); % unit: um
    steps = sqrt( sum(( trajLoc( 2:end, 1:2)- trajLoc( 1:end-1, 1:2)).^2, 2)); % unit: um
    
    tracksFinal(i).steps3D = steps3D'; % unit: um
    tracksFinal(i).steps = steps'; % unit: um
end

fprintf( '\n~~~~~~  Simulation Done  ~~~~~~\n')
toc

tracksLength = cellfun( @length, {tracksFinal.traj}');
steps  = [ tracksFinal.steps]'; % unit: um
steps3D = [ tracksFinal.steps3D]'; % unit: um


% Plot Setting
if plotTrackFlag
    figure( gcf)
    set( gcf, 'Position', [450 450 720 400])
    set( gca, 'FontSize', 12)
    xlabel( 'X (\mum)', 'FontSize', 14)
    ylabel( 'Y (\mum)', 'FontSize', 14)
    zlabel( 'Z (\mum)', 'FontSize', 14)
    % legend( 'FontSize', 14)
    title( 'Simulated Tracks In Cell', 'FontSize', 16)
    axis image
    view(-60, 20) % view(az,el)  az=-37, el = 30
end


%% 3. Diffusion Analysis


maxT = nFrames; % MSD truncation length ( min( nFrames, maxT))

timeStep = frameT; % unit: s    
fprintf( '\n   ~~~~ frame time is  %d ms ~~~~\n', timeStep*1e3)

% initialzation
nTracks = size( tracksFinal, 1);
EnsMSD = nan( nTracks, maxT-1);     EnsTAMSD = nan( nTracks, maxT-1);
alpha = nan( nTracks, 1);           Dalpha = nan( nTracks, 1); 

tic
for i = 1: nTracks
    
    traj = tracksFinal(i).traj(:, 1:2); % unit: um
    
    % ~~~~ Ensemble-Averaged MSD ~~~~~
    dr = traj( 2: min( nFrames, maxT), :) - traj( 1, :);
    EnsMSD( i, 1: size( dr, 1)) = sum( dr.^2, 2); % ensemble MSD for EA-MSD plotting        

    % ~~~~ Time-Averaged MSD ~~~~~
    for tau = 1: min( nFrames, maxT)-1
        dr = traj( tau+1:end, :) - traj( 1:end-tau, :);
        dr2 = sum( dr.^2, 2);
        EnsTAMSD(i, tau) = mean( dr2, 'omitnan'); % ensemble TA-MSD for EATA-MSD plotting
    end    

    % ~~~~ Linear fitting of TA-MSD ~~~~  MSD = 4Dt^alpha
    fitRange = 1: floor( (nFrames-1) / 4); fitTxt = '25% fit';
    % fitRange = 1:5; fitTxt = '1:5 fit';
        tFit = fitRange* timeStep;
        f = polyfit( log( tFit), log( EnsTAMSD( i, fitRange)), 1);
        alpha(i) = f(1);    Dalpha(i) = exp( f(2))/4;
end

fprintf( '\n~~~~~~  Analysis Done  ~~~~~~\n')
toc

% get matrix variables from the structure
steps  = [ tracksFinal.steps]'; % unit: um


%% 4. MSD Plotting & Fitting

plotSim_MSDFit

% simJDPlot

fprintf( '\n~~~~~~ !!!  All Complete  !!! ~~~~~~\n')




%% Function

function posCorr = cellShape( pos, r, l)
% this function corrects position if it goes outside cell

    [ x, y, z] = deal( pos(1), pos(2), pos(3));
    
    inCyl = abs( y) < l & ( x^2 + z^2) < r^2; % in the cylinder region
    inCap = abs( y) >= l & ( x^2 + z^2 + ( y- l*sign(y))^2) < r^2; % in the cap region
    
    if inCyl || inCap % in Cell
        posCorr = pos;
        return
    elseif abs( y) < l % out of cylinder region
        xz = [ x, z]; % vector with respect to the axis of the cylinder
        xzCorr = xz* (2*r- norm( xz))/ norm( xz);
        posCorr = [ xzCorr(1) y xzCorr(2)];
        return
    else % out of the cap region
        radPos = [ x, y- l*sign(y), z]; % vector with respect to sphere center
        radPosCorr = radPos* (2*r - norm( radPos))/ norm( radPos);
        posCorr = radPosCorr + [ 0, l*sign(y), 0];
    end
end
