%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 12/5/2023
    last updated date: 7/7/2025

    ------- this script is adapted from diffSim3D.m -------

Description: this script simulates the diffusion of a free particle in a
confined cell shape in 3D with locError

======Input=======
nTracks, nFrames, frameT, dt, D, locErr

======Output=======
tracksFinal stucture with fields: traj & MSD
    EnsMSD, EnsTAMSD

MSD fitting parameter: 
    diff, locErr, dalpha, alpha (unit adjusted for plotting)

======Random Function=======
rand: Uniformly distributed random numbers in (0,1)
randn: Normally distributed random numbers: mean 0, variance 1
        f(x) = exp( -x^2/2)/ sqrt( 2*pi)

normrnd: normal random numbers (Statistics and Machine Learning Toolbox)
    normrnd( mu, sigma) = mu + sigma* randn()
-------------------------------------------------------------
%}

clear, clc, close all

% simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simParth)
strain = 'sim';


% Tracking parameters
nTracks = 100;  nFrames = 50;   frameT = 10e-3;
plotTrackFlag = 1;

% Diffusion parameters
D = 0.2;        % um^2/s
dt = frameT;    % simulation time step: 10 ms
locErr = 0e-3;  % 40 nm, unit: um

% Cell geometry
cellWid = 0.6;    cellLength = 2; % cell long axis: y & short axis: x-z
r = cellWid/2;    l = ( cellLength - cellWid)/ 2;

fprintf( '~~~~ Simulation Starts ~~~~\n')
fprintf( '   D = %.2f, dt = %d ms, step = %.2f um, locErr = %.2f um\n', D, dt*1e3, sqrt( 4*D*dt), locErr)
fprintf( '   nFrame = %d, total step = %.2f um\n', nFrames, sqrt( 4*D*dt*nFrames))


%% 1. generate randomly distributed origin points on the membrane
    
rng(0) % to make the result repeatable

    % generate randomly distributed origin points
    ori = ( rand( nTracks*2, 3)- 0.5).* [ cellWid, cellLength, cellWid];    
    
    [ x, y, z] = deal( ori(:,1), ori(:,2), ori(:,3));
    
    % exclude points falling outside cells
    inCyl = abs( y) < l & ( x.^2 + z.^2) < r^2; % in the cylinder region
    inCap = abs( y) >= l & ( x.^2 + z.^2 + ( y- l*sign(y)).^2) < r^2; % in the cap region

    ori = ori( inCyl | inCap, :);
    if length( ori) < nTracks
        warning( ' ~~~~~ Not Enough Origin Points !!! ~~~~~')
    end
    
    
%% 2. Diffusion Simulation

clearvars tracksFinal
tracksFinal( nTracks, 1).traj = [];

tic
for i = 1: nTracks

    pos = ori( i, :); % initial position, should be homogeneously distributed in cell
    traj = nan( nFrames, 3);    traj( 1, :) = pos;

    % generate random steps, std = sqrt(2Ddt)
    % jumps = normrnd( 0, sqrt( 2*D*dt), [nFrames-1, 3]); % unit: um
    jumps = sqrt( 2*D*dt)* randn( nFrames-1, 3); % unit: um

    % after each step, check whether it's outside the cell
    for j = 1: nFrames-1
        pos = cellShape( pos+ jumps(j, :), r, l); % correct for cell geometry
        traj( j+1, :) = pos;
    end

    locE = locErr* randn( nFrames, 3); % unit: um
    traj = traj + locE;
    
    % save track X & Y coordinates into tracksFinal structure
    tracksFinal(i).traj = traj; % x & y coordiate, unit: um

    if plotTrackFlag % plot tracks
        plot3( traj(:,1), traj(:,2), traj(:,3), 'LineWidth', 1.5),        hold on        
    end
    
    traj = traj( :, 1:2);
    
    % ~~~~ Calculate Single Steps ~~~~
    steps3D = sqrt( sum(( traj( 2:end, :)- traj( 1:end-1, :)).^2, 2)); % unit: um
    steps = sqrt( sum(( traj( 2:end, 1:2)- traj( 1:end-1, 1:2)).^2, 2)); % unit: um    
    
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
    set( gcf, 'Position', [900, 450, 720, 400])
    set( gca, 'FontSize', 12)
    xlabel( 'X (\mum)', 'FontSize', 14)
    ylabel( 'Y (\mum)', 'FontSize', 14)
    zlabel( 'Z (\mum)', 'FontSize', 14)
    % legend( 'FontSize', 14)
    title( 'Simulated Tracks In Cell', 'FontSize', 16)
    axis image
    view(-60, 20) % view(az,el)  az=-37, el = 30
end


%% Diffusion Analysis

% diffSimAnalysis
% 
% simJDPlot

%     simName = sprintf( 'Sim %s_%d frames_D %.2f_LocE %d.mat', Date, nFrames, D, locError*1e3);
%     save( [simPath 'Data\' simName])

fprintf( '\n~~~~~~ !!!  All Complete  !!! ~~~~~~\n')



%% Function

function posCorr = cellShape( pos, r, l)
% this function corrects the position if the particle diffuse outside cell
% boundary by mirror reflection, note it only works for steps much smaller
% compared to cell geometry, it fails for large jumps

    % global r l
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
