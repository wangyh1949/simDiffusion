%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 12/5/2023
    last updated date: 12/6/2023

Description: this script 'diffSimCell' is adapted from 'diffSim3D'
1) simulates the diffusion of a free particle in a confined cell shape in
3D

======Input=======
nTracks, nFrames, D, dt, locError 

======Output=======
tracksFinal stucture with fields: traj & MSD
    EnsMSD, EnsTAMSD

MSD fitting parameter: 
    diff, locErr, dalpha, alpha (unit adjusted for plotting)

-------------------------------------------------------------
%}

clear, clc, close all

simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)
strain = 'sim';

plotTrackFlag = 1;

nTracks = 100;     nFrames = 50;   maxTau = nFrames - 1;

D = 0.2; % um^2/s
dt = 1e-2; % 10 ms, unit: s

locError = 0e-2; % 40 nm, unit: um

cellWid = 0.6;    cellLength = 2; % cell long axis: y & short axis: x-z
global r l
r = cellWid/2;    l = ( cellLength - cellWid)/ 2;

fprintf( '~~~~ Simulation Starts ~~~~\n')
fprintf( '   D = %.2f, dt = %d ms, step = %.2f um, locErr = %.2f um\n', D, dt*1e3, sqrt( 4*D*dt), locError)
fprintf( '   nFrame = %d, total step = %.2f um\n', nFrames, sqrt( 4*D*dt*nFrames))


%% 1. generate randomly distributed origin points on the membrane
    
rng(0) % to make the result repeatable

    % generate randomly distributed origin points
    ori = (rand( nTracks*2, 3)- 0.5).* [ cellWid, cellLength, cellWid];    
    
    [ x, y, z] = deal( ori(:,1), ori(:,2), ori(:,3));

    inCyl = abs( y) < l & ( x.^2 + z.^2) < r^2; % in the cylinder region
    inCap = abs( y) >= l & ( x.^2 + z.^2 + ( y- l*sign(y)).^2) < r^2; % in the cap region

    ori = ori( inCyl | inCap, :);
    if length( ori) < nTracks
        warning( ' ~~~~~ Not Enough Origin Points !!! ~~~~~')
    end
    
    
%% 2. Diffusion Simulation

tracksFinal( nTracks, 1).traj = [];

tic
for i = 1: nTracks

    pos = ori( i, :); % initial position, should be homogeneously generated in cell
    traj = nan( nFrames, 3);    traj( 1, :) = pos;

    % generate random steps
    % std of the step distribution should be sqrt( 2Dt)
    jumps = normrnd( 0, sqrt( 2*D*dt), [nFrames-1, 3]); % unit: um

    % iterations, after each step judge whether it's outside the cell
    for j = 1: nFrames-1
        pos = cellShape( pos + jumps(j, :));
        traj( j+1, :) = pos;
    end

    locE = normrnd( 0, sqrt( 2*locError^2), [nFrames, 3]); % unit: um
    traj = traj + locE;
    
    % save track X & Y coordinates into tracksFinal structure
    tracksFinal(i).traj = traj; % x & y coordiate, unit: um

    if logical( plotTrackFlag) % plot tracks
        plot3( traj(:,1), traj(:,2), traj(:,3), 'LineWidth', 1.5),        hold on        
    end
    
    traj = traj( :, 1:2);
    
    % ~~~~ Calculate Single Steps ~~~~
    singleSteps = sqrt( sum(( traj( 2:end, :)- traj( 1:end-1, :)).^2, 2)); % unit: um    
    tracksFinal(i).steps = singleSteps'; % unit: um        
    
end

fprintf( '\n~~~~~~  Simulation Done  ~~~~~~\n')
toc

tracksLength = cellfun( @length, {tracksFinal.traj}');
steps  = [ tracksFinal.steps]'; % unit: um


%% Plot Setting
if logical( plotTrackFlag)
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

diffSimAnalysis

msd = EnsTAMSD(:,1); % MSD(1)

simJDPlot

%     simName = sprintf( 'Sim %s_%d frames_D %.2f_LocE %d.mat', Date, nFrames, D, locError*1e3);
%     save( [simPath 'Data\' simName])

fprintf( '\n~~~~~~ !!!  All Complete  !!! ~~~~~~\n')




%%

function posCorr = cellShape( pos)

    global r l
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
