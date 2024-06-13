%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 4/15/2024
    last updated date: 

Description: this script calls the diffusion simulation code and gets an
tracksFinal as output

-------------------------------------------------------------
%}

clear, clc, close all

% track parameters
nTracks = 100; 
nFrames = 100;
frameT = 10e-3; % frame interval
plotTrackFlag = 1;


% diffusion parameters
D = 0.05; % um^2/s
dt = 10e-3; % 10 ms, unit: s (simulation time)
locError = 0e-2; % 40 nm, unit: um


cellWid = 0.6;    cellLength = 2; % cell long axis: y & short axis: x-z
% cellWid = 1;    cellLength = 1; % cell long axis: y & short axis: x-z, unit: um

tf = simCell( nTracks, nFrames, frameT, D, dt, cellWid, cellLength);
figure( 'Position', [100, 400, 580, 400])

for i = 1: nTracks
    
    traj = tf(i).traj;
    
    if logical( plotTrackFlag) % plot tracks
        plot3( traj(:,1), traj(:,2), traj(:,3), 'LineWidth', 1.5), hold on        
    end
        
    % ~~~~ Calculate Single Steps ~~~~
    steps3D = sqrt( sum(( traj( 2:end, :)- traj( 1:end-1, :)).^2, 2)); % unit: um
    steps = sqrt( sum(( traj( 2:end, 1:2)- traj( 1:end-1, 1:2)).^2, 2)); % unit: um
    
    tf(i).steps3D = steps3D'; % unit: um
    tf(i).steps = steps'; % unit: um
    
end

steps  = [ tf.steps]'; % unit: um
steps3D = [ tf.steps3D]'; % unit: um


% Plot Setting
if logical( plotTrackFlag)
    figure( gcf)
%     set( gcf, 'Position', [100, 400, 600, 400])
    set( gca, 'FontSize', 12)
    xlabel( 'X (\mum)', 'FontSize', 14)
    ylabel( 'Y (\mum)', 'FontSize', 14)
    zlabel( 'Z (\mum)', 'FontSize', 14)
    title( 'Simulated Tracks In Cell', 'FontSize', 16)
    axis image
    view(-60, 20) % view(az,el)  az=-37, el = 30
end


%% Diffusion Analysis

% diffSimAnalysis

% msd = EnsTAMSD(:,1); % MSD(1)

% simJDPlot

%     simName = sprintf( 'Sim %s_%d frames_D %.2f_LocE %d.mat', Date, nFrames, D, locError*1e3);
%     save( [simPath 'Data\' simName])

% fprintf( '\n~~~~~~ !!!  All Complete  !!! ~~~~~~\n')

