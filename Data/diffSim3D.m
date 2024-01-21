%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 7/27/2023
    last updated @7/23/2023

Description: this script 'diffSim3D'
1) simulates the diffusion of a free particle in 3D space

======Input=======
nTracks, nFrames, D, timeStep, locError 

======Output=======
tracksFinal stucture with fields: traj & MSD
    EnsMSD, EnsTAMSD

MSD fitting parameter: 
    diff, locErr, dalpha, alpha (unit adjusted for plotting)

-------------------------------------------------------------
%}

clear
close all
clc

tic
nTracks = 1000; 
nFrames = 20;
maxTau = nFrames - 1;
Date = '230727';
simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)

D = 0.5; % um^2/s
timeStep = 20e-3; % 10 ms, unit: s
locError = 50e-3; % 50 nm, unit: um

for locError = (0:20:80)*1e-3
    
    tracksFinal( nTracks, 1).traj = [];
    tracksFinal( nTracks).MSD = []; % TA-MSD for each track

    % initialization of the quantities with the correct size
    EnsMSD = nan( nTracks, maxTau); EnsTAMSD = nan( nTracks, maxTau);
    Diff = nan( nTracks, 1);        LocErr = nan( nTracks, 1);
    Dalpha = nan( nTracks, 1);      alpha = nan( nTracks, 1);

    for i = 1: nTracks

        % generate a random step 
        % std of the step distribution should be sqrt( 2Dt)
        steps = normrnd( 0, sqrt( 2*D*timeStep), [nFrames-1, 2]); % unit: um
        locE = normrnd( 0, sqrt( 2*locError^2), [nFrames, 2]); % unit: um
        traj = ( [0 0; cumsum( steps)]+ locE)* 1e-6; % m

        % save track X & Y coordinates into tracksFinal structure
        tracksFinal(i).traj = traj; % x & y coordiate in the unit of m

        % ~~~~ Ensemble-Averaged MSD ~~~~~
        % calculate first maxTau steps, assign NaN for shorter tracks that
        % doesn't have these steps, use mean(MSD,'omitnan') when plotting EAMSD
        MSD = nan( 1, maxTau);
        cutT = min( nFrames, maxTau+1); % cutoff of time maxTau+1, faster computation
        dr = traj( 2:cutT,:) - traj( 1,:); % displacement from origin with respect to time
        MSD( 1: cutT-1) = sum( dr.^2, 2); % squared displacement (t-1 timelags)
        EnsMSD( i,:) = MSD; % store it in ensemble MSD matrix (fixed length) for EA-MSD plotting

        % ~~~~ Time-Averaged MSD ~~~~~
        % each track has a TA-MSD with 'nFrame-1' data points
        taMSD = nan( 1, maxTau);
        for tau = 1: nFrames-1 % time gap (step size) from 1 to nFrames-1
            dr = traj( tau+1: end, :) - traj( 1: end-tau, :); % displacment of pairs with this timelag
            dr2 = sum( dr.^2, 2);
            taMSD( tau) = mean( dr2, 'omitnan'); % sum over all pairs with this timelag tau
        end

        EnsTAMSD( i,:) = taMSD( 1: maxTau); % store it in ensemble time-averaged MSD (fixed length) for EATA_MSD plotting
        tracksFinal(i).MSD = taMSD( 1: nFrames-1); % stored time-averaged MSD data into tracksFinal.MSD

        % fit the TA-MSD of individual tracks MSD = 4Dt + 4*locErr^2
        fitRange = 1: 3; % points of MSD used for fitting
        p = polyfit( timeStep* fitRange, taMSD( fitRange), 1); % linear fit, y = ax + b
        Diff(i) = p(1)/ 4; % may be negative
        LocErr(i) = sqrt( p(2))/ 2; % unit: m, may be imaginary (should exclude later)    

        % fit the log(MSD) to get Dalpha & alpha
        fitRange = 1: 3; % max(3, floor( nFrames/4)); % 1/4 of the track
        pLog = polyfit( log( timeStep* fitRange), log( taMSD( fitRange)), 1); % log(MSD) = alpha* log(t) + log(4D)
        alpha(i) = pLog(1);
        Dalpha(i) = exp( pLog(2))/4;
    end

    LocErr( LocErr ~= real( LocErr)) = nan; % delete the negative intercepts (imaginary)

    EnsMSD = EnsMSD* 1e12;      % unit: um^2
    EnsTAMSD = EnsTAMSD* 1e12;  % unit: um^2
    diff = Diff* 1e12;          % unit: um^2/s
    dalpha = Dalpha* 1e12;      % unit: um^2/s
    locErr = LocErr* 1e9;       % unit: nm

    tracksLength = cellfun( @length, {tracksFinal.traj}');
    strain = 'Sim';
    extraName = sprintf( 'D: %.2f, Loc: %d', D, locError*1e3);

    simName = sprintf( 'Sim %s_%d frames_D %.2f_LocE %d.mat', Date, nFrames, D, locError*1e3);
    save( [simPath simName])

end

disp( '~~~~~~ !!!  All Complete  !!! ~~~~~~')
toc