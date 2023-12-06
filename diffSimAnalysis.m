%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 12/6/2023
    last updated date: 12/6/2023

Description: this script 'diffSimAnalysis.m' is adapted from 'diffAnalysis.m'
it calculates diffusion quantities of the simulated tracks

======Output=======
tracksFinal stucture with fields: traj & MSD
    EnsMSD, EnsTAMSD

MSD fitting parameter: 
    diff, locErr, dalpha, alpha (unit adjusted for plotting)

-------------------------------------------------------------
%}

tic
% initialization of the quantities with the correct size
EnsMSD = nan( nTracks, maxTau); EnsTAMSD = nan( nTracks, maxTau);
Diff = nan( nTracks, 1);        LocErr = nan( nTracks, 1);
Dalpha = nan( nTracks, 1);      alpha = nan( nTracks, 1);

tracksFinal( nTracks).MSD = [];

for i = 1: nTracks
    
    traj = tracksFinal(i).traj; % x & y coordiate, unit: um
    
    
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
    p = polyfit( dt* fitRange, taMSD( fitRange), 1); % linear fit, y = ax + b
    
%     Diff(i) = p(1)/ 6; % may be negative, 3D case
%     LocErr(i) = sqrt( p(2)/ 6); % unit: um, may be imaginary (should exclude later)
    
    Diff(i) = p(1)/ 4; % may be negative, 2D case
    LocErr(i) = sqrt( p(2)/ 4); % unit: um, may be imaginary (should exclude later)
    
    % fit the log(MSD) to get Dalpha & alpha
    fitRange = 1: 3; % max(3, floor( nFrames/4)); % 1/4 of the track
    pLog = polyfit( log( dt* fitRange), log( taMSD( fitRange)), 1); % log(MSD) = alpha* log(t) + log(4D)
    alpha(i) = pLog(1);
    Dalpha(i) = exp( pLog(2))/4;
%     Dalpha(i) = exp( pLog(2))/6;
end

LocErr( LocErr ~= real( LocErr)) = nan; % delete the negative intercepts (imaginary)

locErr = LocErr* 1e3;       % unit: nm

fprintf( '\n~~~~~~  Analysis Done  ~~~~~~\n')
toc




