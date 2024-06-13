%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 12/6/2023
    last updated date: 1/3/2024

Description: this script 'diffSimAnalysis.m' is adapted from 'diffAnalysis.m'
it calculates diffusion quantities of the simulated tracks

======Output=======
    EnsMSD & EnsTAMSD, both in 3D and 2D

MSD fitting parameter: 
    diff, locErr

-------------------------------------------------------------
%}

tic



%% MSD Analysis in 3D

nTracks = length( tracksFinal);

% initialization of the quantities with the correct size
EnsMSD3D = nan( nTracks, maxTau); 
% EnsTAMSD3D = nan( nTracks, maxTau);
% tracksFinal( nTracks).MSD = [];

for i = 1: nTracks
    
    traj = tracksFinal(i).traj; % x & y coordiate, unit: um    
    
    % ~~~~ Ensemble-Averaged MSD ~~~~~
    % calculate first maxTau steps, assign NaN for shorter tracks that
    % doesn't have these steps, use mean(MSD,'omitnan') when plotting EAMSD
    MSD = nan( 1, maxTau);
    cutT = min( nFrames, maxTau+1); % cutoff of time maxTau+1, faster computation
    dr = traj( 2:cutT,:) - traj( 1,:); % displacement from origin with respect to time
    MSD( 1: cutT-1) = sum( dr.^2, 2); % squared displacement (t-1 timelags)
    EnsMSD3D( i,:) = MSD; % store it in ensemble MSD matrix (fixed length) for EA-MSD plotting    
    
%     % ~~~~ Time-Averaged MSD ~~~~~
%     % each track has a TA-MSD with 'nFrame-1' data points
%     taMSD = nan( 1, maxTau);
%     for tau = 1: nFrames-1 % time gap (step size) from 1 to nFrames-1
%         dr = traj( tau+1: end, :) - traj( 1: end-tau, :); % displacment of pairs with this timelag
%         dr2 = sum( dr.^2, 2);
%         taMSD( tau) = mean( dr2, 'omitnan'); % sum over all pairs with this timelag tau
%     end
%     
%     EnsTAMSD3D( i,:) = taMSD( 1: maxTau); % store it in ensemble time-averaged MSD (fixed length) for EATA_MSD plotting
%     tracksFinal(i).MSD = taMSD( 1: nFrames-1); % stored time-averaged MSD data into tracksFinal.MSD
end

%% MSD Analysis in 2D

tic
% initialization of the quantities with the correct size
EnsMSD = nan( nTracks, maxTau); 
% EnsTAMSD = nan( nTracks, maxTau);
% tracksFinal( nTracks).MSD = [];

for i = 1: nTracks
    
    traj = tracksFinal(i).traj(:, 1:2); % x & y coordiate, unit: um    
    
    % ~~~~ Ensemble-Averaged MSD ~~~~~
    % calculate first maxTau steps, assign NaN for shorter tracks that
    % doesn't have these steps, use mean(MSD,'omitnan') when plotting EAMSD
    MSD = nan( 1, maxTau);
    cutT = min( nFrames, maxTau+1); % cutoff of time maxTau+1, faster computation
    dr = traj( 2:cutT,:) - traj( 1,:); % displacement from origin with respect to time
    MSD( 1: cutT-1) = sum( dr.^2, 2); % squared displacement (t-1 timelags)
    EnsMSD( i,:) = MSD; % store it in ensemble MSD matrix (fixed length) for EA-MSD plotting
    
%     % ~~~~ Time-Averaged MSD ~~~~~
%     % each track has a TA-MSD with 'nFrame-1' data points
%     taMSD = nan( 1, maxTau);
%     for tau = 1: nFrames-1 % time gap (step size) from 1 to nFrames-1
%         dr = traj( tau+1: end, :) - traj( 1: end-tau, :); % displacment of pairs with this timelag
%         dr2 = sum( dr.^2, 2);
%         taMSD( tau) = mean( dr2, 'omitnan'); % sum over all pairs with this timelag tau
%     end
%     
%     EnsTAMSD( i,:) = taMSD( 1: maxTau); % store it in ensemble time-averaged MSD (fixed length) for EATA_MSD plotting
end

fprintf( '\n~~~~~~  Analysis Done  ~~~~~~\n')
toc

%% fitting the MSD & plot

% close all
% 
% eaMSD3D = mean( EnsMSD3D); % EA-MSD, unit: um^2
% % eataMSD3D = mean( EnsTAMSD3D); % EATA-MSD, unit: um^2
% 
%     % fit the MSD of individual tracks MSD = 4Dt + 4*locErr^2
%     fitRange = 1: 3; % points of MSD used for fitting
%     p = polyfit( dt* fitRange, eaMSD3D( fitRange), 1); % linear fit, y = ax + b
% 
%     dif3D = p(1)/ 6; % may be negative, 3D case
%     locErr3D = sqrt( p(2)/ 6); % unit: um, may be imaginary (should exclude later)
% 
%     time = (1: maxTau) * dt;
%     t = (1:20)* dt; MSDFit = p(1)*t + p(2);
% 
%     figure
%     plot( time, eaMSD3D, 'o', 'LineWidth', 1, 'DisplayName', sprintf( 'EA-MSD, D=%.2f', D)), hold on
%     % plot( time, eataMSD3D, 'o', 'LineWidth', 1, 'DisplayName', 'EATA-MSD')
%     plot( t, MSDFit, 'r', 'LineWidth', 2, 'DisplayName',...
%         sprintf( 'D*=%.3f, \\sigma=%.0fnm', dif3D, locErr3D*1e3)) % , [strain extraName]
% 
%     set( figure(gcf), 'Position', [800 400 350 320])
%     set( gca, 'FontSize', 12)
%     xlabel( 'Time (s)', 'FontSize', 14)
%     ylabel( 'MSD (µm^2)', 'FontSize', 14)
%     title( 'Simulation Data 3D', 'FontSize', 14) %
%     legend( 'Location', 'northwest', 'FontSize', 11)

%% 2D projected MSD

eaMSD = mean( EnsMSD); % EA-MSD, unit: um^2
% eataMSD = mean( EnsTAMSD); % EATA-MSD, unit: um^2

    % fit the MSD of individual tracks MSD = 4Dt + 4*locErr^2
    fitRange = 1: 3; % points of MSD used for fitting
    p = polyfit( dt* fitRange, eaMSD( fitRange), 1); % linear fit, y = ax + b

    dif = p(1)/ 4; % may be negative, 3D case
    locErr = sqrt( p(2)/ 4); % unit: um, may be imaginary (should exclude later)

    time = (1: maxTau) * frameT;
    t = (1:20)* dt; MSDFit = p(1)*t + p(2);

%     figure
    
    scatter( time, eaMSD, 30, 'LineWidth', 2, 'MarkerEdgeAlpha', 0.5,... % 'MarkerEdgeColor', colorList(c,:),...
        'DisplayName', sprintf( 'EA-MSD, D=%.2f', D)), hold on
    
    fitRange = 1:3; fitTxt = '1:3 fit'; 
    f = polyfit( log( time( fitRange)), log( eaMSD( fitRange)), 1);
    alphaFit = f(1);    Dapp = exp( f(2))/4;
    t = time( 1: round( end/6));    MSDFit =  4* Dapp* t.^ alphaFit;
    
    % plot( time, eataMSD, 'o', 'LineWidth', 1, 'DisplayName', 'EATA-MSD')
%     plot( t, MSDFit, 'r', 'LineWidth', 2, 'DisplayName',...
%         sprintf( 'D*=%.3f, \\sigma=%.0fnm', dif, locErr*1e3)) % , [strain extraName]

    set( figure(gcf), 'Position', [1160 400 350 320])
    set( gca, 'FontSize', 12)
    xlabel( 'Time (s)', 'FontSize', 14)
    ylabel( 'MSD (µm^2)', 'FontSize', 14)
    title( 'Simulation Data 2D', 'FontSize', 14) %
    legend( 'Location', 'northwest', 'FontSize', 11)

    set( gca, 'Xscale', 'log', 'YScale', 'log')
