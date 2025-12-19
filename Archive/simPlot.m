%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 7/9/2023
    Last updated at 7/21/2023

Description: this script plots the result of tracking data analysis
    1) track length distribution
    2) Diffusion Coefficient
    3) apparent D* coefficient & Dalpha & alpha & localization error (sigma)
       from individual track MSD fitting
    4) EATA-MSD & Fit
    5) EA-MSD vs EATA-MSD & Scaling
    6) localization error (sigma) distribution from     

7/16: plot multiple strains, add D comparison Exp vs Theory
---------------------------------------------------------------------------
%}

%% Setting up the parameters

clear
clc
close all

simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)

list = dir( 'Data\Sim*.mat');
plotNum = getPlotNum( list); % find which files you want to plot 

% colorList = get(gca,'colororder'); close % get the color code for plotting
colorList = get(gca,'colororder');  colorList = repmat( colorList, [2, 1]);  close

lineStyle = {'-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '--'};

count = 0;    locErrStat = [];    diffStat = [];    diffAll = [];   nameAll = {};

for j = plotNum
    
    count = count +1;
    
    % load each file for this strain
    matrixName = list( j).name;
    load( matrixName);   disp(['~~~ ' matrixName ' loaded ~~~'])        
        
    strainName = [strain extraName];

    condDiff = tracksLength >= 12; % long tracks criteria

    tracksLeft = sum( condDiff); 
    gap = round( tracksLeft/ 100); % prepare for CDF plotting

    %% tracksLength 1-CDF

    figure(1)
    [f,xgroup] = ecdf( tracksLength); xgroup = xgroup(2:end); f = f(1:end-1);
    plot( xgroup, 1-f, lineStyle{ count}, 'LineWidth', 2, 'DisplayName',...
        sprintf( '%s, mean=%2.0f, %d tracks', extraName, mean( tracksLength), nTracks))
    hold on

    %% D Coefficient CDF

    figure(2)
    dDiff = diff( condDiff); % unit: um^2/s
    [f,xgroup] = ecdf( dDiff);  xgroup = xgroup([ 2: gap: end-1 end]); f = f([ 2: gap: end-1 end]);
    plot(xgroup, f, lineStyle{ count}, 'LineWidth', 2.5, 'DisplayName',...
        sprintf( '\\langleD\\rangle=%.2f  [%s]', mean( dDiff), strainName)) % strainName
    hold on
    
    diffSet = dDiff( ~isnan( dDiff));
    m = bootstrp( length( diffSet), @mean, diffSet);
    diffSEM = std( diffSet)./ sqrt( length( diffSet));
    
    diffAll = [ diffAll; [ mean(m), std(m), diffSEM]];
    nameAll = [nameAll; extraName( 2:end)];
    
    % Record D statistics info: [meanD, std, SEM, 2*sqrt(4Dt)~pixel]
    diffStat = [diffStat; mean( dDiff), std( dDiff)/ sqrt( length(dDiff)), 2* sqrt( 4* mean(dDiff)*0.0217)/0.16];

    % Histogram
%         [N, edges] = histcounts( dDiff, 'BinWidth', 0.2, 'Normalization', 'probability');
%         tmp = movmean( edges, 2);   centers = tmp( 2:end);
%         plot( centers, N, lineStyle{ count}, 'LineWidth', 2, 'DisplayName',...
%             sprintf( '%s, mean=%.3f', [strain extraName], mean( dDiff)))
%         hold on
        
    %% D_alpha Coefficient CDF

    figure(3)
    dAlpha= dalpha( condDiff);
    [f,xgroup] = ecdf( dAlpha);  xgroup = xgroup([ 2: gap: end-1 end]); f = f([ 2: gap: end-1 end]);
    plot(xgroup, f, lineStyle{ count}, 'LineWidth', 2, 'DisplayName',...
        sprintf( '\\langleD\\alpha\\rangle=%.2f, [%s-%s]', mean( dAlpha), Date, strainName)) % strainName
    hold on

    %% alpha

    figure(4)
    alphaA = alpha( condDiff);
    [f,xgroup] = ecdf( alphaA);  xgroup = xgroup( [ 2: gap: end-1 end]); f = f( [ 2: gap: end-1 end]);
    plot(xgroup, f, lineStyle{ count}, 'LineWidth', 2, 'DisplayName',...
        sprintf( '%s, <\\alpha>=%.2f %s %.2f', [strain extraName], mean( alphaA), char(177), std( alphaA)))
    hold on

    %% EATA-MSD & fitting

    figure(5)
    time = (1: maxTau)* timeStep; % timelag for MSD
    condMSD = condDiff;

    ensTAMSD = EnsTAMSD( condMSD, :);
    eataMSD = mean( ensTAMSD, 1, 'omitnan'); % unit: um^2

    % apparent D linear fitting for EATA-MSD = 4Dt + 4*sigma^2
    fitRange = 1:3; % range for fitting
    p = polyfit( time( fitRange), eataMSD( fitRange), 1); % y=ax+b, may have negative y-intercept
    DFit = p(1)/ 4;     locErrFit = sqrt( p(2))* 1000/2; % unit: nm

    tMSD = 0: 0.02: 0.9;    MSDFit = p(1)*tMSD + p(2);

    plot( time, eataMSD, 'o', 'MarkerSize', 7, 'LineWidth', 0.5, 'color', colorList(count,:),'HandleVisibility','off');  hold on 
    plot( tMSD, MSDFit, lineStyle{ count}, 'LineWidth', 2.5,  'color', colorList(count,:), 'DisplayName',... % 'color', colorList(1,:),
        sprintf( 'D*=%.2f, Loc*=%.0f', DFit, locErrFit))
%         sprintf( 'D*=%.2f, \\sigma=%.0f+%.0fi nm,%s', DFit, real( locErrFit), imag( locErrFit), extraName))

%     % another condition: exclude negative y-intercept (imaginary locErr)
%     condMSD = condDiff & ~isnan( locErr);
%     ensTAMSD = EnsTAMSD( condMSD, :);
%     eataMSD = mean( ensTAMSD, 1, 'omitnan'); % unit: um^2
% 
%     p = polyfit( time( fitRange), eataMSD( fitRange), 1); % y=ax+b, may have negative y-intercept
%     DFit = p(1)/ 4;     locErrFit2 = sqrt( p(2))* 1000/2; % unit: nm
%     tMSD = 0: 0.02: 0.9;    MSDFit = p(1)*tMSD + p(2);
%     plot( time, eataMSD, 'o', 'LineWidth', 0.5, 'color', colorList(2,:), 'HandleVisibility','off');  hold on
%     plot( tMSD, MSDFit, 'LineWidth', 2, 'color', colorList(2,:), 'DisplayName',...
%         sprintf( 'D*=%.2f, \\sigma=%.0f nm (\\sigma>0)', DFit, locErrFit2))


    %% EA-MSD & EATA-MSD

    figure(6)    
    ensMSD = EnsMSD( condDiff,:);    ensTAMSD = EnsTAMSD( condDiff,:);    
    eaMSD = mean( ensMSD, 1, 'omitnan'); % 1*M, ensemble avg MSD(t)
    eataMSD = mean( ensTAMSD, 1, 'omitnan');
    
    visFlag = 'off';
%     if count == length( plotNum)
%         visFlag = 'on';
%     end
    
    plot( time, eaMSD, 'o', 'MarkerSize', 8, 'LineWidth', 1, 'MarkerEdgeColor', colorList(count,:), 'MarkerFaceColor', colorList(count,:), 'DisplayName',...
        sprintf( '%s', extraName)); hold on
    plot( time, eataMSD, 'o', 'MarkerSize', 6, 'LineWidth', 1, 'MarkerEdgeColor', colorList(count,:), 'MarkerFaceColor', 'w', 'DisplayName',...
        sprintf( '%s', extraName), 'HandleVisibility', visFlag)
    
%     fitRange = 1:3;
%     % log fitting for dalpha & alpha
%     pLog = polyfit( log( time( fitRange)), log( eaMSD( fitRange)), 1); % log(MSD) = alpha* log(t) + log(4D)
%     alphaFit = pLog(1);    DalphaFit = exp( pLog(2))/4;
% 
%     tMSD = 0.02: 0.005: 0.08;    MSDFit = 4* DalphaFit* tMSD.^alphaFit*2;
%     plot( tMSD, MSDFit, 'k--', 'LineWidth', 2, 'DisplayName',...
%         sprintf( 'D\\alpha=%.2f, \\alpha=%.2f (1:3 fit)', DalphaFit, alphaFit))
% 
%     fitRange = 6:12;
%     % log fitting for dalpha & alpha
%     pLog = polyfit( log( time( fitRange)), log( eaMSD( fitRange)), 1); % log(MSD) = alpha* log(t) + log(4D)
%     alphaFit = pLog(1);    DalphaFit = exp( pLog(2))/4;
% 
%     tMSD = 0.08: 0.01: 0.5;    MSDFit = 4* DalphaFit* tMSD.^alphaFit/2;
%     plot( tMSD, MSDFit, 'k--', 'LineWidth', 2, 'DisplayName',...
%         sprintf( 'D\\alpha=%.2f, \\alpha=%.2f (6:12 fit)', DalphaFit, alphaFit))


    %% Localization Error Histogram

    figure(7)
    meanLocErr = mean( locErr, 'omitnan'); % unit: nm
    histogram( locErr( condDiff), 'BinWidth', 8, 'LineWidth', 2.5, 'Normalization', 'probability',... % 'DisplayStyle', 'stairs', 'LineWidth', 2,
        'DisplayStyle', 'stairs', 'LineStyle', lineStyle{ count}, 'DisplayName',...
        sprintf( '\\langle\\sigma\\rangle=%.0fnm,%s', meanLocErr, extraName))
    hold on

    % locErrStat = [locErrStat; [locErrFit, meanLocErr]]; % Record locErr Info
    
end
%%

% Diffusion Comparison
figure(8)

% diffTheo = [ 0.55 2.27 3.28];
% plot( [4 5 6], diffTheo, 'o', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerEdgeColor', 'r',...
%     'MarkerFaceColor', 'r', 'DisplayName', 'SE theory')
% hold on

xx = 1: size( diffAll, 1);
for x = 1: size( diffAll, 1)    
    errorbar( x, diffAll(x,1), diffAll(x,2), 'o', 'LineWidth', 2,...
        'MarkerSize', 18, 'CapSize', 12, 'color', colorList( x,:), ...
        'MarkerEdgeColor', colorList( x,:)) % , 'DisplayName', 'Experiment'
    hold on
end
% nameAll{end} = [ 'dex ' nameAll{end}];

set( figure(8), 'Position', [1430 50 450 400])
set(gca, 'FontSize', 14)
ylabel( 'D (\mum^2/s)', 'FontSize', 14)
title( 'Diffusion Coefficient (12+ frames)', 'FontSize', 14)
legend( {'Simulation'}, 'Location', 'northwest', 'FontSize', 16)
xticks( xx);    xticklabels( nameAll)
xlim( [0.5 x+0.5]);   
% ylim( [0 4])

% errorbar( time, eataMSD, eataSEM2, 'o', 'LineWidth', 1, 'MarkerSize', 2, 'color', colorList(count,:), 'DisplayName',...
%     sprintf( '%s EATA-MSD', extraName))
% ydata = cell2mat( diffAll);
% meanDiff = cellfun( @mean, diffAll);    diffTheo = [ 0.55 2.27 3.28];
% % xgroup = [ repmat( nameAll{1}, size( diffAll{1})); repmat( nameAll{2}, size( diffAll{2})); repmat( nameAll{3}, size( diffAll{3}))];
% xgroup = [ ones( size( diffAll{1})); 2*ones( size( diffAll{2})); 3*ones( size( diffAll{3}))];
% b = boxplot( ydata, xgroup, 'Colors', colorList(1,:), 'Symbol', '');
% set( b, 'LineWidth', 1.5)
% hold on
% % plot( diffTheo, 'o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', colorList(5,:), 'DisplayName', 'SE theory')
% plot( meanDiff, '*', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerEdgeColor', colorList(2,:), 'DisplayName', 'mean D')
% hold off
%%
viscosity = 4040./( 6*pi* diffAll(:,1)*50);

fprintf( '\n   <D> = %.4f,  SEM = %.4f, pixel = %.4f, vis = %.3f', [diffStat viscosity]')
fprintf( '\n\n')


%% Plot Setting

% tracksLength
set( figure(1), 'Position', [50 550 450 400])
set( gca, 'FontSize', 12)
xlim( [4 51])
legend( 'FontSize', 12)
xlabel( 'Track Length (frames)', 'FontSize', 14)
ylabel( '1-CDF (portion)', 'FontSize', 14)
title( 'Track Length Distribution', 'FontSize', 14)
xline(12, '--', 'LineWidth', 1.5, 'FontSize', 13, 'HandleVisibility', 'off')

% D coefficient
set( figure(2), 'Position', [510 550 450 400])
set( gca, 'FontSize', 12)
legend( 'Location', 'southeast', 'FontSize', 13)
xlabel( 'D (\mum^2/s)', 'FontSize', 14)
ylabel( 'Cumulative Distribution Function', 'FontSize', 14)
title( 'Diffusion Coefficient (12+ frames)', 'FontSize', 14)
xlim( [-0.5 inf])
xlim( [-0.2 6])

% D_alpha
set( figure(3), 'Position', [970 550 450 400])
set( gca, 'FontSize', 12)
legend( 'Location', 'southeast', 'FontSize', 11)
xlabel( 'D\alpha (\mum^2/s)', 'FontSize', 14)
ylabel( 'Probability', 'FontSize', 14)
title( 'D\alpha Coefficient (12+ frames, 1:3 Fit)', 'FontSize', 14)

% alpha
set( figure(4), 'Position', [1430 550 450 400])
set( gca, 'FontSize', 12)
xlim( [-0.5 2.5])
xlabel( 'Dynamic Exponent \alpha', 'FontSize', 14)
ylabel( 'Probability', 'FontSize', 14)
legend( 'Location', 'northwest', 'FontSize', 13)
title( 'Scaling Factor \alpha (12+ frames), 1:3 Fit', 'FontSize', 14)

% EATA-MSD & fitting
set( figure(5), 'Position', [50 50 450 400])
set( gca, 'FontSize', 12)
xlabel( 'time Lag \tau (s)', 'FontSize', 14)
ylabel( 'EATA-MSD (µm^2)', 'FontSize', 14)
title( 'EATA-MSD (12+ frames, 1:3 fit)', 'FontSize', 14) %  (tl>=12 & D*>0)  (xNorm tracks)
legend( 'Location', 'northwest', 'FontSize', 13)
% set( gca, 'Xscale', 'linear', 'YScale', 'linear')
% ylim( [0 16]) % xlim( [0 1.2]);
ylim ([0 inf])
% set( gca, 'Xscale', 'log', 'YScale', 'log')
% xlim( [0.01 1]);  ylim( [0.01 2]

% EATA vs EA-MSD & fitting
set( figure(6), 'Position', [510 50 450 400])
set( gca, 'FontSize', 12)
xlabel('Time Lag \tau (s)', 'FontSize', 14)
ylabel('MSD (µm^2)', 'FontSize', 14)
title( 'EA vs EATA-MSD of Beads (12+ frames)', 'FontSize', 14)
legend( 'Location', 'northwest', 'FontSize', 13)
% legend( {'EA-MSD' 'EATA-MSD'}, 'Location', 'northwest', 'FontSize', 12)
set( gca, 'Xscale', 'log', 'YScale', 'log')
% xlim( [0.01 1.5]);   ylim( [0.05 100])
% ylim auto
    tMSD = 0.08: 0.01: 0.75;    MSDFit = 4* min( diffAll(:,1))* tMSD/2;
    plot( tMSD, MSDFit, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off')
    text( 0.24, MSDFit( round( end/3))/1.9, 'MSD ~ t', 'FontSize', 14)

% locErr distribution
set( figure(7), 'Position', [970 50 450 400])
set( gca, 'FontSize', 12)
legend( 'FontSize', 13)
xlabel( 'Localization Error \sigma (nm)', 'FontSize', 14)
ylabel( 'Probability', 'FontSize', 14)
title( 'locErr Distribution (12+ frames, \sigma>0)', 'FontSize', 14)
