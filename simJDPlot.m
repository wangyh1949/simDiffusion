%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 12/6/2023
    last updated date: 1/3/2024

Description: this script 'simJDPlot.m' is adapted from 'jdPlot.m'
it plot the 1-CDF of jump distance and do fitting to extract D coefficient

-------------------------------------------------------------
%}

timeStep = frameT; meanD = D; extraName = '';

count = 1;


%% histogram of displacement in 3D vs 2D
figure
colorList = get( gca,'colororder');

h = histogram( steps3D, 40, 'Normalization', 'probability', ... % , 'BinWidth', 0.002, 
    'FaceAlpha', 0.7, 'EdgeColor', 'w', 'DisplayName', '3D'); hold on

histogram( steps, 'BinWidth', h.BinWidth, 'Normalization', 'probability', ... % , 'BinWidth', 0.002, 
    'FaceAlpha', 0.4, 'EdgeColor', 'w', 'DisplayName', '2D')

    set( figure( gcf), 'Position', [200 350 430 380])
    set( gca, 'FontSize', 12)
    xlabel( 'Step Size (\mum)', 'FontSize', 14)
    ylabel( 'Probability', 'FontSize', 14)
    legend( 'FontSize', 14)
    title( 'Simulation 3D', 'FontSize', 16)
    

%% Jump Distance 1-CDF analysis

sqSteps = steps.^2; % r^2 [um^2]
sqSteps3D = steps3D.^2;

gap = round( length( sqSteps)/ 150); % prepare for CDF plotting

% 3D
    [f, x] = ecdf( sqSteps3D); x = x( [ 2: gap: end-1 end]); f = f( [ 2: gap: end-1 end]);                     
    
    % linear population fit
    g = fittype( 'exp( -x/b)');
    cdfFit = fit( x, 1-f, g, 'StartPoint', 4* timeStep*meanD);
    
    figure
    plot( x, 1-f, 'o', 'LineWidth', 0.5, 'color', [colorList( count,:) 0.5],...
        'DisplayName', sprintf( '1-CDF, D=%.2f', D))
    hold on

    DFit = (cdfFit.b - 6* locErr^ 2)/ (6* timeStep);
    %         DFit = cdfFit.b/ (4* timeStep);
    plot( x, cdfFit(x), 'color', [colorList( count,:) 0.5],...
        'LineWidth', 2, 'DisplayName', sprintf( 'D=%.3f, %s', DFit, [strain extraName]))

    % 1-CDF Jump Distance
    set( figure( gcf), 'Position', [1000 250 350 320])
    set( gca, 'FontSize', 12)
    set(gca, 'YScale', 'log')
    ylim( [0.01, 1]);
    legend( 'FontSize', 12)
    xlabel('Step Size r^2 (\mum^2)', 'FontSize', 14)
    ylabel('1 - P(r^2, \tau)', 'FontSize', 14)
    title( 'Simulation Data 3D', 'FontSize', 14)

% 2D    
    [f, x] = ecdf( sqSteps); x = x( [ 2: gap: end-1 end]); f = f( [ 2: gap: end-1 end]);  
    
    g = fittype( 'exp( -x/b)');
    cdfFit = fit( x, 1-f, g, 'StartPoint', 4* timeStep*meanD);
    
    figure
    plot( x, 1-f, 'o', 'LineWidth', 0.5, 'color', [colorList( count,:) 0.5],...
        'DisplayName', sprintf( '1-CDF, D=%.2f', D))
    hold on

    DFit = (cdfFit.b - 4* locErr^ 2)/ (4* timeStep);
    %         DFit = cdfFit.b/ (4* timeStep);
    plot( x, cdfFit(x), 'color', [colorList( count,:) 0.5],...
        'LineWidth', 2, 'DisplayName', sprintf( 'D=%.3f, %s', DFit, [strain extraName]))

    % 1-CDF Jump Distance
    set( figure( gcf), 'Position', [1360 250 350 320])
    set( gca, 'FontSize', 12)
    set(gca, 'YScale', 'log')
    ylim( [0.01, 1]);
    legend( 'FontSize', 12)
    xlabel('Step Size r^2 (\mum^2)', 'FontSize', 14)
    ylabel('1 - P(r^2, \tau)', 'FontSize', 14)
    title( 'Simulation Data 2D', 'FontSize', 14)


    
%% 2-pop Fit function
% 
%     g2 = fittype( 'a*exp( -x/b)+ (1-a)*exp( -x/d)');        
% 
%     % meanD is actually severely underestimating the D
%     cdfFit2 = fit( x, 1-f, g2, 'StartPoint', [0.9, 4*timeStep*meanD, 4*timeStep*meanD*5],...
%         'Lower', [0 0 0], 'Upper', [1 4*timeStep*meanD*5 4*timeStep*meanD*20]);
% 
%     % 2 population fit
%     figure(2)
%     plot( x, 1-f, 'o', 'LineWidth', 0.5, 'Color', [colorList( ceil( count/ colorOrder),:) 0.5],...
%         'HandleVisibility','off')
%     hold on
% 
%     getD = @(x) ( x- 4* locError^ 2)/ (4* timeStep);
% %         getD = @(x) x/ (4* timeStep);
% 
%     if cdfFit2.b < cdfFit2.d
%         fitResult = [ cdfFit2.a*100, getD( cdfFit2.b), getD( cdfFit2.d)];
%     else
%         fitResult = [ (1-cdfFit2.a)*100, getD( cdfFit2.d), getD( cdfFit2.b)];
%     end
%     plot( x, cdfFit2(x), lineStyle{ order( count)}, 'color', colorList( ceil( count/ colorOrder),:),...
%         'LineWidth', 2, 'DisplayName',...
%         sprintf( '[%.0f%%] D1=%.2f, D2=%.2f, %s', fitResult, [strain extraName]))
% 
% 
%     set( figure(2), 'Position', [820 300 400 380])
%     set( gca, 'FontSize', 11)
%     set(gca, 'YScale', 'log')
%     % xlim( [0,0.8]);
%     ylim( [0.01, 1])
%     legend( 'Location', 'southwest', 'FontSize', 10)
%     xlabel('Step Size r^2 (\mum^2)', 'FontSize', 12)
%     ylabel('1 - P(r^2, \tau)', 'FontSize', 12)
%     title( 'JD 2-Population Fit', 'FontSize', 14)        
        