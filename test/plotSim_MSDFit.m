%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 9/4/2025
    last updated date: 9/5/2025

    ------- this script is adapted from plotLociMSDFit.m -------

Description: This code plots the MSD of simulated diffusion data and does
fitting to the MSD (multiple methods)

1. EA-MSD   & linear fitting
2. EATA-MSD & linear fitting
3. EATA-MSD & non-linear fitting
---------------------------------------------------------------------------
%}


%% 
close all

c = 1;

expT = timeStep; % continuous exposure
plotRange = 1: maxT-1; % all tracks have same length, maxT = nFrames
time = (1: maxT-1)* timeStep;
legtxt = sprintf( '%s', strain); % 'sim'

simInfo = sprintf( 'D=%.2f, locE=%dnm\nexpT=%dms, simT=%dms\ncell wid=%.1fum, length=%.1fum\n%d tracks, %d frames', ...
    D, locErr*1e3, frameT*1e3, dt*1e3, cellWid, cellLength, nTracks, nFrames);

fitR = 1:4;    fitR2 = fitR;   fitR3 = 1:12;
% set up fitting parameters
fitTxt1 = sprintf( '%d:%d fit', min( fitR), max( fitR));
fitTxt2 = sprintf( '%d:%d fit', min( fitR2), max( fitR2));
fitTxt3 = sprintf( '%d:%d fit', min( fitR3), max( fitR3));


% set up figures
f1 = figure( 'Position', [400 400 400 380]);
f2 = figure( 'Position', [810 400 400 380]);
f3 = figure( 'Position', [1220 400 400 380]);
colorList = get( gca,'colororder'); colorList = repmat( colorList, [2, 1]);


eaMSD = mean( EnsMSD, 1, 'omitnan'); % unit: um^2
eataMSD = mean( EnsTAMSD, 1, 'omitnan');            

        
    %% 0. linear fit: MSD = 4Dt + 4*loc^2 - 4Ddt/3
    f = polyfit( time( fitR), eataMSD( fitR), 1); % may be negative
    Dapp = f(1)/ 4;     locErrFit = sqrt( f(2)+ timeStep*f(1)/3)* 1000/2; % unit: nm
    

    %% 1. EA-MSD & Linear fitting    
    f = polyfit( log( time( fitR)), log( eaMSD( fitR)), 1);
    alphaFit2 = f(1);    Dapp2 = exp( f(2))/4;
    
    figure( f1)    
    scatter( time( plotRange), eaMSD( plotRange), 30, 'LineWidth', 2, ...% 'MarkerEdgeColor', colorList(c,:),...
        'MarkerEdgeAlpha', 0.5, 'HandleVisibility','off'), hold on

    t =  linspace( time(1), time( 3*fitR(end)), 1000);  MSDFit =  4* Dapp2* t.^ alphaFit2;
    plot( t, MSDFit, 'LineWidth', 1.5, ... %'color', [ colorList(c,:) 0.8], ...
        'DisplayName', sprintf( 'D\\alpha=%.2f, \\alpha=%.2f [%s]', Dapp2, alphaFit2, fitTxt1)) % , [strain extraName]
    
        
    %% 2. EATA-MSD & Linear fitting 
    % log-log fit, MSD = 4Dt^alpha, log(MSD) = alpha*log(t)+ log(4D)        
    f = polyfit( log( time( fitR2)), log( eataMSD( fitR2)), 1);
    alphaFit = f(1);    DFit = exp( f(2))/4;
    
    figure( f2)
    scatter( time( plotRange), eataMSD( plotRange), 30, 'LineWidth', 1.75, ...% 'MarkerEdgeColor', colorList(c,:),...
        'MarkerEdgeAlpha', 0.5, 'HandleVisibility','off'), hold on
    
    t =  linspace( time(1), time( 3*fitR2(end)), 1000);     
    MSDFit =  4* DFit* t.^ alphaFit;
    plot( t, MSDFit, 'LineWidth', 1.5, 'DisplayName', ... % 'color', [ colorList(c,:) 0.8],
        sprintf( 'D\\alpha=%.2f, \\alpha=%.2f [%s]', DFit, alphaFit, fitTxt2))  
    
    text( 0.95, 0.05, simInfo, 'FontSize', 11, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom' , 'Units', 'normalized')


    %% 3. EATA-MSD & Non-linear fitting
    % non-linear fit using a comprehensive form, MSD = 4Dt^a + b
    fun = fittype( 'log( 4*a*(x)^b+4*c)');
    x0 = [ DFit, alphaFit, (locErrFit/1e3)^2];
    xmin = [ 0, 0, -inf];   xmax = [ inf, 2, inf];
    
    f = fit( ( time( fitR3))', log( eataMSD( fitR3))', fun, 'StartPoint', x0, 'Lower', xmin, 'Upper', xmax);
    DFit = f.a;  alphaFit = f.b;  %locErrFit = sign( f.c)* sqrt( abs( f.c)); % unit: nm
    
    motionBlur = 4*2* DFit* expT^alphaFit/ (( 1+alphaFit)* (2+alphaFit));
    locErrFit = sqrt( ( 4*f.c + motionBlur)/ 4); % unit: um
    
    figure( f3)
    scatter( time( plotRange), eataMSD( plotRange), 30, 'LineWidth', 1.5, ...%'MarkerEdgeColor', colorList(c,:),...
        'MarkerEdgeAlpha', 0.5, 'HandleVisibility','off'), hold on
    
    t =  linspace( time(1), time( fitR3(end)), 1000);    MSDFit =  4* DFit* t.^ alphaFit + 4*f.c;
    plot( t, MSDFit, 'LineWidth', 1.5, 'DisplayName', ... % 'color', [ colorList(c,:) 0.8], 
        sprintf( 'D=%.2f, \\alpha=%.2f, \\sigma=%.0f [%s]', DFit, alphaFit, locErrFit*1e3, fitTxt3))

    text( 0.95, 0.05, simInfo, 'FontSize', 11, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom' , 'Units', 'normalized')


%% Figure Setting

% limX = [0.01 50]; limY = [1e-3 0.1]; % 20-20ms
limX = [1e-2 1]; limY = [1e-2 1]; % SK407: LacZ, 20ms
% limX = 'auto'; limY = 'auto';

% EA-MSD & fitting
figure( f1)
set( gca, 'FontSize', 13)
xlabel( 'Time (s)', 'FontSize', 14)
ylabel( 'EA-MSD (µm2)', 'FontSize', 14)
legend( 'Location', 'northwest', 'FontSize', 11)
title( sprintf( 'EA-MSD (linear fit)'), 'FontSize', 14)
set( gca, 'Xscale', 'log', 'YScale', 'log')
box on

xlim( limX) 
ylim( limY)

% EATA-MSD & fitting
figure( f2)
set( gca, 'FontSize', 14)
xlabel( 'Time (s)', 'FontSize', 14)
ylabel( 'EATA-MSD (µm^2)', 'FontSize', 14)
legend( 'Location', 'northwest', 'FontSize', 11)
title( sprintf( 'sim MSD (linear fit)'), 'FontSize', 14)
set( gca, 'Xscale', 'log', 'YScale', 'log')
box on

xlim( limX)
ylim( limY)

% EATA-MSD & nonlinear fitting
figure( f3)
set( gca, 'FontSize', 13)
xlabel( 'Time (s)', 'FontSize', 14)
ylabel( 'EATA-MSD (µm^2)', 'FontSize', 14)
legend( 'Location', 'northwest', 'FontSize', 11)
title( sprintf( 'sim MSD (non-linear fit)'), 'FontSize', 14)
set( gca, 'Xscale', 'log', 'YScale', 'log')
box on

xlim( limX)
ylim( limY)