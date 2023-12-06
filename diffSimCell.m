%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 12/5/2023
    last updated date: 7/23/2023

Description: this script 'diffSimCell' is adapted from 'diffSim3D'
1) simulates the diffusion of a free particle in a confined cell shape in
3D

======Input=======
nTracks, nFrames, D, timeStep, locError 

======Output=======
tracksFinal stucture with fields: traj & MSD
    EnsMSD, EnsTAMSD

MSD fitting parameter: 
    diff, locErr, dalpha, alpha (unit adjusted for plotting)

-------------------------------------------------------------
%}


clear, clc, close all


simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Simulation\'; cd( simPath)

nTracks = 100;     nFrames = 50;   maxTau = nFrames - 1;

D = 0.2; % um^2/s
dt = 1e-2; % 10 ms, unit: s

locError = 5e-2; % 50 nm, unit: um

cellWid = 0.6;    cellLength = 2; % cell long axis: y & short axis: x-z
global r l
r = cellWid/2;    l = ( cellLength - cellWid)/ 2;

fprintf( '~~~~ Simulation Starts ~~~~\n')
fprintf( '   D = %.2f, dt = %d ms, step = %.2f um   ||', D, dt*1e3, sqrt( 4*D*dt))
fprintf( '   nFrame = %d, total step = %.2f um\n', nFrames, sqrt( 4*D*dt*nFrames))

%%
rng(0) % to make the result repeatable
plotTrackFlag = 0;

tic

%     nTracks = 10000; 
    % generate randomly distributed origin points
    ori = (rand( nTracks*2, 3)- 0.5).* [ cellWid, cellLength, cellWid];    
    [ x, y, z] = deal( ori(:,1), ori(:,2), ori(:,3));

    inCyl = abs( y) < l & ( x.^2 + z.^2) < r^2; % in the cylinder region
    inCap = abs( y) >= l & ( x.^2 + z.^2 + ( y- l*sign(y)).^2) < r^2; % in the cap region

    ori = ori( inCyl | inCap, :);
    if length( ori) < nTracks
        warning( ' ~~~~~ Not Enough Origin Points !!! ~~~~~')
    end
        
    %%

    
tracksFinal( nTracks, 1).traj = [];

% initialization of the quantities with the correct size
EnsMSD = nan( nTracks, maxTau); EnsTAMSD = nan( nTracks, maxTau);
Diff = nan( nTracks, 1);        LocErr = nan( nTracks, 1);
Dalpha = nan( nTracks, 1);      alpha = nan( nTracks, 1);


% figure( 'Position', [900, 450, 720, 400])
    
    for i = 1: nTracks
        
%         pos = [0 0 0]; % initial position, should be homogeneously generated in cell
        pos = ori( i, :); % initial position, should be homogeneously generated in cell
        traj = nan( nFrames, 3);    traj( 1, :) = pos;
        
        % generate random steps
        % std of the step distribution should be sqrt( 2Dt)
        jumps = normrnd( 0, sqrt( 2*D*dt), [nFrames-1, 3]); % unit: um
        
        % iterations, each step judge whether it's outside the cell
        for j = 1: nFrames-1
            pos = cellShape( pos + jumps(j, :));
            traj( j+1, :) = pos;
        end
        
%         locE = normrnd( 0, sqrt( 2*locError^2), [nFrames, 2]); % unit: um

        % save track X & Y coordinates into tracksFinal structure
        tracksFinal(i).traj = traj; % x & y coordiate, unit: um
        
        if logical( plotTrackFlag)
            plot3( traj(:,1), traj(:,2), traj(:,3), 'LineWidth', 1.5),        hold on        
        end
        
        
        % ~~~~ Calculate Single Steps ~~~~
        singleSteps = sqrt( sum(( traj( 2:end, :)- traj( 1:end-1, :)).^2, 2)); % unit: um    
        tracksFinal(i).steps = singleSteps'; % unit: um        
        
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
        
        Diff(i) = p(1)/ 6; % may be negative, 3D case
        LocErr(i) = sqrt( p(2))/ 2; % unit: um, may be imaginary (should exclude later)    

        % fit the log(MSD) to get Dalpha & alpha
        fitRange = 1: 3; % max(3, floor( nFrames/4)); % 1/4 of the track
        pLog = polyfit( log( dt* fitRange), log( taMSD( fitRange)), 1); % log(MSD) = alpha* log(t) + log(4D)
        alpha(i) = pLog(1);
        Dalpha(i) = exp( pLog(2))/4;
    end

    LocErr( LocErr ~= real( LocErr)) = nan; % delete the negative intercepts (imaginary)

    EnsMSD = EnsMSD* 1e12;      % unit: um^2
    EnsTAMSD = EnsTAMSD* 1e12;  % unit: um^2
    diff   = Diff* 1e12;          % unit: um^2/s
    dalpha = Dalpha* 1e12;      % unit: um^2/s
    locErr = LocErr* 1e9;       % unit: nm
    steps  = [ tracksFinal.steps]'; % unit: um
    
    tracksLength = cellfun( @length, {tracksFinal.traj}');
    strain = 'Sim';
    extraName = sprintf( 'D: %.2f, Loc: %d', D, locError*1e3);

%     simName = sprintf( 'Sim %s_%d frames_D %.2f_LocE %d.mat', Date, nFrames, D, locError*1e3);
%     save( [simPath 'Data\' simName])

disp( '~~~~~~ !!!  All Complete  !!! ~~~~~~')
toc




% Plot Setting
if logical( plotTrackFlag)
    figure( gcf)
    set( gcf, 'Position', [900, 450, 720, 400])
    set( gca, 'FontSize', 12)
    xlabel( 'X', 'FontSize', 14)
    ylabel( 'Y', 'FontSize', 14)
    zlabel( 'Z', 'FontSize', 14)
    % legend( 'FontSize', 14)
    title( 'Simulated Tracks In Cell', 'FontSize', 16)
    axis image
    view(-60, 20) % view(az,el)  az=-37, el = 30
end


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
