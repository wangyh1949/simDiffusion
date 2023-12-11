%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 12/6/2023
    last updated date: 12/9/2023

Description: this script 'diffSimMemb' is adapted from 'diffSimCell'
1) simulates the diffusion of a free particle on cell membrane in 3D 

cell shape is a spherocylinder (capsule)

~~~~~~ How to Ratate a Vector in 3D ~~~~~~
~~~~~~  Rodriguez rotation formula  ~~~~~~

v′=vcosθ+(r×v)sinθ+r(v*r)(1−cosθ)   where

 v is the original vector, 
 r is the vector about which the rotation on angle θ is performed, 
 v′ is the vector after rotation.


~~~~~~~~~~~~~~ Steps ~~~~~~~~~~~~~

1. Find an orthogonal vector k1 to v at the base
2. Use cross( k1, v) to find the third orthogonal vector k2
        k1-k2 defines the orthogonal plane to vector v
3. Find a random vector in k1-k2 plane as the rotatioin axis
        ratation radius is always r 
        rotation direction would be random
4. Calculate rotation angle by 'step size/r'
5. With the rotation axis and angle, I can get the new pts 
6. Fix the position if it goes out of the hemisphere


12.7: determine rotation axis by the k1-k2 plane orthogonal to v
-------------------------------------------------------------
%}

clear, clc, close all

simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)
strain = 'sim';

plotTrackFlag = 1;

nTracks = 100;     nFrames = 30;   maxTau = nFrames - 1; % for TA-MSD

D = 0.2;    % um^2/s 
dt = 1e-2;  % 10 ms, unit: s

locError = 0e-2; % 40 nm, unit: um

cellWid = 0.6;    cellLength = 2; % cell long axis: y & short axis: x-z
global r l
r = cellWid/2;    l = ( cellLength - cellWid)/ 2;

fprintf( '~~~~ Simulation Starts ~~~~\n')
fprintf( '   D = %.2f, dt = %d ms, step = %.2f um, locErr = %.2f um\n', D, dt*1e3, sqrt( 4*D*dt), locError)
fprintf( '   nFrame = %d, total step = %.2f um\n', nFrames, sqrt( 4*D*dt*nFrames))


%% 1. generate randomly distributed origin points on the membrane
    
rng(0) % to make the result repeatable

    % equivalent to: unifrnd( -cellLength/2, cellLength/2, nTracks, 1);
    y = ( rand( nTracks, 1)- 0.5)* cellLength; % random y from -L/2 ~ L
    phi = rand( nTracks, 1)* 2*pi;
    
    radi = ones( nTracks, 1)* r; % radius at x-z plane, in cylinder is r, if cap is smaller
    inCap = abs( y) >= l; % flag for Cap region
    radi( inCap) = sqrt( r^2 - (y( inCap)- l*sign( y( inCap))).^2); % radius in cap region
    x = radi.* cos( phi);   z = radi.* sin( phi);
    ori = [x, y, z]; % origin points for the diffusion
    
%     scatter3( x, y, z)
%     axis image

%% 2. Diffusion Simulation

close all

% nTracks = 100;     nFrames = 100;

tracksFinal( nTracks, 1).traj = [];

tic
for i = 1: nTracks
    
    traj = nan( nFrames, 3);
    pos = ori( i,:);  traj( 1,:) = pos;
    
    randAngle = 2*pi* rand( nFrames-1, 1); % random angle from 0~2pi
    jumps = raylrnd( sqrt( 2*D*dt), nFrames-1, 1); % step magnitude, unit: um 
%     jumps = ones( nFrames-1, 2)*0.1;
    
    for j = 1: nFrames-1
        
        y = pos(2);
        
        % in the Cap region
        if abs( y) >= l 
            
            v = pos - [0, l* sign(y), 0]; % relative vector on a sphere            
            
            % find the basis k1 & k2 of the orthogonal plane of v
            phiOrtho = atan( v(3)./ v(1))+ pi/2; % orthogonal phi angle at the base of hemisphere
            k1 = [ cos( phiOrtho), 0, sin( phiOrtho)]; % orthogonal vector to v in the plane of base
            k2 = cross( v, k1)/ r; % orthogonal vector to v & k1
            
            % find a random vector in k1-k2 plane as the rotation axis
            k = k1* cos( randAngle(j)) + k2* sin( randAngle(j));
            
            % jump size should have variance = 4Dt, so times sqrt(2)
            rotTheta = jumps(j)* sqrt(2)/ r; % angle of rotation, what should it be 6Dt or 4Dt?

            % v′=vcosθ+(r×v)sinθ+r(v*r)(1−cosθ), dot(v,k)=0
            vRot = v*cos( rotTheta) + cross(k, v)* sin( rotTheta);
            
            % if goes out of the cap region, fix the position to the cylinder
            pos = cellShape( vRot + [0, l* sign(y), 0], 'onCap'); % add back the shift
            
        else
        % in the cylinder region
            
            phi = atan2( pos(3), pos(1)); % phi angle at the base of hemisphere            
            phiRot = jumps(j)* cos( randAngle(j))/ r;            
            pos = cellShape( [ r*cos( phi+ phiRot), pos(2)+jumps(j)*sin( randAngle(j)), r*sin( phi+ phiRot)], 'onCyl');
                        
        end
        
        traj( j+1,:) = pos;
    end

    locE = normrnd( 0, sqrt( 2*locError^2), [nFrames, 3]); % unit: um
    traj = traj + locE;
    
    % save track X & Y coordinates into tracksFinal structure
    tracksFinal(i).traj = traj; % x & y coordiate, unit: um
%     traj = traj( :, 1:2);   
    
    % ~~~~ Calculate Single Steps ~~~~
    singleSteps = sqrt( sum(( traj( 2:end, :)- traj( 1:end-1, :)).^2, 2)); % unit: um    
    tracksFinal(i).steps = singleSteps'; % unit: um        
    
    if logical( plotTrackFlag)
        plot3( traj(:,1), traj(:,2), traj(:,3), 'LineWidth', 1), hold on
    end
end

fprintf( '\n~~~~~~  Simulation Done  ~~~~~~\n')
toc

steps  = [ tracksFinal.steps]'; % unit: um

% Plot Setting

if logical( plotTrackFlag)
    figure( gcf)
    hold off
    set( gcf, 'Position', [900, 450, 720, 400])
    set( gca, 'FontSize', 12)
    xlabel( 'X (\mum)', 'FontSize', 14)
    ylabel( 'Y (\mum)', 'FontSize', 14)
    zlabel( 'Z (\mum)', 'FontSize', 14)
    title( 'Simulation on Cell Membrane', 'FontSize', 16)
    view(-60, 20) % view(az,el)  az=-37, el = 30
    axis image
end

    

%% Diffusion Analysis

% diffSimAnalysis
% 
% msd = EnsTAMSD(:,1); % MSD(1)
% 
% simJDPlot

%     simName = sprintf( 'Sim %s_%d frames_D %.2f_LocE %d.mat', Date, nFrames, D, locError*1e3);
%     save( [simPath 'Data\' simName])

fprintf( '\n~~~~~~ !!!  All Complete  !!! ~~~~~~\n')



%%

function posCorr = cellShape( pos, posFlag)

    global r l
    [ x, y, z] = deal( pos(1), pos(2), pos(3));
    
    posCorr = pos;
    
    if strcmp( posFlag, 'onCap') && abs( y) < l
        % should fix the pos to Cylinder            
        xz = [ x, z]; % vector with respect to the axis of the cylinder
        xzCorr = r* xz/ norm( xz);
        posCorr = [ xzCorr(1) y xzCorr(2)];
    elseif strcmp( posFlag, 'onCyl') && abs( y) >= l
        % should fix the pos to Cap
        xz = [ x, z]; % vector with respect to the axis of the cylinder
        rCorr = sqrt( ( r^2 - (y- l*sign(y))^2));
        xzCorr =  rCorr* xz/ norm( xz);
        posCorr = [ xzCorr(1) y xzCorr(2)];            
    end
    
end

