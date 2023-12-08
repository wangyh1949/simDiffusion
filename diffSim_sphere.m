%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 12/8/2023
    last updated date: 12/8/2023

Description: this script 'diffSim_sphere' is adapted from 'diffSimMemb'
1) simulates the diffusion of a free particle on a 3D sphere 
2) I also added Maggie's simulation code as a comparison

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


-------------------------------------------------------------
%}

clear, clc, close all

simPath = 'C:\Users\yuhuanw2\Documents\MATLAB\simDiffusion\'; cd( simPath)

    nTracks = 100;     nFrames = 30;
    plotTrackFlag = 1;

    D = 0.2; % um^2/s 
    dt = 1e-2; % 10 ms, unit: s
    locError = 0e-2; % unit: um

global r l
r = 1; % radius of the sphere

fprintf( '~~~~ Simulation Starts ~~~~\n')
fprintf( '   D = %.2f, dt = %d ms, step = %.2f um, locErr = %.2f um\n', D, dt*1e3, sqrt( 4*D*dt), locError)
fprintf( '   nFrame = %d, total step = %.2f um\n', nFrames, sqrt( 4*D*dt*nFrames))


%% 1. generate randomly distributed origin points on the sphere

rng(0) % to make the result repeatable
    
    % equivalent to: unifrnd( -cellLength/2, cellLength/2, nTracks, 1);
    y = ( rand( nTracks, 1)- 0.5)* 2*r; % random y from -r ~ r
    phi = rand( nTracks, 1)* 2*pi;    
    radi = sqrt( r^2 - y.^2);    
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
%     jumps = 0.1* ones( nFrames-1, 1); % fixed step size
    
    for j = 1: nFrames-1
        
        v = pos;
            
        % find the basis k1 & k2 of the orthogonal plane of v
        phiOrtho = atan( v(3)./ v(1))+ pi/2; % orthogonal phi angle at the base of hemisphere
        k1 = [ cos( phiOrtho), 0, sin( phiOrtho)]; % orthogonal vector to v in the plane of base
        k2 = cross( v, k1)/ r; % orthogonal vector to v & k1

        % find a random vector in k1-k2 plane as the rotation axis
        k = k1* cos( randAngle(j)) + k2* sin( randAngle(j));

        % jump size should have variance = 2Dt in each dimension
        rotTheta = jumps(j, 1)/ r; % angle of rotation, what should it be 6Dt or 4Dt?

        % v′=vcosθ+(r×v)sinθ+r(v*r)(1−cosθ), dot(v,k)=0
        pos = v*cos( rotTheta) + cross(k, v)* sin( rotTheta); % new position after rotation
             
        traj( j+1,:) = pos;
    end
    
    locE = normrnd( 0, sqrt( 2*locError^2), [nFrames, 3]); % unit: um
    traj = traj + locE;
    
    % save track X & Y coordinates into tracksFinal structure
    tracksFinal(i).traj = traj; % x & y coordiate, unit: um
    
    % ~~~~ Calculate Single Steps ~~~~
    singleSteps = sqrt( sum(( traj( 2:end, :)- traj( 1:end-1, :)).^2, 2)); % unit: um    
    tracksFinal(i).steps = singleSteps'; % unit: um        
    
    plot3( traj(:,1), traj(:,2), traj(:,3)), hold on % , 'LineWidth', 1.5
end

% Plot Setting

if logical( plotTrackFlag)
    figure( gcf)
    hold off
    set( gcf, 'Position', [600, 450, 420, 400])
    set( gca, 'FontSize', 12)
    xlabel( 'X (\mum)', 'FontSize', 14)
    ylabel( 'Y (\mum)', 'FontSize', 14)
    zlabel( 'Z (\mum)', 'FontSize', 14)
%     title( 'Maggie''s Membrane Simulation', 'FontSize', 16)
    title( 'Simulated Tracks on Sphere', 'FontSize', 16)
    axis image
    view(-60, 20) % view(az,el)  az=-37, el = 30
end


%% Maggie's Code of Simulating Membrane Diffusion
% she used a different way to rotate the vector, a little bit confusing,
% but it works well!

figure
for i = 1: nTracks
    
    traj = nan( nFrames, 3);
    pos = ori( i,:);  traj( 1,:) = pos;
     
    jumps = raylrnd( sqrt( 2*D*dt), nFrames-1, 1); % step magnitude, unit: um 
%     jumps = 0.1* ones( nFrames-1, 1); % fixed step size
    
    for j = 1: nFrames-1
        
        R = r; step_size = jumps(j);
        
        %get initial point position info
        z0 = pos(3);
        phi0 = atan2(pos(2),pos(1));
        
        rand_angle = 2*pi*rand;
        theta0 = atan2( sqrt(R.^2-z0.^2), z0);
        
        %see Rodrigues' rotation formula. the point of below is just to
        %locate a point with set distance to initial point
        phi_perp = phi0+pi./2;
        rand_pos_sph = [R,rand_angle,step_size/R];
        
        phi = rand_pos_sph(2);
        theta = rand_pos_sph(3);
        x = R.*cos(phi).*sin(theta);
        y = R.*sin(phi).*sin(theta);
        z = R.*cos(theta);
        
        v = [x,y,z];
        k = [cos(phi_perp),sin(phi_perp),0];
        pos = v.*cos(theta0) +cross(k,v).*sin(theta0) + k.*(sum(k.*v)).*(1-cos(theta0));
        
        traj(j+1,:) = pos;
    end
    
    locE = normrnd( 0, sqrt( 2*locError^2), [nFrames, 3]); % unit: um
    traj = traj + locE;
        
    % ~~~~ Calculate Single Steps ~~~~
    singleSteps = sqrt( sum(( traj( 2:end, :)- traj( 1:end-1, :)).^2, 2)); % unit: um    
    tracksFinal(i).stepsMag = singleSteps'; % unit: um        
    
    plot3( traj(:,1), traj(:,2), traj(:,3)), hold on % , 'LineWidth', 1.5
end


fprintf( '\n~~~~~~  Simulation Done  ~~~~~~\n')
toc

steps  = [ tracksFinal.steps]'; % unit: um
stepsMag  = [ tracksFinal.stepsMag]'; % unit: um


% Plot Setting
if logical( plotTrackFlag)
    hold off
    figure( gcf)
    set( gcf, 'Position', [1030, 450, 420, 400])
    set( gca, 'FontSize', 12)
    xlabel( 'X (\mum)', 'FontSize', 14)
    ylabel( 'Y (\mum)', 'FontSize', 14)
    zlabel( 'Z (\mum)', 'FontSize', 14)
%     title( 'Maggie''s Membrane Simulation', 'FontSize', 16)
    title( 'Maggie''s Simulated Tracks', 'FontSize', 16)
    axis image
    view(-60, 20) % view(az,el)  az=-37, el = 30
end


fprintf( '\n~~~~~~ !!!  All Complete  !!! ~~~~~~\n')



