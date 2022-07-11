% test oriented epipolar constraints on synthetic data
demoFilePath = fileparts(matlab.desktop.editor.getActiveFilename);
cd(demoFilePath); % move to the folder of the demo
addpath(genpath(demoFilePath));
%addpath(genpath('../BrewerMap/'))
addpath(genpath('../ComputerVisionToolkit/'))

clear;
close all;

%\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\

do_debugPlot = true;
mrkSzBig = 100;
mrkSzSmall = 20;
N = 11;
%cmap = brewermap(N,'Set3');
cmap = colormap("lines");
cmap = cmap(1:N,:);
close(gcf);
%\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\

% Instantiate camera pairs (reference frame of the first camera = world reference frame)
I = eye(3);                 % rotation cam 1
R = eul([0,0.2,0]);         % rotation cam 2
t = [-1;0;-0.5];            % translation cam 2
f = 2*rand(1);              % focal
K =diag([f,f,1]);           % shared intrinsic
P1 = K*[I,zeros(3,1)];      % first camera matrix
P2 = K*[R, t];              % second camera matrix
C1 = null(P1);
C1 = C1./C1(end);           % center cam 1
C2 = null(P2);
C2 = C2./C2(end);           % center cam 2
G = [R,t; zeros(1,3),1];    % 3D transformations
G_inv = [R', -R'*t; zeros(1,3),1] ;


% Instantiate world points \\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.
N = 10; % number of points in front of the camera
pScene = [ 2*rand(1,N); 1*rand(1,N); rand(1,N)+7; ones(1,N)];


% Project points to image \\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\

x1h = P1*pScene;                   % homogeneous coordinates
x1 = x1h./ repmat(x1h(3,:),3,1);    % cartesian coordinates
x1N = inv(K)*x1;
x1N = x1N./ repmat(x1N(3,:),3,1);   % normalized coordinates

x2h = P2*pScene;
x2 = x2h./ repmat(x2h(3,:),3,1);    % homogeneous coordinates

x2N = inv(K)*x2;
x2N = x2N./ repmat(x2N(3,:),3,1);   % normalized coordinates
x2N_cam1 =R'*x2N -R'*t;             % normalized coordinates in the reference frame of the first camera


% Compute fundamental matrix and epipoles \\°//.\\°//.\\°//.\\°//.\\°//.\\°
F_corr = fund_lin(x2N(1:2,:), x1N(1:2,:));
F_cam = fund(P1,P2);
F = F_cam;
e1 = epipole(F);
e1 = e1./e1(3);
e1N = inv(K)*e1;
e2 = epipole(F');
e2 = e2./e2(3);
e2N = inv(K)*e2;
e2N_cam1 =R'*e2N -R'*t;

%% Test 1: all points in front of the cameras

ok = isConsistentlyOriented(F,x1,x2);
assert(ok,'isConsistentlyOriented fails for points in front of the cameras. Ok should be true!');

if(do_debugPlot)
    figure;
    displayCamera(I,C1,'b');
    hold all;
    displayCamera(R,C2,'r');
    axis equal;
    scatter3(pScene(1,:),pScene(2,:),pScene(3,:),mrkSzBig,cmap(1:N,:),'filled');
    scatter3(x1N(1,:),x1N(2,:),x1N(3,:),mrkSzSmall,cmap(1:N,:),'filled');
    scatter3(x2N_cam1(1,:),x2N_cam1(2,:),x2N_cam1(3,:),mrkSzSmall,cmap(1:N,:),'filled');

    id = 1;
    W = pScene(:,id);
    line([x1N(1,id),W(1)],[x1N(2,id),W(2)], [x1N(3,id),W(3)],'Color','b','LineStyle','--');
    line([C1(1),W(1)],[C1(2),W(2)], [C1(3),W(3)],'Color','b','LineStyle','--');
    line([x2N_cam1(1,id),W(1)],[x2N_cam1(2,id),W(2)], [x2N_cam1(3,id),W(3)],'Color','r','LineStyle','--');
    line([C2(1),W(1)],[C2(2),W(2)], [C2(3),W(3)],'Color','r','LineStyle','--');
    plot3(e2N_cam1(1),e2N_cam1(2),e2N_cam1(3),'r*');
    line([C2(1),e1N(1)],[C2(2),e1N(2)], [C2(3),e1N(3)],'Color','c','LineStyle','--');
    line([C1(1),e2N_cam1(1)],[C1(2),e2N_cam1(2)], [C1(3),e2N_cam1(3)],'Color','c','LineStyle','--');
    plot3(e1N(1),e1N(2),e1N(3),'r*');
    title('Point in front of the camera')
end


%% Test 2: a points behind the second camera

N_back = 1;
fprintf('Added %d points behind the cameras\n', N_back);
pScene(:,end+1) = [ 2*rand(1,N_back); 1*rand(1,N_back); 0.3*rand(1,N_back); 1];
N =size(pScene, 2);


% Project points to image \\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\

x1h = P1*pScene;                   % homogeneous coordinates
x1 = x1h./ repmat(x1h(3,:),3,1);    % cartesian coordinates
x1N = inv(K)*x1;
x1N = x1N./ repmat(x1N(3,:),3,1);   % normalized coordinates

x2h = P2*pScene;
x2 = x2h./ repmat(x2h(3,:),3,1);    % homogeneous coordinates

x2N = inv(K)*x2;
x2N = x2N./ repmat(x2N(3,:),3,1);   % normalized coordinates
x2N_cam1 =R'*x2N -R'*t;             % normalized coordinates in the reference frame of the first camera




ok = isConsistentlyOriented(F,x1,x2);
assert(~ok,'isConsistentlyOriented pass for configuration with a point behind the second camera and in front to the first. Ok should be false');

if(do_debugPlot)
    figure;
    displayCamera(I,C1,'b');
    hold all;
    displayCamera(R,C2,'r');
    axis equal;
    scatter3(pScene(1,:),pScene(2,:),pScene(3,:),mrkSzBig,cmap(1:N,:),'filled');
    scatter3(x1N(1,:),x1N(2,:),x1N(3,:),mrkSzSmall,cmap(1:N,:),'filled');
    scatter3(x2N_cam1(1,:),x2N_cam1(2,:),x2N_cam1(3,:),mrkSzSmall,cmap(1:N,:),'filled');

    id = N;
    W = pScene(:,id);
    line([x1N(1,id),W(1)],[x1N(2,id),W(2)], [x1N(3,id),W(3)],'Color','b','LineStyle','--');
    line([C1(1),W(1)],[C1(2),W(2)], [C1(3),W(3)],'Color','b','LineStyle','--');
    line([x2N_cam1(1,id),W(1)],[x2N_cam1(2,id),W(2)], [x2N_cam1(3,id),W(3)],'Color','r','LineStyle','--');
    line([C2(1),W(1)],[C2(2),W(2)], [C2(3),W(3)],'Color','r','LineStyle','--');
    plot3(e2N_cam1(1),e2N_cam1(2),e2N_cam1(3),'r*');
    line([C2(1),e1N(1)],[C2(2),e1N(2)], [C2(3),e1N(3)],'Color','c','LineStyle','--');
    line([C1(1),e2N_cam1(1)],[C1(2),e2N_cam1(2)], [C1(3),e2N_cam1(3)],'Color','c','LineStyle','--');
    plot3(e1N(1),e1N(2),e1N(3),'r*');
    title('A single point behind the camera')
end
