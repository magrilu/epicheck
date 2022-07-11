function [ok,F] = isConsistentlyOriented(F,m1,m2,k_tol,e2)
% ISCONSISTENTLYORIENTED
% checks if the Fundamental matrix F satisfy the orientation constraints for
% all the pairs of corresponding points m1, m2
% Correspondences must be given in homogeneous coordinates as 3xn matrix
% If not provided in input the epipole in the second image is computed
%
% INPUT:
%        F: 3x3 fundamental matrix
%        m1: 3xn homogeneous coords in the first image
%        m2: 3xn homogeneous coords in the second image
%        k_tol: tolerance to test that correspondences satisfy the oriented
%        epipolar constraints. Default value 1e-6 (line 52). The lower the
%        more restrictive is the test.
%        e2: 3x1 epipole in the second image (optional)
% OUTPUT:
%        ok: A boolean that is true if all the points are on the same side of the
%        cameras, false otherwise
%        F: output with fixed sign according to the first correspondence
%
% PARAMS:
%        do_verbose: log information
%
% Reference:
% Chum, Werner, Matas,
% Epipolar Geometry Estimation via RANSAC Benefits from the Oriented Epipolar Constraint
% https://core.ac.uk/download/pdf/47168857.pdf
% see Eq.2 
%
% June 2022
% luca.magri@polimi.it
%
% //.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°//.\\°

% INTERNAL PARAMS
do_verbose = false;


% check arguments
if(size(m1,1)~=3) % m1 in homogenous coords
    error('Error in checking the orientation of epipolar constraints:\n m1 is a %dx%d, but shoud be a 3xn, points must be columns of homogeneous coords',size(m1,1),size(m1,2));
end
if(size(m2,1)~=3) % m2 in homogeneous coords
    error('Error in checking the orientation of epipolar constraints:\n m2 is a %dx%d, but shoud be a 3xn, points must be columns of homogeneous coords',size(m2,1),size(m2,2));
end
if(size(m1,2)~=size(m2,2)) % same number of points
    error('Error in checking the orientation of epipolar constraints:\n m1 has %d points, while m2 has %d points, the number of points in m1 and m2 should be equal',size(m1,2),size(m2,2));
end
N  = size(m1,2); % number of correspondences

if(~exist('k_tol','var') || isempty(k_tol))
    % compute the epipole in the second image
    k_tol = 1e-6;
end

if(~exist('e2','var') || isempty(e2))
    % compute the epipole in the second image
    [~,~,V]= svd(F');
    e2 = V(:,3);
    e2 = e2./e2(3);
end

jj = 1;
% fix the sign of F, so that the first point is in front of the camera
A = cross(e2,m2(:,jj));
B = F*m1(:,jj);
if(any(sign(A)~=sign(B)))
    F = -F;
end
B = F*m1(:,jj);
assert(all(sign(B)==sign(A)),'After fixing the sign of F, A and B should have the same signs');

ok = true;

while((jj < N) && ok)
    jj = jj+1;
    A = cross(e2,m2(:,jj));
    B = F*m1(:,jj);
    % fix the scale of A and B
    A = A./norm(A);
    B = B./norm(B);
    ok = all(norm(A-B) < k_tol);
    if(do_verbose)
        if(~ok)
            g1 = sprintf("%.3f ",A);
            g2 = sprintf("%.3f ",B);
            fprintf('oreintation test failed for the %d-th corr:\n A = [%s] and\n B = [%s]\n',jj,g1,g2);
        end
    end
end


end

