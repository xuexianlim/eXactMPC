function rotMat = rotMat(theta)
%{
Generate rotation matrix from an angle
Input: Angle in radians
Output: 2D rotation matrix
%}
    rotMat = [cos(theta), -sin(theta);
              sin(theta), cos(theta)];
end