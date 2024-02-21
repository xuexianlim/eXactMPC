function transMat = transMat(theta, pos)
%{
Generate transformation matrix from an angle and position vector
Input: Angle in radians and position vector
Output: 2D transformation matrix
%}
    transMat = [rotMat(theta), pos;
                0, 0, 1];
end