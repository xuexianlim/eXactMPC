classdef articulatedArmModel < handle
    properties
        len1;
        len2;
        len3;
        ang1;
        ang2;
        ang3;
    end

    methods
        % Constructor
        function obj = articulatedArmModel(len1, len2, len3, ang1, ang2, ang3)
            obj.len1 = len1;
            obj.len2 = len2;
            obj.len3 = len3;
            obj.ang1 = ang1;
            obj.ang2 = ang2;
            obj.ang3 = ang3;
        end

        %{
        Perform forward kinematics to calculate end effector pose
        Input: Obj
        Output: End effector pose
        %}
        function poseEE = forwardKinematics(obj)
            posEE = [obj.len1*cos(obj.ang1) + ...
                obj.len2*cos(obj.ang1 + obj.ang2) + ....
                obj.len3*cos(obj.ang1 + obj.ang2 + obj.ang3);
                obj.len1*sin(obj.ang1) + ...
                obj.len2*sin(obj.ang1 + obj.ang2) + ....
                obj.len3*sin(obj.ang1 + obj.ang2 + obj.ang3)];
            %{
            posEE = transMat(obj.ang1, [0; 0])*...
                transMat(obj.ang2, [obj.len1; 0])*...
                transMat(obj.ang3, [obj.len2; 0])*...
                [obj.len3; 0; 1];
            %}
            angEE = obj.ang1 + obj.ang2 + obj.ang3;
            poseEE = [posEE(1:2); angEE];
        end

        %{
        Perform inverse kinematics to obtain joint angles
        [Robotics: Modelling, Planning and Control (Siciliano et al., 2008,
        pp. 92-93)]
        Input: End effector pose
        Output: Vectors of possible joint angles
        %}
        function [anglesI, anglesII] = inverseKinematics(obj, pose)
            posEEx = pose(1);
            posEEy = pose(2);
            angEE = pose(3);
            
            posJoint3x = posEEx - obj.len3*cos(angEE);
            posJoint3y = posEEy - obj.len3*sin(angEE);
            
            % First possible option (I)
            cosAng2 = (posJoint3x^2 + posJoint3y^2 - obj.len1^2 - obj.len2^2)/...
                (2*obj.len1*obj.len2);
            sinAng2 = sqrt(1 - cosAng2^2);
            
            sinAng1 = ((obj.len1 + obj.len2*cosAng2)*posJoint3y - ...
                obj.len2*sinAng2*posJoint3x)/...
                (posJoint3x^2 + posJoint3y^2);
            cosAng1 = ((obj.len1 + obj.len2*cosAng2)*posJoint3x + ...
                obj.len2*sinAng2*posJoint3y)/...
                (posJoint3x^2 + posJoint3y^2);

            anglesI(1) = atan2(sinAng1, cosAng1);
            anglesI(2) = atan2(sinAng2, cosAng2);
            anglesI(3) = angEE - anglesI(1) - anglesI(2);

            % Second possible option (II)
            sinAng2 = -sinAng2;

            sinAng1 = ((obj.len1 + obj.len2*cosAng2)*posJoint3y - ...
                obj.len2*sinAng2*posJoint3x)/...
                (posJoint3x^2 + posJoint3y^2);
            cosAng1 = ((obj.len1 + obj.len2*cosAng2)*posJoint3x + ...
                obj.len2*sinAng2*posJoint3y)/...
                (posJoint3x^2 + posJoint3y^2);

            anglesII(1) = atan2(sinAng1, cosAng1);            
            anglesII(2) = atan2(sinAng2, cosAng2);
            anglesII(3) = angEE - anglesII(1) - anglesII(2);
        end
    end
end