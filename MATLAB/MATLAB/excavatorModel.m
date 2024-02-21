classdef excavatorModel < articulatedArmModel
    properties
        c; %excavatorConstants
        lenBoom;
        lenArm;
        lenBucket;
    end

    methods
        % Constructor
        function obj = excavatorModel(excavatorConstants, alpha, beta, gamma)
            obj = obj@articulatedArmModel(excavatorConstants.lenBA, ...
                excavatorConstants.lenAL, ...
                excavatorConstants.lenLM, ...
                alpha, beta, gamma);
            obj.c = excavatorConstants;
        end
        
        %{
        Set the boom, arm and bucket actuator lengths and update the
        corresponding angles
        Input: Boom, arm and bucket actuator lengths
        Output: Boom, arm and bucket angles
        %}
        function [alpha, beta, gamma] = setLengths(obj, lenBoom, lenArm, lenBucket)
            obj.lenBoom = lenBoom;
            obj.lenArm = lenArm;
            obj.lenBucket = lenBucket;
            
            % Boom (0.905 - 1.305 m)
            R = obj.c.lenBC; % R formula (cos version)
            theta = atan2(obj.c.iBC(2), obj.c.iBC(1));
            obj.ang1 = acos((-lenBoom^2 + obj.c.lenBD^2 + obj.c.lenBC^2) ...
                /(2*R*obj.c.lenBD)) + theta - obj.c.angABD;
            alpha = obj.ang1;

            % Arm (1.046 - 1.556 m)
            R = obj.c.lenAE; % R formula (sin version)
            theta = atan2(obj.c.bAE(1), obj.c.bAE(2));
            obj.ang2 = asin((-lenArm^2 + obj.c.lenAF^2 + obj.c.lenAE^2) ...
                /(2*R*obj.c.lenAF)) - theta - obj.c.angFAL;
            beta = obj.ang2;

            % Bucket (0.84 - 1.26 m)
            R = obj.c.lenJG; % R formula (sin version)
            theta = atan2(obj.c.aJG(1), obj.c.aJG(2));
            angLJH = asin((-lenBucket^2 + obj.c.lenHJ^2 + obj.c.lenJG^2) ...
                /(2*R*obj.c.lenHJ)) - theta;

            R = sqrt((obj.c.lenJL - obj.c.lenHJ*cos(angLJH))^2 + (obj.c.lenHJ*sin(angLJH))^2); % R formula (cos version)
            theta = atan2(obj.c.lenHJ*sin(angLJH), obj.c.lenJL - obj.c.lenHJ*cos(angLJH));
            obj.ang3 = acos((obj.c.lenHK^2 - obj.c.lenJL^2 - obj.c.lenLK^2 - obj.c.lenHJ^2 + ...
                2*obj.c.lenJL*obj.c.lenHJ*cos(angLJH))/(2*R*obj.c.lenLK)) - ...
                theta - obj.c.angKLM;
            gamma = obj.ang3;
        end
        
        %{
        Set the boom, arm and bucket angles and update the
        corresponding actuator lengths
        Input: Boom, arm and bucket angles
        Output: Boom, arm and bucket actuator lengths
        %}        
        function [lenBoom, lenArm, lenBucket] = setAngles(obj, alpha, beta, gamma)
            obj.ang1 = alpha;
            obj.ang2 = beta;
            obj.ang3 = gamma;

            % Boom
            R = obj.c.lenBC; % R formula (cos version)
            theta = atan2(obj.c.iBC(2), obj.c.iBC(1));
            obj.lenBoom = sqrt(-2*R*obj.c.lenBD*cos(alpha - theta + obj.c.angABD) + obj.c.lenBD^2 + obj.c.lenBC^2);
            lenBoom = obj.lenBoom;

            % Arm
            R = obj.c.lenAE; % R formula (sin version)
            theta = atan2(obj.c.bAE(1), obj.c.bAE(2)); 
            obj.lenArm = sqrt(-2*R*obj.c.lenAF*sin(beta + theta + obj.c.angFAL) + obj.c.lenAF^2 + obj.c.lenAE^2);
            lenArm = obj.lenArm;

            %Bucket
            aAK = obj.c.aAL + obj.c.lenLK*[cos(gamma + obj.c.angKLM); sin(gamma + obj.c.angKLM)];
            aJK = aAK - obj.c.aAJ;
            lenJK = norm(aJK);
            cosAngHJK = (lenJK^2 + obj.c.lenHJ^2 - obj.c.lenHK^2)/(2*lenJK*obj.c.lenHJ);
            sinAngHJK = sqrt(1 - cosAngHJK^2);
            aAH = obj.c.aAJ + (obj.c.lenHJ/lenJK)*[cosAngHJK -sinAngHJK; sinAngHJK cosAngHJK]*aJK;
            obj.lenBucket = sqrt((aAH(1) - obj.c.aAG(1))^2 + (aAH(2) - obj.c.aAG(2))^2);
            lenBucket = obj.lenBucket;
        end

        %{
        Calculate boom, arm and bucket virtual joint torques
        Input: Boom, arm and bucket angular accelerations and velocities
        Output: Boom, arm and bucket virtual joint torques
        %}
        function [TBoom, TArm, TBucket] = inverseDynamics(obj, ang1Dot, ang2Dot, ang3Dot, ang1DDot, ang2DDot, ang3DDot, F)
            qDot = [ang1Dot; ang2Dot; ang3Dot];
            qDDot = [ang1DDot; ang2DDot; ang3DDot];
            
            % Jacobians of the centre of masses
            iJacBoom = [-obj.c.lenBCoMBoom*sin(obj.ang1 + obj.c.angABCoMBoom), 0, 0; ...
                obj.c.lenBCoMBoom*cos(obj.ang1 + obj.c.angABCoMBoom), 0, 0; ...
                1, 0, 0];
            iJacArm = [-obj.c.lenBA*sin(obj.ang1) - obj.c.lenACoMArm*sin(obj.ang1 + obj.ang2 + obj.c.angLACoMArm), ...
                -obj.c.lenACoMArm*sin(obj.ang1 + obj.ang2 + obj.c.angLACoMArm), ...
                0; ...
                obj.c.lenBA*cos(obj.ang1) + obj.c.lenACoMArm*cos(obj.ang1 + obj.ang2 + obj.c.angLACoMArm), ...
                obj.c.lenACoMArm*cos(obj.ang1 + obj.ang2 + obj.c.angLACoMArm), ...
                0; ...
                1, 1, 0];
            iJacBucket = [-obj.c.lenBA*sin(obj.ang1) - obj.c.lenAL*sin(obj.ang1 + obj.ang2) - obj.c.lenLCoMBucket*sin(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket), ...
                -obj.c.lenAL*sin(obj.ang1 + obj.ang2) - obj.c.lenLCoMBucket*sin(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket), ...
                -obj.c.lenLCoMBucket*sin(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket); ...
                obj.c.lenBA*cos(obj.ang1) + obj.c.lenAL*cos(obj.ang1 + obj.ang2) + obj.c.lenLCoMBucket*cos(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket), ...
                obj.c.lenAL*cos(obj.ang1 + obj.ang2) + obj.c.lenLCoMBucket*cos(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket), ...
                obj.c.lenLCoMBucket*cos(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket); ...
                1, 1, 1];
            iJacBoomDot = [-obj.c.lenBCoMBoom*cos(obj.ang1 + obj.c.angABCoMBoom)*ang1Dot, 0, 0; ...
                -obj.c.lenBCoMBoom*sin(obj.ang1 + obj.c.angABCoMBoom)*ang1Dot, 0, 0; ...
                0, 0, 0];
            iJacArmDot = [-obj.c.lenBA*cos(obj.ang1)*ang1Dot - obj.c.lenACoMArm*cos(obj.ang1 + obj.ang2 + obj.c.angLACoMArm)*(ang1Dot + ang2Dot), ...
                -obj.c.lenACoMArm*cos(obj.ang1 + obj.ang2 + obj.c.angLACoMArm)*(ang1Dot + ang2Dot), ...
                0; ...
                -obj.c.lenBA*sin(obj.ang1)*ang1Dot - obj.c.lenACoMArm*sin(obj.ang1 + obj.ang2 + obj.c.angLACoMArm)*(ang1Dot + ang2Dot), ...
                -obj.c.lenACoMArm*sin(obj.ang1 + obj.ang2 + obj.c.angLACoMArm)*(ang1Dot + ang2Dot), ...
                0; ...
                0, 0, 0];
            iJacBucketDot = [-obj.c.lenBA*cos(obj.ang1)*ang1Dot - obj.c.lenAL*cos(obj.ang1 + obj.ang2)*(ang1Dot + ang2Dot) - obj.c.lenLCoMBucket*cos(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket)*(ang1Dot + ang2Dot + ang3Dot), ...
                -obj.c.lenAL*cos(obj.ang1 + obj.ang2)*(ang1Dot + ang2Dot) - obj.c.lenLCoMBucket*cos(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket)*(ang1Dot + ang2Dot + ang3Dot), ...
                -obj.c.lenLCoMBucket*cos(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket)*(ang1Dot + ang2Dot + ang3Dot); ...
                -obj.c.lenBA*sin(obj.ang1)*ang1Dot - obj.c.lenAL*sin(obj.ang1 + obj.ang2)*(ang1Dot + ang2Dot) - obj.c.lenLCoMBucket*sin(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket)*(ang1Dot + ang2Dot + ang3Dot), ...
                -obj.c.lenAL*sin(obj.ang1 + obj.ang2)*(ang1Dot + ang2Dot) - obj.c.lenLCoMBucket*sin(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket)*(ang1Dot + ang2Dot + ang3Dot), ...
                -obj.c.lenLCoMBucket*sin(obj.ang1 + obj.ang2 + obj.ang3 + obj.c.angMLCoMBucket)*(ang1Dot + ang2Dot + ang3Dot); ...
                0, 0, 0];
            iJacM = [-obj.c.lenBA*sin(obj.ang1) - obj.c.lenAL*sin(obj.ang1 + obj.ang2) - obj.c.lenLM*sin(obj.ang1 + obj.ang2 + obj.ang3), ...
                -obj.c.lenAL*sin(obj.ang1 + obj.ang2) - obj.c.lenLM*sin(obj.ang1 + obj.ang2 + obj.ang3), ...
                -obj.c.lenLM*sin(obj.ang1 + obj.ang2 + obj.ang3); ...
                obj.c.lenBA*cos(obj.ang1) + obj.c.lenAL*cos(obj.ang1 + obj.ang2) + obj.c.lenLM*cos(obj.ang1 + obj.ang2 + obj.ang3), ...
                obj.c.lenAL*cos(obj.ang1 + obj.ang2) + obj.c.lenLM*cos(obj.ang1 + obj.ang2 + obj.ang3), ...
                obj.c.lenLM*cos(obj.ang1 + obj.ang2 + obj.ang3); ...
                1, 1, 1];
            
            % Mass matrix
            M = iJacBoom(1:2,:)'*obj.c.massBoom*iJacBoom(1:2,:) + iJacBoom(3,:)'*obj.c.moiBoom*iJacBoom(3,:) + ...
                iJacArm(1:2,:)'*obj.c.massArm*iJacArm(1:2,:) + iJacArm(3,:)'*obj.c.moiArm*iJacArm(3,:) + ...
                iJacBucket(1:2,:)'*obj.c.massBucket*iJacBucket(1:2,:) + iJacBucket(3,:)'*obj.c.moiBucket*iJacBucket(3,:);

            % Coriolis and centrifugal terms
            b = iJacBoom(1:2,:)'*obj.c.massBoom*iJacBoomDot(1:2,:)*qDot + ...
                iJacBoom(3,:)'*(obj.c.moiBoom*iJacBoomDot(3,:)*qDot) + ...
                iJacArm(1:2,:)'*obj.c.massArm*iJacArmDot(1:2,:)*qDot + ...
                iJacArm(3,:)'*(obj.c.moiArm*iJacArmDot(3,:)*qDot) + ...
                iJacBucket(1:2,:)'*obj.c.massBucket*iJacBucketDot(1:2,:)*qDot + ...
                iJacBucket(3,:)'*(obj.c.moiBucket*iJacBucketDot(3,:)*qDot);

            % Gravitational terms
            g = -iJacBoom(1:2,:)'*obj.c.massBoom*obj.c.g(1:2) - ...
                iJacArm(1:2,:)'*obj.c.massArm*obj.c.g(1:2) - ...
                iJacBucket(1:2,:)'*obj.c.massBucket*obj.c.g(1:2);

            % External force
            extF = iJacM(1:2,:)'*F;
            
            T = M*qDDot + b + g - extF;
            TBoom = T(1);
            TArm = T(2);
            TBucket = T(3);
        end
        
        %{
        Calculate boom, arm and bucket actuator forces
        Input: Boom, arm and bucket virtual joint torques
        Output: Boom, arm and bucket actuator forces
        %}
        function [FBoom, FArm, FBucket] = calcForces(obj, TBoom, TArm, TBucket)
            % Boom
            angBDC = acos((obj.lenBoom^2 + obj.c.lenBD^2 - obj.c.lenBC^2)/...
                (2*obj.lenBoom*obj.c.lenBD));
            theta = pi/2 - angBDC; % Angle between tangential velocity vector of D and linear velocity of boom actuator
            FBoom = TBoom/(obj.c.lenBD*cos(theta));

            % Arm
            angAFE = acos((obj.lenArm^2 + obj.c.lenAF^2 - obj.c.lenAE^2)/...
                (2*obj.lenArm*obj.c.lenAF));
            theta = angAFE - pi/2; % Angle between tangential velocity vector of F and linear velocity of arm actuator
            FArm = TArm/(obj.c.lenAF*cos(theta));

            % Bucket
            FBucket = TBucket/-(0.0192*obj.ang3^3 + 0.0864*obj.ang3^2 + 0.045*obj.ang3 - 0.1695);
        end

        function visualise(obj)
            hold on;
            axis equal;
            
            % Boom (0.77 - 1.35 m)
            iBA = rotMat(obj.ang1)*obj.c.bBA;
            iBD = rotMat(obj.ang1)*obj.c.bBD;
            iBE = rotMat(obj.ang1)*obj.c.bBE;
            boom = polyshape([[0; 0], iBD, iBA, iBE]');
            
            plot(boom, FaceColor="k");
            plot([0, obj.c.iBC(1)], [0, obj.c.iBC(2)], Color="k"); % BC
            plot([obj.c.iBC(1), iBD(1)], [obj.c.iBC(2), iBD(2)], Color="r"); % Boom actuator
            
            % Arm (1.01 - 1.59 m)
            iBF = transMat(obj.ang1 + obj.ang2, iBA)*[obj.c.aAF; 1];
            iBF = iBF(1:2);
            iBG = transMat(obj.ang1 + obj.ang2, iBA)*[obj.c.aAG; 1];
            iBG = iBG(1:2);
            iBJ = transMat(obj.ang1 + obj.ang2, iBA)*[obj.c.aAJ; 1];
            iBJ = iBJ(1:2);
            iBL = transMat(obj.ang1 + obj.ang2, iBA)*[obj.c.aAL; 1];
            iBL = iBL(1:2);
            
            
            R = 2*obj.c.lenJG; % R formula (sin version)
            theta = atan2(obj.c.aJG(1), obj.c.aJG(2));
            angLJH = asin((-obj.lenBucket^2 + obj.c.lenHJ^2 + obj.c.lenJG^2) ...
                /(R*obj.c.lenHJ)) - theta;
            aAH = transMat(angLJH, obj.c.aAJ)*[obj.c.lenHJ; 0; 1];
            aAH = aAH(1:2);
            iBH = transMat(obj.ang1 + obj.ang2, iBA)*[aAH; 1];
            iBH = iBH(1:2);
            

            arm = polyshape([iBA, iBJ, iBL, iBG, iBF]');

            plot(arm, FaceColor="k");
            plot([iBE(1), iBF(1)], [iBE(2), iBF(2)], Color="r"); % Arm actuator
            plot([iBG(1), iBH(1)], [iBG(2), iBH(2)], Color="r"); % Bucket actuator
            plot([iBJ(1), iBH(1)], [iBJ(2), iBH(2)], Color="k");
            
            % Bucket (0.83 - 1.26 m)
            iBK = transMat(obj.ang1 + obj.ang2 + obj.ang3, iBL)*[obj.c.lLK; 1];
            iBK = iBK(1:2);
            iBM = transMat(obj.ang1 + obj.ang2 + obj.ang3, iBL)*[obj.c.lLM; 1];
            iBM = iBM(1:2);
            bucket = polyshape([iBL, iBM, iBK]');

            plot(bucket, FaceColor="k");
            plot([iBH(1), iBK(1)], [iBH(2), iBK(2)], Color="k");

            % Ground
            yline(obj.c.yGround, Color="g");

            hold off;
        end
    end
end