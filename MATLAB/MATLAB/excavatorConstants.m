classdef excavatorConstants
    properties (Constant)
        % Boom
        iBC = [0.135; -0.264];
        bBA = [2.050; 0];
        bBD = [1.025; 0.278];
        bBE = [0.977; 0.737];
        bDA = excavatorConstants.bBA - excavatorConstants.bBD;
        bAE = excavatorConstants.bBE - excavatorConstants.bBA;
        bCoMBoom = [1.025; 0.384];
        
        lenBC = norm(excavatorConstants.iBC);
        lenBA = norm(excavatorConstants.bBA);
        lenBD = norm(excavatorConstants.bBD);
        lenBE = norm(excavatorConstants.bBE);
        lenDA = norm(excavatorConstants.bDA);
        lenAE = norm(excavatorConstants.bAE);
        lenBCoMBoom = norm(excavatorConstants.bCoMBoom);

        angABD = acos((excavatorConstants.lenBA^2 + ...
            excavatorConstants.lenBD^2 - ...
            excavatorConstants.lenDA^2)/...
            (2*excavatorConstants.lenBA*excavatorConstants.lenBD));
        angABCoMBoom = atan2(excavatorConstants.bCoMBoom(2), ...
            excavatorConstants.bCoMBoom(1));

        massBoom = 227.343;
        moiBoom = 67.768;
        
        % Arm
        aAF = [-0.251; 0.158];
        aAG = [-0.134; 0.320];
        aAJ = [0.880; 0];
        aAL = [1.050; 0];
        aFL = excavatorConstants.aAL - excavatorConstants.aAF;
        aJG = excavatorConstants.aAG - excavatorConstants.aAJ;
        aJL = excavatorConstants.aAL - excavatorConstants.aAJ;
        aCoMArm = [0.225; 0.227];

        lenHJ = 0.240;
        lenHK = 0.240;
        lenAF = norm(excavatorConstants.aAF);
        lenAG = norm(excavatorConstants.aAG);
        lenAJ = norm(excavatorConstants.aAJ);
        lenAL = norm(excavatorConstants.aAL);
        lenFL = norm(excavatorConstants.aFL);
        lenJG = norm(excavatorConstants.aJG);
        lenJL = norm(excavatorConstants.aJL);
        lenACoMArm = norm(excavatorConstants.aCoMArm);

        angFAL = acos((excavatorConstants.lenAL^2 + ...
            excavatorConstants.lenAF^2 - ...
            excavatorConstants.lenFL^2)/...
            (2*excavatorConstants.lenAF*excavatorConstants.lenAL));
        angLACoMArm = atan2(excavatorConstants.aCoMArm(2), ...
            excavatorConstants.aCoMArm(1));

        massArm = 130.123;
        moiArm = 30.258;

        % Bucket
        lLK = [-0.014; 0.164];
        lLM = [0.567; 0];
        %lLN = [0.248; 0.433];
        lKM = excavatorConstants.lLM - excavatorConstants.lLK;
        lCoMBucket = [0.289; 0.166];

        lenLK = norm(excavatorConstants.lLK);
        lenLM = norm(excavatorConstants.lLM);
        lenKM = norm(excavatorConstants.lKM);
        lenLCoMBucket = norm(excavatorConstants.lCoMBucket);

        angKLM = acos((excavatorConstants.lenLM^2 + ...
            excavatorConstants.lenLK^2 - ...
            excavatorConstants.lenKM^2)/...
            (2*excavatorConstants.lenLK*excavatorConstants.lenLM));
        angMLCoMBucket = atan2(excavatorConstants.lCoMBucket(2), ...
            excavatorConstants.lCoMBucket(1));

        massBucket = 53.000;
        moiBucket = 3.021;

        % Environment
        yGround = -0.95737;
        g = [0; -9.81; 0];
    end
end