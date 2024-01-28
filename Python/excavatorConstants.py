import casadi as csd

# Boom
iBC = csd.vertcat(0.135, -0.264)
bBA = csd.vertcat(2.050, 0)
bBD = csd.vertcat(1.025, 0.278)
bBE = csd.vertcat(0.977, 0.737)
bDA = bBA - bBD
bAE = bBE - bBA
bCoMBoom = csd.vertcat(1.025, 0.384)
        
lenBC = csd.norm_2(iBC)
lenBA = csd.norm_2(bBA)
lenBD = csd.norm_2(bBD)
lenBE = csd.norm_2(bBE)
lenDA = csd.norm_2(bDA)
lenAE = csd.norm_2(bAE)
lenBCoMBoom = csd.norm_2(bCoMBoom)

angABD = csd.acos((lenBA**2 + lenBD**2 - lenDA**2)/
              (2*lenBA*lenBD))
angABCoMBoom = csd.atan2(bCoMBoom[1], bCoMBoom[0])

massBoom = 227.343
moiBoom = 67.768
        
# Arm
aAF = csd.vertcat(-0.251, 0.158)
aAG = csd.vertcat(-0.134, 0.320)
aAJ = csd.vertcat(0.880, 0)
aAL = csd.vertcat(1.050, 0)
aFL = aAL - aAF
aJG = aAG - aAJ
aJL = aAL - aAJ
aCoMArm = csd.vertcat(0.225, 0.227)

lenHJ = 0.240
lenHK = 0.240
lenAF = csd.norm_2(aAF)
lenAG = csd.norm_2(aAG)
lenAJ = csd.norm_2(aAJ)
lenAL = csd.norm_2(aAL)
lenFL = csd.norm_2(aFL)
lenJG = csd.norm_2(aJG)
lenJL = csd.norm_2(aJL)
lenACoMArm = csd.norm_2(aCoMArm)

angFAL = csd.acos((lenAL**2 + lenAF**2 - lenFL**2)/
              (2*lenAF*lenAL))
angLACoMArm = csd.atan2(aCoMArm[1], aCoMArm[0])

massArm = 130.123
moiArm = 30.258

# Bucket
lLK = csd.vertcat(-0.014, 0.164)
lLM = csd.vertcat(0.567, 0)
lKM = lLM - lLK
lCoMBucket = csd.vertcat(0.289, 0.166)

lenLK = csd.norm_2(lLK)
lenLM = csd.norm_2(lLM)
lenKM = csd.norm_2(lKM)
lenLCoMBucket = csd.norm_2(lCoMBucket)

angKLM = csd.acos((lenLM**2 + lenLK**2 - lenKM**2)/
              (2*lenLK*lenLM))
angMLCoMBucket = csd.atan2(lCoMBucket[1], lCoMBucket[0])

massBucket = 53.000
moiBucket = 3.021

# Environment
yGround = -0.95737
g = csd.vertcat(0, -9.81)