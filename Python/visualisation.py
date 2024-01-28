import casadi as csd
import excavatorConstants as C
import matplotlib.pyplot as plt
import os
import cv2

plt.rcParams["figure.dpi"] = 300
visFolder = './Python/Plots/'

# Utils
def rotMat(theta):
    return csd.vertcat(csd.horzcat(csd.cos(theta), -csd.sin(theta)),
                       csd.horzcat(csd.sin(theta), csd.cos(theta)))

def transMat(theta, pos):
    return csd.vertcat(csd.horzcat(rotMat(theta), pos),
                       csd.horzcat(0, 0, 1))

# Plotting
def plotExcavator(ax, q, **kwargs):
    alpha = q[0]
    beta = q[1]
    gamma = q[2]
    
    color = 'default'
    alphaTransparency = 0
    if 'color' in kwargs.keys():
        color = kwargs['color']
    if 'alpha' in kwargs.keys():
        alphaTransparency = kwargs['alpha']

    # Boom
    iBA = rotMat(alpha)@C.bBA
    iBD = rotMat(alpha)@C.bBD
    iBE = rotMat(alpha)@C.bBE

    temp = [iBD, iBA, iBE]
    xBoom = [0]
    yBoom = [0]
    for pos in temp:
        xBoom += [pos[0].__float__()]
        yBoom += [pos[1].__float__()]
    if color == 'default':
        ax.fill(xBoom, yBoom, facecolor='yellow', edgecolor='orange', linewidth=1)
        plt.plot([C.iBC[0].__float__(), iBD[0].__float__()], [C.iBC[1].__float__(), iBD[1].__float__()], color='r', linewidth=1)
    else:
        ax.fill(xBoom, yBoom, facecolor='none', edgecolor=color, linewidth=1, alpha=alphaTransparency)
        plt.plot([C.iBC[0].__float__(), iBD[0].__float__()], [C.iBC[1].__float__(), iBD[1].__float__()], color=color, linewidth=1, alpha=alphaTransparency)

    # Arm
    iBF = transMat(alpha + beta, iBA)@csd.vertcat(C.aAF, 1)
    iBF = iBF[0:2]
    iBG = transMat(alpha + beta, iBA)@csd.vertcat(C.aAG, 1)
    iBG = iBG[0:2]
    iBJ = transMat(alpha + beta, iBA)@csd.vertcat(C.aAJ, 1)
    iBJ = iBJ[0:2]
    iBL = transMat(alpha + beta, iBA)@csd.vertcat(C.aAL, 1)
    iBL = iBL[0:2]

    temp = [iBA, iBJ, iBL, iBG, iBF, iBA]
    xBoom = []
    yBoom = []
    for pos in temp:
        xBoom += [pos[0].__float__()]
        yBoom += [pos[1].__float__()]
    if color == 'default':
        ax.fill(xBoom, yBoom, facecolor='yellow', edgecolor='orange', linewidth=1)
    else:
        ax.fill(xBoom, yBoom, facecolor='none', edgecolor=color, linewidth=1, alpha=alphaTransparency)

    R = 2*C.lenJG # R formula (sin version)
    theta = csd.atan2(C.aJG[0], C.aJG[1])
    lenBucket = 0.0048*gamma**4 + 0.0288*gamma**3 + 0.0225*gamma**2 - 0.1695*gamma + 0.9434
    angLJH = csd.asin((-lenBucket**2 + C.lenHJ**2 + C.lenJG**2)/(R*C.lenHJ)) - theta
    aAH = transMat(angLJH, C.aAJ)@csd.vertcat(C.lenHJ, 0, 1)
    aAH = aAH[0:2]
    iBH = transMat(alpha + beta, iBA)@csd.vertcat(aAH, 1)
    iBH = iBH[0:2]
    if color == 'default':
        plt.plot([iBE[0].__float__(), iBF[0].__float__()], [iBE[1].__float__(), iBF[1].__float__()], color='r', linewidth=1)
        plt.plot([iBG[0].__float__(), iBH[0].__float__()], [iBG[1].__float__(), iBH[1].__float__()], color='r', linewidth=1)
        plt.plot([iBJ[0].__float__(), iBH[0].__float__()], [iBJ[1].__float__(), iBH[1].__float__()], color='k', linewidth=1)
    else:
        plt.plot([iBE[0].__float__(), iBF[0].__float__()], [iBE[1].__float__(), iBF[1].__float__()], color=color, linewidth=1, alpha=alphaTransparency)
        plt.plot([iBG[0].__float__(), iBH[0].__float__()], [iBG[1].__float__(), iBH[1].__float__()], color=color, linewidth=1, alpha=alphaTransparency)
        plt.plot([iBJ[0].__float__(), iBH[0].__float__()], [iBJ[1].__float__(), iBH[1].__float__()], color=color, linewidth=1, alpha=alphaTransparency)        

    # Bucket
    iBK = transMat(alpha + beta + gamma, iBL)@csd.vertcat(C.lLK, 1)
    iBK = iBK[0:2]
    iBM = transMat(alpha + beta + gamma, iBL)@csd.vertcat(C.lLM, 1)
    iBM = iBM[0:2]

    temp = [iBL, iBM, iBK, iBL]
    xBoom = []
    yBoom = []
    for pos in temp:
        xBoom += [pos[0].__float__()]
        yBoom += [pos[1].__float__()]
    if color == 'default':
        ax.fill(xBoom, yBoom, facecolor='yellow', edgecolor='orange', linewidth=1)
        plt.plot([iBH[0].__float__(), iBK[0].__float__()], [iBH[1].__float__(), iBK[1].__float__()], color='k', linewidth=1)
    else:
        ax.fill(xBoom, yBoom, facecolor='none', edgecolor=color, linewidth=1, alpha=alphaTransparency)
        plt.plot([iBH[0].__float__(), iBK[0].__float__()], [iBH[1].__float__(), iBK[1].__float__()], color=color, linewidth=1, alpha=alphaTransparency)        

def visualise(q, qOld, qDes, t, k):
    fig, ax = plt.subplots()

    if qOld is not None:
        plotExcavator(ax, qOld, color='k', alpha=0.2)
    if qDes is not None:
        plotExcavator(ax, qDes, color='lime', alpha=0.2)
    plotExcavator(ax, q)

    # Ground
    plt.axhline(y=C.yGround, color='g', linewidth=1)

    plt.xlim([-1, 4])
    plt.ylim([-1.5, 3])
    ax.set_aspect('equal')
    ax.set_title("t = {t} s, k = {k}".format(t=t, k=k))
    plt.savefig(visFolder + "Excavator_{y}.jpg".format(y=k))
    plt.close()

def graph(tStart, tEnd, interval, name, **kwargs):
    n = (tEnd - tStart)/interval
    x = csd.linspace(tStart, tEnd, int(n + 1))
    
    fig, ax = plt.subplots()
    for label, y in kwargs.items():
        if y.ndim == 1:
            if y.shape[0] == n + 1:
                plt.plot(x, y, label="{label}".format(label=label) , linewidth=1)
            elif y.shape[0] == n:
                plt.plot(x[1:], y, label="{label}".format(label=label) , linewidth=1)
        elif y.ndim == 2:
            for row in range(y.shape[0]):
                if y.shape[1] == n + 1:
                    plt.plot(x, y[row, :], label="{label}{row}".format(label=label, row=row) , linewidth=1)
                elif y.shape[1] == n:
                    plt.plot(x[1:], y[row, :], label="{label}{row}".format(label=label, row=row) , linewidth=1)
    ax.legend()
    ax.set_title(name)
    plt.legend(bbox_to_anchor=(1.04, 1), loc='upper left')
    plt.savefig(visFolder + "Graph_{name}.jpg".format(name=name), bbox_inches='tight')
    plt.close()

def createVideo(kStart, kEnd, name, fps):
    videoName = visFolder + "{name}.mp4".format(name=name)
    imgs = []

    k = kStart
    while "Excavator_{k}.jpg".format(k=k) in os.listdir(visFolder) and k <= kEnd:
        imgs += ["Excavator_{k}.jpg".format(k=k)]
        k += 1
    
    frame = cv2.imread(os.path.join(visFolder, imgs[0])) 
    height, width, layers = frame.shape   
  
    video = cv2.VideoWriter(videoName, cv2.VideoWriter_fourcc(*'avc1'), fps, (width, height))  
  
    # Appending the images to the video one by one 
    for img in imgs:  
        video.write(cv2.imread(os.path.join(visFolder, img)))
      
    # Deallocating memories taken for window creation 
    cv2.destroyAllWindows()  
    video.release()  # releasing the video generated 