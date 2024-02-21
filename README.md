# eXactMPC

https://github.com/xuexianlim/eXactMPC

# C++
This program implements kinematic MPC for the eXact electric excavator arm. Casadi and Ipopt (already included in Casadi) are required for it to run. The installation steps can be found here: https://github.com/casadi/casadi/wiki/InstallationLinux

To build this repository, run the following commands:
```
$ cd build
$ cmake ..
$ make
```

To run the program, run the following command:
```
$ ./eXactMPC
```

# Python
The file to run is `MPCSimulation.py`. At least Python 3.10 is required to run the program due to the presence of match-case syntax.

In order to generate videos from the visualisation module, OpenCV should be installed via the following command:
```
$ sudo apt-get install python3-opencv
```

Using pip3 to install opencv-python might not work due to the absence of the necessary codec for encoding the video.

# MATLAB
The code runs as it is. The file to run is `simulation.m`.