# eXactMPC
# C++
This program implements kinematic MPC for the eXact electric excavator arm. Casadi and Ipopt are required for it to run. The installation steps can be found here: https://github.com/casadi/casadi/wiki/InstallationLinux

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
For the visualisation video generation, OpenCV should be installed via the following command:
```
$ sudo apt-get install python3-opencv
```

Using pip3 to install opencv-python might not work due to the absence of the necessary codec for encoding the video.