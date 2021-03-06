import numpy as np
import RJplots as RJ

### Loads the 4 example shapes which are stored in the Data folder
shape1 = np.load("Data/Shape1.npy")
shape2 = np.load("Data/Shape2.npy")
shape3 = np.load("Data/Shape3.npy")
shape4 = np.load("Data/Shape4.npy")

### Sets array to store RJ-values, as well as structure classifications
rj1 = np.zeros(4,dtype=float)
rj2 = np.zeros(4,dtype=float)
cla = np.zeros(4,dtype=int)

### Calls get_rj function on each shape to obtain the RJ-values and structure classification
rj1[0],rj2[0],cla[0] = RJ.get_rj(shape1)
rj1[1],rj2[1],cla[1] = RJ.get_rj(shape2)
rj1[2],rj2[2],cla[2] = RJ.get_rj(shape3)
rj1[3],rj2[3],cla[3] = RJ.get_rj(shape4)

### Makes an RJ-plot for these 4 shapes and prints out results
RJ.make_rj_plot(rj1,rj2,"example.png",markercolor="C2",markersize=10)
print("R1 values:", rj1)
print("R2 values:", rj2)
print("Classifications:", cla)

