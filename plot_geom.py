# Licensed under a 3-clause BSD style license - see LICENSE.rst
import taco
from Quaternion import Quat

q = Quat([30,20,0])

planes = taco.make_taco()
radiator = taco.make_radiator()

T = q.transform

def plotline(p0, p1):
    plot([p0[1], p1[1]], [p0[2], p1[2]], '-b')

figure(1)
clf()

ioff()
for plane in planes:
    p0 = np.dot(T, plane.p0)
    p1 = np.dot(T, plane.p1)
    p2 = np.dot(T, plane.p2)
    plotline(p0, p1)
    plotline(p0, p2)
    plotline(p2, p1)
    
for p in radiator:
    p = np.dot(T, p)
    plot([p[1]], [p[2]], '.r')

xlim(-1,1)
ylim(-.1,1)
ion()
show()

