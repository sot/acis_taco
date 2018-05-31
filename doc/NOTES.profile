# These two lines are standard and probably already there.
import IPython.ipapi
ip = IPython.ipapi.get()

# These two are the important ones.
import line_profiler
ip.expose_magic('lprun', line_profiler.magic_lprun)

import numpy as np
p_earth_body = np.array([-11913349.37481491,   1600513.79810546,   6787847.04879577])

orbit_xyz = -p_earth_body
att = [0,0,0]
import taco
lprun -f taco.calc_earth_vis -f taco.plane_line_intersect vis,illum,rays = taco.calc_earth_vis(orbit_xyz, att) 

import taco2
lprun -f taco2.calc_earth_vis execfile('test_refl2.py')


OR NOT??
---+++ Code timing performance

NRAW has written a very fast code that runs about 100 times faster than the TLA code.  This was initially perplexing because both codes are fundamentally ray-trace codes that determine if points on the radiator surface have visibility to a set (typically 100) of points on the Earth surface.  The TLA code is written in Python but spends most of its time within a highly-optimized C subroutine to determine the ray intersections.

It would appear that the difference stems from an approximation within the NRAW code below in the main routine to calculate the illumination:
<verbatim>
void calcEarthIllum(...)
{ 
.. 
 //determine the max angle and create a grid of points on the Earth
  max_ang=calcMaxPitch(acisangle_r,distance);

 // Some setup and then loop over the 100 rays to the Earth
}
</verbatim>

Here the max pitch angle is calculated using the single roll angle for the center of the Earth.  This max pitch angle is applied as a filter within the rays loop to determine which rays have visibility to the Earth.  However, near perigee the effective roll angle for rays to the Earth varies significantly so this approximation breaks down.  The <verb>calcMaxPitch()</verb> is (I believe) computationally expensive and corresponds effectively to the C ray-trace code in the TLA code.
