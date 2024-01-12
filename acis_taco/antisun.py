# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function

import numpy as np
import ska_quatutil
from numpy import arccos, arctan2, cos, degrees, radians, sin, sqrt


class AntiSun(object):
    def __init__(self, phys_x0, phys_y0, img_pix_scale):
        """
        :param phys_x0: center of physical coord system in image coords
        :param phys_y0: center of physical coord system in image coords
        :param img_pix_scale: pixel scale of image coord system in degrees/pixel
        """
        self.phys_x0 = phys_x0  
        self.phys_y0 = phys_y0
        self.img_pix_scale = img_pix_scale 

    def img2phys(self, x, y):
        # rotate physical by 90 deg clockwise to get to image
        phys_y = (x - self.phys_x0) * self.img_pix_scale
        phys_x = -(y - self.phys_y0) * self.img_pix_scale
        return phys_x, phys_y

    def img2polar(self, x, y):
        phys_x, phys_y = self.img2phys(x, y)
        r, phi = self.phys2polar(phys_x, phys_y)
        return r, phi

    def phys2polar(self, phys_x, phys_y):
        # phi = 0 is along the +y_phys axis
        r = sqrt(phys_x**2 + phys_y**2)
        phi = arctan2(-phys_x, phys_y)
        return r, phi

    def img2eci(self, x, y, sun_eci):
        phys_x, phys_y = self.img2phys(x, y)
        eci = self.phys2eci(phys_x, phys_y, sun_eci)
        return eci

    def phys2eci(self, phys_x, phys_y, sun_eci):
        r, phi = self.phys2polar(phys_x, phys_y)
        theta = radians(r)
        sin_theta = sin(theta)
        eci_img = np.array([cos(theta),
                            sin(phi) * sin_theta,
                            cos(phi) * sin_theta])
        self.q_x_to_antisun = ska_quatutil.quat_x_to_vec(-sun_eci)
        eci = np.dot(self.q_x_to_antisun.transform, eci_img)
        return eci

    def img2sky(self, x, y, sun_eci):
        phys_x, phys_y = self.img2phys(x, y)
        ra, dec = self.phys2sky(phys_x, phys_y, sun_eci)
        return ra, dec

    def phys2sky(self, phys_x, phys_y, sun_eci):
        eci = self.phys2eci(phys_x, phys_y, sun_eci)
        ra, dec = ska_quatutil.eci2radec(eci)
        return ra, dec

    def eci2polar(self, eci, sun_eci):
        """
        phi = arctan2(sin(phi), cos(phi))
            = arctan2(sin(theta) * sin(phi), sin(theta) * cos(phi))
            = arctan2(y, z)
        """
        q_x_to_antisun = ska_quatutil.quat_x_to_vec(-sun_eci)
        eci_img = np.dot(q_x_to_antisun.transform.transpose(), eci)
        theta = arccos(eci_img[0])
        phi = arctan2(eci_img[1], eci_img[2])
        r = degrees(theta)
        return r, phi
    
    def eci2phys(self, eci, sun_eci):
        q_x_to_antisun = ska_quatutil.quat_x_to_vec(-sun_eci)
        eci_img = np.dot(q_x_to_antisun.transform.transpose(), eci)
        r = degrees(arccos(eci_img[0]))
        r12 = np.sqrt(eci_img[1]**2 + eci_img[2]**2)
        phys_x = -r * eci_img[1] / r12
        phys_y = r * eci_img[2] / r12
        return phys_x, phys_y

    def phys2img(self, phys_x, phys_y):
        x = phys_y / self.img_pix_scale + self.phys_x0
        y = -phys_x / self.img_pix_scale + self.phys_y0
        return x, y

    def eci2img(self, eci, sun_eci):
        phys_x, phys_y = self.eci2phys(eci, sun_eci)
        x, y = self.phys2img(phys_x, phys_y)
        return x, y

if __name__ == '__main__':
    a = AntiSun(25., 25., (180-45) / 27.)
    print(a.img2phys(25., 25.), a.img2phys(25., 30.), a.img2phys(30., 25.))
    print(a.img2polar(25., 25.), a.img2polar(25., 30.), a.img2polar(30., 25.))
    sun_eci = np.array([-1.0, .0, 0])

    print(a.phys2polar(0., 0.), a.phys2polar(0., 0.1), a.phys2polar(0.1, 0.))
    print(a.phys2eci(0., 0., sun_eci), a.phys2eci(0., 0.1, sun_eci), a.phys2eci(0.1, 0., sun_eci))
    print(a.phys2sky(0., 0., sun_eci), a.phys2sky(0., 0.1, sun_eci), a.phys2sky(0.1, 0., sun_eci))
    print()

    sun_eci = np.array([0.0, -1.0, 0.000])
    print(a.phys2eci(0., 0., sun_eci), a.phys2eci(0., 0.1, sun_eci), a.phys2eci(0.1, 0., sun_eci))
    print(a.phys2sky(0., 0., sun_eci), a.phys2sky(0., 0.1, sun_eci), a.phys2sky(0.1, 0., sun_eci))
    print()

    sun_eci = np.array([0.0, -1.0, -1.0])
    print(a.phys2eci(0., 0., sun_eci), a.phys2eci(0., 0.1, sun_eci), a.phys2eci(0.1, 0., sun_eci))
    print(a.phys2sky(0., 0., sun_eci), a.phys2sky(0., 0.1, sun_eci), a.phys2sky(0.1, 0., sun_eci))
    print()

    phys_x, phys_y = a.img2phys(12.0, 19.0)
    print(phys_x, phys_y)
    print(a.img2polar(12.0, 19.0))
    print(a.img2eci(12.0, 19.0, sun_eci))
    eci = a.phys2eci(phys_x, phys_y, sun_eci)
    print(eci)
    print(a.eci2phys(eci, sun_eci))
    print(a.eci2polar(eci, sun_eci))
    print(a.eci2img(eci, sun_eci))
    print()



