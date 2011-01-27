import numpy as np
from numpy import sin, cos, radians, arctan2, sqrt
import Ska.quatutil

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
        phys_x = (x - self.phys_x0) * self.img_pix_scale
        phys_y = (y - self.phys_y0) * self.img_pix_scale
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
        eci = self.phys2eci(phys_x, phys_y)
        return eci

    def phys2eci(self, phys_x, phys_y, sun_eci):
        r, phi = self.phys2polar(phys_x, phys_y)
        theta = radians(r)
        sin_theta = sin(theta)
        eci_img = np.array([cos(theta),
                            sin(phi) * sin_theta,
                            cos(phi) * sin_theta])  # OFLS uses -cos(phi)
        self.q_x_to_antisun = Ska.quatutil.quat_x_to_vec(-sun_eci)
        eci = np.dot(self.q_x_to_antisun.transform, eci_img)
        return eci

    def img2sky(self, x, y, sun_eci):
        phys_x, phys_y = self.img2phys(x, y)
        ra, dec = self.phys2sky(phys_x, phys_y, sun_eci)
        return ra, dec

    def phys2sky(self, phys_x, phys_y, sun_eci):
        eci = self.phys2eci(phys_x, phys_y, sun_eci)
        ra, dec = Ska.quatutil.eci2radec(eci)
        return ra, dec

if __name__ == '__main__':
    a = AntiSun(25., 25., (180-45) / 27.)
    print a.img2phys(25., 25.), a.img2phys(25., 30.), a.img2phys(30., 25.)
    print a.img2polar(25., 25.), a.img2polar(25., 30.), a.img2polar(30., 25.)
    sun_eci = np.array([-1.0, .0, 0])

    print a.phys2polar(0., 0.), a.phys2polar(0., 0.1), a.phys2polar(0.1, 0.)
    print a.phys2eci(0., 0., sun_eci), a.phys2eci(0., 0.1, sun_eci), a.phys2eci(0.1, 0., sun_eci)
    print a.phys2sky(0., 0., sun_eci), a.phys2sky(0., 0.1, sun_eci), a.phys2sky(0.1, 0., sun_eci)
    print

    sun_eci = np.array([0.0, -1.0, 0.000])
    print a.phys2eci(0., 0., sun_eci), a.phys2eci(0., 0.1, sun_eci), a.phys2eci(0.1, 0., sun_eci)
    print a.phys2sky(0., 0., sun_eci), a.phys2sky(0., 0.1, sun_eci), a.phys2sky(0.1, 0., sun_eci)
    print

    sun_eci = np.array([0.0, -1.0, -1.0])
    print a.phys2eci(0., 0., sun_eci), a.phys2eci(0., 0.1, sun_eci), a.phys2eci(0.1, 0., sun_eci)
    print a.phys2sky(0., 0., sun_eci), a.phys2sky(0., 0.1, sun_eci), a.phys2sky(0.1, 0., sun_eci)
    print

    sun_eci = np.array([-1.0, .0, -1.0])
    if 0:
        print a.phys2eci(0., 0., sun_eci), a.phys2eci(0., 0.1, sun_eci), a.phys2eci(0.1, 0., sun_eci)
        print a.phys2sky(0., 0., sun_eci), a.phys2sky(0., 0.1, sun_eci), a.phys2sky(0.1, 0., sun_eci)

