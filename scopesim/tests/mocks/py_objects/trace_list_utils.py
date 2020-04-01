import numpy as np


def rotate(x, y, x0, y0, angle):
    """ Rotate a line by ``angle`` [deg] around the point (x0, y0)"""
    angle /= 57.29
    xnew = x0 + (x - x0) * np.cos(angle) - (y - y0) * np.sin(angle)
    ynew = y0 + (x - x0) * np.sin(angle) + (y - y0) * np.cos(angle)

    return xnew, ynew


def shear(x, y, x0, y0, mx=0., my=0.):
    """ Shear a line by factors (mx, my) relative to point (x0, y0)"""
    xnew = x + mx * (y - y0)
    ynew = y + my * (x - x0)

    return xnew, ynew


def stretch(x, y, nx=0., ny=0.):
    """ Stretch a line by factors (nx, ny)"""
    xnew = x * nx
    ynew = y * ny

    return xnew, ynew


def translate(x, y, dx=0., dy=0.):
    """ Offset a line by (dx, dy)"""
    xnew = x + dx
    ynew = y + dy

    return xnew, ynew


def curvify(x, y, x0, y0):
    """
    Curves a straight line onto the circumference of a circle

    The curved line is the arc of a circle centred on the point (x0, y0)
    with a radius equal to the minimum distance between the (x0, y0) and the
    straight line.

    Note: The length of the arc is not equal to the length of the straight line

    Parameters
    ----------
    x, y : array
        The x and y coordinates to curvify
    x0, y0 : float
        The central point of the circle which defines the radius of curvature

    Returns
    -------
    xnew, ynew : array
        The altered coordinates

    """
    dx = x - x0
    dy = y - y0
    dist = (dx**2 + dy**2)**0.5
    radius = np.min(dist)
    scale = radius / dist
    xnew = x0 + dx * scale
    ynew = y0 + dy * scale

    return xnew, ynew
