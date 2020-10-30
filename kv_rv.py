import numpy as np
from math import pi, acos, sin, cos, sqrt
from numpy import cross, dot, matmul, transpose
from numpy.linalg import norm


def rv2kp(rr: np.array or list, vv: np.array or list, mu: float) -> (float, float, float, float, float, float):
    x = np.array([1, 0, 0])
    y = np.array([0, 1, 0])
    z = np.array([0, 0, 1])

    # Necessary quantities
    r: float = norm(rr)
    v: float = norm(vv)
    E: float = (v ** 2) / 2 - mu / r
    hh: np.array = cross(rr, vv)
    h: float = norm(hh)
    ee: np.array = cross(vv, hh) / mu - rr / r
    N: np.array = cross(z, hh) / norm(cross(z, hh))

    # 1 - Semimajor Axis
    a: float = -mu / (2 * E)

    # 2 - Eccentricity
    e: float = norm(ee)

    # 3 - Inclination
    i: float = acos(dot(z, hh) / h)

    # 4 - Right ascension of the ascending node
    if dot(y, N) >= 0:
        OM: float = acos(dot(x, N) / norm(x) * norm(N))
    else:
        OM: float = 2 * pi - acos(dot(x, N) / norm(x) * norm(N))

    # 5 - Argument of periapsis
    if dot(z, ee) >= 0:
        om = acos(dot(N, ee) / e)
    else:
        om = 2 * pi - acos(dot(N, ee) / e)

    # 6 - True Anomaly
    if dot(rr, vv) >= 0:
        theta: float = acos(dot(rr, ee) / (r * e))
    else:
        theta: float = 2 * pi - acos(dot(rr, ee) / (r * e))

    return a, e, i, OM, om, theta


def kp2rv(a: float, e: float, i: float, OM: float, om: float, theta: float, mu: float) -> (np.array, np.array):
    # 1 - Position and velocity in perifocal
    p: float = a * (1 - (e ** 2))
    r: float = p / (1 + e * cos(theta))

    rr_PF: np.array = np.array([r * cos(theta), r * sin(theta), 0])
    vv_PF: np.array = sqrt(mu / p) * np.array([-sin(theta), e + cos(theta), 0])

    # 2 - Rotation Matrix
    R_OM: np.array = np.array([[cos(OM), sin(OM), 0], [-sin(OM), cos(OM), 0], [0, 0, 1]])
    R_i: np.array = np.array([[1, 0, 0], [0, cos(i), sin(i)], [0, -sin(i), cos(i)]])
    R_om: np.array = np.array([[cos(om), sin(om), 0], [-sin(om), cos(om), 0], [0, 0, 1]])

    R_GE: np.array = transpose(matmul(matmul(R_om, R_i), R_OM))

    # 3 - Position and velocity in geocentric equatorial
    rr_GE: np.array = matmul(R_GE, rr_PF)
    vv_GE: np.array = matmul(R_GE, vv_PF)

    return rr_GE, vv_GE
