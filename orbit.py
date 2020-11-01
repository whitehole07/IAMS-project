__all__ = ["OrbitPosition"]

from kv_rv import rv2kp, kp2rv
from math import pi, cos, sin, sqrt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from utilities import wrap_to_2pi


class OrbitPosition(object):
    def __init__(self, mu: float = 398600.4, **kwargs) -> None:
        """
        :argument
            a, e, i, OM, om, theta, /mu/   -   1x1 Orbital parameters (float)

            or

            rr, vv, /mu/                   -   1x3 Position and velocity (np.array or list), 1x1 Planetary constant (float)

        :raises
            ValueError, if passed wrong arguments
        """

        for key in kwargs.keys():
            if key not in ['a', 'e', 'i', 'OM', 'om', 'theta', 'rr', 'vv', 'mu']:
                raise ValueError("Arguments must be (rr, vv, /mu/) or (a, e, i, OM, om, theta, /mu/)")

        self.__dict__, self.mu = {**self.__dict__, **{f"{key}": value for key, value in kwargs.items()}}, mu

        if len(kwargs) == 6:
            self.rr, self.vv, = kp2rv(self.a, self.e, self.i, self.OM, self.om, self.theta, self.mu)
        elif len(kwargs) == 2:
            self.rr = np.array(self.rr) if type(self.rr) is list else self.rr
            self.vv = np.array(self.vv) if type(self.rr) is list else self.vv

            self.a, self.e, self.i, self.OM, self.om, self.theta = rv2kp(self.rr, self.vv, self.mu)
        else:
            raise ValueError("Arguments must be (rr, vv, /mu/) or (a, e, i, OM, om, theta, /mu/)")

    def get_hh(self, normed: bool = True) -> np.array:
        if normed:
            return np.cross(self.rr, self.vv) / np.linalg.norm(np.cross(self.rr, self.vv))
        else:
            return np.cross(self.rr, self.vv)

    def get_r(self, theta: float) -> float:
        return (self.a * (1 - self.e**2))/(1 + self.e*cos(theta))

    def get_rr(self, theta: float, dim: int = 3) -> np.array:
        if dim == 3:
            return kp2rv(self.a, self.e, self.i, self.OM, self.om, theta, self.mu)[0]
        elif dim == 2:
            r: float = (self.a * (1 - self.e**2))/(1 + self.e*cos(theta))
            return np.array([r*cos(theta), r*sin(theta), 0])
        else:
            raise ValueError("'dim' should be equals to 2 or 3")

    def get_vv(self, theta: float, dim: int = 3, cart: bool = False) -> np.array:
        if dim == 3:
            return kp2rv(self.a, self.e, self.i, self.OM, self.om, theta, self.mu)[1]
        elif dim == 2:
            p: float = self.a * (1 - self.e**2)
            v_theta: float = sqrt(self.mu/p) * (1 + self.e*cos(theta))
            v_r: float = sqrt(self.mu/p) * (self.e*sin(theta))
            if not cart:
                return np.array([v_r, v_theta, 0])
            else:
                return np.array([v_r*cos(theta) - v_theta*sin(theta), v_r*sin(theta) + v_theta*cos(theta), 0])
        else:
            raise ValueError("'dim' should be equals to 3")

    def get_orbit_arc(self, dim: int = 3, theta_i: float = 0.0, theta_f: float = 2*pi, resolution: int = 1000) -> np.array:
        thetas: np.array = wrap_to_2pi(np.linspace(theta_i, theta_f + 2*pi, resolution)) if theta_f < theta_i else np.linspace(theta_i, theta_f, resolution)
        return np.transpose(np.array([self.get_rr(theta, dim=dim) for theta in thetas]))

    @staticmethod
    def get_fig_ax(title: str, dim: int = 3, figsize: tuple = (8, 8)) -> tuple:
        params: dict = {
            'legend.fontsize': 'x-large',
            'figure.figsize': figsize,
            'axes.labelsize': 'x-large',
            'axes.titlesize': 22,
            'axes.labelpad': 10,
            'axes.titleweight': 'normal',
            'xtick.labelsize': 'large',
            'ytick.labelsize': 'large',
            'lines.linewidth': 2.5
        }

        pylab.rcParams.update(params)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d') if dim == 3 else fig.add_subplot(111)

        ax.axis('equal') if dim == 2 else None
        ax.set_xlabel('X [km]')
        ax.set_ylabel('Y [km]')
        ax.set_zlabel('Z [km]') if dim == 3 else None

        ax.view_init(elev=20., azim=125) if dim == 3 else None

        ax.set_title(title)

        ax.plot([0], [0], [0], 'sk', label='Focus', markersize=10) if dim == 3 else ax.plot([0], [0], 'sk', label='Focus', markersize=10)

        return fig, ax

    def plot_orbit(self, orbit_label: str = "Orbit", mode: str = 'show', dim: int = 3, theta_i: float = 0.0, theta_f: float = 2*pi,
                   resolution: int = 1000, figsize: tuple = (8, 8), fig=None, ax=None, **kwargs) -> None:
        if dim not in (2, 3):
            raise ValueError("'dim' should be equals to 2 or 3")

        if fig is None or ax is None:
            fig, ax = self.get_fig_ax("Orbit", dim=dim, figsize=figsize)

        rr_arc: np.array = self.get_orbit_arc(dim=dim, theta_i=theta_i, theta_f=theta_f, resolution=resolution)

        ax.plot(rr_arc[0, :], rr_arc[1, :], rr_arc[2, :],  zdir='z', label=orbit_label, zorder=1) if dim == 3 else ax.plot(rr_arc[0, :], rr_arc[1, :], label=orbit_label, zorder=1)

        if not kwargs["only_orbit"]:
            ap_arc: np.array = np.array(
                [[0, self.get_rr(pi, dim=dim)[0]],
                 [0, self.get_rr(pi, dim=dim)[1]],
                 [0, self.get_rr(pi, dim=dim)[2]]])
            pe_arc: np.array = np.array(
                [[0, self.get_rr(0, dim=dim)[0]],
                 [0, self.get_rr(0, dim=dim)[1]],
                 [0, self.get_rr(0, dim=dim)[2]]])
            n_arc: np.array = np.array(
                [[0, self.get_rr(2 * pi - self.om, dim=dim)[0]],
                 [0, self.get_rr(2 * pi - self.om, dim=dim)[1]],
                 [0, self.get_rr(2 * pi - self.om, dim=dim)[2]]])

            ax.plot(ap_arc[0, :], ap_arc[1, :], ap_arc[2, :], '--', zdir='z', label='Apoapsis') if dim == 3 else ax.plot(ap_arc[0, :], ap_arc[1, :], '--', label='Apoapsis')
            ax.plot(pe_arc[0, :], pe_arc[1, :], pe_arc[2, :], '--', zdir='z', label='Periapsis') if dim == 3 else ax.plot(pe_arc[0, :], pe_arc[1, :], '--', label='Periapsis')
            ax.plot(n_arc[0, :], n_arc[1, :], n_arc[2, :], '--', zdir='z', label='Ascending right') if dim == 3 else ax.plot(n_arc[0, :], n_arc[1, :], '--', label='Ascending right')

        if kwargs["velocity"] is True:
            velocity_kwargs: dict = {key.replace("velocity_", ""): value for key, value in kwargs.items() if key.startswith("velocity_")}

            rr, vv = self.get_rr(self.theta, dim=dim)[:dim], self.get_vv(self.theta, dim=dim, cart=True)[:dim]
            ax.plot(*[[x] for x in rr], 'ok')
            ax.quiver(*rr, *vv, color='black', zorder=2, **velocity_kwargs) if dim == 3 else ax.quiver(*rr, *vv, color='black', zorder=2, **velocity_kwargs)

        if kwargs["angular_momentum"] is True and dim == 3:
            angular_momentum_kwargs: dict = {key.replace("angular_momentum_", ""): value for key, value in kwargs.items() if key.startswith("angular_momentum_")}

            hh: np.array = self.get_hh()
            ax.plot([0], [0], [0], 'ok')
            ax.quiver(0, 0, 0, *hh, color='black', zorder=3, **angular_momentum_kwargs)

        self.plot(fig, ax, mode=mode)

    @staticmethod
    def plot(fig, ax, mode: str = "show") -> None:
        if mode == "show":
            ax.legend(loc='upper right')
            plt.tight_layout()
            plt.show()
        elif mode == "save":
            ax.legend(loc='upper right')
            plt.tight_layout()
            plt.savefig("orbit.png", dpi=fig.dpi)
        elif mode == 'add':
            return
        else:
            raise ValueError(f"'{mode}' mode does not exist: available are 'show', 'save' or 'add'")

    def update_rv(self) -> None:
        self.rr, self.vv, = kp2rv(self.a, self.e, self.i, self.OM, self.om, self.theta, self.mu)

    def update_kp(self) -> None:
        self.a, self.e, self.i, self.OM, self.om, self.theta = rv2kp(self.rr, self.vv, self.mu)
