__all__ = ["TransferAbstract"]

from orbit import OrbitPosition
from math import cos, sin, acos, sqrt, pi, atan, tan, atan2
from numpy import sign
from utilities import wrap_to_2pi
from copy import deepcopy


class NonCoaxialOrbitsError(Exception):
    pass


class NonOrbitObjectPassed(Exception):
    pass


class TransferAbstract(object):
    def __init__(self, orbit_i: OrbitPosition, orbit_f: OrbitPosition) -> None:
        """
        :argument
            orbit_i, orbit_f   -   1x1 Orbit object (OrbitPosition)

        :raises
            NonOrbitObjectPassed, if at least one of the objects passed is not an instance of OrbitPosition
        """

        if type(orbit_i) is not OrbitPosition or type(orbit_f) is not OrbitPosition:
            raise NonOrbitObjectPassed(f"Orbits must be OrbitPosition instances: type orbit_i: {type(orbit_i)}, type orbit_f: {type(orbit_f)}")

        self.orbit_i, self.orbit_f = orbit_i, orbit_f
        self.actual: OrbitPosition = deepcopy(self.orbit_i)

        self.states: list = []

    def get_full_states(self) -> list:
        return [*[{"state": self.orbit_i, "description": "IO"}], *self.states, *[{"state": self.orbit_f, "description": "FO"}]]

    # noinspection PyTypeChecker
    def plot_orbits(self, states: tuple = (), dim: int = 3, figsize: tuple = (8, 8), only_orbit: bool = True) -> None:
        full_states: list = [x for x in self.get_full_states() if "from_to" not in x.keys()]
        states = states if states else tuple(range(len(full_states)))

        if max(states) >= len(full_states):
            raise IndexError("Too many elements requested")

        fig, ax = OrbitPosition.get_fig_ax(f"Orbits ({min(states) + 1} -> {max(states) + 1})", dim=dim, figsize=figsize)

        for state in states:
            additional: dict = {"orbit_label": full_states[state]["description"]}
            full_states[state]["state"].plot_orbit(mode='add', dim=dim, figsize=figsize, fig=fig, ax=ax, only_orbit=only_orbit, **additional)
        else:
            OrbitPosition.plot(fig, ax, mode="show")

    # noinspection PyTypeChecker
    def plot_transfer(self, states: tuple = (), dim: int = 3, figsize: tuple = (8, 8)) -> None:
        full_states: list = [x for x in self.get_full_states() if "from_to" in x.keys()]
        states = states if states else tuple(range(len(full_states)))

        if max(states) >= len(full_states):
            raise IndexError("Too many elements requested")

        fig, ax = OrbitPosition.get_fig_ax(f"Transfer Plot ({min(states) + 1} -> {max(states) + 1})", dim=dim, figsize=figsize)

        for state in states:
            additional: dict = {"orbit_label": full_states[state]["description"],
                                "theta_i": full_states[state]["from_to"][0],
                                "theta_f": full_states[state]["from_to"][1]}

            full_states[state]["state"].plot_orbit(mode='add', dim=dim, figsize=figsize, fig=fig, ax=ax, only_orbit=True, **additional)
        else:
            OrbitPosition.plot(fig, ax, mode="show")

    def waiting(self, maneuver: str, actual: OrbitPosition, to_: float, from_: float = None) -> None:
        theta_i: float = actual.theta if from_ is None else from_
        theta_f: float = to_

        actual.theta = theta_f
        actual.update_rv()
        self.states.append({
            "state": actual,
            "description": f"{maneuver}-W",
            "from_to": (theta_i, theta_f),
            "dt": [self.compute_time(actual.a, actual.e, theta_i, theta_f, actual.mu)]
        })

    def wait(self, theta_f: float = None) -> None:
        theta_f = self.orbit_f.theta if theta_f is None else theta_f

        self.actual = deepcopy(self.actual)

        self.waiting("OR", self.actual, theta_f)

    @staticmethod
    def compute_time(a: float, e: float, theta_i: float, theta_f: float, mu: float) -> float:
        T: float = sqrt((a ** 3) / mu) * 2*pi

        def compute_time_from_periapsis(theta: float) -> float:
            E: float = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(theta / 2))
            if theta <= pi:
                return (sqrt((a ** 3) / mu)) * (E - (e * sin(E)))
            else:
                return T + (sqrt((a ** 3) / mu)) * (E - (e * sin(E)))

        t_i: float = compute_time_from_periapsis(theta_i)
        t_f: float = compute_time_from_periapsis(theta_f)

        if theta_f > theta_i:
            return t_f - t_i
        else:
            return t_f - t_i + T

    def del_last(self, num: int = 1) -> None:
        self.states = self.states[:-num]
        self.actual = self.states[-1]["state"] if self.states else self.orbit_i

    def del_all(self) -> None:
        self.states = []
        self.actual = deepcopy(self.orbit_i)

    def change_plane(self, where: str, i_f: float = None, OM_f: float = None) -> None:
        if where not in ["first", "second"]:
            raise ValueError("'where' must be 'first' or 'second'")

        if not i_f or OM_f:
            i_f: float = self.orbit_f.i
            OM_f: float = self.orbit_f.OM

        self.actual = deepcopy(self.actual)

        # 1 - Rotations
        dOM: float = OM_f - self.actual.OM
        di: float = i_f - self.actual.i

        # 2 - Direction change
        alpha: float = sign(dOM) * acos(sin(self.actual.i) * sin(i_f) * cos(dOM) + cos(self.actual.i) * cos(i_f))

        # 3 - True Anomaly
        cos_u_i: float = (cos(alpha)*cos(self.actual.i)-cos(i_f))/(sin(alpha)*sin(self.actual.i)) if dOM*di > 0 else (-cos(alpha)*cos(self.actual.i)+cos(i_f))/(sin(alpha)*sin(self.actual.i))
        sin_u_i: float = (sin(dOM) / sin(alpha)) * sin(i_f)
        u_i: float = atan2(sin_u_i, cos_u_i)

        cos_u_f: float = (-cos(alpha)*cos(i_f)+cos(self.actual.i))/(sin(alpha)*sin(i_f)) if dOM*di > 0 else (cos(alpha)*cos(i_f)-cos(self.actual.i))/(sin(alpha)*sin(i_f))
        sin_u_f: float = (sin(dOM) / sin(alpha)) * sin(self.actual.i)
        u_f: float = atan2(sin_u_f, cos_u_f)

        theta_i: float = u_i - self.actual.om if dOM*di > 0 else 2*pi - u_i - self.actual.om
        theta_f: float = wrap_to_2pi(theta_i) if where == "second" else wrap_to_2pi(theta_i + pi)

        # >> Move to theta
        self.waiting("PC", self.actual, theta_f)

        self.actual = deepcopy(self.actual)

        # 4 - New Argument of periapsis
        self.actual.om = wrap_to_2pi(u_f - theta_i) if dOM*di > 0 else wrap_to_2pi(2*pi - u_f - theta_i)
        self.actual.OM = OM_f
        self.actual.i = i_f

        # >> Perform maneuver
        self.actual.update_rv()
        self.states.append({
            "state": self.actual,
            "description": "PC",
            "dv": [abs(2*(sqrt(self.actual.mu/(self.actual.a*(1-self.actual.e**2)))*(1+self.actual.e*cos(self.actual.theta)))*sin(alpha/2))]
        })

    def rotate_periapsis(self, where: str, how: str, om_f: float = None) -> None:
        if where not in ["first", "second"]:
            raise ValueError("'where' must be 'first' or 'second'")

        if how not in ["same", "opposite"]:
            raise ValueError("'how' must be 'first' or 'second'")

        if not om_f:
            om_f: float = self.orbit_f.om if how == "same" else wrap_to_2pi(pi + self.orbit_f.om)

        self.actual = deepcopy(self.actual)

        # 1 - Rotation entity
        dom: float = om_f - self.actual.om

        # 2 - Where perform the maneuver
        theta_f: float = wrap_to_2pi(dom/2) if where == "first" else wrap_to_2pi((dom/2) + pi)

        # >> Move to theta
        self.waiting("PR", self.actual, theta_f)

        self.actual = deepcopy(self.actual)

        # 3 - Calculate new orbit's theta
        self.actual.theta = 2*pi - self.states[-1]["state"].theta
        self.actual.om = om_f

        # >> Perform maneuver
        self.actual.update_rv()
        self.states.append({
            "state": self.actual,
            "description": "PR",
            "dv": [abs(2*(sqrt(self.actual.mu/(self.actual.a*(1-self.actual.e**2)))*self.actual.e*sin(dom/2)))]
        })

    def reshape_orbit(self, where: str, a_f: float = None, e_f: float = None) -> None:
        if where not in ["pericentre", "apocentre"]:
            raise ValueError("'where' must be 'pericentre' or 'apocentre'")

        def useful_quantities(a: float, e: float) -> (float, float):
            p: float = a * (1 - e**2)
            r_p: float = p/(1 + e)
            r_a: float = p/(1 - e)
            return r_p, r_a

        def get_v_theta(a: float, e: float, theta: float):
            p: float = a * (1 - e**2)
            return sqrt(self.actual.mu/p) * (1 + e*cos(theta))

        def discriminate(rp_i, ra_i, rp_f, ra_f, om_i, conc_: int, man_: int) -> dict:
            if conc_ not in range(1, 3) or man_ not in range(1, 3):
                raise AttributeError("'conc' and 'man' must be between 1 and 2 included")

            map_dict: dict = {
                1: {  # Vettori eccentricità concordi
                    1: {  # P -> A
                        "rp_t": min([rp_i, ra_f]),
                        "ra_t": max([rp_i, ra_f]),
                        "om_t": om_i if rp_i < ra_f else wrap_to_2pi(om_i + pi),
                        "theta_i": 0,
                        "theta_t": (0, pi) if rp_i < ra_f else (pi, 0),
                        "theta_f": pi
                    },
                    2: {  # A -> P
                        "rp_t": min([ra_i, rp_f]),
                        "ra_t": max([ra_i, rp_f]),
                        "om_t": wrap_to_2pi(om_i + pi) if ra_i < rp_f else om_i,
                        "theta_i": pi,
                        "theta_t": (0, pi) if ra_i < rp_f else (pi, 0),
                        "theta_f": 0
                    }
                },
                2: {  # Vettori eccentricità discordi
                    1: {  # P -> P
                        "rp_t": min([rp_i, rp_f]),
                        "ra_t": max([rp_i, rp_f]),
                        "om_t": om_i if rp_i < rp_f else wrap_to_2pi(om_i + pi),
                        "theta_i": 0,
                        "theta_t": (0, pi) if rp_i < rp_f else (pi, 0),
                        "theta_f": 0
                    },
                    2: {  # A -> A
                        "rp_t": min([ra_i, ra_f]),
                        "ra_t": max([ra_i, ra_f]),
                        "om_t": wrap_to_2pi(om_i + pi) if ra_i < ra_f else om_i,
                        "theta_i": pi,
                        "theta_t": (0, pi) if ra_i < ra_f else (pi, 0),
                        "theta_f": pi
                    }
                }
            }

            return map_dict[conc_][man_].values()

        if not a_f or not e_f:
            orbit_f = deepcopy(self.orbit_f)
        else:
            orbit_f = deepcopy(self.actual)
            orbit_f.a = a_f
            orbit_f.e = e_f

        self.actual = deepcopy(self.actual)

        if self.actual.om == orbit_f.om:
            conc, man = (1, 1) if where == "pericentre" else (1, 2)
        elif self.actual.om == wrap_to_2pi(orbit_f.om + pi):
            conc, man = (2, 1) if where == "pericentre" else (2, 2)
        else:
            raise NonCoaxialOrbitsError(f"om_i={self.actual.om} while om_f={orbit_f.om}")

        rp_t, ra_t, om_t, theta_f, from_to, new_theta = discriminate(*useful_quantities(self.actual.a, self.actual.e),
                                                                     *useful_quantities(orbit_f.a, orbit_f.e),
                                                                     self.actual.om, conc, man)
        self.waiting("RE", self.actual, theta_f)

        a_t: float = (ra_t + rp_t) / 2
        e_t: float = (ra_t - rp_t) / (ra_t + rp_t)

        dv1: float = abs(get_v_theta(a_t, e_t, from_to[0]) - get_v_theta(self.actual.a, self.actual.e, self.actual.theta))

        self.actual = deepcopy(self.actual)
        self.actual.a = a_t
        self.actual.e = e_t
        self.actual.om = om_t

        self.waiting("RE-T", self.actual, to_=from_to[1], from_=from_to[0])

        self.actual = deepcopy(self.actual)

        self.actual.a = orbit_f.a
        self.actual.e = orbit_f.e
        self.actual.theta = new_theta
        self.actual.om = wrap_to_2pi(om_t + pi) if conc == 2 and man == 1 else om_t

        dv2: float = abs(get_v_theta(self.actual.a, self.actual.e, self.actual.theta) - get_v_theta(a_t, e_t, from_to[1]))

        self.actual.update_rv()
        self.states.append({
            "state": self.actual,
            "description": "RE",
            "dv": [dv1, dv2]
        })

    def get_combination_array(self, where_CP: str, where_RP: str, type_RP: str, where_RE: str) -> list:
        # 1 - Clear all
        self.del_all()

        # 2 - Change Plane
        self.change_plane(where_CP)

        # 3 - Rotate Periapsis
        self.rotate_periapsis(where_RP, type_RP)

        # 4 - Reshape Orbit
        self.reshape_orbit(where_RE)

        # 5 - Wait Maneuver
        self.wait()

        # 6 - Get Data
        dv_list: list = [x["dv"] for x in self.get_full_states() if "dv" in x]
        dt_list: list = [x["dt"] for x in self.get_full_states() if "dt" in x]

        dv: float = sum([dv for x in dv_list for dv in x])
        dt: float = sum([dt for x in dt_list for dt in x])

        return [where_CP, where_RP, type_RP, where_RE, dv, dt, dv_list, dt_list]

    def _get_all_combinations_plot(self):
        for where_CP in ["first", "second"]:
            self.change_plane(where_CP)
            for RP in [("first", "same"), ("second", "opposite"), ("first", "opposite"), ("second", "same")]:
                self.rotate_periapsis(RP[0], RP[1])
                for where_RE in ["pericentre", "apocentre"]:
                    self.reshape_orbit(where_RE)
                    self.wait()
                    self.plot_transfer()
                    self.del_last(4)
                self.del_last(2)
            self.del_last(2)
