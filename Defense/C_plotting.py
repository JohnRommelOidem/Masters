import matplotlib.pyplot as plt
import matplotlib
import os
from mpl_toolkits import mplot3d

from B_numerics import *

matplotlib.rcParams.update(
    {'font.size': 12,
     'lines.linewidth': 1,
     'figure.dpi': 300,
     'markers.fillstyle': 'full',
     'image.cmap': 'turbo',
     "figure.figsize": (4.5, 3)}
)
turbo_cmap = matplotlib.cm.get_cmap('turbo')


figsize1 = (3.5, 5.6)
figsize2 = (7, 4)
figsize3 = (4.5, 3)
figsize4 = (4.5, 4.5)
figsize5 = (8, 8)
figsize6 = (20, 20)
figsize7 = (5.5, 4.5)


folder = os.path.expanduser("~/Desktop/MS Thesis") + '/Figures/'


def savefig(name: str, path_name: str = '', pdf=True):
    os.makedirs(os.path.dirname(f'{folder}{path_name}{name}.pdf'), exist_ok=True)
    if pdf:
        plt.savefig(f'{folder}{path_name}{name}.pdf', bbox_inches="tight")
    else:
        plt.savefig(f'{folder}{path_name}{name}.png', bbox_inches="tight")


def save(name: str, array: ndarray, path_name: str = ''):
    os.makedirs(os.path.dirname(f'{folder}{path_name}{name}.pdf'), exist_ok=True)
    np.save(f'{folder}{path_name}{name}', array)


def load(name: str, path_name: str = ''):
    return np.load(f'{folder}{path_name}{name}.npy')


class Plot:
    """
    Generic class for plots
    """

    def __init__(self, path_name=''):
        self.path_name = path_name

    def save(self, name: str, array: ndarray):
        save(name, array, self.path_name+'arrays/')

    def savefig(self, name: str, pdf=True):
        savefig(name, self.path_name+'figures/', pdf=pdf)

    def load(self, name: str):
        return load(name, self.path_name+'arrays/')


class MNCTrajectoryPlot(MNCValues):
    def __init__(
            self,
            L: float, E: float, T: float,
            r: float, theta: float = 0, z: float = 0, pr: float = 0, pz: float = 0,
            dt: float = 0.01, unknown: str = 'pz',
            is_E: bool = True
    ):
        super().__init__(L, E, is_E)

        if unknown == 'pz':
            pz = np.sqrt(2*self.get_kinetic_energy(r, z)-pr**2)
        elif unknown == 'pr':
            pr = np.sqrt(2*self.get_kinetic_energy(r, z)-pz**2)
        else:
            raise ValueError('invalid unknown coordinate')
        self.dt = dt
        self.qp = np.array((
            (r, theta, z),
            (pr, self.L, pz)
        ), float)
        self.t_points = np.arange(0, T, dt)
        self.q_points = np.empty((self.qp.shape[1], self.t_points.shape[0]))
        self.p_points = np.empty_like(self.q_points)

    def get_trajectory(self):
        self.q_points, self.p_points = MNCintegrate(
            self.qp, self.q_points, self.p_points, self.t_points, mnc_de, self.dt
        )
        return self.q_points, self.p_points


class MBHTrajectoryPlot(MBHValues):
    def __init__(
            self,
            r_esc: float, p: float, E: float, T: float,
            r: float, th: float = PI/2, pr: float = 0,
            dt: float = 0.01
    ):
        super().__init__(r_esc, p, E)
        r_bounds = self.get_r_bounds()
        if len(r_bounds) == 1:
            if not (1 < r <= r_bounds[0]):
                print(r_bounds)
                raise ValueError(f"r is invalid")
        else:
            if not (r_bounds[-2] < r <= r_bounds[-1]):
                print(r_bounds)
                raise ValueError(f"r is invalid")

        pr = np.sqrt(MBHescape_energy(self.B, self.L)) * pr / (1 - 1 / r)
        pth = -r*np.sqrt((self.E-self.get_potential_energy(r, th)-((1-1/r)*pr)**2)/(1-1/r))
        self.qp = np.array((
            (r, th),
            (pr, pth)
        ), float)
        self.T = T
        self.dt = dt
        self.t_points = np.arange(0, T, dt)
        self.q_points = np.empty((self.qp.shape[1], self.t_points.shape[0]))
        self.p_points = np.empty_like(self.q_points)

    def get_trajectory(self):
        self.q_points, self.p_points = MBHintegrate(
            self.qp, self.q_points, self.p_points, self.t_points, self.B, self.L, self.E, self.dt
        )
        return self.q_points, self.p_points, self.t_points


class NewtonTrajectoryPlot(NewtonValues):
    def __init__(
            self,
            L: float, E: float, T: float,
            r: float, theta: float = 0, z: float = 0, pr: float = 0, pz: float = 0,
            dt: float = 0.01, unknown: str = 'pz'
    ):
        super().__init__(L, E)

        if unknown == 'pz':
            pz = np.sqrt(2*self.get_kinetic_energy(r, z)-pr**2)
        elif unknown == 'pr':
            pr = np.sqrt(2*self.get_kinetic_energy(r, z)-pz**2)
        else:
            raise ValueError('invalid unknown coordinate')
        self.dt = dt
        self.qp = np.array((
            (r, theta, z),
            (pr, self.L, pz)
        ), float)
        self.t_points = np.arange(0, T, dt)
        self.q_points = np.empty((self.qp.shape[1], self.t_points.shape[0]))
        self.p_points = np.empty_like(self.q_points)

    def get_trajectory(self):
        self.q_points, self.p_points = MNCintegrate(
            self.qp, self.q_points, self.p_points, self.t_points, newton_de, self.dt
        )
        return self.q_points, self.p_points


class MagnetTrajectoryPlot(MagnetValues):
    def __init__(
            self,
            L: float, E: float, T: float,
            r: float, theta: float = 0, z: float = 0, pr: float = 0, pz: float = 0,
            dt: float = 0.01, unknown: str = 'pz'
    ):
        super().__init__(L, E)

        if unknown == 'pz':
            pz = np.sqrt(2*self.get_kinetic_energy(r, z)-pr**2)
        elif unknown == 'pr':
            pr = np.sqrt(2*self.get_kinetic_energy(r, z)-pz**2)
        else:
            raise ValueError('invalid unknown coordinate')
        self.dt = dt
        self.qp = np.array((
            (r, theta, z),
            (pr, self.L, pz)
        ), float)
        self.t_points = np.arange(0, T, dt)
        self.q_points = np.empty((self.qp.shape[1], self.t_points.shape[0]))
        self.p_points = np.empty_like(self.q_points)

    def get_trajectory(self):
        self.q_points, self.p_points = MNCintegrate(
            self.qp, self.q_points, self.p_points, self.t_points, magnet_de, self.dt
        )
        return self.q_points, self.p_points