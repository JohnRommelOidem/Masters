from A_integrator import *
import matplotlib.pyplot as plt

PLUS_MINUS = np.array([1, -1])
PI = np.pi


def cylindrical_to_cartesian(array):
    r, theta, z = array
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return np.array((x, y, z), float)


def cartesian_to_cylindrical(array):
    x, y, z = array
    r = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan(y / x)
    return np.array((r, theta, z), float)


def sphere_to_cyl(R, th):
    return R*np.sin(th), R*np.cos(th)


def cyl_to_sphere(r, z):
    return np.sqrt(r**2+z**2), np.arctan(r/z)


def mnc_potential_energy(L, r, z):
    R = np.sqrt(r * r + z * z)
    return (L / r - r / 2) ** 2 / 2 - 1 / R


def newton_potential_energy(L, r, z):
    R = np.sqrt(r * r + z * z)
    return (L / r) ** 2 / 2 - 1 / R


def magnet_potential_energy(L, r):
    return (L / r - r / 2) ** 2 / 2


def get_h_eq(L: float):
    """
    Minimum Hamiltonian energy for a given L.
    """
    r_eq = np.roots([1, 0, 0, 4, -4 * L ** 2])
    r_eq = np.abs(
        r_eq[
            np.logical_and(r_eq >= 0, np.isreal(r_eq))
        ]
    )[0]
    return mnc_potential_energy(L, r_eq, 0)


def energy_h_to_E(L: float, h: float):
    r_eq = np.roots([1, 0, 0, 4, -4 * L ** 2])
    r_eq = np.abs(r_eq[np.logical_and(r_eq >= 0, np.isreal(r_eq))])[0]
    h_eq = mnc_potential_energy(L, r_eq, 0)
    h_esc = (np.abs(L) - L) / 2
    E = (h - h_eq) / (h_esc - h_eq)
    return E, h_eq, h_esc


def energy_E_to_h(L: float, E: float):
    r_eq = np.roots([1, 0, 0, 4, -4 * L ** 2])
    r_eq = np.abs(r_eq[np.logical_and(r_eq >= 0, np.isreal(r_eq))])[0]
    h_eq = mnc_potential_energy(L, r_eq, 0)
    h_esc = (np.abs(L) - L) / 2
    h = h_eq + E * (h_esc - h_eq)
    return h, h_eq, h_esc


class MNCValues:
    def __init__(self, L: float, energy: float, is_E: bool = True):
        if is_E:
            self.E = energy
            self.h, self.h_eq, self.h_esc = energy_E_to_h(L, energy)
        else:
            self.h = energy
            self.E, self.h_eq, self.h_esc = energy_h_to_E(L, energy)
        if self.E < 0:
            raise ValueError('energy given is below the minimum')
        self.L = L

    def get_potential_energy(self, r: float, z: float):
        return mnc_potential_energy(self.L, r, z)

    def get_kinetic_energy(self, r: float, z: float):
        return self.h-mnc_potential_energy(self.L, r, z)

    def get_r_eq(self):
        r_eq = np.roots([1, 0, 0, 4, -4 * self.L ** 2])
        r_eq = np.abs(r_eq[np.logical_and(r_eq >= 0, np.isreal(r_eq))])[0]
        return r_eq

    def get_r_bounds(self):
        r_bounds = np.roots([1, 0, -4 * (self.L + 2 * self.h), -8, 4 * self.L ** 2])
        r_bounds = np.sort(np.abs(r_bounds[np.logical_and(r_bounds > 0, np.isreal(r_bounds))]))
        return r_bounds

    def plot_zvc(self, z_bounds: ndarray, r_buff: float = 0.1, resolution: int = 100):
        r_bounds = self.get_r_bounds()
        r_mesh, z_mesh = np.meshgrid(
            np.linspace(
                r_bounds[0] - r_buff * np.diff(r_bounds),
                r_bounds[1] + r_buff * np.diff(r_bounds),
                resolution
            ),
            np.linspace(np.min(z_bounds), np.max(z_bounds), resolution)
        )
        H = self.get_potential_energy(r_mesh, z_mesh)
        cs = plt.contour(r_mesh, z_mesh, H, np.array([self.h]), colors='k', linestyles='--')
        return cs.allsegs[0]

    def get_asymptotes(self):
        if self.E > 1:
            return np.sqrt(2 * (2 * self.h + self.L) + PLUS_MINUS * 4 * np.sqrt(self.h ** 2 + self.L * self.h))
        # there is only a single asymptote if the energy is h_esc itself
        elif self.E == 1.:
            return np.array([np.sqrt(2 * (2 * self.h + self.L))])
        elif self.E < 1:
            raise ValueError('Energy should be above the minimum energy of escape.')


def MBHpotential_energy(B, L, r, th):
    return (1-1/r)*(1+(L/(r*np.sin(th))-B*(r*np.sin(th)))**2)


def MBHpotential_energy_cyl(B, L, r, z):
    return (1-1/np.sqrt(r**2+z**2))*(1+(L/r-B*r)**2)


def MBHescape_energy(B, L):
    return 1+2*B*(abs(L)-L)


def MBHreparametrization(r_esc: float, p: float, E: float):
    L_plus = r_esc*((3-p)*(3*p-1))**(1/4)/(np.sqrt(2*(4*p**2-9*p+3+np.sqrt((3*p-1)*(3-p)))))
    B_plus = L_plus/r_esc**2
    E_plus = E
    if p > (5+np.sqrt(13))/4:
        L_minus = -r_esc*((3-p)*(3*p-1))**(1/4)/(np.sqrt(2*(4*p**2-9*p+3-np.sqrt((3*p-1)*(3-p)))))
        B_minus = -L_minus/r_esc**2
        E_minus = E*(1-4*B_minus*L_minus)
        return B_plus, L_plus, E_plus, B_minus, L_minus, E_minus
    else:
        return B_plus, L_plus, E_plus


class MBHValues:
    def __init__(self, r_esc: float, p: float, E: float):
        self.L = r_esc*((3-p)*(3*p-1))**(1/4)/np.sqrt(2*(4*p**2-9*p+3+np.sqrt((3*p-1)*(3-p))))
        self.B = self.L / r_esc ** 2
        self.r_esc = r_esc
        if E < 0:
            raise ValueError('energy given is below the minimum')
        else:
            self.E = E * MBHescape_energy(self.B, self.L)

    def get_potential_energy(self, r: float, th: float):
        return MBHpotential_energy(self.B, self.L, r, th)

    def get_potential_energy_cyl(self, r: float, z: float):
        return MBHpotential_energy_cyl(self.B, self.L, r, z)

    def get_kinetic_energy(self, r: float, th: float):
        return MBHpotential_energy(self.B, self.L, r, th) - self.E

    def get_asymptotes(self):
        if self.E < MBHescape_energy(self.B, self.L):
            raise ValueError('Energy should be above the minimum energy of escape.')
        else:
            asym = np.roots(
                [self.B**2, 0, 1-self.E-2*self.L*self.B, 0, self.L**2]
            )
            return np.sort(np.abs(asym[np.logical_and(asym > 0, np.isreal(asym))]))

    def get_r_bounds(self):
        r_bounds = np.roots(
            [self.B ** 2, -self.B ** 2, 1 - self.E - 2 * self.B * self.L, 2 * self.B * self.L - 1, self.L ** 2,
             -self.L ** 2]
        )
        r_bounds = np.sort(np.abs(r_bounds[np.logical_and(r_bounds > 0, np.isreal(r_bounds))]))
        return r_bounds

    def get_r_eq(self):
        r_eq = np.roots([2*self.B**2, -self.B**2, 0, 1-2*self.B*self.L, -2*self.L**2, 3*self.L**2])
        r_eq = np.sort(np.abs(r_eq[np.logical_and(r_eq > 0, np.isreal(r_eq))]))
        return r_eq

    def get_U_eq(self):
        r_eq = self.get_r_eq()
        if r_eq.shape[0] == 1:
            return self.get_potential_energy(r_eq[0], PI/2)
        elif r_eq.shape[0] == 2:
            return np.array([self.get_potential_energy(r_eq[0], PI/2), self.get_potential_energy(r_eq[1], PI/2)])
        else:
            #print("No equilibrium points found")
            return np.array([])

    def get_pr_bounds(self):
        r_eq = self.get_r_eq()
        U_eq = self.get_U_eq()
        if U_eq.shape[0] == 2 and U_eq[0] > self.E:
            return PLUS_MINUS*np.sqrt(self.E-self.get_potential_energy(r_eq[1], PI/2))
        else:
            #print("No equilibrium points found")
            return PLUS_MINUS*np.sqrt(self.E)

    def plot_zvc(self, z_bounds: ndarray, r_buff: float = 0.1, resolution: int = 100, r_bounds=None, color='k'):
        """
        Plots the (r,z) curve with zero velocity for a given energy. Known as the zero velocity curve or zvc.
        :param z_bounds: any array whose minimum and maximum values gives up to where the curve is plotted along z.
        :param r_buff: amount of extra space along r outside the curve.
        :param resolution: resolution
        :return: all the points of the zvc; it is also plotted.
        """
        if not r_bounds:
            r_bounds = self.get_r_bounds()
        if len(r_bounds) == 1:
            bounds = np.array([0, r_bounds[0]])
        else:
            bounds = r_bounds[-2:]

        r_mesh, z_mesh = np.meshgrid(
            np.linspace(
                bounds[0] - r_buff * np.diff(bounds),
                bounds[1] + r_buff * np.diff(bounds),
                resolution
            ),
            np.linspace(np.min(z_bounds), np.max(z_bounds), resolution)
        )
        E = self.get_potential_energy(*cyl_to_sphere(r_mesh, z_mesh))
        cs = plt.contour(r_mesh, z_mesh, E, np.array([self.E]), colors=color, linestyles='--')
        if len(r_bounds) == 1:
            plt.contour(r_mesh, z_mesh, r_mesh**2+z_mesh**2, np.array([1]), colors="k")
        return cs.allsegs[0]


class NewtonValues:
    def __init__(self, L: float, E: float):
        self.E = E
        self.L = L

    def get_potential_energy(self, r: float, z: float):
        return newton_potential_energy(self.L, r, z)

    def get_kinetic_energy(self, r: float, z: float):
        return self.E-newton_potential_energy(self.L, r, z)

    def get_r_bounds(self):
        r_bounds = (-1+PLUS_MINUS*np.sqrt(1+2*self.E*self.L**2))/(2*self.E)
        return r_bounds


class MagnetValues:
    def __init__(self, L: float, E: float):
        self.E = E
        self.L = L

    def get_potential_energy(self, r: float, z: float):
        return newton_potential_energy(self.L, r, z)

    def get_kinetic_energy(self, r: float, z: float):
        return self.E-newton_potential_energy(self.L, r, z)

    def get_r_bounds(self):
        r_bounds = np.sqrt(2*self.E+2*self.L**2)+PLUS_MINUS*np.sqrt(2*self.E)
        return r_bounds
