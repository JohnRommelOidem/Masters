import matplotlib.pyplot as plt
import numpy as np

from C_plotting import *
from D1_poincare import *
from D2_escape_time import *
from D3_SALI import *
import timeit


def PlotNewtonian():
    def U(r):
        return 3/r**2-6/r

    def r_bounds(U):
        return (-6+PLUS_MINUS*np.sqrt(36+12*U))/(2*U)

    dr = 0.01
    r_points = np.arange(0.4, 10, dr)
    plt.figure(figsize=figsize3)
    plt.axhline(linestyle="dashed", c="k")

    ell_r = np.arange(*r_bounds(-1), dr)
    plt.plot(ell_r, np.full_like(ell_r, -1), label="elliptic")

    par_r = np.arange(r_bounds(0.2)[0], r_points[-1], dr)
    plt.plot(par_r, np.full_like(par_r, 0.2), label="hyperbolic")

    plt.plot(r_points, U(r_points), label="potential", c="k")
    plt.xlim(right=r_points[-1])
    plt.ylim(top=(U(r_points[0])))
    plt.legend()
    savefig("Newtonian", "figures/")


def plot3D_1():
    L = 1
    E = 0.5
    r = 1.7844903
    pr = 0
    T = 1000
    plt.figure(figsize=figsize3)
    plot = PlotTrajectory(L, E, T, r, pr=pr)
    plot = PlotTrajectory(L, get_h_eq(L), T, plot.get_r_eq(), pr=pr, is_E=False)
    plot.plot_3d_trajectory()
    plot.savefig(f'{plot3D_1.__name__}')


def plot3D_2():
    L = 1
    E = 0.5
    r = 0.494551951
    pr = 0
    T = 500
    plt.figure(figsize=figsize3)
    plot = PlotTrajectory(L, E, T, r, pr=pr)
    plot.qp[-1, -1] = 0
    plot.plot_3d_trajectory()
    plot.savefig(f'{plot3D_2.__name__}')


def plot3D_3():
    L = 1
    E = 0.5
    r = 1.8
    pr = 0
    T = 1000
    plt.figure(figsize=figsize3)
    plot = PlotTrajectory(L, E, T, r, pr=pr)
    plot.plot_3d_trajectory()
    plot.savefig(f'{plot3D_3.__name__}')


def plot3D_4():
    L = 1
    E = 0.5
    r = 1.85
    pr = 0
    T = 1000
    plt.figure(figsize=figsize3)
    plot = PlotTrajectory(L, E, T, r, pr=pr)
    plot.plot_3d_trajectory()
    plot.savefig(f'{plot3D_4.__name__}')


def plot3D_5():
    L = 1.
    E = 1.1
    r = np.sqrt(2 * np.abs(L))
    T = 100
    plt.figure(figsize=figsize3)
    plot = PlotTrajectory(L, E, T, r)
    plot.plot_3d_trajectory()
    plot.savefig(f'{plot3D_5.__name__}')


def plot3D_6():
    L = 1.
    E = 1.1
    r = np.sqrt(2 * np.abs(L)) + 0.2
    T = 550
    plt.figure(figsize=figsize3)
    plot = PlotTrajectory(L, E, T, r)
    plot.plot_3d_trajectory()
    plot.savefig(f'{plot3D_6.__name__}')


def plot3D_7():
    L = 5.0
    E = 1.2
    r = 3.85
    T = 1_000
    plt.figure(figsize=figsize3)
    plot = PlotTrajectory(L, E, T, r)
    plot.plot_3d_trajectory()
    plot.savefig(f"trapped_3d")


def poster_escaping():
    L = 1.
    P = 1.1
    r = np.sqrt(2 * np.abs(L))-0.1
    T = 40
    plt.figure(figsize=figsize8)
    plot = PlotTrajectory(L, P, T, r)
    plot.qp[1, 2] = -plot.qp[1, 2]
    plot.plot_asym()
    plot.plot_trajectory()
    r = plot.q_points[0]
    r = r[len(r)//2:]
    plt.axvline(np.max(r), linestyle='--', c='g', zorder=2)
    plt.axvline(np.min(r), linestyle='--', c='g', zorder=2)
    plt.ylim(top=0.)
    ax = plt.gca()

    # Hide X and Y axes label marks
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_tick_params(labelleft=False)

    # Hide X and Y axes tick marks
    ax.set_xticks([])
    ax.set_yticks([])
    plot.savefig(f'poster_escaping', pdf=False)
    plt.close("all")


def poster_inescaping():
    L = 1.
    P = 1.1
    r = np.sqrt(2 * np.abs(L)) + 0.3
    T = 100
    plt.figure(figsize=figsize8)
    plot = PlotTrajectory(L, P, T, r)
    plot.qp[1, 2] = -plot.qp[1, 2]
    plot.plot_asym()
    plot.plot_trajectory()
    plt.ylim(top=0.)
    r = plot.q_points[0]
    r = r[len(r)//2:]
    plt.axvline(np.max(r), linestyle='--', c='g', zorder=2)
    plt.axvline(np.min(r), linestyle='--', c='g', zorder=2)
    ax = plt.gca()

    # Hide X and Y axes label marks
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_tick_params(labelleft=False)

    # Hide X and Y axes tick marks
    ax.set_xticks([])
    ax.set_yticks([])
    plot.savefig(f'poster_inescaping', pdf=False)
    plt.close("all")


def poster_sample_chaotic():
    L = 1
    E = 0.5
    r = 1.85
    T = 20
    plt.clf()
    plt.figure(figsize=figsize8)
    plt.xlabel('$r$')
    plt.ylabel('$z$')
    plot = PlotTrajectory(L, E, T, r-0.02)
    plot.plot_trajectory()
    plot = PlotTrajectory(L, E, T, r)
    plot.plot_trajectory(color="red")

    ax = plt.gca()
    # Hide X and Y axes label marks
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_tick_params(labelleft=False)
    # Hide X and Y axes tick marks
    ax.set_xticks([])
    ax.set_yticks([])

    plot.savefig("poster_chaotic_orbit", pdf=False)
    plt.clf()


def poster_sample_ordered():
    L = 1
    E = 0.5
    r = 1.7844903
    T = 20

    plt.figure(figsize=figsize8)
    plt.xlabel('$r$')
    plt.ylabel('$z$')
    plot = PlotTrajectory(L, E, T, 1.8)
    plot.plot_trajectory()
    plot = PlotTrajectory(L, E, T, r)
    plot.plot_trajectory(color="red")
    ax = plt.gca()

    # Hide X and Y axes label marks
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_tick_params(labelleft=False)

    # Hide X and Y axes tick marks
    ax.set_xticks([])
    ax.set_yticks([])
    plot.savefig("poster_periodic_orbit", pdf=False)
    plt.clf()


def test():
    L = 1.0
    E = 1.05
    resolution = 300
    T = 1_000
    plot = EscapePlot(
        L,
        E,
        resolution,
        T
    )
    time_mesh = plot.load(f"escape_time/_sample")
    pass_mesh = plot.load(f"escape_pass/_sample")
    r_mesh = np.copy(plot.r_mesh)
    pr_mesh = plot.pr_mesh
    print(f"{np.min(r_mesh)}, {np.max(r_mesh)}, {np.min(pr_mesh)}, {np.max(pr_mesh)}")
    time_mesh[pass_mesh != 2] = np.nan
    time_mesh[time_mesh > 200] = np.nan
    r_mesh[time_mesh != time_mesh] = np.nan
    r_center = 1.1
    pr_center = 0.25
    #r_inesc = 0.8475439238215814
    #pr_inesc = 0.09497932029522116
    index = np.unravel_index(np.nanargmin(time_mesh), time_mesh.shape)
    r = r_mesh[index]
    pr = pr_mesh[index]
    T = time_mesh[index]
    print(f"r_inesc = {r}\n" f"pr_inesc = {pr}\n" f"T = {T}")
    PlotTrajectory(L, E, T, r, pr=pr).plot_trajectory()
    plt.show()
    plt.clf()


def sample_SALI():
    L = 1.0
    E = 0.42
    sali = PlotSALI(L, E, resolution=100)
    sali.plot_SALI(name="hires", load_array=True)


def sample_SALI500():
    L = 1.0
    E = 0.42
    sali = PlotSALI(L, E, resolution=500)
    sali.plot_SALI(name="hires_500", load_array=True)

def sample_LCN():
    L = 1.0
    E = 0.42
    sali = PlotSALI(L, E, resolution=100)
    sali.plot_LCN(name="hires_LCN", load_array=True)


def LCN_threshold():
    L = 5.0
    E = 1.0
    sali = PlotSALI(L, E, resolution=100)
    sali.plot_LCN(name="LCN_threshold")


def SALI_threshold():
    L = 5.0
    E = 1.0
    sali = PlotSALI(L, E, resolution=100)
    sali.plot_SALI(name="SALI_threshold")


def fanxion():
    L = 1.0
    E = 0.42
    sali = PlotSALI(L, E, resolution=100)
    SALI_MESH = sali.load(name="hires_LCN")

    size = np.shape(SALI_MESH)[0]
    plt.plot(sali.r_mesh[size//2, :], SALI_MESH[size//2, :])
    plt.show()

    """
    plt.xlabel("t")
    plt.ylabel("LCN")

    L = 1
    E = 0.5
    r = 1.7844903
    pr = 0
    plot = PlotTrajectory(L, E, 1000, r, pr)
    q_points, p_points, var_sum, var_diff, SALI, tpoints, LCN = calculate_SALI_trajectory(qp=plot.qp, T=plot.T)
    plt.plot(tpoints, LCN, label="periodic", alpha=0.8)

    r = 1.8
    plot = PlotTrajectory(L, E, 1000, r, pr)
    q_points, p_points, var_sum, var_diff, SALI, tpoints, LCN = calculate_SALI_trajectory(qp=plot.qp, T=plot.T)
    plt.plot(tpoints, LCN, label="quasi periodic", alpha=0.8)

    r = 1.85
    plot = PlotTrajectory(L, E, 1000, r, pr)
    q_points, p_points, var_sum, var_diff, SALI, tpoints, LCN = calculate_SALI_trajectory(qp=plot.qp, T=plot.T)
    plt.plot(tpoints, LCN, label="chaotic", alpha=0.8)

    plt.legend()
    plot.savefig("LCN_trajectory")
    plt.close("all")
    """


def main():
    sample_SALI500()


if __name__ == '__main__':
    main()
