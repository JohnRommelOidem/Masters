from abc import ABC
from manim import *
from C_plotting import *

PLOT_LABEL_FONT: float = 36


class MNCenter(Dot3D, ABC):
    def __init__(self):
        super().__init__(color=ORANGE, radius=1.0)


class MNCTestParticle2D(Dot, ABC):
    config = {
        'stroke_opacity': 1,
        'color': WHITE
    }

    def __init__(
            self,
            L: float, E: float,
            r: float, z: float = 0, theta: float = 0,
            pr: float = 0,
            fill_color=BLUE, axis=Axes(), speed=ValueTracker(1.),
            **kwargs
    ):
        trajectory_class = MNCTrajectoryPlot(
            L, E, 0, r, theta, z, pr
        )
        self.qp = trajectory_class.qp
        self.axis = axis
        self.speed = speed
        super().__init__(
            point=self.axis.coords_to_point(self.qp[0, 0], self.qp[0, 2]),
            fill_color=fill_color, **self.config, **kwargs
        )

    @staticmethod
    def updater(mob, dt):
        dt *= mob.speed.get_value()
        mob.previous_velocity = [mob.qp[1, 0], mob.qp[1, 2]]
        mob.qp = MNCintegrator(mob.qp, mnc_de, dt)
        mob.move_to(
            mob.axis.coords_to_point(
                mob.qp[0, 0], mob.qp[0, 2]
            )
        )

    def move(self):
        self.add_updater(self.updater)

    def stop(self):
        self.remove_updater(self.updater)


class MNCTestParticle3D(Dot3D, ABC):
    config = {
        'radius': 0.1,
        'fill_opacity': 1
    }

    def __init__(
            self,
            L: float, E: float,
            r: float, z: float = 0, theta: float = 0,
            pr: float = 0,
            axis=ThreeDAxes(), speed: int = 1, fill_color=BLUE_E, color=BLUE, **kwargs
    ):
        trajectory_class = MNCTrajectoryPlot(
            L, E, 0, r, theta, z, pr, dt=0.01
        )
        self.qp = trajectory_class.qp
        self.axis = axis
        self.speed = ValueTracker(speed)
        super().__init__(
            point=self.axis.coords_to_point(*cylindrical_to_cartesian(self.qp[0])), fill_color=fill_color, color=color,
            **self.config, **kwargs
        )
        #self.move()

    @staticmethod
    def updater(mob, dt):
        dt *= mob.speed.get_value()
        mob.qp = MNCintegrator(mob.qp, mnc_de, dt)
        mob.move_to(
            mob.axis.coords_to_point(
                *cylindrical_to_cartesian(mob.qp[0])
            )
        )

    def move(self):
        self.add_updater(self.updater)

    def stop(self):
        self.remove_updater(self.updater)


class MNCDotParticle3D(Dot, ABC):
    config = {
        'radius': 0.1,
        'fill_opacity': 1
    }

    def __init__(
            self,
            L: float, E: float,
            r: float, z: float = 0, theta: float = 0,
            pr: float = 0,
            axis=ThreeDAxes(), speed: int = 1, fill_color=BLUE_E, color=BLUE, **kwargs
    ):
        trajectory_class = MNCTrajectoryPlot(
            L, E, 0, r, theta, z, pr, dt=0.01
        )
        self.qp = trajectory_class.qp
        self.axis = axis
        self.speed = ValueTracker(speed)
        super().__init__(
            point=self.axis.coords_to_point(*cylindrical_to_cartesian(self.qp[0])), fill_color=fill_color, color=color,
            **self.config, **kwargs
        )
        #self.move()

    @staticmethod
    def updater(mob, dt):
        dt *= mob.speed.get_value()
        mob.qp = MNCintegrator(mob.qp, mnc_de, dt)
        mob.move_to(
            mob.axis.coords_to_point(
                *cylindrical_to_cartesian(mob.qp[0])
            )
        )

    def move(self):
        self.add_updater(self.updater)

    def stop(self):
        self.remove_updater(self.updater)


class NewtonTestParticle3D(Dot3D, ABC):
    config = {
        'radius': 0.1,
        'fill_opacity': 1
    }

    def __init__(
            self,
            L: float, E: float,
            r: float, z: float = 0, theta: float = 0,
            pr: float = 0,
            axis=ThreeDAxes(), speed: int = 1, fill_color=YELLOW_E, color=YELLOW, **kwargs
    ):
        trajectory_class = NewtonTrajectoryPlot(
            L, E, 0, r, theta, z, pr
        )
        self.qp = trajectory_class.qp
        self.axis = axis
        self.speed = ValueTracker(speed)
        super().__init__(
            point=self.axis.coords_to_point(*cylindrical_to_cartesian(self.qp[0])), fill_color=fill_color, color=color,
            **self.config, **kwargs
        )
        #self.move()

    @staticmethod
    def updater(mob, dt):
        dt *= mob.speed.get_value()
        mob.qp = MNCintegrator(mob.qp, newton_de, dt)
        mob.move_to(
            mob.axis.coords_to_point(
                *cylindrical_to_cartesian(mob.qp[0])
            )
        )

    def move(self):
        self.add_updater(self.updater)


class MagnetTestParticle3D(Dot3D, ABC):
    config = {
        'radius': 0.1,
        'fill_opacity': 1
    }

    def __init__(
            self,
            L: float, E: float,
            r: float, z: float = 0, theta: float = 0,
            pr: float = 0, pz: float = 0,
            axis=ThreeDAxes(), speed: float = 1, fill_color=YELLOW_E, color=YELLOW, **kwargs
    ):
        trajectory_class = MagnetTrajectoryPlot(
            L, E, 0, r, theta, z, pr, pz
        )
        self.qp = trajectory_class.qp
        self.axis = axis
        self.speed = ValueTracker(speed)
        super().__init__(
            point=self.axis.coords_to_point(*cylindrical_to_cartesian(self.qp[0])), fill_color=fill_color, color=color,
            **self.config, **kwargs
        )
        #self.move()

    @staticmethod
    def updater(mob, dt):
        dt *= mob.speed.get_value()
        mob.qp = MNCintegrator(mob.qp, magnet_de, dt)
        mob.move_to(
            mob.axis.coords_to_point(
                *cylindrical_to_cartesian(mob.qp[0])
            )
        )

    def move(self):
        self.add_updater(self.updater)


class NewtonTrajectory3D(VMobject, ABC):

    def __init__(self, L, E, r, T, axis=ThreeDAxes()):
        self.L = L
        self.E = E
        self.r = r
        self.T = T
        self.axis = axis
        super().__init__()
        self.start_new_path(axis.c2p(r.get_value(), 0., 0.))
        self.set_color(YELLOW)
        self.updater(self)
        self.move()

    @staticmethod
    def updater(mob):
        traj = NewtonTrajectoryPlot(mob.L.get_value(), mob.E.get_value(), mob.T.get_value(), mob.r.get_value())
        traj.get_trajectory()
        traj_points = np.empty((traj.q_points.shape[1], traj.q_points.shape[0]))
        for i, _ in enumerate(traj_points):
            traj_points[i] = mob.axis.c2p(*cylindrical_to_cartesian(traj.q_points[:, i]))
        mob.set_points_smoothly(traj_points)

    def move(self):
        self.add_updater(self.updater)


class MagnetTrajectory3D(VMobject, ABC):

    def __init__(self, L, E, r, axis, T):
        self.L = L
        self.E = E
        self.r = r
        self.T = T
        self.axis = axis
        super().__init__()
        self.start_new_path(axis.c2p(r.get_value(), 0., 0.))
        self.set_color(YELLOW)
        self.updater(self)
        self.move()

    @staticmethod
    def updater(mob):
        traj = MagnetTrajectory3D(mob.L.get_value(), mob.E.get_value(), mob.T.get_value(), mob.r.get_value())
        traj.get_trajectory()
        traj_points = np.empty((traj.q_points.shape[1], traj.q_points.shape[0]))
        for i, _ in enumerate(traj_points):
            traj_points[i] = mob.axis.c2p(*cylindrical_to_cartesian(traj.q_points[:, i]))
        mob.set_points_smoothly(traj_points)

    def move(self):
        self.add_updater(self.updater)


class MNCTrajectory2D(VMobject, ABC):

    def __init__(self, L, E, r, axis, T, pr=ValueTracker(0.)):
        self.L = L
        self.E = E
        self.r = r
        self.T = T
        self.pr = pr
        self.axis = axis
        super().__init__()
        self.start_new_path(axis.c2p(r.get_value(), 0., 0.))
        self.set_color(BLUE)
        self.updater(self)
        self.move()

    @staticmethod
    def updater(mob):
        traj = MNCTrajectoryPlot(
            mob.L.get_value(), mob.E.get_value(), mob.T.get_value(), mob.r.get_value(), pr=mob.pr.get_value()
        )
        traj.get_trajectory()
        traj_points = np.empty((traj.q_points.shape[1], traj.q_points.shape[0]))
        for i, _ in enumerate(traj_points):
            traj_points[i] = mob.axis.c2p(traj.q_points[0, i], traj.q_points[2, i])
        mob.set_points_smoothly(traj_points)

    def move(self):
        self.add_updater(self.updater)


class MBHTrajectory2D(VMobject, ABC):

    def __init__(self, r_esc, p, E, r, axis, T, pr=ValueTracker(0.)):
        self.r_esc = r_esc
        self.p = p
        self.E = E
        self.r = r
        self.T = T
        self.pr = pr
        self.axis = axis
        super().__init__()
        self.start_new_path(axis.c2p(r.get_value(), 0., 0.))
        self.set_color(BLUE)
        self.updater(self)
        self.move()

    @staticmethod
    def updater(mob):
        traj = MBHTrajectoryPlot(
            mob.r_esc.get_value(), mob.p.get_value(), mob.E.get_value(),
            mob.T.get_value(), mob.r.get_value(), pr=mob.pr.get_value()
        )
        traj.get_trajectory()
        traj_points = np.empty((traj.q_points.shape[1], 3))
        for i, _ in enumerate(traj_points):
            traj_points[i] = mob.axis.c2p(*sphere_to_cyl(traj.q_points[0, i], traj.q_points[1, i]))
        mob.set_points_smoothly(traj_points)

    def move(self):
        self.add_updater(self.updater)


class MNCTrajectory3D(VMobject, ABC):

    def __init__(self, L, E, r, pr, axis, T):
        self.L = L
        self.E = E
        self.r = r
        self.pr = pr
        self.T = T
        self.axis = axis
        super().__init__()
        self.start_new_path(axis.c2p(r.get_value(), 0., 0.))
        self.set_color(YELLOW)
        self.updater(self)
        self.move()

    @staticmethod
    def updater(mob):
        traj = MNCTrajectoryPlot(
            mob.L.get_value(), mob.E.get_value(), mob.T.get_value(), mob.r.get_value(), pr=mob.pr.get_value()
        )
        traj.get_trajectory()
        traj_points = np.empty((traj.q_points.shape[1], traj.q_points.shape[0]))
        for i, _ in enumerate(traj_points):
            traj_points[i] = mob.axis.c2p(*cylindrical_to_cartesian(traj.q_points[:, i]))
        mob.set_points_smoothly(traj_points)

    def move(self):
        self.add_updater(self.updater)


class BlackHole(VGroup, ABC):
    def __init__(self, **kwargs):
        BH_width = 5
        black_hole = Dot(radius=BH_width * 0.2/2, color=BLACK)
        stream_1 = StreamLines(
            lambda pos: UP*pos[0]+LEFT*pos[1],
            x_range=[-BH_width * 0.2/2, BH_width * 0.2/2, 0.3],
            y_range=[-BH_width * 0.2/2, BH_width * 0.2/2, 0.3],
            color=YELLOW_D,
            noise_factor=0.1
        )
        stream_2 = StreamLines(
            lambda pos: UP*pos[0]+LEFT*pos[1],
            x_range=[-BH_width * 0.2/2, BH_width * 0.2/2, 0.1],
            y_range=[-BH_width * 0.2/2, BH_width * 0.2/2, 0.1],
            color=ORANGE,
            noise_factor=0.1
        )
        #stream_1.start_animation(warm_up=False, flow_speed=1.5, time_width=0.5)
        #stream_2.start_animation(warm_up=False, flow_speed=1.3, time_width=5.0)
        BH_back = Dot(radius=BH_width * 0.25/2, color=ORANGE)
        super().__init__(BH_back, stream_2, stream_1, black_hole)


class MNCZeroVelocityCurve(VGroup, ABC):
    curve_config = {
        'color': GRAY
    }

    def __init__(self, L, E, axis=Axes(), **kwargs):
        self.kwargs = kwargs
        self.L = L
        self.E = E
        self.axis = axis
        self.axis_bounds = self.axis.y_range[:-1]
        super().__init__()

        self.updater(self)
        self.move()

    @staticmethod
    def updater(mob):
        mob.values = MNCValues(mob.L.get_value(), mob.E.get_value())
        mob.h = ValueTracker(mob.values.h)
        mob_copy = VGroup()

        for point_group in mob.values.plot_zvc(mob.axis_bounds):
            curve_points = np.append(point_group, np.zeros((point_group.shape[0], 1)), axis=1)
            curve = VMobject().set_points_smoothly(
                np.array([mob.axis.c2p(*curve_point) for curve_point in curve_points])
            )
            mob_copy.add(curve)
        mob.become(mob_copy)

    def move(self):
        self.add_updater(self.updater)

    def stop(self):
        self.remove_updater(self.updater)


class MNCAsymptotes(VGroup, ABC):
    line_config = {
        'color': RED
    }
    dashed_config = {
        'num_dashes': 15
    }

    def __init__(self, L, E, axis=Axes(), **kwargs):
        self.kwargs = kwargs
        self.L = L
        self.E = E
        self.axis = axis
        self.axis_bounds = self.axis.y_range[:-1]
        super().__init__()
        self.updater(self)
        self.move()

    @staticmethod
    def updater(mob):
        mob.values = MNCValues(mob.L.get_value(), mob.E.get_value())
        mob_copy = VGroup()
        if mob.E.get_value() >= 1:
            for asymptote in mob.values.get_asymptotes():
                asymptote_line = Line(
                    mob.axis.coords_to_point(asymptote, mob.axis_bounds[0], 0),
                    mob.axis.coords_to_point(asymptote, mob.axis_bounds[1], 0),
                    **mob.line_config
                )
                dashed_line = DashedVMobject(asymptote_line, **mob.dashed_config, **mob.kwargs)
                mob_copy.add(dashed_line)
        mob.become(mob_copy)

    def move(self):
        self.add_updater(self.updater)


class MBHZeroVelocityCurve(VGroup, ABC):
    curve_config = {
        'color': GRAY
    }

    def __init__(self, r_esc, p, E, axis=Axes(), **kwargs):
        self.kwargs = kwargs
        self.r_esc = r_esc
        self.p = p
        self.E = E
        self.axis = axis
        self.axis_bounds = self.axis.y_range[:-1]
        super().__init__()

        self.updater(self)
        self.move()

    @staticmethod
    def updater(mob):
        mob.values = MBHValues(mob.r_esc.get_value(), mob.p.get_value(), mob.E.get_value())
        mob_copy = VGroup()

        for point_group in mob.values.plot_zvc(mob.axis_bounds, r_buff=0):
            curve_points = np.append(point_group, np.zeros((point_group.shape[0], 1)), axis=1)
            curve = VMobject().set_points_smoothly(
                np.array([mob.axis.c2p(*curve_point) for curve_point in curve_points])
            )
            mob_copy.add(curve)
        mob.become(mob_copy)

    def move(self):
        self.add_updater(self.updater)

    def stop(self):
        self.remove_updater(self.updater)


class MBHAsymptotes(VGroup, ABC):
    line_config = {
        'color': RED
    }
    dashed_config = {
        'num_dashes': 15
    }

    def __init__(self, r_esc, p, E, axis=Axes(), **kwargs):
        self.kwargs = kwargs
        self.r_esc = r_esc
        self.p = p
        self.E = E
        self.axis = axis
        self.axis_bounds = self.axis.y_range[:-1]
        super().__init__()
        self.updater(self)
        self.move()

    @staticmethod
    def updater(mob):
        mob.values = MBHValues(mob.r_esc.get_value(), mob.p.get_value(), mob.E.get_value())
        mob_copy = VGroup()
        if mob.E.get_value() >= 1:
            for asymptote in mob.values.get_asymptotes():
                asymptote_line = Line(
                    mob.axis.coords_to_point(asymptote, mob.axis_bounds[0], 0),
                    mob.axis.coords_to_point(asymptote, mob.axis_bounds[1], 0),
                    **mob.line_config
                )
                dashed_line = DashedVMobject(asymptote_line, **mob.dashed_config, **mob.kwargs)
                mob_copy.add(dashed_line)
        mob.become(mob_copy)

    def move(self):
        self.add_updater(self.updater)


class Amplitudes(VGroup, ABC):
    def __init__(self, trajectory, axis):
        self.trajectory = trajectory
        self.axis = axis
        super().__init__(
            DashedLine(color=GREEN),
            DashedLine(color=GREEN)
        )
        self.updater(self)
        self.add_updater(self.updater)

    @staticmethod
    def updater(mob):
        points = mob.trajectory.get_all_points()
        points = points[round(points.shape[0]-30_00):]
        min, max, _ = mob.axis.y_range
        mob[0].put_start_and_end_on(
            start=mob.axis.c2p(*(min*UP))[1]*UP + np.max(points[:, 0])*RIGHT,
            end=mob.axis.c2p(*(max*UP))[1]*UP + np.max(points[:, 0])*RIGHT
        )
        mob[1].put_start_and_end_on(
            start=mob.axis.c2p(*(min*UP))[1]*UP + np.min(points[:, 0])*RIGHT,
            end=mob.axis.c2p(*(max*UP))[1]*UP + np.min(points[:, 0])*RIGHT
        )
