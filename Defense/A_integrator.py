import numpy as np
from math import sin, cos, sqrt

ndarray = np.ndarray


def mnc_de(qp: ndarray):
    q, p = qp
    r, t, z = q
    pr, l, pz = p

    return np.array((
        (pr,
         l / r ** 2 - 1 / 2,
         pz),
        (l ** 2 / r ** 3 - r / 4 - r / (r ** 2 + z ** 2) ** 1.5,
         0.,
         -z / (r ** 2 + z ** 2) ** 1.5)
    ))


def newton_de(qp: ndarray):
    q, p = qp
    r, t, z = q
    pr, l, pz = p

    return np.array((
        (pr,
         l / r ** 2,
         pz),
        (l ** 2 / r ** 3 - r / (r ** 2 + z ** 2) ** 1.5,
         0.,
         -z / (r ** 2 + z ** 2) ** 1.5)
    ))


def magnet_de(qp: ndarray):
    q, p = qp
    r, t, z = q
    pr, l, pz = p

    return np.array((
        (pr,
         l / r ** 2 - 1 / 2,
         pz),
        (l ** 2 / r ** 3 - r / 4,
         0.,
         0.)
    ))


def mbh_de(qp: ndarray, B, L, E):
    q, p = qp

    r, th = q
    pr, pth = p
    if r <= 1:
        return np.array([[0., 0.], [0., 0.]])
    else:
        return np.array((
            (
                pr*(1-1/r),
                pth/r**2
            ),
            (
                -(E/(1-1/r)**2+pr**2)/(2*r**2)+pth**2/r**3+L**2/(r**3*sin(th)**2)-B**2*r*sin(th)**2,
                L**2*cos(th)/(r**2*sin(th)**3)-B**2*r**2*sin(th)*cos(th)
            )
        ))


optimal_a = np.array([
    0.5153528374311229364,
    -0.085782019412973646,
    0.4415830236164665242,
    0.1288461583653841854
])
optimal_b = np.array([
    0.1344961992774310892,
    -0.2248198030794208058,
    0.7563200005156682911,
    0.3340036032863214255
])


def MNCintegrator(qp: ndarray, df, dt: float = 0.01):
    """
    4th order numerical integrator step that optimally conserves energy for the 4th order.

    :param qp: Coordinates
    :param df: A function that returns the time derivative of the coordinates
    :param dt: Time step
    :return: Coordinates after integrating one time step
    """
    output = np.copy(qp)
    for i, _ in enumerate(optimal_a):
        output[1] += dt * optimal_b[i] * df(output)[1]
        output[0] += dt * optimal_a[i] * df(output)[0]
    return output


def MNCintegrate(qp: ndarray, q_points, p_points, t_points, diff_eqn, dt: float = 0.01):
    # getting the value of each integration step
    for t, _ in enumerate(t_points):
        qp = MNCintegrator(qp, diff_eqn, dt)
        q_points[:, t], p_points[:, t] = qp[:, :]
    return q_points, p_points


def MBHintegrator(qp: ndarray, B, L, E, dt: float = 0.01):
    """
    4th order numerical integrator step that optimally conserves energy for the 4th order.

    :param qp: Coordinates
    :param df: A function that returns the time derivative of the coordinates
    :param dt: Time step
    :return: Coordinates after integrating one time step
    """
    output = np.copy(qp)
    for i, _ in enumerate(optimal_a):
        output[1] += dt * optimal_b[i] * mbh_de(output, B, L, E)[1]
        output[0] += dt * optimal_a[i] * mbh_de(output, B, L, E)[0]
    return output


def MBHintegrate(qp: ndarray, q_points, p_points, t_points, B, L, E, dt: float = 0.01):
    # getting the value of each integration step
    for t, _ in enumerate(t_points):
        qp = MBHintegrator(qp, B, L, E, dt)
        q_points[:, t], p_points[:, t] = qp[:, :]
    return q_points, p_points
