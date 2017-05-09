import time
from collections import OrderedDict

import numpy as np

from Launch import globals

conn = globals.conn
space_center = globals.space_center
vessel = globals.vessel
mu = globals.mu
g0 = globals.g0
# For ease of use
orbref = vessel.orbit.body.non_rotating_reference_frame


# Trigonometric function in degrees :)


def r2d(x):
    return x * 180 / np.pi


def d2r(x):
    return x * np.pi / 180


def cosd(x):
    return np.cos(d2r(x))


def acosd(x):
    np.clip(x, -1, 1)
    return r2d(np.arccos(x))


def sind(x):
    return np.sin(d2r(x))


def asind(x):
    np.clip(x, -1, 1)
    return r2d(np.arcsin(x))


def tand(x):
    return np.tan(d2r(x))


def atand(x):
    return r2d(np.arctan(x))


def atan2d(x, y):
    return r2d(np.arctan2(x, y))


# Another simplified function, and yeah i am that lazy
# This one is for vector operation
def norm(x):
    return np.linalg.norm(x)


def unit(x):
    x = np.asarray(x)
    if norm(x) == 0:
        return x
    else:
        return x / norm(x)


def cross(x, y):
    return np.cross(x, y)


def dot(x, y):
    return np.vdot(x, y)


def vang(x, y):
    x = unit(x)
    y = unit(y)
    return acosd(np.clip(dot(x, y), -1, 1))


class Struct:
    pass


def launch_targeting(periapsis,
                     apoapsis,
                     inclination,
                     lan,
                     true_anomaly=0,
                     slip=-1.5):
    # surfref = vessel.orbit.body.reference_frame
    local_x = [1, 0, 0]
    # local_y = [0, 1, 0]
    local_z = [0, 0, 1]
    apoapsis *= 1000
    periapsis *= 1000
    apoapsis += vessel.orbit.body.equatorial_radius
    periapsis += vessel.orbit.body.equatorial_radius

    semimajor_axis = (apoapsis + periapsis) / 2
    ecc = (apoapsis - periapsis) / (apoapsis + periapsis)
    velocity_periapsis = np.sqrt(
        (mu * apoapsis) / (periapsis * semimajor_axis))
    radius = (semimajor_axis * (1 - ecc ** 2)) / (1 + ecc * cosd(true_anomaly))
    velocity = np.sqrt((velocity_periapsis ** 2) + 2 * mu *
                       ((1 / radius) - (1 / periapsis)))
    angle = acosd((periapsis * velocity_periapsis) / (radius * velocity))
    descending = False
    if inclination < 0:
        descending = True
        inclination = abs(inclination)
    if inclination < vessel.flight().latitude:
        azimuth = 90
    else:
        beta_inertial = asind(cosd(inclination) /
                              cosd(vessel.flight().latitude))
        if descending:
            if beta_inertial <= 90:
                beta_inertial = 180 - beta_inertial
            elif beta_inertial >= 270:
                beta_inertial = 540 - beta_inertial
        earth_velocity = (vessel.orbit.body.rotational_speed *
                          vessel.orbit.body.equatorial_radius)
        velocity_x = (velocity * sind(beta_inertial) - earth_velocity *
                      cosd(vessel.flight().latitude))
        velocity_y = (velocity *
                      cosd(beta_inertial))
        azimuth = atan2d(velocity_x, velocity_y)

    relative_longitude = asind(tand(vessel.flight().latitude) /
                               tand(inclination))
    if descending:
        relative_longitude = 180 - relative_longitude
    if conn.krpc.get_status().version == '0.3.6':
        prime_meridian = vessel.orbit.body.msl_position(0, 0, orbref)
        rotational_angle = atan2d(dot(local_z, prime_meridian),
                                  dot(local_x, prime_meridian))
        if rotational_angle < 0:
            rotational_angle += 360
    else:
        rotational_angle = r2d(vessel.orbit.body.rotation_angle)
    print(rotational_angle)
    geo_longitude = lan + relative_longitude - rotational_angle
    geo_longitude = np.mod(geo_longitude + 360, 360)
    node_angle = geo_longitude - vessel.flight().longitude
    node_angle = np.mod(node_angle + 360 + slip, 360)
    launch_time = (node_angle / 360) * vessel.orbit.body.rotational_period

    rx = np.array([[1, 0, 0],
                   [0, cosd(inclination), -sind(inclination)],
                   [0, sind(inclination), cosd(inclination)]])
    ry = np.array([[cosd(0), 0, sind(0)],
                   [0, 1, 0],
                   [-sind(0), 0, cosd(0)]])
    rz = np.array([[cosd(lan), -sind(lan), 0],
                   [sind(lan), cosd(lan), 0],
                   [0, 0, 1]])

    m0 = np.matmul(rz, ry)
    m1 = np.matmul(m0, rx)
    m2 = np.transpose([0, 0, -1])

    target_plane_normal = np.matmul(m1, m2).transpose()
    target = Struct()
    target.radius = radius
    target.normal = target_plane_normal
    target.angle = angle
    target.velocity = velocity

    return azimuth, launch_time, target


def angle_from_vec(x, ref, angle):
    east = [0, 0, 1]
    north = [0, 1, 0]
    up = [1, 0, 0]
    surface_frame = vessel.surface_reference_frame
    vector = space_center.transform_direction(x, ref, surface_frame)
    if angle == 'pitch':
        return 90 - vang(up, vector)
    elif angle == 'yaw':
        out = atan2d(dot(east, vector), dot(north, vector))
        if out < 0:
            out += 360
        return out


def upfg(vehicle, target, previous):
    gamma = target.angle
    iy = np.asarray(target.normal)
    rdval = target.radius
    vdval = target.velocity
    t = globals.universal_time()
    m = globals.state_mass()
    r = globals.orbital_position()
    r = np.array([r[0], r[2], r[1]])
    v = globals.orbital_velocity()
    v = np.array([v[0], v[2], v[1]])
    cser = previous.cser
    rbias = np.asarray(previous.rbias)
    rd = np.asarray(previous.rd)
    rgrav = np.asarray(previous.rgrav)
    tp = previous.time
    vprev = np.asarray(previous.v)
    vgo = np.asarray(previous.vgo)

    n = len(vehicle)
    md = list()
    ve = list()
    f_t = list()
    a_t = list()
    tu = list()
    tb = list()
    if n > 1 and globals.state_thrust() == 0:
        stage_controller(vehicle)
        return upfg(vehicle, target, previous)

    for i in range(n):
        f_t.append(vehicle[i].fT)
        ve.append(vehicle[i].ve)
        md.append(f_t[i] / ve[i])
        a_t.append(f_t[i] / vehicle[i].m0)
        tu.append(ve[i] / a_t[i])
        tb.append(vehicle[i].maxT)

    dt = t - tp
    dvsensed = v - vprev
    vgo = vgo - dvsensed
    # vgo1 = vgo
    tb[0] = tb[0] - previous.tb
    f_t[0] = globals.state_thrust()
    a_t[0] = f_t[0] / m
    tu[0] = ve[0] / a_t[0]
    l_ = 0
    li_ = list()

    for i in range(n - 1):
        li_.append(ve[i] * np.log(tu[i] / (tu[i] - tb[i])))
        l_ += li_[i]
        if l_ > norm(vgo):
            vehicle.remove(vehicle[-1])
            print('We have more than what we need')
            return upfg(vehicle, target, previous)
    li_.append(norm(vgo) - l_)

    tgoi = list()

    for i in range(n):
        tb[i] = tu[i] * (1 - np.exp(-li_[i] / ve[i]))
        if i == 0:
            tgoi.append(tb[i])
        else:
            tgoi.append(tgoi[i - 1] + tb[i])

    # l1 = li[0]
    tgo = tgoi[n - 1]
    if tgo > 5:
        theta_max = d2r(60)
    else:
        theta_max = d2r(1)

    l_ = 0
    s_ = 0
    j_ = 0
    q_ = 0
    p_ = 0
    h_ = 0
    ji_ = list()
    si_ = list()
    qi_ = list()
    pi_ = list()
    tgoi1 = 0

    for i in range(n):
        if i > 0:
            tgoi1 = tgoi[i - 1]
        ji_.append(tu[i] * li_[i] - ve[i] * tb[i])
        si_.append(tb[i] * li_[i] - ji_[i])
        qi_.append(si_[i] * (tu[i] + tgoi1) - 0.5 * ve[i] * tb[i] ** 2)
        pi_.append(qi_[i] * (tu[i] + tgoi1) - 0.5 * ve[i] * tb[i] ** 2 * (tb[i] / 3 + tgoi1))

        ji_[i] += li_[i] * tgoi1
        si_[i] += l_ * tb[i]
        qi_[i] += j_ * tb[i]
        pi_[i] += h_ * tb[i]

        l_ += li_[i]
        j_ += ji_[i]
        s_ += si_[i]
        q_ += qi_[i]
        p_ += pi_[i]
        h_ = j_ * tgoi[i] - q_

    lamb = unit(vgo)
    # rgrav1 = rgrav
    if previous.tgo != 0:
        rgrav = (tgo / previous.tgo) ** 2 * rgrav
    rgo = rd - (r + v * tgo + rgrav)
    # rgo1 = rgo
    iz = unit(cross(rd, iy))
    # iz1 = iz
    rgoxy = rgo - dot(iz, rgo) * iz
    rgoz = (s_ - dot(lamb, rgoxy)) / dot(lamb, iz)
    rgo = rgoxy + rgoz * iz + rbias
    lambdade = q_ - s_ * j_ / l_
    lambdadot = (rgo - s_ * lamb) / lambdade
    if (norm(lambdadot) * j_ / l_) > theta_max:
        lambdadotmag = theta_max / (j_ / l_)
        lambdadot = unit(lambdadot) * lambdadotmag
        rgo = s_ * lamb + lambdade * lambdadot
    i_f = unit(lamb - lambdadot * j_ / l_)
    phi = np.arccos(dot(i_f, lamb))
    phidot = -phi * l_ / j_
    vthrust = (l_ - 0.5 * l_ * phi ** 2 - j_ * phi *
               phidot - 0.5 * h_ * phidot ** 2) * lamb
    vthrust = vthrust - (l_ * phi + j_ * phidot) * unit(lambdadot)
    rthrust = (s_ - 0.5 * s_ * phi ** 2 - q_ * phi *
               phidot - 0.5 * p_ * phidot ** 2) * lamb
    rthrust = rthrust - (s_ * phi + q_ * phidot) * unit(lambdadot)
    vbias = vgo - vthrust
    rbias = rgo - rthrust

    i_f = [i_f[0], i_f[2], i_f[1]]
    pitch = angle_from_vec(i_f, orbref, 'pitch')
    yaw = angle_from_vec(i_f, orbref, 'yaw')

    rc1 = r - 0.1 * rthrust - vthrust * tgo / 30
    vc1 = v + rthrust * 1.2 / tgo - 0.1 * vthrust
    [rc2, vc2, cser] = cse_routine(rc1, vc1, tgo, cser)
    vgrav = vc2 - vc1
    rgrav = rc2 - rc1 - vc1 * tgo

    rp = r + v * tgo + rgrav + rthrust
    rp = rp - dot(rp, iy) * iy
    rd = rdval * unit(rp)
    ix = unit(rd)
    iz = cross(ix, iy)
    m1 = np.transpose([ix, iy, iz])
    m2 = np.transpose([sind(gamma), 0, cosd(gamma)])
    mx = np.matmul(m1, m2).transpose()
    vd = vdval * mx
    vgo = vd - v - vgrav + vbias

    previous.cser = cser
    previous.rbias = rbias
    previous.rd = rd
    previous.rgrav = rgrav
    previous.tb = previous.tb + dt
    previous.time = t
    previous.tgo = tgo
    previous.v = v
    previous.vgo = vgo

    guidance = Struct()
    guidance.pitch = pitch
    guidance.yaw = yaw
    guidance.tgo = tgo

    return [previous, guidance]


def cse_routine(r0, v0, dt, last):
    r0 = np.asarray(r0)
    v0 = np.asarray(v0)
    if last.dtcp == 0:
        dtcp = dt
    else:
        dtcp = last.dtcp

    xcp = last.xcp
    x = xcp
    a = last.a
    d = last.d
    e = last.e

    kmax = 10
    imax = 10

    if dt > 0:
        f0 = 1
    else:
        f0 = -1

    n = 0
    r0m = norm(r0)

    f1 = f0 * np.sqrt(r0m / mu)
    f2 = 1 / f1
    f3 = f2 / r0m
    f4 = f1 * r0m
    f5 = f0 / np.sqrt(r0m)
    f6 = f0 * np.sqrt(r0m)

    ir0 = r0 / r0m
    v0s = f1 * v0
    sigma0s = dot(ir0, v0s)
    b0 = dot(v0s, v0s) - 1
    alphas = 1 - b0

    xguess = f5 * x
    xlast = f5 * xcp
    xmin = 0
    dts = f3 * dt
    dtlast = f3 * dtcp
    dtmin = 0

    xmax = np.divide(2 * np.pi, np.sqrt(np.abs(alphas)))

    if alphas > 0:
        dtmax = xmax / alphas
        xp = xmax
        ps = dtmax
        while dts >= ps:
            n = n + 1
            dts = dts - ps
            dtlast = dtlast - ps
            xguess = xguess - xp
            xlast = xlast - xp
    else:
        [dtmax, _, _, _] = k_t_t_i(xmax, sigma0s, alphas, kmax)
        if dtmax < dts:
            while dtmax >= dts:
                dtmin = dtmax
                xmin = xmax
                xmax = np.multiply(2, xmax)
                [dtmax, _, _, _] = k_t_t_i(xmax, sigma0s, alphas, kmax)

    if xmin >= xguess or xguess >= xmax:
        xguess = 0.5 * (xmin + xmax)

    [dtguess, _, _, _] = k_t_t_i(xguess, sigma0s, alphas, kmax)

    if dts < dtguess:
        if xguess < xlast < xmax:
            if dtguess < dtlast < dtmax:
                xmax = xlast
                dtmax = dtlast
    else:
        if xmin < xlast < xguess:
            if dtmin < dtlast < dtguess:
                xmin = xlast
                dtmin = dtlast

    [xguess, dtguess, a, d, e] = k_i_l(
        imax, dts, xguess, dtguess, xmin, dtmin, xmax, dtmax,
        sigma0s, alphas, kmax, a, d, e)

    rs = 1 + 2 * (b0 * a + sigma0s * d * e)
    b4 = 1 / rs

    if n > 0:
        # noinspection PyUnboundLocalVariable
        xc = f6 * (xguess + n * xp)
        # noinspection PyUnboundLocalVariable
        dtc = f4 * (dtguess + n * ps)
    else:
        xc = f6 * xguess
        dtc = f4 * dtguess

    last.dtcp = dtc
    last.xcp = xc
    last.a = a
    last.d = d
    last.e = e

    f = 1 - 2 * a
    gs = 2 * (d * e + sigma0s * a)
    fts = -2 * b4 * d * e
    gt = 1 - 2 * b4 * a

    r = r0m * (f * ir0 + gs * v0s)
    v = f2 * (fts * ir0 + gt * v0s)

    return [r, v, last]


def k_t_t_i(xarg, s0s, a, kmax):
    u1 = uss(xarg, a, kmax)
    zs = 2 * u1
    e = 1 - 0.5 * a * zs ** 2
    w = np.sqrt(max(0.5 + e / 2, 0))
    d = w * zs
    a = d ** 2
    b = 2 * (e + s0s * d)
    q = qcf(w)
    t = d * (b + a * q)

    return [t, a, d, e]


def uss(xarg, a, kmax):
    du1 = np.divide(xarg, 4)
    u1 = du1
    f7 = -a * du1 ** 2
    k = 3
    while k < kmax:
        du1 = f7 * du1 / (k * (k - 1))
        u1old = u1
        u1 = u1 + du1
        if u1 == u1old:
            break
        k += 2
    return u1


def qcf(w):
    if w < 1:
        xq = 21.04 - 13.04 * w
    elif w < 4.625:
        xq = (5 / 3) * (2 * w + 5)
    elif w < 13.846:
        xq = (10 / 7) * (w + 12)
    elif w < 44:
        xq = 0.5 * (w + 60)
    elif w < 100:
        xq = 0.25 * (w + 164)
    else:
        xq = 70

    b = 0
    y = (w - 1) / (w + 1)
    j = np.floor(xq)
    b = y / (1 + (j - 1) / (j + 2) * (1 - b))
    while j > 2:
        j = j - 1
        b = y / (1 + (j - 1) / (j + 2) * (1 - b))

    q = 1 / w ** 2 * (1 + (2 - b / 2) / (3 * w * (w + 1)))
    return q


def k_i_l(imax, dts, xguess, dtguess, xmin, dtmin, xmax, dtmax,
          s0s, alphas, kmax, a, d, e):
    i = 1
    while i < imax:
        dterror = dts - dtguess

        if abs(dterror) < 0.0000001:
            break

        [dxs, xmin, dtmin, xmax, dtmax] = si(
            dterror, xguess, dtguess, xmin, dtmin, xmax, dtmax)
        xold = xguess
        xguess = xguess + dxs

        if xguess == xold:
            break

        dtold = dtguess
        [dtguess, a, d, e] = k_t_t_i(xguess, s0s, alphas, kmax)

        if dtguess == dtold:
            break

        i += 1

    return [xguess, dtguess, a, d, e]


def si(dterror, xguess, dtguess, xmin, dtmin, xmax, dtmax):
    etp = 0.0000001
    dtminp = dtguess - dtmin
    dtmaxp = dtguess - dtmax
    if abs(dtminp) < etp or abs(dtmaxp) < etp:
        dxs = 0
    else:
        if dterror < 0:
            dxs = (xguess - xmax) * (dterror / dtmaxp)
            if (xguess + dxs) <= xmin:
                dxs = (xguess - xmin) * (dterror / dtminp)
            xmax = xguess
            dtmax = dtguess
        else:
            dxs = (xguess - xmin) * (dterror / dtminp)
            if (xguess + dxs) >= xmax:
                dxs = (xguess - xmax) * (dterror / dtmaxp)
            xmin = xguess
            dtmin = dtguess

    return [dxs, xmin, dtmin, xmax, dtmax]


def rodrigues(vector, axis, angle):
    vector = np.asarray(vector)
    axis = unit(axis)
    rotated = vector * cosd(angle)
    rotated += cross(axis, vector) * sind(angle)
    rotated += axis * dot(axis, vector) * (1 - cosd(angle))

    return rotated


def throttle_control(vehicle, g_limit, q_limit):
    min_thrust = vessel.available_thrust * vehicle[0].minThrottle
    max_thrust = vessel.available_thrust * vehicle[0].maxThrottle
    if max_thrust == 0 or vehicle[0].minThrottle == 1:
        return 1
    g_thrust = g_limit * globals.state_mass() * g0
    g_throttle = (g_thrust - min_thrust) / (max_thrust - min_thrust)

    q_ratio = globals.state_q() / q_limit
    q_throttle = 1 - 15 * (q_ratio - 1)
    the_throttle = np.minimum(q_throttle, g_throttle)
    the_throttle = np.clip(the_throttle, 0.01, 1)

    return the_throttle


def analyze_vehicle():
    stage = list()
    for part in vessel.parts.all:
        if part.engine is not None:
            stage.append(part.engine.part.decouple_stage)

    stage = list(OrderedDict.fromkeys(stage))
    print(stage)
    m0 = list()
    m1 = list()
    f_t = list()
    ve = list()
    a_t = list()
    tu = list()
    l1 = list()
    tb = list()
    max_throttle = list()
    min_throttle = list()
    mass = 0

    for i in range(len(stage)):
        part_list = vessel.parts.in_decouple_stage(stage[i])
        fuel_name = list()
        thrust = 0
        max_thrust = 0
        min_thrust = 0
        flow_rate = 0
        isp = 0
        for part in part_list:
            mass += part.mass
            if part.engine is not None:
                thrust += part.engine.max_vacuum_thrust
                isp += part.engine.vacuum_specific_impulse * g0
                flow_rate += thrust / isp
                for fuel in part.engine.propellants:
                    fuel_name.append(fuel.name)
                prev_limit = part.engine.thrust_limit
                part.engine.thrust_limit = 0
                min_thrust += part.engine.available_thrust
                part.engine.thrust_limit = prev_limit
                max_thrust += part.engine.available_thrust

        resources_list = vessel.resources_in_decouple_stage(
            stage[i] + 1, False)
        fuel_mass = 0
        fuel_name = list(OrderedDict.fromkeys(fuel_name))
        for fuel in fuel_name:
            fuel_amount = resources_list.amount(fuel)
            fuel_density = resources_list.density(fuel)
            fuel_mass += fuel_amount * fuel_density

        m0.append(mass)
        m1.append(fuel_mass)
        f_t.append(thrust)
        ve.append(thrust / flow_rate)
        a_t.append(thrust / mass)
        tb.append(fuel_mass / flow_rate)
        tu.append(ve[i] * mass / thrust)
        l1.append(ve[i] * np.log(tu[i] / (tu[i] - tb[i])))
        min_throttle.append(min_thrust / max_thrust)
        max_throttle.append(max_thrust / max_thrust)

    m0.reverse()
    m1.reverse()
    f_t.reverse()
    ve.reverse()
    a_t.reverse()
    tu.reverse()
    l1.reverse()
    tb.reverse()
    min_throttle.reverse()
    max_throttle.reverse()
    vehicle = list()
    for i in range(len(stage)):
        stages = Struct()
        stages.m0 = m0[i]
        stages.fT = f_t[i]
        stages.ve = ve[i]
        stages.l1 = l1[i]
        stages.maxThrottle = max_throttle[i]
        stages.minThrottle = min_throttle[i]
        stages.maxT = tb[i]
        vehicle.append(stages)

    return vehicle


def stage_controller(vehicle, delay=2, ullage=True, booster=False):
    if not booster:
        vessel.control.throttle = 0
        time.sleep(delay)
        vessel.control.activate_next_stage()
        if ullage:
            vessel.control.rcs = True
            vessel.control.forward = 1
            time.sleep(2 * delay)
            vessel.control.forward = 0
        vessel.control.activate_next_stage()
        vehicle.remove(vehicle[0])
        if vessel.control.throttle == 0:
            vessel.control.throttle = 1
    if booster:
        vessel.control.activate_next_stage()
    while globals.state_thrust() < vessel.available_thrust:
        time.sleep(0.01)
