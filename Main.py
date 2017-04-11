import Global
import numpy as np
import sys
import time
import upfg
# from multiprocessing.pool import ThreadPool

conn = Global.conn
space_center = Global.space_center
vessel = Global.vessel
vehicle = upfg.analyze_vehicle()

position = Global.orbital_position
velocity = Global.orbital_velocity
surface_speed = Global.surface_speed

g0 = Global.g0
mu = Global.mu

g_lim = 4
q_lim = 27000

target_apoapsis = float(sys.argv[1])
target_periapsis = float(sys.argv[2])
target_inclination = float(sys.argv[3])
target_lan = float(sys.argv[4])
target_true_anomaly = float(sys.argv[5])
if sys.argv[6] == 'RTLS':
    meco_speed = 1700
    turn_speed = 30
elif sys.argv[6] == 'ASDS':
    meco_speed = 2300
    turn_speed = 30
elif sys.argv[6] == 'EXP':
    meco_speed = 2700
    turn_speed = 30
else:
    meco_speed = vehicle[0].l1 - 1000
    turn_speed = 30
print(meco_speed)
[azimuth, launch_time, target] = upfg.launchTargeting(target_periapsis,
                                                      target_apoapsis,
                                                      target_inclination,
                                                      target_lan,
                                                      target_true_anomaly)

game_launch_time = space_center.ut + launch_time
space_center.warp_to(game_launch_time - 10)

while (space_center.ut - game_launch_time) < 0:
    print('Time to launch %f' % (space_center.ut - game_launch_time))
    time.sleep(1)

vessel.control.throttle = 1
vessel.control.activate_next_stage()
while Global.state_thrust() < vessel.available_thrust:
    time.sleep(0.2)

vessel.auto_pilot.engage()
vessel.auto_pilot.target_heading = 90
vessel.auto_pilot.target_roll = 0
vessel.auto_pilot.target_pitch = 90
vessel.control.activate_next_stage()
vessel.auto_pilot.wait()

print('Proceeding Launch..')

while surface_speed() < turn_speed:
    pass

print('Clear from launch tower..')
print('Begin Pitch Program..')

vessel.auto_pilot.target_heading = azimuth

while True:
    vessel.control.throttle = upfg.throttle_control(vehicle, g_lim, q_lim)

    pitch1 = upfg.atand((900 - 2 * turn_speed) /
                        (surface_speed() - turn_speed))
    pitch2 = upfg.angle_from_vec(Global.surface_velocity(),
                                 Global.body_reference_frame,
                                 'pitch')
    vessel.auto_pilot.target_pitch = min(pitch1, pitch2)
    if surface_speed() > meco_speed or vessel.available_thrust == 0:
        break
    time.sleep(0.01)

print('Main Engine Cutoff')

upfg.stageController(vehicle)

fairing_jettison = False
cser = upfg.struct()
cser.dtcp = 0
cser.xcp = 0
cser.a = 0
cser.d = 0
cser.e = 0

rdinit = upfg.rodrigues(upfg.unit(position()), -target.normal, 20)
rdinit = np.multiply(rdinit, target.radius)
vdinit = np.multiply(target.velocity, upfg.unit(
    upfg.cross(-target.normal, rdinit)))
vdinit = vdinit - velocity()
upfg_internal = upfg.struct()
upfg_internal.cser = cser
upfg_internal.rbias = [0, 0, 0]
upfg_internal.rd = rdinit
upfg_internal.rgrav = np.multiply(
    np.multiply(-(mu / 2), position()), 1 / upfg.norm(position())**3)
upfg_internal.tb = 0
upfg_internal.time = Global.universal_time()
upfg_internal.tgo = 0
upfg_internal.v = velocity()
upfg_internal.vgo = vdinit
converged = False
upfg_guided = upfg.struct()
iteration = 0

while converged is False:
    [upfg_internal, upfg_guided] = upfg.upfg(vehicle, target, upfg_internal)
    t1 = upfg_internal.tgo
    [upfg_internal, upfg_guided] = upfg.upfg(vehicle, target, upfg_internal)
    t2 = upfg_internal.tgo
    if abs(t1 - t2) / t2 < 0.01:
        converged = True
    iteration += 1

print('Guidance converged after %f iteration' % iteration)

while True:
    vessel.control.throttle = upfg.throttle_control(vehicle, g_lim, q_lim)
    [upfg_internal, upfg_guided] = upfg.upfg(vehicle, target, upfg_internal)
    t = upfg_guided.tgo
    if t > 1:
        vessel.auto_pilot.target_heading = upfg_guided.yaw
        vessel.auto_pilot.target_pitch = upfg_guided.pitch
    if t < 0.1:
        vessel.control.throttle = 0
        break
    if Global.orbital_speed() > target.velocity:
        vessel.control.throttle = 0
        break
    if Global.surface_altitude() > 110000 and not fairing_jettison:
        for fairing in vessel.parts.fairings:
            for module in fairing.part.modules:
                if module.has_event('Jettison'):
                    module.trigger_event('Jettison')
                    fairing_jettison = True
print('Mission Success')
