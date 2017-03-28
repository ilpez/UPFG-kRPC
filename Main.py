import Global
import numpy as np
import sys
import time
import upfg
# from multiprocessing.pool import ThreadPool

conn = Global.conn
space_center = Global.space_center
vessel = Global.vessel

position = Global.orbital_position
velocity = Global.orbital_velocity
surface_speed = Global.surface_speed

g0 = Global.g0
mu = Global.mu

target_apoapsis = float(sys.argv[1])
target_periapsis = float(sys.argv[2])
target_inclination = float(sys.argv[3])
target_lan = float(sys.argv[4])
if sys.argv[5] == 'RTLS':
    meco_speed = 1700
    turn_speed = 50
elif sys.argv[5] == 'ASDS':
    meco_speed = 2300
    turn_speed = 30
elif sys.argv[5] == 'EXP':
    meco_speed = 2700
    turn_speed = 30
else:
    meco_speed = 3000
    turn_speed = 30
print(meco_speed)
[azimuth, launch_time, target] = upfg.launchTargeting(target_periapsis,
                                                      target_apoapsis,
                                                      target_inclination,
                                                      target_lan,
                                                      -1.0)

game_launch_time = space_center.ut + launch_time
space_center.warp_to(game_launch_time - 10)

while (space_center.ut - game_launch_time) < 0:
    print('Time to launch %f' % (space_center.ut - game_launch_time))
    time.sleep(1)

vessel.control.throttle = 1
for engine in vessel.parts.engines:
    if not engine.active:
        print('There is no active engine, checking Propellant condition')
        for engine in engine.part.modules:
            if engine.has_field('Propellant'):
                if engine.get_field('Propellant') == 'Very Stable':
                    print('Engine is ready')
                    vessel.control.activate_next_stage()
while vessel.thrust < vessel.available_thrust:
    time.sleep(0.01)

vessel.auto_pilot.engage()
vessel.auto_pilot.target_heading = 0
vessel.auto_pilot.target_pitch = 90
vessel.control.activate_next_stage()

print('Proceeding Launch..')

while surface_speed() < turn_speed:
    time.sleep(0.1)

print('Clear from launch tower..')
print('Begin Pitch and Roll Program..')

vessel.auto_pilot.target_heading = azimuth
vessel.auto_pilot.target_roll = azimuth
vessel.auto_pilot.attenuation_angle = (1, 1, 0.2)
while True:
    if vessel.auto_pilot.target_roll > 0:
        vessel.auto_pilot.target_roll -= 0.1
    else:
        vessel.auto_pilot.target_roll = 0
    pitch = upfg.atand((900 - turn_speed) / (surface_speed() - turn_speed))
    vessel.auto_pilot.target_pitch = pitch
    if surface_speed() > meco_speed or vessel.available_thrust == 0:
        vessel.control.throttle = 0
        time.sleep(2)
        break
    time.sleep(0.01)

print('Main Engine Cutoff')

vessel.control.activate_next_stage()
vessel.control.rcs = True
vessel.control.forward = 1
time.sleep(2)

for engine in vessel.parts.engines:
    if not engine.active:
        print('There is no active engine, checking Propellant condition')
        for engine in engine.part.modules:
            if engine.has_field('Propellant'):
                if engine.get_field('Propellant') == 'Very Stable':
                    print('Engine is ready')
                    vessel.control.forward = 0
                    vessel.control.throttle = 1
                    vessel.control.activate_next_stage()

while vessel.thrust < vessel.available_thrust:
    time.sleep(0.01)

vessel.auto_pilot.target_roll = 0
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
upfg_internal.time = time.time()
upfg_internal.tgo = 0
upfg_internal.v = velocity()
upfg_internal.vgo = vdinit

upfg_guided = upfg.struct()

[upfg_internal, upfg_guided] = upfg.upfg(
    time.time(), position(), velocity(), target, upfg_internal)

while True:
    [upfg_internal, upfg_guided] = upfg.upfg(
        time.time(), position(), velocity(), target, upfg_internal)
    t = upfg_internal.tgo
    [upfg_internal, upfg_guided] = upfg.upfg(
        time.time(), position(), velocity(), target, upfg_internal)
    t1 = upfg_internal.tgo
    if abs(t1 - t) / t1 < 0.01 and t1 > 1:
        vessel.auto_pilot.target_heading = upfg_guided.yaw
        vessel.auto_pilot.target_pitch = upfg_guided.pitch
    if upfg_guided.tgo < 0.1:
        vessel.control.throttle = 0
        break
    if Global.surface_altitude() > 100000 and not fairing_jettison:
        for part in vessel.parts.all:
            for module in part.modules:
                if module.has_event('Jettison'):
                    module.trigger_event('Jettison')
                    fairing_jettison = True

print('Mission Success')
