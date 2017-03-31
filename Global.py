import krpc

conn = krpc.connect(name='Launch to orbit')
space_center = conn.space_center
vessel = space_center.active_vessel
g0 = vessel.orbit.body.surface_gravity
mu = vessel.orbit.body.gravitational_parameter

inertial_reference_frame = vessel.orbit.body.non_rotating_reference_frame
body_reference_frame = vessel.orbit.body.reference_frame

orbital_velocity = conn.add_stream(vessel.velocity,
                                   inertial_reference_frame)
orbital_position = conn.add_stream(vessel.position,
                                   inertial_reference_frame)
orbital_speed = conn.add_stream(getattr,
                                vessel.flight(inertial_reference_frame),
                                'speed')

surface_velocity = conn.add_stream(vessel.velocity,
                                   body_reference_frame)
surface_speed = conn.add_stream(getattr,
                                vessel.flight(body_reference_frame),
                                'speed')
surface_altitude = conn.add_stream(getattr,
                                   vessel.flight(),
                                   'mean_altitude')
state_q = conn.add_stream(getattr, vessel.flight(), 'dynamic_pressure')
state_mass = conn.add_stream(getattr, vessel, 'mass')
state_thrust = conn.add_stream(getattr, vessel, 'thrust')
universal_time = conn.add_stream(getattr, space_center, 'ut')
