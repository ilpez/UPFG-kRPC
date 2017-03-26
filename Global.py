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
orbital_hspeed = conn.add_stream(getattr,
                                 vessel.flight(inertial_reference_frame),
                                 'horizontal_speed')
orbital_vspeed = conn.add_stream(getattr,
                                 vessel.flight(inertial_reference_frame),
                                 'vertical_speed')
orbital_altitude = conn.add_stream(getattr,
                                   vessel.orbit,
                                   'radius')

surface_velocity = conn.add_stream(vessel.velocity,
                                   body_reference_frame)
surface_position = conn.add_stream(vessel.position,
                                   body_reference_frame)
surface_speed = conn.add_stream(getattr,
                                vessel.flight(body_reference_frame),
                                'speed')
surface_hspeed = conn.add_stream(getattr,
                                 vessel.flight(body_reference_frame),
                                 'horizontal_speed')
surface_vspeed = conn.add_stream(getattr,
                                 vessel.flight(body_reference_frame),
                                 'vertical_speed')
surface_altitude = conn.add_stream(getattr,
                                   vessel.flight(),
                                   'mean_altitude')
