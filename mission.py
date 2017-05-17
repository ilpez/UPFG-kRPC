import os


def run_mission(what_mission):
    mission = str(what_mission).capitalize()
    if mission == 'Launch':
        our_mission = list()
        our_mission.append(mission)
        print("Orbital Parameter")
        our_mission.append(input("Orbit apogee : "))
        our_mission.append(input("Orbit perigee : "))
        our_mission.append(input("Orbit inclination: "))
        our_mission.append(input("Orbit LAN : "))
        our_mission.append(input("Orbit true anomaly : "))
        our_mission.append(input("Landing mode(RTLS,ASDS,EXP) : ").upper())
        our_mission = ' '.join(our_mission)
        os.system("python -m " + our_mission)
    else:
        raise NotImplementedError


if __name__ == '__main__':
    run_mission(input("Choose Mission(Launch, Land, Rendezvous) : "))
