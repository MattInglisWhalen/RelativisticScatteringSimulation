
# external libs
from vpython import *

# internal classes
from rss.src.RSSvectors import Vec3
from rss.src.RSSsystem import System


def lepton_lepton_collision():

    scene.forward = vector(1,0,0)
    # Visualization settings
    scene.width = 1366
    scene.height = 768
    # scene.range = (24, 24, 8)
    scene.autoscale = 0
    # scene.center = vector(0, 8, 0)
    scene.background = color.black
    beam_line = cylinder( pos = vector(0,0,-15), axis=vector(0,0,30), radius=0.005, color = color.yellow )

    # Simulation settings
    # thermalization_energy = 2.2
    # thermalization_time = 10
    CoM_energy = 10
    collision_time = 30

    scattering_system = System( mode = "lepton-lepton" , initial_energy = CoM_energy )
    scattering_system.simulate_for( seconds = collision_time )

    print("Done!")
    sleep(1)
    scene.delete()
    print("Exiting...")


def hide_vpython_residuals( my_scene ):
    for obj in my_scene.objects:
        if obj.pos == Vec3(0,0,0):
            obj.visible = False


if __name__ == "__main__" :

    lepton_lepton_collision()
    raise SystemExit
