
# external libs
from vpython import *

# internal classes
from rss.src.RSSvectors import *
from rss.src.RSSsystem import System
from rss.src.RSSfundamental import *


def proton_proton_collision():

    scene.forward = vector(1,0,0)
    # Visualization settings
    scene.width = 1366
    scene.height = 768
    # scene.range = vector(30, 30, 30)
    scene.autoscale = 0
    # scene.center = Vec(0, 8, 0)
    scene.background = color.black
    beam_line = cylinder( pos = vector(0,0,-15), axis=vector(0,0,30), radius=0.005, color = color.yellow )

    # Simulation settings
    thermalization_energy = 2.01
    thermalization_time = 10
    CoM_energy = 100
    collision_time = 30

    scattering_system = System( mode = "composite-composite" , initial_energy = thermalization_energy )
    hide_vpython_residuals(scene)
    scattering_system.simulate_for( seconds = thermalization_time )

    scattering_system.boost_to_collision_energy(CoM_energy)
    scattering_system.print()
    scattering_system.simulate_for( seconds = collision_time )
    print("\n\nAt end of simulation: ")
    scattering_system.print()

    print("Done!")
    sleep(1)
    # scene.delete()
    print("Exiting...")


def hide_vpython_residuals( my_scene ):
    for obj in my_scene.objects:
        if obj.pos == vector(0,0,0):
            obj.visible = False


if __name__ == "__main__" :

    proton_proton_collision()
    raise SystemExit


