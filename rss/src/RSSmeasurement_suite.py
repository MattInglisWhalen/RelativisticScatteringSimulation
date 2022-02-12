
from vpython import *

from os.path import dirname, abspath
import sqlite3

# external libraries
import matplotlib.pyplot as plt
import numpy as np

# internal classes
from rss.src.RSSsql import Squirrel
from rss.src.RSSvectors import *
from rss.src.RSScomposite import Composite

"""
No need for these functions to be a class
"""


class MeasurementSuite:

    def __init__(self, data_dir_path = ""):

        self._data_location = data_dir_path
        if data_dir_path == "" :
            # keep stepping back from the current directory until we are in the directory /rss
            loc = dirname(abspath(__file__))
            while loc[-3:] != "rss" :
                loc = dirname(loc)
                if loc == dirname(loc):
                    print(f"""MeasurementSuite init: python script {__file__} is not in the RSS package's directory. 
                              You will need to provide a path to the data-storage directory in 
                              the "data_dir_path" variable""")

            self._data_location = loc + f"/data"
        self._plots_location = self._data_location + "/plots"


    def measure_truth_level_parton_momenta_at_energy(self, energy, show_new_data = True, updates_database=True):

        sql = Squirrel()
        # sql.delete_table_named("truths")
        sql.create_table_named("truths")

        scene.forward = vector(0, 0, 1)
        scene.center = vector(5, 5, 0)
        scene.width = 1366
        scene.height = 768

        # setup options
        num = 100

        # visualization options
        timestep = 0.01
        frames_per_second = 30
        thermalization_time = 5

        protons = []
        # generate 100 protons
        for i in range(num):
            new_proton = Composite(name="Proton",
                                   pos=Vec3(i % 10, (i - i % 10) / 10, 0))  # grid with 1/GeV (hbar*c/GeV) spacings
            protons.append(new_proton)

        # let them thermalize
        for iframe in range(frames_per_second * thermalization_time):
            for iproton in protons:
                iproton.update_comp(dt=timestep)

        # boost them to the desired energy
        for iproton in protons:
            print(f"\nBefore {iproton.net_pos}")
            iproton.translate_constituents_to_rest_frame()
            iproton.E = energy
            iproton.translate_constituents_to_lab_frame()
            print(iproton)
            print(f"After {iproton.net_pos} with {iproton.net_mom} {iproton.net_mass} ")

        # let them thermalize again in the boosted frame
        for iframe in range(frames_per_second * thermalization_time):
            for iproton in protons:
                iproton.update_comp(dt=timestep)

        proton_energies = []
        energies = []
        px= []
        py = []
        pz = []
        lc_fraction_minus = []
        lc_fraction_plus = []
        for iproton in protons:
            for ifund in iproton.constituents:
                proton_energies.append(iproton.E)
                energies.append(ifund.E)
                px.append(ifund.px)
                py.append(ifund.py)
                pz.append(ifund.pz)
                lc_fraction_minus.append(ifund.p.pm / iproton.p.pm)
                lc_fraction_plus.append(ifund.p.pp / iproton.p.pp)

        if updates_database:
            sql.add_rows_to_database(proton_energies, energies, px, py, pz,
                                     lc_fraction_minus, lc_fraction_plus)

        # show in GUI
        if show_new_data:
            num_bins = 20
            plt.hist(energies)
            plt.show()
            plt.hist(px)
            plt.show()
            plt.hist(py)
            plt.show()
            plt.hist(pz)
            plt.show()
            plt.hist(lc_fraction_minus)
            plt.show()
            plt.hist(lc_fraction_plus)
            plt.show()


    def update_all_plots(self):

        sql = Squirrel()
        print(sql.database_location)

        rows = np.array(sql.get_all_rows_from_truths())
        for irow in rows[0:5] :
            print(irow)

        num_bins = 30

        names_list = [ f"energy_distribution_at_E={rows[0][1]}.png",
                        f"px_distribution_at_E={rows[0][1]}.png",
                        f"py_distribution_at_E={rows[0][1]}.png",
                        f"pz_distribution_at_E={rows[0][1]}.png",
                        f"large_frac_distribution_at_E={rows[0][1]}.png",
                        f"small_frac_distribution_at_E={rows[0][1]}.png" ]

        for idx, name in enumerate(names_list) :
            plt.hist( rows[:,idx+2], num_bins )
            plt.savefig(f"{self._plots_location}/{name}")
            plt.close()

    def show_all_plots(self):

        sql = Squirrel()
