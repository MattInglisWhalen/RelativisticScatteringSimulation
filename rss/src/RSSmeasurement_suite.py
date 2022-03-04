
# default libs
from math import sqrt
from os.path import dirname, abspath

# external libraries
import matplotlib.pyplot as plt
import numpy as np
from vpython import scene, vector, color, cylinder
from vpython import *

# internal classes
from rss.src.RSSsql import Squirrel
from rss.src.RSSvectors import *
from rss.src.RSScomposite import Composite
from rss.src.RSSsystem import System
import rss.src.RSSconstants as const

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

        self._sql = Squirrel()

    @property
    def sql(self):
        return self._sql

    def create_table_named(self, table):
        if table == "truths" :
            self.sql.create_truths_table()
        elif table == "dis" :
            self.sql.create_dis_table()
        else :
            raise NotImplementedError

    def measure_parton_momenta_from_DIS_at_energy(self, energy, updates_database = True) :

        scene.forward = vector(1, 0, 0)
        # Visualization settings
        scene.width = 1366
        scene.height = 768
        # scene.range = vector(30, 30, 30)
        scene.autoscale = 0
        # scene.center = Vec(0, 8, 0)
        scene.background = color.black
        beam_line = cylinder(pos=vector(0, 0, -15), axis=vector(0, 0, 30), radius=0.005, color=color.yellow)

        # Simulation settings
        thermalization_energy = 2.01
        timestep_per_frame = 0.01

        # visualization options
        fps = 30
        thermalization_frames = fps*10
        thermalization_time = thermalization_frames*timestep_per_frame

        CoM_energy = energy
        collision_time = 3

        scattering_system = System(mode="lepton-proton", initial_energy=thermalization_energy)
        hide_vpython_residuals(scene)
        scattering_system.simulate_for(seconds=thermalization_frames/fps, frames_per_second=fps)

        scattering_system.boost_to_collision_energy(CoM_energy)
        lepton_p_before = scattering_system.fundamentals[0].p
        proton_p_before = scattering_system.composites[0].p
        scattering_system.print()
        scattering_system.simulate_for(seconds=collision_time)
        print("\n\nAt end of simulation: ")
        scattering_system.print()

        lepton_p_after = scattering_system.fundamentals[0].p
        # remnant_p = scattering_system.net_mom() - lepton_p_after
        mom_transfer = lepton_p_before-lepton_p_after

        """
        Measured variables
        """

        measured_nu = (mom_transfer*proton_p_before)/const.PROTON_MASS
        if measured_nu == 0 :
            return
        measured_Qsqr = -mom_transfer*mom_transfer
        measured_x = measured_Qsqr/(2*const.PROTON_MASS*measured_nu)
        measured_y = (mom_transfer*proton_p_before)/(lepton_p_before*proton_p_before)
        measured_Wsqr = const.PROTON_MASS**2 + 2*const.PROTON_MASS*measured_nu - measured_Qsqr

        measured_final_proton_mass = scattering_system.composites[0].M
        measured_final_jet_mass = (scattering_system.net_mom() - scattering_system.composites[0].p
                                   - scattering_system.fundamentals[0].p).M

        print(f"Q^2: {measured_Qsqr}\nnu: {measured_nu}\nx: {measured_x}\ny: {measured_y}\n"
              + f"W^2: {measured_Wsqr}\nM_P^2: {measured_final_proton_mass}\nM_J^2: {measured_final_jet_mass}")

        if updates_database:
            self.sql.add_row_to_dis_table(CoM_energy,thermalization_time,measured_Qsqr,measured_nu,
                                          measured_x,measured_y,measured_final_proton_mass,measured_final_jet_mass)





    def measure_truth_level_parton_momenta_at_energy(self, energy, show_new_data = True, updates_database=True):

        scene.forward = vector(0, 0, 1)
        scene.center = vector(5, 5, 0)
        scene.width = 1366
        scene.height = 768

        # setup options
        num = 100
        timestep_per_frame = 0.01

        # visualization options
        frames_per_second = 30
        thermalization_frames = frames_per_second*2
        thermalization_time = thermalization_frames*timestep_per_frame

        protons = []
        # generate 100 protons
        for i in range(num):
            new_proton = Composite(name="Proton",
                                   pos=Vec3(i % 10, (i - i % 10) / 10, 0))  # grid with 1/GeV (hbar*c/GeV) spacings
            protons.append(new_proton)

        # scene.delete()

        # let them thermalize
        for iframe in range(thermalization_frames):
            for iproton in protons:
                iproton.update_comp(dt=timestep_per_frame)

        # boost them to the desired energy
        for iproton in protons:
            print(f"\nBefore {iproton.net_pos}")
            iproton.translate_constituents_to_rest_frame()
            iproton.E = energy
            iproton.translate_constituents_to_lab_frame()
            print(iproton)
            print(f"After {iproton.net_pos} with {iproton.net_mom} {iproton.net_mass} ")

        # let them thermalize again in the boosted frame
        for iframe in range(thermalization_frames):
            for iproton in protons:
                iproton.update_comp(dt=timestep_per_frame)

        proton_energies = []
        thermalization_times = []
        energies = []
        px = []
        py = []
        pz = []
        lc_fraction_minus = []
        lc_fraction_plus = []
        for iproton in protons:
            for ifund in iproton.constituents:
                proton_energies.append(iproton.E)
                thermalization_times.append(thermalization_time)
                energies.append(ifund.E)
                px.append(ifund.px)
                py.append(ifund.py)
                pz.append(ifund.pz)
                lc_fraction_minus.append( ifund.p.pm/iproton.p.pm )
                lc_fraction_plus.append(  ifund.p.pp/iproton.p.pp )

        if updates_database:
            self.sql.add_rows_to_truths_table(proton_energies, thermalization_times,
                                              energies, px, py, pz,
                                              lc_fraction_minus, lc_fraction_plus)

        # show in GUI
        if show_new_data:
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

        rows = np.array(self.sql.get_all_rows_from_truths_table())
        NTOT = len(rows)
        print(f"{NTOT} rows collected")
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
            plt.hist( rows[:,idx+3], num_bins )
            plt.savefig(f"{self._plots_location}/{name}")
            plt.close()

        """
        More detailed plots
        """

        #
        # The lightcone PDF
        #

        zetas = rows[:,7]  # p- / P-
        MIN, MAX, NBINS = 0.001, 1.0, 50

        counts, boundaries = np.histogram(zetas, bins = np.geomspace(MIN,MAX,NBINS+1) )

        bin_widths = [ boundaries[i+1]-boundaries[i] for i in range(NBINS) ]
        normed_masses = [ counts[i]/NTOT for i in range(NBINS) ]  # sum over all normed (probability) masses gives 1

        densities = [ normed_masses[i]/bin_widths[i] for i in range(NBINS) ]  # integrate over bin to get bin mass
        sigma_densities = [ sqrt( densities[i]/(NTOT*bin_widths[i]) ) for i in range(NBINS) ]
        # uncertainty of count N_i in unnormalized histogram is sqrt(N_i), so when converted to a density p_i
        # its uncertainty will be sqrt(Ni)/(Ntot*bin_width_i)

        xvals = [ sqrt(boundaries[i]*boundaries[i+1]) for i in range(NBINS)]
        # xpx = [ densities[i] for i in range(NBINS)]
        xpx = [ xvals[i]*densities[i] for i in range(NBINS) ]  # the usual way to display PDFs
        sigma_xpx = [ xvals[i]*sigma_densities[i] for i in range(NBINS)]  # usual propagation of uncertainty

        plt.scatter(xvals,xpx)
        plt.xlim([0.001,1])
        plt.xscale("log")
        plt.errorbar(xvals,xpx,yerr=sigma_xpx, fmt='o')
        plt.xlabel("x = p.n / P.n")
        plt.ylabel("x f_i/N(x)")
        plt.title("Parton Lightcone Momentum Fraction Probabilities")
        plt.savefig(f"{self._plots_location}/PartonLightconeFractionDistribution.png")
        plt.show()
        plt.close()

        #
        # The radial PDF
        #

        radii = [ sqrt( rows[i,4]**2 + rows[i,5]**2 ) for i in range(NTOT) ]
        MIN, MAX, NBINS = 0.001, 1.0, 50

        counts, boundaries = np.histogram(radii, bins=np.geomspace(MIN, MAX, NBINS + 1))

        bin_widths = [boundaries[i + 1] - boundaries[i] for i in range(NBINS)]
        normed_masses = [counts[i] / NTOT for i in range(NBINS)]  # sum over all normed (probability) masses gives 1

        densities = [normed_masses[i] / bin_widths[i] for i in range(NBINS)]  # integrate over bin to get bin mass
        sigma_densities = [sqrt(densities[i] / (NTOT * bin_widths[i])) for i in range(NBINS)]
        # uncertainty of count N_i in unnormalized histogram is sqrt(N_i), so when converted to a density p_i
        # its uncertainty will be sqrt(Ni)/(Ntot*bin_width_i)

        xvals = [sqrt(boundaries[i] * boundaries[i + 1]) for i in range(NBINS)]
        xpx = [densities[i] for i in range(NBINS)]
        # xpx = [ xvals[i]*densities[i] for i in range(NBINS) ]  # the usual way to display PDFs
        sigma_px = [ sigma_densities[i] for i in range(NBINS)]  # usual propagation of uncertainty

        plt.scatter(xvals, xpx)
        plt.xlim([0.001, 1])
        plt.xscale("log")
        plt.errorbar(xvals, xpx, yerr=sigma_xpx, fmt='o')
        plt.xlabel("p_T")
        plt.ylabel("f_i/N(p_T)")
        plt.title("Parton Transverse Momentum Probabilities")
        plt.savefig(f"{self._plots_location}/RadialMomentumDistribution.png")
        plt.show()
        plt.close()


def hide_vpython_residuals( my_scene ):
    for obj in my_scene.objects:
        if obj.pos == vector(0,0,0):
            obj.visible = False
