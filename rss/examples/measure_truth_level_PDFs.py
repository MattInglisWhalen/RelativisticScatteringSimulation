
# local classes
from rss.src.RSSmeasurement_suite import *


def measure_truth_PDFs():

    suite = MeasurementSuite()
    suite.measure_truth_level_parton_momenta_at_energy(100, show_new_data=True)
    suite.update_all_plots()

    print("Exiting...")


if __name__ == "__main__" :

    measure_truth_PDFs()
    raise SystemExit
