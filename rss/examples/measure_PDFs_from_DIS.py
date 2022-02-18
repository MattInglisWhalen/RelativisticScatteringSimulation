
# local classes
from rss.src.RSSmeasurement_suite import *


def measure_DIS_PDFs():

    suite = MeasurementSuite()

    num_trials = 1
    show_plots = True

    for i in range(num_trials) :
        if num_trials > 1 :
            show_plots = False
        suite.measure_parton_momenta_from_DIS_at_energy(100, show_new_data=show_plots)
        print(f"Done {i}")
    # suite.update_all_plots()

    print("Exiting...")


if __name__ == "__main__" :

    measure_DIS_PDFs()
    raise SystemExit
