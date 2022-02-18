
# local classes
from rss.src.RSSmeasurement_suite import *


def measure_truth_PDFs():

    suite = MeasurementSuite()
    suite.create_table_named("truths")

    num_trials = 1
    show_plots = True

    for i in range(num_trials) :
        if num_trials > 1 :
            show_plots = False
        suite.measure_truth_level_parton_momenta_at_energy(100, show_new_data=show_plots)
        print(f"Done {i}")
    suite.update_all_plots()

    print("Exiting...")


if __name__ == "__main__" :

    measure_truth_PDFs()
    raise SystemExit
