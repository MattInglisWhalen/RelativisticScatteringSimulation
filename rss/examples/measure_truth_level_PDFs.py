
# external libs
import numpy as np

# internal classes
from rss.src.RSSmeasurement_suite import *


def measure_truth_PDFs():

    suite = MeasurementSuite()
    suite.create_table_named("truths")

    rows = np.array(suite.sql.get_all_rows_from_truths_table())
    NTOT = len(rows)
    print(f"{NTOT} rows so far in DIS table")
    for idx in reversed(range(NTOT)) :
        print(rows[idx])
        if idx < NTOT - 5 :
            break

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
