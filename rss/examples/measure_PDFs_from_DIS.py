
# external libs
import numpy as np

# internal classes
from rss.src.RSSmeasurement_suite import *


def measure_DIS_PDFs():

    suite = MeasurementSuite()
    # suite.sql.delete_table_named("dis")
    suite.create_table_named("dis")

    rows = np.array(suite.sql.get_all_rows_from_dis_table())
    NTOT = len(rows)
    print(f"{NTOT} rows so far in DIS table")
    for idx in reversed(range(NTOT)) :
        print(rows[idx])
        if idx < NTOT - 5 :
            break

    num_trials = 1000

    for i in range(num_trials) :
        if num_trials > 1 :
            show_plots = False
        suite.measure_parton_momenta_from_DIS_at_energy(100)
        print(f"Done {i}")
    # suite.update_all_plots()

    print("Exiting...")


if __name__ == "__main__" :

    measure_DIS_PDFs()
    raise SystemExit
