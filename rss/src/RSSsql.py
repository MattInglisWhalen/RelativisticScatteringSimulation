from os.path import dirname, abspath
import sqlite3


class Squirrel:
    """
    For SQL queries and manips of a local (serverless) database
    """
    def __init__(self, filename = "PDF_raw_data.db", loc = "" ):

        self._db_location = loc + f"/data/{filename}"
        if loc == "" :
            # keep stepping back from the current directory until we are in the directory /rss
            loc = dirname(abspath(__file__))

            while loc[-3:] != "rss" :
                loc = dirname(loc)
                if loc == dirname(loc):
                    print(f"""Squirrel init: python script {__file__} is not in the RSS package's directory. 
                              You will need to provide a path to the database's directory in the "loc" variable""")

            self._db_location = loc + f"/data/{filename}"

        self._connection = sqlite3.connect(self._db_location)
        self._cursor = self._connection.cursor()

    def __repr__(self):
        return f"{self._db_location}"

    def __del__(self):
        self._connection.commit()  # Kinda silly to write this after every execution

    @property
    def cur(self):
        return self._cursor
    @property
    def database_location(self):
        return self._db_location


    def add_row_to_database(self, proton_energy, thermalization_time,E, px, py, pz, pnb_frac, pn_frac):
        sql_statement_string = f"""
            INSERT INTO truths 
                (proton_energy, energy, px, py, pz, frac_large, frac_small)
            VALUES
                ({proton_energy},{thermalization_time},{E},{px},{py},{pz},{pnb_frac},{pn_frac});
        """
        self.cur.execute(sql_statement_string)

    def add_rows_to_database(self, proton_energy_list = [], thermalization_time_list = [],
                                E_list = [], px_list = [], py_list = [], pz_list = [],
                                pnb_frac_list = [], pn_frac_list = []):
        for idx, _ in enumerate( E_list ) :
            sql_statement_string = f"""
                INSERT INTO truths 
                    (proton_energy, therm_time, energy, px, py, pz, frac_large, frac_small)
                VALUES
                    ({proton_energy_list[idx]},{thermalization_time_list[idx]},
                        {E_list[idx]},{px_list[idx]},{py_list[idx]},{pz_list[idx]},
                        {pnb_frac_list[idx]},{pn_frac_list[idx]});
            """
            self.cur.execute(sql_statement_string)
            self._connection.commit()

    def get_row_from_truths_with_id(self, ID):
        sql_statement_string = f"""
            SELECT
                *
            FROM
                truths
            WHERE
                id = {ID} ;
        """
        self.cur.execute(sql_statement_string)
        row = self.cur.fetchall()
        return row

    def get_all_rows_from_truths(self):
        sql_statement_string = f"""
            SELECT
                *
            FROM
                truths;
        """
        self.cur.execute(sql_statement_string)
        rows = self.cur.fetchall()
        return rows


    def make_string_primary_key_for_table(self, key="id", table = "truths"):
        sql_statement_string = f"""
            ALTER TABLE {table}
            ADD PRIMARY KEY({key});
        """
        self.cur.execute(sql_statement_string)

    def create_table_named(self, table = "truths"):
        sql_statement_string = f"""
                CREATE TABLE IF NOT EXISTS {table}(id INTEGER PRIMARY KEY, 
                                                    proton_energy INTEGER, 
                                                    therm_time FLOAT,                                                    
                                                    energy FLOAT, 
                                                    px FLOAT, 
                                                    py FLOAT, 
                                                    pz FLOAT, 
                                                    frac_large FLOAT, 
                                                    frac_small FLOAT);
        """
        self.cur.execute(sql_statement_string)

    def delete_table_named(self, table = "truths"):
        sql_statement_string = f"""
                DROP TABLE {table} ;
        """
        self.cur.execute(sql_statement_string)

