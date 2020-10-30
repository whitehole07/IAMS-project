__all__ = ["OrbitsHandler"]

from transfer import TransferAbstract
import pandas as pd
import os


class IncompatibleObjectError(Exception):
    pass


class OrbitsHandler(object):
    def __init__(self, **kwargs: list or TransferAbstract) -> None:
        """
        :argument
            [orbit_i, orbit_f]     -   1x2 List of Orbit objects ([OrbitPosition, OrbitPosition])

            or

            transfer               -   1x1 Transfer object (TransferAbstract)

        :raises
            IncompatibleObjectError, if the object passed is not a list of OrbitPosition instances or TransferAbstract object
        """

        self.transfers: dict = {}
        for key, value in kwargs.items():
            if type(value) is list and len(value) == 2:
                self.transfers = {**self.transfers, **{key: {"transfer": TransferAbstract(value[0], value[1])}}}
            elif isinstance(value, TransferAbstract):
                self.transfers = {**self.transfers, **{key: {"transfer": value}}}
            else:
                raise IncompatibleObjectError(f"Incompatible object passed: kwarg = '{key}': {type(value)}'")

    def get_comparison_matrix(self, transfer_name: str) -> pd.DataFrame:
        try:
            tr: TransferAbstract = self.transfers[transfer_name]["transfer"]
        except KeyError:
            raise KeyError(f"'{transfer_name}' not found in transfers' list")

        df: pd.DataFrame = pd.DataFrame(columns=["where_CP", "where_RP", "type_RP", "where_RE", "dv_tot", "dt_tot", "dv_list", "dt_list"])
        for where_CP in ["first", "second"]:
            for RP in [("first", "same"), ("second", "opposite"), ("first", "opposite"), ("second", "same")]:
                for where_RE in ["pericentre", "apocentre"]:
                    df.loc[len(df)] = tr.get_combination_array(where_CP, RP[0], RP[1], where_RE)

        self.transfers[transfer_name]["comparison_matrix"] = df
        return df

    def comparison_matrix_to_file(self, transfer_name: str, *, to_excel: bool = True, single_file: bool = True) -> None:
        try:
            tr: TransferAbstract = self.transfers[transfer_name]["transfer"]
            df: pd.DataFrame = self.transfers[transfer_name]["comparison_matrix"]
        except KeyError:
            raise KeyError("Comparison Matrix or TransferAbstract object not founded, you should run get_comparison_matrix function before")

        if not single_file and to_excel:
            if not os.path.exists("csvs"):
                os.mkdir("csvs")

            name: str = "%s_init(a=%.2f, e=%.2f, i=%.2f)_to_final(a=%.2f, e=%.2f, i=%.2f)" \
                        % (transfer_name,
                           tr.orbit_i.a, tr.orbit_i.e, tr.orbit_i.i,
                           tr.orbit_f.a, tr.orbit_f.e, tr.orbit_f.i)

            df.to_csv(f"csvs/{name}.csv")

        elif single_file and to_excel:
            if not os.path.exists("csvs"):
                os.mkdir("csvs")

            name: str = "%s_init(a=%.2f, e=%.2f, i=%.2f)_to_final(a=%.2f, e=%.2f, i=%.2f)" \
                        % (transfer_name,
                           tr.orbit_i.a, tr.orbit_i.e, tr.orbit_i.i,
                           tr.orbit_f.a, tr.orbit_f.e, tr.orbit_f.i)

            if not os.path.exists(f"csvs/{name}"):
                os.mkdir(f"csvs/{name}")

            df[:]["dv_tot"].to_excel(f"csvs/{name}/dv_tot.xlsx")
            df[:]["dt_tot"].to_excel(f"csvs/{name}/dt_tot.xlsx")

            for index, row in df.iterrows():
                dir_name: str = f"where_CP='{row['where_CP']}'_where_RP='{row['where_RP']}'_type_RP='{row['type_RP']}'_where_RE='{row['where_RE']}'"

                if not os.path.exists(f"csvs/{name}/{dir_name}"):
                    os.mkdir(f"csvs/{name}/{dir_name}")

                for data in ("dv", "dt"):
                    pd.DataFrame([d for x in row[f"{data}_list"] for d in x]).to_excel(f"csvs/{name}/{dir_name}/{data}.xlsx")
