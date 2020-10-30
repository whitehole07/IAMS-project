__all__ = ["OrbitsHandler"]

from transfer import TransferAbstract
import pandas as pd


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
                self.transfers = {**self.transfers, **{key: TransferAbstract(value[0], value[1])}}
            elif isinstance(value, TransferAbstract):
                self.transfers = {**self.transfers, **{key: value}}
            else:
                raise IncompatibleObjectError(f"Incompatible object passed: kwarg = '{key}': {type(value)}'")

    def get_comparison_matrix(self, transfer_name: str, *, to_csv: bool = True) -> pd.DataFrame:
        try:
            tr: TransferAbstract = self.transfers[transfer_name]
        except KeyError:
            raise KeyError(f"'{transfer_name}' not found in transfers' list")

        df: pd.DataFrame = pd.DataFrame(columns=["where_CP", "where_RP", "type_RP", "where_RE", "dv_tot", "dt_tot", "dv_list", "dt_list"])
        for where_CP in ["first", "second"]:
            for RP in [("first", "same"), ("second", "opposite"), ("first", "opposite"), ("second", "same")]:
                for where_RE in ["pericentre", "apocentre"]:
                    df.loc[len(df)] = tr.get_combination_array(where_CP, RP[0], RP[1], where_RE)

        if to_csv:
            df.to_csv("csvs/%s_init(a=%.2f, e=%.2f, i=%.2f)_to_final(a=%.2f, e=%.2f, i=%.2f).csv"
                      % (transfer_name,
                         tr.orbit_i.a, tr.orbit_i.e, tr.orbit_i.i,
                         tr.orbit_f.a, tr.orbit_f.e, tr.orbit_f.i))

        return df
