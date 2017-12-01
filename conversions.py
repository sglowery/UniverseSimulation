class Conversions:

    meters_in_mpc = 3.086e22
    kilometers_in_mpc = meters_in_mpc / 1e3

    def meters_to_mpc(self, meters):
        return meters / Conversions.meters_in_mpc

    def mpc_to_meters(self, mpc):
        return mpc * Conversions.meters_in_mpc
