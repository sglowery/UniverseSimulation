class Conversions:

    metersInMpc = 3.086e22

    def meters_to_mpc(self, meters):
        return meters / Conversions.metersInMpc

    def mpc_to_meters(self, mpc):
        return mpc * Conversions.metersInMpc
