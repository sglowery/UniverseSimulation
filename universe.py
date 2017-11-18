import numpy as np
from conversions import conversions

class Universe:

    def __init__(self, matter_density, lambda_density, rad_density, hubble, scale=1, curvature=0):

        self.matter_density = matter_density    # unitless
        self.lambda_density = lambda_density    # unitless
        self.rad_density = rad_density          # unitless
        self.hubble = hubble                    # km s^1 Mpc^-1
        self.scale = scale
        self.curvature = curvature

    def num_components(self):
        return int(self.matter_density > 0) + int(self.rad_density > 0) + int(self.lambda_density > 0)

    def deceleration_constant(self, scale):
        dm, dl, dr = self.density_at_scale(scale)
        return dr + .5 * dm - dl

    def cosmic_time(self, scale, steps, startpoint=1e-10):
        # H_0*t = integral from 0 to a: da * [Ω_r0/a^2 + Ω_m0/a + Ω_l0/a^2 + (1-Ω0)]^(-1/2)
        time = 0
        old = 0
        for da in np.logspace(np.log10(startpoint), np.log10(scale), num=steps, endpoint=True):
            diff = da - old
            time += diff / np.sqrt(self.rad_density / da**2 + self.matter_density / da +
                                            self.lambda_density / da**2 +
                                            (1 - self.lambda_density+self.matter_density+self.rad_density))
            old = da
        #time /= (self.hubble/conversions.metersInMpc)
        #time /= (3600 * 24 * 365 * 1e9)

        return time

    def density_at_scale(self, scale):
        comps = self.num_components()
        if(comps == 3):
            if(scale == 1):
                return self.matter_density, self.lambda_density, self.rad_density
            else:
                # For a universe with matter, radiation and lambda components, the following equations are used to
                # determine their density parameters as a function of scale factor:

                # ---------- EQUATION 1: Ω_m + Ω_l + Ω_r = 1 ----------
                # with Ω_i = ε_i / ε_c
                # where ε_i is the energy density (J m^-3) of the ith component and ε_c is the critical density of
                # the universe (which determines curvature).
                # Since Ω_i is a fractional, dimensionless number, the sum of all the density parameters of each
                # component should equal 1
                #
                # ---------- EQUATION 2: Ω_i / Ω_j = (Ω_i0 / a^(3(1+w_i))) / (Ω_j0 / a^(3(1+w_j)))) ----------
                # where Ω_i0 = Ω_i at the current time and w_i refers to the pressure value of the ith component, with:
                # w_m = 0, w_r = 1/3 and w_l = -1
                # A simpler way of putting Equation 2: Ω_i / Ω_j = (Ω_i0 / Ω_j0) * a^(3(w_j - w_i))
                #
                # By using equation 1 and combining two different ratios in equation 2, we get three equations for three
                # unknowns. I choose to solve for everything in terms of the density parameter of dark energy, but it's
                # no less valid to choose any other density parameter. I merely have to solve for two using ratios,
                # plug in and solve for the rest:
                # Ω_m/Ω_l = (Ω_m0/a^3)/Ω_l0 => Ω_m = [(Ω_m0/a^3)/Ω_l0]*Ω_l
                # Ω_r/Ω_l = (Ω_r0/a^4)/Ω_l0 => Ω_r = [(Ω_r0/a^4)/Ω_l0]*Ω_l
                # Then, substitute into Equation 1 to solve for lambda density:
                # [(Ω_m0/a^3)/Ω_l0]*Ω_l + [(Ω_r0/a^4)/Ω_l0]*Ω_l + Ω_l = 1
                # => Ω_l*[1 + (Ω_m0/a^3)/Ω_l0 + (Ω_r0/a^4)/Ω_l0] = 1
                # => Ω_l = 1/[1 + (Ω_m0/a^3)/Ω_l0 + (Ω_r0/a^4)/Ω_l0]

                # first: ratio of matter to lambda
                density_m = self.matter_density / self.lambda_density / scale**3 # * density_l
                density_r = self.rad_density / self.lambda_density / scale**4 # * density_l
                # (note: the "* density_l" part at the end of the previous two lines is just meant to illustrate that
                # I only need to multiply these values by that number to obtain the correct value, which will be done
                # after density_l is determined

                # Plug in density_m and density_r into the Ω_m + Ω_l + Ω_r = 1 equation and solve for density_l
                density_l = 1/(1+density_m+density_r)

                # As mentioned earlier, just multiply the density_m and density_r values by density_l to get the
                # correct value!
                density_m *= density_l
                density_r *= density_l

                return density_m, density_l, density_r

        elif(comps == 2):
            # TODO: verify these calculations
            if(self.matter_density == 0):
                density_r = (self.rad_density / self.lambda_density) / self.scale**4
                density_l = 1/(1+density_r)
                density_r *= density_l
                return density_l, density_r
            elif(self.lambda_density == 0):
                density_r = (self.rad_density / self.matter_density) / self.scale
                density_m = 1/(1+density_r)
                density_r *= density_m
                return density_m, density_r
            elif(self.rad_density == 0):
                density_m = (self.matter_density / self.lambda_density) / self.scale**3
                density_l = 1/(1+density_m)
                density_m *= density_l
                return density_m, density_l

