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
        return int(self.matter_density != 0) + int(self.rad_density != 0) + int(self.lambda_density != 0)

    def deceleration_constant(self, scale):
        dm, dl, dr = self.density_at_scale(scale)
        return dr + .5 * dm - dl

    def cosmic_time(self, target_scale, steps=1e5, start_point=1e-12, scales=-1):
        # This method takes an input scale factor to correspond to a time after the Big Bang, in gigayears
        #
        # ---------- EQUATION 1: H_0*t = integral from 0 to a: da * [Ω_r0/a^2 + Ω_m0/a + Ω_l0/a^2 + (1-Ω0)]^(-1/2) -----
        # This equation was derived from the general Friedmann equation (Ryden, p. 85)
        time = 0
        old = 0

        # allows for some default coordinates to be defined if the user wants to override them
        if scales == -1:
            scales = np.logspace(np.log10(start_point), np.log10(target_scale), num=steps)
        for scale_factor in scales:
            diff = scale_factor - old
            time += diff / np.sqrt(self.rad_density / scale_factor**2 + self.matter_density / scale_factor +
                                            self.lambda_density * scale_factor**2 +
                                            (1 - (self.lambda_density+self.matter_density+self.rad_density)))
            old = scale_factor
        time *= (conversions.kilometers_in_mpc / self.hubble)
        time /= (3600 * 24 * 365 * 1e9)

        return time

    def scale_at_time(self, time):
        scale_rm = self.rad_density/self.matter_density
        scale_ml = np.float_power(self.matter_density/self.lambda_density, 1.0/3.0)
        time_rm = (4.0/3.0)*(1-1/np.sqrt(2)) * np.power(scale_rm, 2)/np.sqrt(self.rad_density) \
                  * conversions.kilometers_in_mpc / self.hubble
        time_ml = 2 / self.hubble * conversions.kilometers_in_mpc / (3*np.sqrt(1-self.matter_density)) * np.log(1 + np.sqrt(2))
       # print("time rm:", time_rm/(3600*24*365*1e6), "\ntime ml:",time_ml/(3600*24*365*1e9))
        if time > 0 and time < .75*time_rm:
            return np.sqrt((2*np.sqrt(self.rad_density)*time*self.hubble/conversions.kilometers_in_mpc))
        elif time >= .75*time_rm and time < time_rm + .25*time_ml:
            return np.float_power(time/(2/(3*self.hubble/conversions.kilometers_in_mpc)), 2.0/3.0)
        elif time >= time_rm + .25*time_ml and time < 1.25*time_ml:
            return np.float_power(3.0 / 2.0 * np.sqrt(self.matter_density) * self.hubble/conversions.kilometers_in_mpc
                                  * time, 1.0/2.0)
        elif time >= 1.25*time_ml:
            return scale_ml*np.exp(np.sqrt(1-self.matter_density)*time*self.hubble/conversions.kilometers_in_mpc/2.975)
        else:
            print("Time must be greater than 0!")

    def density_at_scale(self, scale):
        comps = self.num_components()
        if comps == 3:
            if scale == 1:
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
                density_m = self.matter_density / self.lambda_density / scale**3  # * density_l
                density_r = self.rad_density / self.lambda_density / scale**4  # * density_l
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

        elif comps == 2:
            if self.matter_density == 0:
                density_r = (self.rad_density / self.lambda_density) / scale**4
                density_l = 1/(1+density_r)
                density_r *= density_l
                return density_l, density_r
            elif self.lambda_density == 0:
                density_r = (self.rad_density / self.matter_density) / scale
                density_m = 1/(1+density_r)
                density_r *= density_m
                return density_m, density_r
            elif self.rad_density == 0:
                density_m = (self.matter_density / self.lambda_density) / scale**3
                density_l = 1/(1+density_m)
                density_m *= density_l
                return density_m, density_l
