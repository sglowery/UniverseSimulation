import numpy as np

class Universe:

    def __init__(self, matter_density, lambda_density, rad_density, hubble, scale=1):

        self.matter_density = matter_density    # unitless
        self.lambda_density = lambda_density    # unitless
        self.rad_density = rad_density          # unitless
        self.hubble = hubble                    # km s^1 Mpc^-1
        self.scale = scale

    def num_components(self):
        return int(self.matter_density > 0) + int(self.rad_density > 0) + int(self.lambda_density > 0)

    def deceleration_constant(self, scale):
        dm, dl, dr = self.density_at_scale(scale)
        return dr + .5 * dm - dl

    def density_at_scale(self, scale):
        comps = self.num_components()
        if(comps == 3):
            if(scale == 1):
                return self.matter_density, self.lambda_density, self.rad_density
            else:
                # Ω_m + Ω_l + Ω_r = 1
                # solve for two using ratios, plug in, solve for the rest
                # Ω_m/Ω_l = (Ω_m0/a^3)/Ω_l0 => Ω_m = [(Ω_m0/a^3)/Ω_l0]*Ω_l
                # Ω_r/Ω_l = (Ω_r0/a^4)/Ω_l0 => Ω_r = [(Ω_r0/a^4)/Ω_l0]*Ω_l
                # then, substitute into first equation to solve for lambda density:
                # [(Ω_m0/a^3)/Ω_l0]*Ω_l + [(Ω_r0/a^4)/Ω_l0]*Ω_l + Ω_l = 1
                # => Ω_l*[1 + (Ω_m0/a^3)/Ω_l0 + (Ω_r0/a^4)/Ω_l0] = 1
                # => Ω_l = 1/[1 + (Ω_m0/a^3)/Ω_l0 + (Ω_r0/a^4)/Ω_l0]
                # this particular equation is only

                # first: ratio of matter to lambda
                density_m = self.matter_density / self.lambda_density / scale**3 # * density_l
                density_r = self.rad_density / self.lambda_density / scale**4 # * density_l

                # plug in density_m and density_r into the Ω_m + Ω_l + Ω_r = 1 equation and solve for density_l
                density_l = 1/(1+density_m+density_r)

                # since we initially defined density_m and density_r in terms of density_l, we just need to multiply it
                # through to find their actual values
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

