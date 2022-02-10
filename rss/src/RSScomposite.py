
# external libs
# from vpython import *
import numpy as np


# internal classes
from rss.src.RSSvectors import *
from rss.src.RSSfundamental import Fundamental
from rss.src.RSSmediator import Mediator
import rss.src.RSSconstants as const


class Composite:
    """ Class variables """
    #
    # name [string]
    # constituents[]  [array of Fundamentals]
    # mediators[] [array of Mediators]
    # dict_constituents_mediators{} [dictionary with constituents as a key for the mediators value]
    # _r3 -- position [Vec3]
    # _p -- momentum [Mom4]
    # N -- num_constituents  [int]
    # M -- mass   [float]

    # billiards[] [secondary list of constituents which have been scattered]


    def __init__(self,
                 name="CoreCompositeModel",
                 pos=Vec3(0, 0, 0),
                 energy=1,
                 direction=Vec3(0, 0, 1),
                 mass=1):

        self.name = name

        # Names
        if name[0:6] == "Proton" :
            self._M = const.PROTON_MASS
        elif name[0:5] == "Meson" :
            self._M = const.MESON_MASS
        else:
            self._M = mass

        if energy < self._M:
            print("Tachyonic composite created because of energy " , energy , " and mass " , self._M )
            self.name = "Tachyonic Composite"

        self.constituents = []
        self.mediators = []
        self.dict_constituents_mediators = dict()
        """A dictionary is important here because a particle can be removed from the
            constituents list, but their mediator holds their potential energy
            which needs to remain after their associated particle is removed """
        self.billiards = []

        self._r3 = pos
        # self._p = Mom4()
        # self.E = energy
        # self.dir = direction
        p3_mag = sqrt( energy**2 - self._M**2 )
        self._p = Mom4( energy, p3_mag*direction.hat.x, p3_mag*direction.hat.y, p3_mag*direction.hat.z)
        self.N = 0
        self.offshell_counter = 0

        self.create_constituents()
        self.translate_constituents_to_lab_frame()
        self.boost_constituents_to_lab_frame()
        if self.name == "Meson::Uncollidable" :
            for ifund in self.constituents :
                ifund.name += "::Uncollidable"

    def __repr__(self):
        return self.name + " : Mass " + str(self.M) + " with " + str(self.p) + " at " + str(self.r3)  + "\n" \
               + "\t\t\tNumber of constituents is " + str(self.N)

    def __del__(self):
        """
        Default destructor isn't sufficient since the core's image stays in the scene
        """

    def safe_remove(self, other):
        # can probably speed this up later with class type-checks
        for ifund in self.constituents[:] :
            if ifund == other :
                ifund.image.visible = 0
                self.constituents.remove(ifund)
                break
        for ifund in self.billiards[:] :
            if ifund == other :
                ifund.image.visible = 0
                self.billiards.remove(ifund)
                break

    def create_constituents(self):
        # create the constituents with random energy fractions of the total
        # and set them in motion in random radial directions in composite frame

        cutoffFraction = 0.05
        energyLeft = totalEnergy = self.M

        fail_counter = 0

        # do the creation of fundamentals in the composite's rest frame
        while energyLeft > cutoffFraction * totalEnergy:

            # generate energetic emission from uniform distribution of energy left
            # energy_draw = energyLeft * random()

            if fail_counter > 9998 :
                break

            # generate energetic emission from Sudakov (lognormal) distribution of energy left
            peak_energy = energyLeft / totalEnergy * np.e / 20  # Mode = 1/20, mean ~ 1/5
            energy_draw = energyLeft * np.random.lognormal(np.log(peak_energy), 1)
            if 2 * energy_draw > energyLeft:
                fail_counter += 1
                continue

            random_dir = Vec3.random_dir()

            fundamental_to_add1 = Fundamental( name="Quark",
                                               energy=energy_draw,
                                               direction=random_dir)
            if fundamental_to_add1.name == "Tachyon":
                fail_counter += 1
                continue
            fundamental_to_add2 = Fundamental( name="Quark",
                                               energy=energy_draw,
                                               direction=-random_dir)

            self.constituents.append(fundamental_to_add1)
            self.constituents.append(fundamental_to_add2)

            # add the mediators
            mediator_to_add1 = Mediator()
            mediator_to_add2 = Mediator()
            self.mediators.extend( [ mediator_to_add1 , mediator_to_add2 ] )
            self.dict_constituents_mediators[ fundamental_to_add1 ] = mediator_to_add1
            self.dict_constituents_mediators[ fundamental_to_add2 ] = mediator_to_add2

            energyLeft -= 2 * energy_draw

        print(fail_counter+1 , " tries to create " , self.name , "'s constituents ")

        # finish by creating a fundamental with the remaining energy
        final_fundamental_to_add = Fundamental( name="Core", energy=energyLeft, mass=energyLeft)
        self.constituents.append( final_fundamental_to_add )

        # add its mediator
        final_mediator_to_add = Mediator()
        self.mediators.append( final_mediator_to_add )
        self.dict_constituents_mediators[final_fundamental_to_add] = final_mediator_to_add

        self.N = len(self.constituents)

    def translate_constituents_to_lab_frame(self):
        for particle in self.constituents:
            particle.r3 = particle.r3.as_seen_in_frame_with( CoM_pos3=self.r3 , CoM_mom4=self.p )

    def translate_constituents_to_rest_frame(self):
        for particle in self.constituents:
            particle.r3 = particle.r3.as_seen_in_rest_frame_of( CoM_pos3=self.r3 , CoM_mom4=self.p )

    def boost_constituents_to_lab_frame(self):
        boost_mom = self._p
        for particle in self.constituents:
            particle.boost_away_with_momentum(boost_mom)
        for mediator in self.mediators:
            mediator.p = mediator.p.boosted_away_with_momentum(boost_mom)

    def boost_constituents_to_rest_frame(self):
        boost_mom = self._p
        for particle in self.constituents:
            particle.boost_to_rest_frame_of(boost_mom)
        for mediator in self.mediators:
            mediator.p = mediator.p.boosted_to_rest_frame_of(boost_mom)

    def hide_all_constituents(self):
        for ifund in self.constituents :
            ifund.image.visible = 0

    # primary getters and setters
    @property
    def p(self):
        return self._p
    @p.setter
    def p(self, vec4):
        # print("Inside Mom4 setter for ",self,". Careful, this changes the composite's mass!")
        self._M = sqrt(vec4*vec4)
        self._p = vec4

    @property
    def r3(self):
        return self._r3
    @r3.setter
    def r3(self, vec3):
        self._r3 = vec3

    @property
    def M(self):
        return self._M

    # secondary getters and setters
    @property
    def x(self):
        return self.r3.x
    @x.setter
    def x(self, posx):
        self.r3.x = posx
    @property
    def y(self):
        return self.r3.y
    @y.setter
    def y(self, posy):
        self.r3.y = posy
    @property
    def z(self):
        return self.r3.z
    @z.setter
    def z(self, posz):
        self.r3.z = posz

    @property
    def E(self):
        return self.p.E
    @E.setter
    def E(self, momE):
        print("You're manually setting energy=",momE,"for",self.name,". This changes the 3-momentum from",self.p3)
        self.boost_constituents_to_rest_frame()
        new_p3mag = momE**2 - self.M**2
        try:
            new_p3mag = sqrt( momE**2 - self.M**2 )
        except ValueError :
            print("{}^2 - {}^2 < 0 : Cannot find new 3-momentum for {}".format( momE, self.M, self))
        self._p = Mom4(momE, self.dir.x*new_p3mag, self.dir.y*new_p3mag, self.dir.z*new_p3mag)
        self.boost_constituents_to_lab_frame()
        print("\t\t\t\t\t to",self.p3)
    @property
    def px(self):
        return self.p.px
    @px.setter
    def px(self, momx):
        new_energy = sqrt( momx**2 + self.py**2 + self.pz**2 + self.M**2 )
        self._p = Mom4(new_energy, momx, self.py, self.pz)
    @property
    def py(self):
        return self.p.py
    @py.setter
    def py(self, momy):
        new_energy = sqrt( self.px**2 + momy**2 + self.pz**2 + self.M**2 )
        self._p = Mom4(new_energy, self.px, momy, self.pz)
    @property
    def pz(self):
        return self.p.pz
    @pz.setter
    def pz(self, momz):
        new_energy = sqrt( self.px**2 + self.py**2 + momz**2 + self.M**2 )
        self._p = Mom4(new_energy, self.px, self.py, momz)
    @property
    def p3(self):
        return Vec3(self.px, self.py, self.pz)
    @p3.setter
    def p3(self, vec3):
        new_energy = sqrt( vec3.mag2 + self.M**2 )
        self._p = Mom4(new_energy, vec3.x, vec3.y, vec3.z)
    @property
    def p3mag(self):
        return self.p3.mag
    @property
    def M2(self):
        return self.M*self.M

    @property
    def dir(self):
        return self.p.dir
    @dir.setter
    def dir(self, vec3):
        self._p.dir = vec3

    @property
    def vx(self):
        return self._p.v3.x
    @property
    def vy(self):
        return self._p.v3.y
    @property
    def vz(self):
        return self._p.v3.z
    @property
    def v3(self):
        return self._p.v3
    @property
    def v3mag(self):
        return self._p.v3.mag
    @property
    def v3mag2(self):
        return self._p.v3.mag2
    @property
    def gamma(self):
        return self._p.gamma

    def core_spring_force( self , mass , future_timestep , r3 = Vec3(0,0,0) , p3=Vec3(0,0,0)  ):

        # S = int dt ( -m/gamma_rvec - k rvec^2/2) --> relativistic harmonic oscillator
        # then in boosted frame, use t=gamma_beta(t' - beta r_par'),
        #                            r_perp = rperp', r_par = gamma_beta(r_par' - beta*t')
        # the use ELE to find eq of motion

        delta = r3.as_seen_in_rest_frame_of( CoM_pos3 = self.r3 + self.v3*future_timestep, CoM_mom4 = self.p)

        delta_par = (delta*self.v3.hat)*self.v3.hat
        delta_perp = delta - delta_par

        v3_par = self.v3.hat*(p3*self.v3.hat)/sqrt( p3.mag2 + mass**2 )
        v3_perp = p3/sqrt( p3.mag2 + mass**2 ) - v3_par

        rest_frame_force_par  = -const.FORCE_CONSTANT*delta_par
        rest_frame_force_perp = -const.FORCE_CONSTANT*delta_perp

        lab_frame_force_par = rest_frame_force_par + self.gamma*self.v3*(rest_frame_force_perp*v3_perp)
        lab_frame_force_perp = self.gamma*rest_frame_force_perp*(1-self.v3*v3_par)

        # forceX = delta.x * self.gamma * ( 1 - p3.z * self.vz / sqrt( p3.mag2 + mass**2 ) )
        # forceY = delta.y * self.gamma * ( 1 - p3.z * self.vz / sqrt( p3.mag2 + mass**2 ) )
        # forceZ = delta.z + ( self.vz * self.gamma
        #                        * ( delta.x * p3.x + delta.y * p3.y ) / sqrt( p3.mag2 + mass**2 ) )
        # Above not symmetric in x,y and z, since assumed composite travelled entirely in z direction

        return lab_frame_force_par + lab_frame_force_perp

    def phasespace_force(self, mass , timestep_into_future , r3=Vec3(0,0,0) , p3 = Vec3(0,0,0) ):

        # returns a ( Vec3 , Vec3 ) pair representing the instantaneous change in ( r3, p3 )
        rDot = p3 / sqrt( p3.mag2 + mass**2 )
        pDot = self.core_spring_force( r3 = r3 , p3 = p3 , mass = mass , future_timestep = timestep_into_future )

        return ( rDot , pDot )

    def phasespace_force_RK4(self, time_step , fundamental):

        pos_now = fundamental.r3
        mom_now = fundamental.p3
        mass_now = fundamental.M

        K1 = self.phasespace_force( timestep_into_future = 0 , r3 = pos_now , p3 = mom_now , mass = mass_now )
        K2 = self.phasespace_force( timestep_into_future = time_step/2 , r3 = pos_now + (time_step/2)*K1[0] ,
                                    p3 = mom_now + (time_step/2)*K1[1] , mass = mass_now   )
        K3 = self.phasespace_force( timestep_into_future = time_step/2 , r3 = pos_now + (time_step/2)*K2[0] ,
                                    p3 = mom_now + (time_step/2)*K2[1] , mass = mass_now  )
        K4 = self.phasespace_force( timestep_into_future = time_step , r3 = pos_now + time_step*K3[0] ,
                                    p3 = mom_now + time_step*K3[1] , mass = mass_now )

        return ( (K1[0] + 2*K2[0] + 2*K3[0] + K4[0])/6 , (K1[1] + 2*K2[1] + 2*K3[1] + K4[1])/6 )

    def update_comp(self, dt):
        self.r3 += self.v3 * dt
        for particle in self.constituents:
            phase_space_forces = self.phasespace_force_RK4( time_step = dt , fundamental = particle )
            particle.update(dt=dt,
                            ps_forces = phase_space_forces,
                            mediators = [self.dict_constituents_mediators[particle]]
                           )

    def get_collision_pairs(self, other_fundamentals = []):

        annihilation_pairs = []
        billiard_pairs = []

        for particle in self.constituents:
            annihilation_pair, billiard_pair = particle.get_fund_collision_pairs( other_fundamentals )
            annihilation_pairs.extend( annihilation_pair )
            billiard_pairs.extend( billiard_pair)

        return annihilation_pairs, billiard_pairs

    def spring_breaks_and_kinematically_allowed(self, fund):

        mediator = self.dict_constituents_mediators[fund]

        possible_meson_energy = sqrt(fund.p3mag**2 + const.MESON_MASS**2)
        possible_meson_p3 = fund.p3

        # Needs to be energetically favourable to create a meson. Two conditions:
        # 1) From |Emed| + Equark - Emeson = |Emed'| > 0
        min_mediator_energy = possible_meson_energy - fund.E
        # 2) From P -> P' + pmeson, P'^2 = P^2 + pmeson^2 - 2P.pmeson > Lambda^2 > 0
        remnant_virtuality = (self.M**2 + const.MESON_MASS**2
                              - 2*(self.p.E*possible_meson_energy - self.p3*possible_meson_p3)
                              )
        if -mediator.E > min_mediator_energy and remnant_virtuality > 4*const.MESON_MASS**2 :
            print("\nBREAKING STRINGS!")
            print("Remnant mass predicted to be",sqrt(remnant_virtuality))
            return True
        return False

    @property
    def net_mass(self):
        try:
            return sqrt(self.net_mom*self.net_mom)
        except ValueError:
            print(self.net_mom," has negative mass")
    @property
    def net_pos(self):
        return self.center_of_energy
    @property
    def center_of_energy(self):
        sum_weighted_pos = Vec3()
        sumE = 0
        for particle in self.constituents:
            sum_weighted_pos += particle.E * particle.r3
            sumE += particle.E

        return sum_weighted_pos / sumE
    @property
    def net_mom(self):
        sum_4mom = Mom4()
        for particle in self.constituents:
            sum_4mom += particle.p
        for mediator in self.mediators :
            sum_4mom += mediator.p
        return sum_4mom

    def recalculate_composite_properties(self):
        print("Before recalculating, we currently think",self)
        p_prime = self.net_mom
        self._p = p_prime
        try:
            self._M = sqrt(p_prime*p_prime)
        except ValueError:
            print( "{} <= 0: no mass for spacelike {} in {}".format(p_prime*p_prime,p_prime,self) )

        self.N = len(self.constituents)
        print("So after recalculating:",self)

    def print(self):
        print(self, "\n")
