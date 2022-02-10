
# external libs
from vpython import *

# Internal classes
from rss.src.RSSvectors import *
import rss.src.RSSconstants as const


class Fundamental(object):

    # Class variables
    # name
    # _r3 -- position
    # _p -- momentum
    # M -- mass
    # _image

    # billiards[] -- [a required list to work symbiotically with Composite functionality]

    def __init__( self , name="Point" , pos=Vec3(0,0,0) , energy=1 , direction = Vec3(0,0,1) , mass=0.01 ):

        self.name = name

        # Additional tags can be added after the first few letters to designate uncollidables
        if self.name[0:6] == "Lepton":
            self._M = const.LEPTON_MASS
        elif self.name[0:5] == "Quark" :
            self._M = const.QUARK_MASS
        else:
            self._M = mass

        if energy < self._M :  # if is faster than try because I expect it to happen a lot
            self.name = "Tachyon"
            self._M = 0  # won't actually be used but this stops the sqrt from throwing an error

        self._r3 = pos
        abs_p3 = sqrt( energy**2 - self._M**2 )
        self._p = Mom4(energy, abs_p3*direction.x, abs_p3*direction.y, abs_p3*direction.z)

        self._image = sphere( pos = self.r3.to_vpython() , radius = min(0.2,max(self._M,0.001)) )
        if name[0:4] == "Core" :
            self._image.color = color.red
        elif name[0:7] == "Virtual":
            self._image.color = color.yellow
            self._image.radius = 0.05
        elif name[0:6] == "Lepton" :
            self._image.color = color.cyan
            self._image.radius = 0.05
        elif name[0:3] == "Jet" :
            self._image.visible = 0
            self._image = cone( radius = 0.1/self.v3mag , pos = self.r3.to_vpython(), axis=-self.v3.to_vpython()/1000 )
            self._image.color = color.purple

        self.billiards = []

    def __repr__(self):
        """
        Returns the string representation of the particle
        """
        return self.name + " : Mass {} with {} at {}".format(self.M,self.p,self.r3)

    def __del__(self):
        """
        Default destructor isn't sufficient since the image stays in the scene
        """
        self._image.visible = 0

    @property
    def r3(self):
        return self._r3
    @r3.setter
    def r3(self, vec3):
        self._r3 = vec3
        self.image.pos = self.r3.to_vpython()

    @property
    def M(self):
        return self._M
    @M.setter
    def M(self, val):
        print("Setting the mass of",self,". Careful, a fundamental's mass shouldn't change!")
        self._M = val

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
    def p(self):
        return self._p
    @p.setter
    def p(self, vec4):
        try:
            print(f"Inside Mom4 setter for {self}. Careful, this changes the fundamental's mass to {sqrt(vec4 * vec4)}")
            self._M = sqrt(vec4*vec4)
        except ValueError:
            print(vec4*vec4,"<= 0 : Cannot set mass for",self,"with",vec4)
        self._p = vec4
        print(f"So now we have set {self}")
    @property
    def E(self):
        return self.p.E
    @E.setter
    def E(self, momE):
        new_p3mag = sqrt( momE**2 - self.M**2 )
        self._p = Mom4( momE , self.dir.x*new_p3mag , self.dir.y*new_p3mag , self.dir.z*new_p3mag )
    @property
    def px(self):
        return self.p.px
    @px.setter
    def px(self, momx):
        new_energy = sqrt( momx**2 + self.py**2 + self.pz**2 + self.M**2 )
        self._p = Mom4( new_energy , momx , self.py , self.pz )
    @property
    def py(self):
        return self.p.py
    @py.setter
    def py(self, momy):
        new_energy = sqrt( self.px**2 + momy**2 + self.pz**2 + self.M**2 )
        self._p = Mom4( new_energy , self.px , momy , self.pz )
    @property
    def pz(self):
        return self.p.pz
    @pz.setter
    def pz(self, momz):
        new_energy = sqrt( self.px**2 + self.py**2 + momz**2 + self.M**2 )
        self._p = Mom4( new_energy , self.px , self.py , momz )

    @property
    def p3(self):
        return self.p.p3
    @p3.setter
    def p3(self, vec3):
        new_energy = sqrt(vec3.mag2 + self.M**2)
        self._p = Mom4( new_energy , vec3.x , vec3.y , vec3.z )
    @property
    def p3mag(self):
        return self.p3.mag
    @property
    def p3mag2(self):
        return self.p3.mag2
    @property
    def M2(self):
        return self._M**2

    @property
    def dir(self):
        return self.p.dir
    @dir.setter
    def dir(self, vec3):
        self._p = Mom4(self.E, vec3.hat.x*self.p3mag, vec3.hat.y*self.p3mag, vec3.hat.z*self.p3mag )

    @property
    def vx(self):
        return self.p.v3.x
    @property
    def vy(self):
        return self.p.v3.y
    @property
    def vz(self):
        return self.p.v3.z
    @property
    def v3(self):
        return self.p.v3
    @property
    def v3mag(self):
        return self.p.v3.mag
    @property
    def v3mag2(self):
        return self.p.v3.mag2
    @property
    def gamma(self):
        return self.p.gamma


    @property
    def image(self):
        return self._image
    @image.setter
    def image(self, vispy_object):
        self._image = vispy_object

    def print(self):
        print(self,"\n")

    def boost_away_with_momentum(self, boost_mom):
        self._p = self.p.boosted_away_with_momentum(boost_mom)

    def boost_to_rest_frame_of(self, boost_mom):
        self._p = self.p.boosted_to_rest_frame_of(boost_mom)

    """
    Collision handling
    """
    def collidable_with(self,other):
        collidable_list = ["Lepto","Quark"]
        if self.name[0:5] in collidable_list and other.name[0:5] in collidable_list :
            return True
        return False
    def are_collided(self, other ):
        collision_zone = self.image.radius + other.image.radius
        if (self.r3 - other.r3).mag < collision_zone and self.collidable_with(other) :
            return True
        return False
    def which_collided(self , others):
        collided_list = []
        for iother in others:
            if self.are_collided(iother):
                collided_list.append(iother)
        return collided_list

    def do_billiards_collision_with(self,other):

        print("In fundamental billiards method,", self.p ,"+",other.p,"=",self.p+other.p)
        print("Positions are",self.r3,"and",other.r3)
        # scattered fundamentals are coloured blue
        self.image.color = color.blue
        other.image.color = color.blue

        # find CoM momentum, and its conjugate for boosting with boost_particle routine
        CoM_mom = self.p + other.p
        CoM_pos = (self.E*self.r3 + other.E*other.r3)/(self.E+other.E)

        # boost to CoM frame, net 3 momenta there should be zero
        boosted_CoM_p = CoM_mom.boosted_to_rest_frame_of( CoM_mom )
        boosted_self_p = self.p.boosted_to_rest_frame_of( CoM_mom )
        boosted_other_p = other.p.boosted_to_rest_frame_of( CoM_mom )

        print("Boosted p3s are", boosted_self_p.p3 ,"+",boosted_other_p.p3,"=", boosted_self_p.p3+boosted_other_p.p3)

        boosted_self_p3 = boosted_self_p.p3
        boosted_other_p3 = boosted_other_p.p3

        # masses for shorter expressions
        CoM_M2 = boosted_CoM_p.M2
        self_M2 = boosted_self_p.M2
        other_M2 = boosted_other_p.M2
        Com_E = boosted_CoM_p.E

        # (p1 + p2 - p1Prime)^2 = p2Prime^2 gives the following for the new p1Prime energy
        E_self_before = boosted_self_p.E
        E_self_after  = ( CoM_M2 + self_M2 - other_M2 ) / ( 2*Com_E )
        p3mag_self = sqrt( E_self_after**2 - self_M2 )

        E_other_before = boosted_other_p.E
        E_other_after  = ( CoM_M2 + other_M2 - self_M2 ) / ( 2*Com_E )
        p3mag_other = sqrt( E_other_after**2 - other_M2 )

        # Impulse is only applied in the impact's normal direction
        # so change in energy is only due to change in momentum in normal direction, in CoM frame
        pos_self_boosted = self.r3.as_seen_in_rest_frame_of( CoM_pos3=CoM_pos , CoM_mom4=CoM_mom )
        pos_other_boosted = other.r3.as_seen_in_rest_frame_of( CoM_pos3=CoM_pos , CoM_mom4=CoM_mom )
        normal_vec3 = Vec3( pos_self_boosted.x - pos_other_boosted.x ,
                            pos_self_boosted.y - pos_other_boosted.y ,
                            pos_self_boosted.z - pos_other_boosted.z  ).hat
        print("nHat=",normal_vec3)

        normal_self_p3mag = boosted_self_p3*normal_vec3
        normal_self_p3 = normal_self_p3mag*normal_vec3
        tang_self_p3 = boosted_self_p3 - normal_self_p3
        print("so pPerp=", tang_self_p3,"and pPar=",normal_self_p3)
        print(E_self_after,"^2 + ",normal_self_p3mag,"^2 compared with", E_self_before, "^2")
        normal_self_p3mag_after = -normal_self_p3mag
        try :
            normal_self_p3mag_after *= sqrt( 1 + (E_self_after**2 - E_self_before**2) / normal_self_p3mag**2  )
        except ZeroDivisionError:
            pass
        p3dir_self = (tang_self_p3 + normal_self_p3mag_after*normal_vec3).hat

        normal_other_p3mag = boosted_other_p3*normal_vec3
        normal_other_p3 = normal_other_p3mag*normal_vec3
        tang_other_p3 = boosted_other_p3 - normal_other_p3
        print(E_other_after,"^2 + ",normal_other_p3mag,"^2 compared with", E_other_before, "^2")
        normal_other_p3mag_after = -normal_other_p3mag
        try :
            normal_other_p3mag_after *= sqrt( 1 + (E_other_after**2 - E_other_before**2) / normal_other_p3mag**2  )
        except ZeroDivisionError:
            pass

        p3dir_other = (tang_other_p3 + normal_other_p3mag_after*normal_vec3).hat

        print("In CoM frame dirs before=",boosted_self_p.dir,boosted_other_p.dir)
        print("               and after=",p3dir_self,p3dir_other)

        # reconstruct the 4-momenta after the collision, in the CoM frame
        boosted_self_p_after = Mom4( E_self_after , p3dir_self.x * p3mag_self
                                                  , p3dir_self.y * p3mag_self
                                                  , p3dir_self.z * p3mag_self )

        boosted_other_p_after = Mom4( E_other_after , p3dir_other.x * p3mag_other
                                                    , p3dir_other.y * p3mag_other
                                                    , p3dir_other.z * p3mag_other )

        print("End fundamental billiards method,",  boosted_self_p_after.boosted_away_with_momentum( CoM_mom ),"+",
                                                    boosted_other_p_after.boosted_away_with_momentum( CoM_mom ),"=",
              boosted_self_p_after.boosted_away_with_momentum( CoM_mom ) + boosted_other_p_after.boosted_away_with_momentum( CoM_mom ))

        # boost back to lab frame
        self._p = boosted_self_p_after.boosted_away_with_momentum( CoM_mom )
        other._p = boosted_other_p_after.boosted_away_with_momentum( CoM_mom )

    def update(self, dt, ps_forces = (Vec3(0, 0, 0) , Vec3(0, 0, 0)), mediators = [] ):

        # need to store the momentum because the 3-force doesn't say anything about the energy change
        mom4_before = self.p

        # Update based on RK4 forces
        self.r3 = self.r3 + dt*ps_forces[0]
        self.p3 = self.p3 + dt*ps_forces[1]

        for imed in mediators :  # using a for loop over length-one list to avoid if statement
            # print(mom4_before - self.p, " ", mom4_before , " " , self.p)
            imed.p = imed.p + ( mom4_before - self.p )

        self.image.pos = self.r3.to_vpython()
        self.image.axis -= dt*ps_forces[0].to_vpython()
        self.billiards = []

    def get_fund_collision_pairs(self, other_fundamentals = []):

        # Now handle collision logic
        assert(self.name != "Virtual")

        elastic_probability = 0.5
        """ Probability any collision will simply be elastic.
            Conversely, the remaining probability is for the pair to become a virtual particle """

        collisions = self.which_collided( other_fundamentals )
        annihilation_pair = []
        billiard_pair = []
        for other in collisions:
            if random() < elastic_probability or (self.p + other.p)*(self.p + other.p) < 4*const.LEPTON_MASS**2 :
                # elastic, or not enough energy to create a virtual, so do a billiards collision
                print("Scattered",self , " and " , other)
                self.name = "Billiard"
                other.name = "Billiard"
                billiard_pair.append( (self,other) )
            else:  # Annihilation reaction
                print("Annihilated", self , " and " , other)
                self.name = "Annihilating"
                other.name = "Annihilating"
                annihilation_pair.append( (self , other) )
            break  # only consider two-body collisions, so just look at the first in the list

        return annihilation_pair, billiard_pair

    def get_decay_pairs(self, dt) :
        assert( self.name == "Virtual" or self.name == "Jet" )

        parent_plus_decay_pair = []
        # decay probability in time deltaT is p0 = 1 - exp( -M*deltaT )
        # which gives a lifetime tau = 1/M
        if random() < ( 1 - exp( -self.M*dt ) ) :

            lepton_branching_fraction = 0.1

            # boost rest frame of the virtual CoM frame, net 3 momenta there should be zero
            boosted_CoM_p = self.p.boosted_to_rest_frame_of( self.p )

            # masses for shorter expressions
            boosted_energy_daughters = boosted_CoM_p.E/2
            rand_dir = Vec3.random_dir()

            from rss.src.RSScomposite import Composite
            if self.name == "Virtual" :

                # can decay leptonically
                if random() < lepton_branching_fraction or self.M < 2*const.MESON_MASS :
                    daughter_mass = const.LEPTON_MASS
                    daughter_type = "Lepton::Uncollidable"
                # or it can decay hadronically
                else :
                    if self.M > 4*const.MESON_MASS :
                        # create jets
                        daughter_mass = random()*(boosted_energy_daughters-2*const.MESON_MASS) + 2*const.MESON_MASS
                        daughter_type = "Jet"
                    else :
                        # create mesons
                        daughter_mass = const.MESON_MASS
                        daughter_type = "Meson::Uncollidable"
            else :  # self.name == "Jet"
                if self.M > 4*const.MESON_MASS:
                    # create jets
                    daughter_mass = random()*(boosted_energy_daughters - 2*const.MESON_MASS) + 2*const.MESON_MASS
                    daughter_type = "Jet"
                else:
                    # create mesons
                    daughter_mass = const.MESON_MASS
                    daughter_type = "Meson::Uncollidable"

            boosted_p3mag_daughters = sqrt( boosted_energy_daughters**2 - daughter_mass**2 )

            # construct decay products in the CoM frame
            boosted_decay1 = Mom4( boosted_energy_daughters , rand_dir.x * boosted_p3mag_daughters
                                                            , rand_dir.y * boosted_p3mag_daughters
                                                            , rand_dir.z * boosted_p3mag_daughters )

            boosted_decay2 = Mom4( boosted_energy_daughters , -rand_dir.x * boosted_p3mag_daughters
                                                            , -rand_dir.y * boosted_p3mag_daughters
                                                            , -rand_dir.z * boosted_p3mag_daughters )

            if daughter_type[0:6] == "Lepton" or daughter_type[0:3] == "Jet" :
                daughter1 = Fundamental(name = daughter_type,
                                        pos = self.r3,
                                        energy = boosted_decay1.boosted_away_with_momentum( self.p ).E,
                                        direction = boosted_decay1.boosted_away_with_momentum( self.p ).dir,
                                        mass = daughter_mass)
                daughter2 = Fundamental(name = daughter_type,
                                        pos = self.r3,
                                        energy = boosted_decay2.boosted_away_with_momentum( self.p ).E,
                                        direction = boosted_decay2.boosted_away_with_momentum( self.p ).dir,
                                        mass = daughter_mass )
            else :  # meson
                daughter1 = Composite(name = daughter_type,
                                      pos = self.r3,
                                      energy = boosted_decay1.boosted_away_with_momentum( self.p ).E,
                                      direction = boosted_decay1.boosted_away_with_momentum( self.p ).dir,
                                      mass = daughter_mass)
                daughter2 = Composite(name = daughter_type,
                                      pos = self.r3,
                                      energy = boosted_decay2.boosted_away_with_momentum( self.p ).E,
                                      direction = boosted_decay2.boosted_away_with_momentum( self.p ).dir,
                                      mass = daughter_mass )

            print(f"\n{self.name} ~Decaying~~ to {daughter1} and {daughter2}" )
            parent_plus_decay_pair.append( (self,daughter1,daughter2) )

        return parent_plus_decay_pair

    def recalculate_composite_properties(self):
        # this is not a composite, so
        pass
