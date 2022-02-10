
# math sqrt is faster than numpy sqrt for single values. Only use numpy sqrt for array manips
from math import sqrt, sin, cos, acos, pi
import vpython as vp

# for random Vec3 directions
from random import random

"""
4-vector representation as a tuple, with convenient functions and syntactic sugar
"""


class Vec4:
    """
    Representation of a Lorentz 4-vector.
    """
    def __init__(self, t=0., x=0., y=0., z=0.):
        self._values = [t,x,y,z]

    def __repr__(self):
        """
        Returns a string representation of the object.
        """
        return "{}[{:.5f},{:.5f},{:.5f},{:.5f}]".format(self.__class__.__name__,self[0],self[1],self[2],self[3])

    @property
    def values(self):
        return self._values
    @values.setter
    def values(self, other):
        self._values = other.values

    def __getitem__(self, i):
        """
        Square-brackets return component
        """
        return self._values[i]
    def __setitem__(self, i, val):
        """
        Square-brackets set component
        """
        self._values[i] = val

    def __eq__(self,other):
        return self[0] == other[0] and self[1] == other[1] and self[2] == other[2] and self[3] == other[3]

    def __add__(self, other):
        """
        Vec4(self) + Vec4(other)
        """
        new_vec = self.__class__(self[0]+other[0],self[1]+other[1],self[2]+other[2],self[3]+other[3])
        return new_vec
    def __iadd__(self, other):
        """
        Vec4(self) += Vec4(other)
        """
        for i in range(0,4) :
            self[i] += other[i]
        return self

    def __sub__(self, other):
        """
        Vec4(self) - Vec4(other)
        """
        new_vec = self.__class__(self[0]-other[0],self[1]-other[1],self[2]-other[2],self[3]-other[3])
        return new_vec
    def __isub__(self, other):
        """
        Vec4(self) -= Vec4(other)
        """
        for i in range(0,4) :
            self[i] -= other[i]
        return self

    def __mul__(self, other):
        """
        Vec4(self) * Vec4(other)
        """
        scalar = self[0]*other[0] - self[1]*other[1] - self[2]*other[2] - self[3]*other[3]
        return scalar

    @property
    def invariant(self):
        return self*self

    # pylint: disable=too-many-arguments
    def lorentz_transformation(self, betax = 0., betay = 0., betaz = 0.):

        gamma = 1 / sqrt( 1 - betax**2 - betay**2 - betaz**2 )
        rat   = gamma**2 / ( 1 + gamma )

        mat = [ [       gamma,    -gamma*betax,     -gamma*betay,    -gamma*betaz],
                [-gamma*betax,  1+rat*betax**2,  rat*betax*betay, rat*betax*betaz],
                [-gamma*betay, rat*betax*betay,   1+rat*betay**2, rat*betay*betaz],
                [-gamma*betaz, rat*betax*betaz,  rat*betay*betaz,  1+rat*betaz**2]  ]

        vec4_prime = [0, 0, 0, 0]
        for i in range(0, 4):
            for j in range(0, 4):
                vec4_prime[i] += mat[i][j]*self[j]

        return self.__class__(vec4_prime[0],vec4_prime[1],vec4_prime[2],vec4_prime[3])


class Mom4(Vec4):
    """
    Extending Lorentz 4-vector functionality for a 4-momentum
    """

    @property
    def E(self):
        return self[0]
    @E.setter
    def E(self, val):
        self[0] = val
    @property
    def px(self):
        return self[1]
    @px.setter
    def px(self, val):
        self[1] = val
    @property
    def py(self):
        return self[2]
    @py.setter
    def py(self, val):
        self[2] = val
    @property
    def pz(self):
        return self[3]
    @pz.setter
    def pz(self, val):
        self[3] = val
    @property
    def p3(self):
        return Vec3(self[1], self[2], self[3])
    @p3.setter
    def p3(self, v3):
        self[1] = v3.x
        self[2] = v3.y
        self[3] = v3.z
    @property
    def v3(self):
        energy = self[0]
        return Vec3( self[1]/energy , self[2]/energy, self[3]/energy )
    @property
    def v3sqr(self):
        return self.v3 * self.v3
    @property
    def gamma(self):
        return sqrt( 1 / ( 1 - self.v3 * self.v3 ) )
    @property
    def dir(self):
        return Vec3(self[1], self[2], self[3]).hat
    @dir.setter
    def dir(self, v3):
        v3norm = v3.hat
        p3abs = self.p3mag
        self[1] = p3abs * v3norm.x
        self[2] = p3abs * v3norm.y
        self[3] = p3abs * v3norm.z

    @property
    def M2(self):
        return self*self
    @property
    def M(self):
        mass2 = self.M2
        if mass2 >= 0 :
            return sqrt(mass2)
        return 1j*sqrt(-mass2)

    @property
    def p3mag(self):
        return sqrt(self.p3mag2)
    @property
    def p3mag2(self):
        return self.px**2 + self.py**2 + self.pz**2

    def boosted_away_with_momentum(self, p4):
        try:
            return self.lorentz_transformation( -p4.v3.x, -p4.v3.y, -p4.v3.z )
        except ZeroDivisionError:
            print(p4*p4,"= 0 : no rest frame for lightlike momenta like",p4)
            raise ValueError
        except ValueError:
            print(p4*p4,"< 0 : no rest frame for spacelike momenta like",p4)
            raise ValueError
    def boosted_to_rest_frame_of(self, p4):
        try:
            return self.lorentz_transformation( p4.v3.x, p4.v3.y, p4.v3.z )
        except ZeroDivisionError:
            print(p4 * p4, "= 0 : no rest frame for lightlike momenta like", p4)
            raise ValueError
        except ValueError:
            print(p4 * p4, "< 0 : no rest frame for spacelike momenta like", p4)
            raise ValueError


class Pos4(Vec4):
    """
    Extending Lorentz 4-vector functionality for a 4-position
    """

    @property
    def t(self):
        return self[0]
    @t.setter
    def t(self,val):
        self[0] = val
    @property
    def x(self):
        return self[1]
    @x.setter
    def x(self,val):
        self[1] = val
    @property
    def y(self):
        return self[2]
    @y.setter
    def y(self,val):
        self[2] = val
    @property
    def z(self):
        return self[3]
    @z.setter
    def z(self,val):
        self[3] = val

    @property
    def r3(self):
        return Vec3(self[1],self[2],self[3])
    @r3.setter
    def r3(self,v3):
        self[1] = v3.x
        self[2] = v3.y
        self[3] = v3.z

    def as_vec3(self):
        return Vec3(*self._values[1:])
    def as_pos4(self):
        return Pos4(*self._values[0:])


class Vec3:

    def __init__(self, x=0., y=0., z=0.):
        self._values = [x,y,z]

    def __repr__(self):
        """
        Returns a string representation of the object.
        """
        return "{}[{:.5f},{:.5f},{:.5f}]".format(self.__class__.__name__,self[0],self[1],self[2])

    @property
    def values(self):
        return self._values
    @values.setter
    def values(self, other):
        self._values = other.values

    def __getitem__(self, i):
        """
        Square-brackets return component
        """
        return self._values[i]
    def __setitem__(self, i, val):
        """
        Square-brackets set component
        """
        self._values[i] = val

    def __eq__(self,other):
        return self[0] == other[0] and self[1] == other[1] and self[2] == other[2]

    def __add__(self, other):
        """
        Vec3(self) + Vec3(other)
        """
        new_vec = Vec3(self[0] + other[0], self[1] + other[1], self[2] + other[2])
        return new_vec
    def __iadd__(self, other):
        """
        Vec3(self) += Vec3(other)
        """
        for i in range(0, 3):
            self[i] += other[i]
        return self

    def __neg__(self):
        """
        unary (-self)
        :return:
        """
        return Vec3(-self[0], -self[1], -self[2])
    def __sub__(self, other):
        """
        Vec3(self) - Vec3(other) = Vec3(self) + ( -Vec3(other) )
        """
        new_vec = self + (-other)
        return new_vec
    def __isub__(self, other):
        """
        Vec3(self) -= Vec3(other)
        """
        for i in range(0, 3):
            self[i] -= other[i]
        return self

    def __mul__(self, other):
        """
        scalar = Vec3(self) * Vec3(other)
        or
        Vec(3) = Vec3(self) * scalar
        I would really like a branch-free way of doing this
        """
        if isinstance(other, (float, int)):
            return Vec3( self.x*other, self.y*other, self.z*other)
        elif isinstance(other,Vec3):
            return self.x*other.x + self.y*other.y + self.z*other.z
        else:
            return NotImplemented
    def __rmul__(self, other):
        return self*other

    def __truediv__(self, other):
        return self * (1/other)

    @property
    def x(self):
        return self[0]
    @x.setter
    def x(self,val):
        self[0] = val
    @property
    def y(self):
        return self[1]
    @y.setter
    def y(self,val):
        self[1] = val
    @property
    def z(self):
        return self[2]
    @z.setter
    def z(self,val):
        self[2] = val

    @property
    def mag2(self):
        return self[0]**2 + self[1]**2 + self[2]**2
    @property
    def mag(self):
        return sqrt( self.mag2 )
    @mag.setter
    def mag(self, val):
        direction = self.hat
        self[0] = val * direction[0]
        self[1] = val * direction[1]
        self[2] = val * direction[2]

    @property
    def hat(self):
        norm = self.mag
        try :
            return self / norm
        except ZeroDivisionError:
            return Vec3(0,0,1)  # still need to return a vector of length 1
    @hat.setter
    def hat(self, v3):
        norm = self.mag
        direction = v3.hat
        self[0] = norm * direction[0]
        self[1] = norm * direction[1]
        self[2] = norm * direction[2]

    def as_vec3(self):
        return Vec3(*self._values[0:])
    def as_pos4(self):
        return Pos4(0,*self._values[0:])

    def as_seen_in_rest_frame_of(self, CoM_pos3 = Pos4() , CoM_mom4 = Mom4(1,0,0,0)):
        CoM_pos3 = CoM_pos3.as_vec3()
        try:
            delta4 = Pos4()
            delta4.r3 = self - CoM_pos3
            return delta4.lorentz_transformation(CoM_mom4.v3.x,CoM_mom4.v3.y,CoM_mom4.v3.z).as_vec3()
        except ZeroDivisionError:
            print(CoM_mom4 * CoM_mom4, " = 0 : no rest frame for lightlike momenta like ", CoM_mom4)
            raise ValueError
        except ValueError:
            print(CoM_mom4 * CoM_mom4, " < 0 : no rest frame for spacelike momenta like ", CoM_mom4)
            raise ValueError

    def as_seen_in_frame_with(self, CoM_pos3 = Pos4() , CoM_mom4 = Mom4(1,0,0,0) ):
        CoM_pos3 = CoM_pos3.as_vec3()
        try:  # have to implement this ourselves so we don't deal with time mixing
            vec_par = (self*CoM_mom4.dir) * CoM_mom4.dir
            vec_perp = self - vec_par
            vec_par_prime = vec_par / CoM_mom4.gamma
            vec_perp_prime = vec_perp
            return CoM_pos3 + vec_perp_prime + vec_par_prime
        except ZeroDivisionError:
            print(CoM_mom4 * CoM_mom4, " = 0 : no rest frame for lightlike momenta like ", CoM_mom4)
            raise ValueError
        except ValueError:
            print(CoM_mom4 * CoM_mom4, " < 0 : no rest frame for spacelike momenta like ", CoM_mom4)
            raise ValueError

    def to_vpython(self):
        return vp.vector(self.x,self.y,self.z)

    @staticmethod
    def random_dir():
        mu  = (2*random() - 1)
        phi = (2*random() - 1)*pi
        theta = acos(mu)
        return Vec3( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) )
    @staticmethod
    def random_dir_2D_around(other) :
        other_hat = other.hat
        rand_dir = Vec3.random_dir()
        parallel_part = (rand_dir*other_hat)*other_hat
        perp_part = rand_dir - parallel_part
        return perp_part.hat


