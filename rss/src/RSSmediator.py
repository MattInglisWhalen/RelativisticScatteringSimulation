
from rss.src.RSSvectors import *

""" The point of this object is to handle the energy and momentum stored in a spring connecting
    one Fundamental with another point in space. Composite should have a list of Mediators, and
    these Mediators should only be associated to a Fundamental through a map/dict """


class Mediator(object):

    # Class variables
    # p -- momentum

    def __init__( self , mom4 = Mom4(0,0,0,0) ):
        self._p = mom4

    def __repr__(self):
        return str(self.p)

    @property
    def p(self):
        return self._p
    @p.setter
    def p(self, vec4):
        self._p = vec4
    @property
    def E(self):
        return self.p.E
    @E.setter
    def E(self, val):
        self.p.E = val
    @property
    def px(self):
        return self.p.px
    @px.setter
    def px(self, val):
        self.p.px = val
    @property
    def py(self):
        return self.p.py
    @py.setter
    def py(self, val):
        self.p.py = val
    @property
    def pz(self):
        return self.p.pz
    @pz.setter
    def pz(self, val):
        self.p.pz = val

    @property
    def M(self):
        return np.sqrt( self.M2 )
    @property
    def M2(self):
        return self.p*self.p

    def print(self):
        print(self,"\n")
