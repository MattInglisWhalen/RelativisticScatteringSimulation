
# external libs
from vpython import *

# internal classes
from rss.src.RSSvectors import *
from rss.src.RSScomposite import Composite
from rss.src.RSSfundamental import Fundamental
import rss.src.RSSconstants as const


"""

--- Types of particles in this model ---

The "lepton": mass 0.3 GeV. Interacts via contact interactions with quarks and other leptons. 
                            Two leptons can annihilate to produce a virtual.

The "quark": mass 0.01 Gev, never "seen" in isolation. Bound to the proton's core with a spring force. Can scatter 
                elastically off leptons or quarks in other hadrons, or can 
                annihilate with another quark to produce a virtual.

The "proton": a composite particle with mass 1 GeV. Made of 
                quarks that are bound to the proton's core through a spring force
                
The "mediator": the spring which holds quarks to the proton's core. Can be severed under 
                energetically favourable (and kinematically allowed) conditions to produce a meson or a proton
                
The "meson": a composite particle with mass 0.5 GeV. You cannot start with a meson in the initial state.
                
The "virtual": a stand-in for an actual force-carrying boson. Has an arbitrary mass which serves to place the
                virtual always on-shell. Always decays to leptons, with lifetime ~ 1/M

"""


class System:

    # Class variables
    #
    # composites[] -- Composite particles which contain a collection of Fundamentals
    # fundamentals[] -- A list of fundamental particles, not including the constituents of the composites
    # virtuals[] -- Also a list of Fundamentals, but it's useful to keep them separate
    # initial collision energy

    # Class Methods
    # __init__ , __repr__
    #
    # boost_to_collision_energy()
    # handle_annihilations()
    #

    def __init__(self, mode="fundamental-fundamental", initial_energy=1):

        self.fundamentals = []
        self.composites = []

        self.virtuals = []
        self.collidable_composites = []

        starting_distance = 2
        energy_each = initial_energy / 2

        if mode == "lepton-lepton":
            self.fundamentals.append(Fundamental(name="Lepton1",
                                                 # ID=self.maxID+1, # TODO think if IDs are useful
                                                 pos=Vec3(const.LEPTON_MASS*random()/2,
                                                          const.LEPTON_MASS*random()/2,
                                                          -starting_distance),
                                                 energy=energy_each,
                                                 direction=Vec3(0, 0, 1),
                                                 )
                                     )
            self.fundamentals.append(Fundamental(name="Lepton2",
                                                 pos=Vec3(const.LEPTON_MASS*random()/2,
                                                          const.LEPTON_MASS*random()/2,
                                                          starting_distance),
                                                 energy=energy_each,
                                                 direction=Vec3(0, 0, -1),
                                                 )
                                     )
        elif mode == "lepton-proton":
            self.fundamentals.append(Fundamental(name="Lepton",
                                                 pos=Vec3(const.LEPTON_MASS*random()/2,
                                                          const.LEPTON_MASS*random()/2,
                                                          -starting_distance),
                                                 energy=energy_each,
                                                 direction=Vec3(0, 0, 1),
                                                 )
                                     )
            self.composites.append(Composite(name="Proton",
                                             pos=Vec3(const.PROTON_MASS*random()/5,
                                                      const.PROTON_MASS*random()/5,
                                                      starting_distance),
                                             energy=energy_each,
                                             direction=Vec3(0, 0, -1),
                                             )
                                   )
        elif mode == "proton-lepton":
            self.composites.append(Composite(name="Proton",
                                             pos=Vec3(const.PROTON_MASS*random()/5,
                                                      const.PROTON_MASS*random()/5,
                                                      -starting_distance),
                                             energy=energy_each,
                                             direction=Vec3(0, 0, 1),
                                             )
                                   )
            self.fundamentals.append(Fundamental(name="Lepton",
                                                 pos=Vec3(const.LEPTON_MASS*random()/2,
                                                          const.LEPTON_MASS*random()/2,
                                                          starting_distance),
                                                 energy=energy_each,
                                                 direction=Vec3(0, 0, -1),
                                                 )
                                     )

        else:  # (mode == "proton-proton")
            self.composites.append(Composite(name="Proton1",
                                             pos=Vec3(const.PROTON_MASS*random()/50,
                                                      const.PROTON_MASS*random()/50,
                                                      -starting_distance),
                                             energy=energy_each,
                                             direction=Vec3(0, 0, 1),
                                             )
                                   )
            self.composites.append(Composite(name="Proton2",
                                             pos=Vec3(const.PROTON_MASS*random()/50,
                                                      const.PROTON_MASS*random()/50,
                                                      starting_distance),
                                             energy=energy_each,
                                             direction=Vec3(0, 0, -1),
                                             )
                                   )
        self.num_fundamentals = len(self.fundamentals)
        self.num_composites = len(self.composites)
        self.collidable_composites.extend( self.composites )

    def __repr__(self):
        """
        Returns the string representation of the system
        """
        returnString = " == System stats ==\n"
        returnString += "Net mass : " + str(self.invariant_mass)
        returnString += "\nNet momentum : " + str(self.net_mom) + "\n"
        for ifund in self.fundamentals:
            returnString += str(ifund) + "\n"
        for icomp in self.composites:
            returnString += str(icomp) + "\n"
            for ifund in icomp.constituents:
                returnString += "\t" + str(ifund) + "\n"
        for ivirt in self.virtuals:
            returnString += str(ivirt) + "\n"
        return returnString

    def safe_remove(self, other):
        # can probably speed this up later with class type-checks
        if isinstance(other, Fundamental) :
            for ifund in self.fundamentals[:] :
                if ifund == other :
                    self.fundamentals.remove(ifund)
                    break
            for icomp in self.composites :
                for ifund in icomp.constituents[:] :
                    if ifund == other :
                        icomp.safe_remove(ifund)
                        break
            for icomp in self.collidable_composites :
                for ifund in icomp.constituents[:]:
                    if ifund == other :
                        icomp.safe_remove(ifund)
                        break
            for ivirt in self.virtuals[:] :
                if ivirt == other :
                    self.virtuals.remove(ivirt)
                    break
        elif isinstance(other, Composite):
            for icomp in self.composites[:]:
                if icomp == other:
                    self.composites.remove(icomp)
                    break
            for icomp in self.collidable_composites :
                if icomp == other:
                    self.collidable_composites.remove(icomp)
                    break

    def print(self):
        print(repr(self))

    @property
    def net_mom(self):
        all_sum = Mom4(0,0,0,0)
        for ifund in self.fundamentals :
            all_sum += ifund.p
        for icomp in self.composites :
            all_sum += icomp.net_mom
        for ivirt in self.virtuals :
            all_sum += ivirt.p
        return all_sum
    @property
    def invariant_mass(self):
        return sqrt(self.net_mom*self.net_mom)

    def boost_to_collision_energy(self, CoM_energy = 10 ):
        print("\n\n  ====> <====  Boosting to collision energy \n\n")
        energy_each = CoM_energy / 2
        for ifund in self.fundamentals :
            ifund.E = energy_each
        for icomp in self.composites :
            icomp.translate_constituents_to_rest_frame()
            icomp.E = energy_each
            icomp.translate_constituents_to_lab_frame()

    def handle_annihilations( self , annihilation_pairs = [] ):

        # creates a new virtual particle from each pair, and deletes
        # the associated parent particles from the relevant lists
        for pair in annihilation_pairs:

            print("\nHandling annihilations")

            virt_energy = pair[0].E + pair[1].E
            virt_p3 = pair[0].p3 + pair[1].p3
            virt_r3 = (pair[0].r3*pair[0].E + pair[1].r3*pair[1].E)/(pair[0].E+pair[1].E)

            new_virtual = Fundamental(name="Virtual",
                                      pos=virt_r3,
                                      energy=virt_energy,
                                      direction=virt_p3.hat,
                                      mass=sqrt(virt_energy**2 - virt_p3.mag2)
                                      )
            self.virtuals.append(new_virtual)
            print(new_virtual , " with mass " , new_virtual.M )

            # delete the annihilating fundamentals from their lists
            self.safe_remove(pair[0])
            self.safe_remove(pair[1])

            # gc.collect()

        for _ in annihilation_pairs[0:1] :  # branchless if len(annihilation_pairs) > 0
            for icomp in self.composites :
                icomp.recalculate_composite_properties()

    def find_owner(self, fund):
        owner = 0
        for ifund in self.fundamentals:
            if ifund == fund:
                owner = ifund
        for icomp in self.composites :
            for ifund in icomp.constituents:
                if ifund == fund:
                    owner = icomp
        return owner

    def handle_billiards_collisions(self, billiards_pairs) :

        # loop to do collisions, and then loop again to do recalculation of composite momentum
        for pair in billiards_pairs :

            print("\nHandling contact scattering")
            # find the owner of each fundamental in the pair -- the owner of a lepton fundamental is the lepton itself,
            # while the owner of a quark is the composite which contains it
            owners = self.find_owner(pair[0]), self.find_owner(pair[1])
            # change the momentum of each fundamental in an elastic collision
            pair[0].do_billiards_collision_with(pair[1])
            try:
                owners[0].billiards.append(pair[0])
                owners[1].billiards.append(pair[1])
            except AttributeError:
                print("How did an int get to be the owner??? Pair is {} and owners are {}/{}".format(pair,
                                                                                                     owners[0],
                                                                                                     owners[1]))
            print(pair)

        for _ in billiards_pairs[0:1] :  # branchless if len(billiards_pairs) > 0
            print("\n\n Recalculations of all composite momenta")
            for icomp in self.collidable_composites :
                print(icomp.name)
                icomp.recalculate_composite_properties()

    def handle_decays(self, parent_plus_decay_pairs = [] ) :
        for parent, daughter1, daughter2 in parent_plus_decay_pairs:
            print("Cleaning up decay...",parent)
            print("\t\t",daughter1)
            print("\t\t",daughter2)

            parent.image.visible = 0  # can't rely on garbage cleanup to get rid of the last in a list
            self.safe_remove(parent)

            if daughter1.name[0:6] == "Lepton" :
                self.fundamentals.append(daughter1)
                self.fundamentals.append(daughter2)
            elif daughter1.name[0:3] == "Jet" :
                self.virtuals.append(daughter1)
                self.virtuals.append(daughter2)
            else :  # Meson
                self.composites.append(daughter1)
                self.composites.append(daughter2)

    def break_energetic_springs(self):
        for icomp in self.composites :
            for ifund in icomp.billiards[:] :
                if icomp.spring_breaks_and_kinematically_allowed( ifund ) :
                    # snap the spring, giving the spring's potential energy to the quark at its end
                    # and transforming the quark into a new core of a meson

                    print("Because of mediator: ", icomp.dict_constituents_mediators[ifund])
                    print("and billiard: ", ifund.p)

                    # take enough energy from the spring to create the meson, with the same momentum as
                    # the quark that created the meson
                    new_meson = Composite( name = "Meson::Uncollidable" ,
                                           pos = ifund.r3,
                                           energy = sqrt(ifund.p3mag**2 + const.MESON_MASS**2),
                                           direction = ifund.dir
                                         )
                    self.composites.append(new_meson)
                    print("Created: " , new_meson )

                    # update the mediator to preserve energy
                    icomp.dict_constituents_mediators[ifund].E = ( icomp.dict_constituents_mediators[ifund].E
                                                                    + ifund.E - new_meson.E )

                    # delete the progenitor quark from the proton's list of constituents
                    print(len(icomp.constituents),len(icomp.billiards))
                    icomp.safe_remove( ifund )
                    print(len(icomp.constituents), len(icomp.billiards))

                    # update self momentum
                    icomp.recalculate_composite_properties()

    def nearest_composite_to(self, composite):
        new_nearest = composite
        best_dist_sqr = 1000
        for icomp in self.composites :
            test_dist_sqr = (icomp.r3-composite.r3).mag2
            if 0.01 < test_dist_sqr < best_dist_sqr :
                new_nearest = icomp
                best_dist_sqr = test_dist_sqr
        return new_nearest

    def return_composites_to_mass_shell(self):

        """
        For composites below their mass-shell:
        -
        The idea is to subtract energy-momentum δ from the off-shell composite P*, which gets added to
        the nearest neighbour K of the off-shell. We want to maintain the mass of the nearest neighbour,
        so K = K' + delta obeys delta^2 = -2 K.δ This means that δ is spacelike. Also means Eδ < |vecδ|
        so we can use R = Eδ/|vecδ| as a perturbative parameter
        -
        If vecδ is orthogonal to vecK, then we get a nice equation Eδ = sqrt(EK^2 + vecδ^2) - EK (exactness required)
        Then since P* - delta = P, P^2 = P*^2 -2 Eδ(E* + EK) + 2 vecP* . vecδ .
        -
        Define pT as vecP* . vecδhat
        Define Delta^2 = P^2 - P*^2
        Define Sigma = E* + EK
         Then |vecδ| = Delta^2/(2pT) * (1- Sigma R / pT )^-1
        -
        Expand perturbatively in R, using also δ^2 = -2 Eδ EK to write R = (1-R^2)|vecδ|/(2EK)
        Find that, up to O(R), |vecδ| = (Delta^2/2pT)[ 1+ Delta^2*Sigma/(4 EK pT^2)]
        """
        """
        For composites above their mass-shell:
        Simply replace the composite with a jet, which has already been handled
        """

        wiggle_room_factor_low = 0.98
        wiggle_room_factor_high = 1.05


        # composites need to end up on-shell at the end of the experiment
        for icomp in self.composites :
            correct_mass = const.MESON_MASS if icomp.name[0:5] == "Meson" else const.PROTON_MASS
            # low mass condition
            if icomp.M < correct_mass*wiggle_room_factor_low:
                icomp.offshell_counter += 1
                if icomp.offshell_counter > 100 :

                    print("\nStealing energy for ",icomp,"since", icomp.M ,"<",correct_mass*wiggle_room_factor_low)

                    nearest = self.nearest_composite_to(icomp)
                    # go to CoM frame of (icomp and nearest)
                    transverse_dir = Vec3.random_dir_2D_around( nearest.p3 )
                    while transverse_dir*icomp.p3 < 0 :
                        transverse_dir = Vec3.random_dir_2D_around(nearest.p3)

                    Delta = correct_mass**2 - icomp.M**2
                    pT = icomp.p3*transverse_dir
                    # this is not an exact method but under repeated iteration it will move towards onshellness
                    mag_transverse = ( Delta**2/(2*pT) )*( 1 + Delta**2*(icomp.E + nearest.E)/(4*pT**2*nearest.E))
                    delta = Mom4( sqrt(nearest.E**2 + mag_transverse**2) - nearest.E,
                                  transverse_dir.x*mag_transverse,
                                  transverse_dir.y * mag_transverse,
                                  transverse_dir.z * mag_transverse,
                                  )

                    if (icomp.p - delta)*(icomp.p - delta) < 0 :
                        print(f"{(icomp.p - delta)*(icomp.p - delta)} < 0 : ABORTING")
                        continue

                    print(f"Giving 4-mom {delta} to {nearest}")
                    nearest.E += delta.E
                    print(f"so now {nearest}")

                    print(f"Subtracting 4-mom {delta} from {icomp}'s {icomp.constituents[-1]}")
                    icomp.constituents[-1].p = icomp.constituents[-1].p - delta

                    icomp.recalculate_composite_properties()
                    nearest.recalculate_composite_properties()

                    # energy_to_steal = sqrt( correct_mass**2 + icomp.p3.mag2 ) - icomp.E
                    # need to steal energy from a nearby composite to boost up the mass
                    # picture a pion exchanging momentum
                    # should be a while(notpossible) loop here

                    # icomp.E += energy_to_steal
                    # nearest.E -= energy_to_steal
                    # print("Stole",energy_to_steal,"so now ",icomp,nearest)
            # awkward stage where a meson above its mass-shell but it can't decay to another meson
            elif icomp.name[0:5] == "Meson" and wiggle_room_factor_high*const.MESON_MASS < icomp.M < 2*const.MESON_MASS :
                print("\nHARD PART: Losing energy since", icomp.M ,">",wiggle_room_factor_high*const.MESON_MASS)
                energy_to_lose = icomp.E - sqrt( correct_mass**2 + icomp.p3.mag2 )
                # need to steal energy from a nearby composite to boost up the mass
                # picture a pion exchanging momentum
                # should be a while(notpossible) loop here
                nearest = self.nearest_composite_to(icomp)
                icomp.E -= energy_to_lose
                nearest.E += energy_to_lose
                print("Lost",energy_to_lose,"so now ",icomp,nearest)
            elif icomp.M > correct_mass*wiggle_room_factor_high :
                # need to begin decaying to lower-mass hadronic states
                icomp.offshell_counter += 1
                if icomp.offshell_counter > 200 :  # time or position would be better
                    print(f"Offshell {icomp} so now decaying")
                    print(f"{len(scene.objects)} and {scene.objects}")
                    # generate a hadron jet from the off-shell composite
                    new_jet = Fundamental( name="Jet",
                                           pos = icomp.r3,
                                           energy = icomp.E,
                                           direction = icomp.dir,
                                           mass = icomp.M
                                         )
                    self.virtuals.append( new_jet )
                    self.safe_remove( icomp )
            # end if
        # end for

    def potential_scatterers_of_fundamental(self,ifund):
        others = []
        for iother in self.fundamentals:
            if iother != ifund:
                others.append(iother)
        for iother in self.collidable_composites:
            others.extend(iother.constituents)
        return others

    def potential_scatterers_of_composite(self, icomp):
        others = []

        for iother in self.fundamentals:
            others.append(iother)
        for iother in self.collidable_composites :
            if iother != icomp:  # only consider inter-composite collisions
                try:
                    others.extend(iother.constituents)
                except AttributeError:
                    print(iother,"is not a composite. How did it get in that list?")
        return others

    def simulate_for(self, seconds=10):

        frames_per_second = 30
        dt_of_frame = 0.01

        for i in range(100) :
            self.composites[0].update_comp(dt=dt_of_frame)
            self.composites[1].update_comp(dt=dt_of_frame)


        for iframe in range(0, floor(frames_per_second * seconds)):

            rate(frames_per_second)

            annihilation_pairs = []
            billiards_pairs = []
            parent_plus_decay_pairs = []

            for ifund in self.fundamentals:
                # fundamentals are free particles so their phase-space forces are simply those which conserve momentum
                phasespace_forces = ( ifund.v3, Vec3() )
                # Forward timestep
                ifund.update(dt=dt_of_frame, ps_forces = phasespace_forces)
                # Need to know the positions of other fundamentals in order to apply collision logic
                others = self.potential_scatterers_of_fundamental(ifund)
                new_annihilation_pair, new_billiards_pair = ifund.get_fund_collision_pairs(other_fundamentals = others)
                annihilation_pairs.extend( new_annihilation_pair )
                billiards_pairs.extend( new_billiards_pair)

            for icomp in self.composites :
                # Forward timestep
                icomp.update_comp(dt=dt_of_frame)
                # Collision detection and extraction
                others = self.potential_scatterers_of_composite(icomp)
                new_annihilation_pair, new_billiards_pair =  icomp.get_collision_pairs(other_fundamentals = others)
                annihilation_pairs.extend( new_annihilation_pair )
                billiards_pairs.extend( new_billiards_pair)

            for ivirt in self.virtuals:
                # Virtuals are free particles so their phase-space forces are simply those which conserve momentum
                phasespace_forces = ( ivirt.v3, Vec3() )
                # Forward timestep
                ivirt.update(dt=dt_of_frame, ps_forces = phasespace_forces)
                # Decay detection and extraction
                new_parent_plus_decay_pair = ivirt.get_decay_pairs(dt=dt_of_frame)
                parent_plus_decay_pairs.extend( new_parent_plus_decay_pair )

            self.handle_annihilations( annihilation_pairs )
            self.handle_billiards_collisions( billiards_pairs )
            self.handle_decays( parent_plus_decay_pairs )

            # Billiards should no longer be part of a composite if binding energy exceeds some threshold
            # In such cases, a new composite is created from the spring's potential energy
            self.break_energetic_springs()
            self.return_composites_to_mass_shell()


