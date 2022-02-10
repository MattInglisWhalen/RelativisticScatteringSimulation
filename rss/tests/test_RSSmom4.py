
# external libraries
import unittest
import numpy as np

# class to be tested
from rss.src.RSSvectors import Mom4


class Mom4TestCase(unittest.TestCase):

    def test_boosted_to_rest_frame_of(self):

        """
        Test that the boosting function works as advertised
        """
        p4 = Mom4(10,0,0,9.9)  # Alon z-axis
        print(p4[-1])
        p4_in_rest_frame = p4.boosted_to_rest_frame_of( p4 )
        self.assertAlmostEqual(p4_in_rest_frame.E, np.sqrt(10**2 - 9.9**2))
        self.assertAlmostEqual(p4_in_rest_frame.px,0)
        self.assertAlmostEqual(p4_in_rest_frame.py,0)
        self.assertAlmostEqual(p4_in_rest_frame.pz,0)

        p4 = Mom4(10,3,3,3)  # Along arbitrary axis
        p4_in_rest_frame = p4.boosted_to_rest_frame_of( p4 )
        self.assertAlmostEqual(p4_in_rest_frame.E, np.sqrt(10**2 - 3 * 3**2))
        self.assertAlmostEqual(p4_in_rest_frame.px,0)
        self.assertAlmostEqual(p4_in_rest_frame.py,0)
        self.assertAlmostEqual(p4_in_rest_frame.pz,0)

        p4 = Mom4(10,0,0,6)
        k4 = Mom4(10,0,1,6)  # Boosting a different momentum
        k4_in_rest_frame = k4.boosted_to_rest_frame_of( p4 )
        self.assertAlmostEqual(k4_in_rest_frame.E, np.sqrt(10**2 - 6**2))
        self.assertAlmostEqual(k4_in_rest_frame.px,0)
        self.assertAlmostEqual(k4_in_rest_frame.py,1)
        self.assertAlmostEqual(k4_in_rest_frame.pz,0)

    def test_boosted_away_with_momentum(self):
        """
        Test that the boosting function works as advertised
        """
        rest_mom = Mom4(1,0,0,0)
        p4 = Mom4(10,0,0,6)  # Along z-axis
        rest_mom_boosted_away = rest_mom.boosted_away_with_momentum( p4 )
        self.assertAlmostEqual(rest_mom_boosted_away.E, 10/np.sqrt(10**2-6**2) )
        self.assertAlmostEqual(rest_mom_boosted_away.px,0)
        self.assertAlmostEqual(rest_mom_boosted_away.py,0)
        self.assertAlmostEqual(rest_mom_boosted_away.pz,10*(6/10)/np.sqrt(10**2-6**2))

        moving = Mom4(10,0,-7,0)
        push = Mom4(10,0,+7,0)  # Along arbitrary axis
        moving_pushed_opposite = moving.boosted_away_with_momentum( push )
        self.assertAlmostEqual(moving_pushed_opposite.E, np.sqrt(10**2 - 7**2))
        self.assertAlmostEqual(moving_pushed_opposite.px,0)
        self.assertAlmostEqual(moving_pushed_opposite.py,0)
        self.assertAlmostEqual(moving_pushed_opposite.pz,0)

