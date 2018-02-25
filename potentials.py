#!/usr/bin/env python3

"""
Plot interaction potentials

   Copyright (C) 2015-2018  Vishnu V. Krishnan : vishnugb@gmail.com

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import matplotlib.pyplot as plt
import numpy as np


def lennard_jones(distance, epsillon, sigma):
    """
    The Lennard Jones or 12-6 potential
    It models the interactions between two neutral,
    spherically symetric 'atoms' that can form neutral dipoles.
    """

    return (4*epsillon) * (((sigma/distance)**12) - ((sigma/distance)**6))


#X = np.delete(np.linspace(0.0, 2.0, num=100, dtype='float64'), 0)
X = np.linspace(0.85, 2.0, num=500, dtype='float64')
Y00 = lennard_jones(X, 1.0, 1.00)
Y01 = lennard_jones(X, 1.5, 0.80)
Y11 = lennard_jones(X, 0.5, 0.88)

plt.plot(X, Y00, label='00')
plt.plot(X, Y01, label='01')
plt.plot(X, Y11, label='11')

plt.title('The Lennard-Jones interaction potential for Kob-Andersen pairs')
plt.xlabel('r')
plt.ylabel('Potential Energy')
plt.legend()
plt.show()
