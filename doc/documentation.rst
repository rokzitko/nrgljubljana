.. _documentation:

Documentation
=============

NRG Ljubljana solves quantum impurity problems such as the single-impurity Anderson model (SIAM), described by the following Hamiltonian:

.. math::

   H  &= \sum_{k\sigma} \epsilon_k c^\dagger_{k\sigma} c_{k\sigma} + \sum_\sigma \epsilon_d d^\dagger_\sigma d_\sigma + U n_\uparrow n_\downarrow \\
      &+ \sum_{k\sigma} \left( V_k c^\dagger_{k\sigma} d_\sigma + \text{H.c.} \right)
              
with :math:`n_\sigma=d^\dagger_\sigma d_\sigma`. Here c are conduction-band operators, while d are impurity operators.

Table of Contents
-----------------

.. toctree::
   :maxdepth: 5

   index
   install
   issues
   changelog
   about

