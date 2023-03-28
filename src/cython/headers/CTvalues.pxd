# File: CTvalues.pxd
# distutils: language = c++
# cython: language_level=3


""" ----------------------------------------------------------------------------
C++ Source Files
---------------------------------------------------------------------------- """

""" ----------------------------------------------------------------------------
End of Source Files
---------------------------------------------------------------------------- """


# ------------------------------------------------------------------------------
# C++ Header files + exported definitions
# ------------------------------------------------------------------------------

cdef extern from "src/cpp/CTvalues_optimization.hpp" namespace "ctv":
  pass

cdef extern from "src/cpp/CTvalues_kde.hpp" namespace "ctv":
  pass
