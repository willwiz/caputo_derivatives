#type: ignore
try:
  from src.py.matlaws._HModels import (
    caputo_initialize,  pyHOG2D, NeoHookean2D, HOGstruc2D, GetScaledParameters_cpp,
    GetScaledParametersScaled_cpp,
    HESimulation, HESimulation_Scaled, CaputoSimulation_C_M, CaputoSimulation_C_M_Scaled,
  )
except ImportError:
  raise ImportError('ERROR<<< The _HModels Module has not been compiled.')