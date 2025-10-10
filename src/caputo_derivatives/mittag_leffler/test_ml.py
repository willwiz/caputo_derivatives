import numpy as np
from scipy.special import erfc

from .mittag_leffler import ml

z = np.linspace(-2.0, 2.0, 50)
assert np.allclose(ml(z, 1.0), np.exp(z))

z = np.linspace(-2.0, 2.0, 50)
assert np.allclose(ml(z**2, 2.0), np.cosh(z))

z = np.linspace(0.0, 2.0, 50)
assert np.allclose(ml(np.sqrt(z), 0.5), np.exp(z) * erfc(-np.sqrt(z)))
