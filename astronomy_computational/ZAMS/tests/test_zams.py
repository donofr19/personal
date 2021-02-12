import pytest
from zams import Teff, surface_luminosity
from astro_const import Rsun, Lsun
from numcheck import within_tolerance

# delete this line when your Teff interpolation routine is ready for testing
@pytest.mark.skip(reason="Teff routine not ready")
def test_Teff():
    # does interpolation respect table values?
    for mass,value in zip([0.1,0.15,0.2,0.3],[2800.0,3150.0,3300.0,3400.0]):
        assert Teff(mass) == value
    
    # is interpolation monotonic?
    assert Teff(0.17) > Teff(0.12)
    assert Teff(0.17) > Teff(0.16)
    assert Teff(0.25) > Teff(0.18)
    assert Teff(0.28) > Teff(0.26)

# delete this line when your surface_luminosity routine is ready for testing
@pytest.mark.skip(reason="test_surface_luminosity routine not ready")
def test_surface_luminosity():
    Teff = 5772.0
    radius = Rsun
    assert within_tolerance(surface_luminosity(Teff,radius),Lsun,tol=1.0e-3)
