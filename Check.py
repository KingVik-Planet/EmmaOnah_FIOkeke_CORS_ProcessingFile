try:
    from pykrige.ok import OrdinaryKriging
    print("pykrige library is installed and accessible.")
except ImportError:
    print("pykrige library is not installed or not accessible.")
