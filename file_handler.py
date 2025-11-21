from astropy.io import fits
import numpy as np
import glob
from astropy.time import Time
import astropy.units as u
from astropy.time import Time
from datetime import datetime, timedelta

def BeSSSpectra(file):
    """
    Load a 1-D BeSS spectrum from FITS.

    Returns:
        wavelengths: np.ndarray (linear solution from header)
        spectrum: np.ndarray
        header: FITS header
    """

    try:
        f = fits.open(file)
        header = f[0].header
        spectrum = f[0].data
        f.close()

        # Linear wavelength solution from header
        crval1 = header.get('CRVAL1', 6563.0)  # default near Hα
        cdelt1 = header.get('CDELT1', 1.0)
        crpix1 = header.get('CRPIX1', 1.0)

        wavelengths = crval1 + (np.arange(len(spectrum)) + 1 - crpix1) * cdelt1

        return wavelengths, spectrum, header

    except Exception as e:
        print(f"⚠️ Failed to read {file}: {e}")
        return None, None, None

# ----------------------------
# sort_files: wavelength filter + DATE-OBS from header
# ----------------------------
def sort_files(file_paths, spec_feature=6562.8):
    """
    Collect FITS files covering Hα, store observation dates from header.

    Returns:
        valid_files: list of paths
        valid_times: list of astropy Time objects (from DATE-OBS)
    """
    files = glob.glob(file_paths)
    valid_files = []
    valid_times = []

    print(f"Found {len(files)} files in directory.")

    for file in files:
        try:
            wavelengths, spectrum, header = BeSSSpectra(file)
            lam_min = np.min(wavelengths)
            lam_max = np.max(wavelengths)
            if wavelengths is None:
                continue

            # Wavelength coverage filter
            if lam_min > spec_feature or lam_max < spec_feature or np.abs(lam_min - spec_feature) > 1500\
                or np.abs(lam_max - spec_feature) > 1500:
                continue

            # Get observation date from header
            date_obs = header.get('DATE-OBS', None)
            if date_obs is None:
                print(f"⚠️ {file} has no DATE-OBS")
                continue

            try:
                date_obj = Time(date_obs, format='isot', scale='utc')
            except Exception:
                # Fallback for non-standard DATE-OBS
                date_obj = Time(datetime.strptime(date_obs[:10], "%Y-%m-%d"))

            valid_files.append(file)
            valid_times.append(date_obj)

        except Exception as e:
            print(f"⚠️ Skipping {file} due to error: {e}")
            continue

    # Sort by date
    if valid_files:
        sorted_pairs = sorted(zip(valid_times, valid_files))
        valid_times, valid_files = zip(*sorted_pairs)
    else:
        valid_times, valid_files = [], []

    print(f"{len(valid_files)} files with spectra on {spec_feature} Å.")

    return list(valid_files), list(valid_times)