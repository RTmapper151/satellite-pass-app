# tle_downloader.py
from skyfield.api import load
import os

def download_tle_files(save_folder, groups, max_days=1.0):
    """
    Download TLE files for specified groups into save_folder.
    Skips downloading if file is fresh (younger than max_days).
    """
    os.makedirs(save_folder, exist_ok=True)
    base_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP={group}&FORMAT=TLE'

    for group in groups:
        filename = os.path.join(save_folder, f'{group}.tle')
        url = base_url.format(group=group)

        if not load.exists(filename) or load.days_old(filename) >= max_days:
            print(f"Downloading fresh TLE data for group: {group}")
            load.download(url, filename=filename)
        else:
            print(f"TLE data for group '{group}' is fresh. Skipping download.")

    print(f"\nTLE files saved in folder: {save_folder}")
