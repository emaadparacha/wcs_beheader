#!/work/mccleary_group/vassilakis.g/miniconda3/envs/astrometry-net/bin/python

import os
import shutil
import argparse
import subprocess
from astropy.io import fits

# Set up argument parser
parser = argparse.ArgumentParser(description='Solve field and update FITS file headers.')
parser.add_argument('files', metavar='F', type=str, nargs='+',
                    help='A list of files to solve field and update headers for')

# Parse command line arguments
args = parser.parse_args()

# Set the directory paths
old_dir = os.getcwd()
new_dir = os.path.join(old_dir, 'astrometry-out')

# Ensure new directory exists
if not os.path.exists(new_dir):
    os.makedirs(new_dir)

# Construct the solve-field command
# Construct the new solve-field command
solve_field_cmd = [
    'solve-field', '--scale-low', '0.1', '--scale-high', '180.0',
    '--scale-units', 'degwidth', '--downsample', '2',
    '--objs', '1000', '--tweak-order', '4', '--overwrite',
    '-D', new_dir
] + args.files

print(f"Running command: {' '.join(solve_field_cmd)}")

# Run solve-field command on all files at once
print(f"Running solve-field on all files: {args.files}")
subprocess.run(solve_field_cmd)

# Iterate over files to update the WCS information
for file in args.files:
    # construct the file paths
    old_file = os.path.join(old_dir, file)
    new_file = os.path.join(new_dir, file.replace('.fits', '.new'))

    # Check if the new file exists before proceeding
    if not os.path.exists(new_file):
        print(f"New FITS file {new_file} does not exist. Skipping header update for {file}.")
        continue

    # open the old fits file
    with fits.open(old_file) as hdul_old:
        # open the new fits file
        with fits.open(new_file) as hdul_new:
            # replace the header of the old file with that of the new file
            hdul_old[0].header = hdul_new[0].header
            # write the changes to a temporary fits file
            hdul_old.writeto(old_file + ".temp", overwrite=True)

    # remove the old fits file
    os.remove(old_file)
    # rename the temporary fits file to the old fits file name
    shutil.move(old_file + ".temp", old_file)

    print(f"Header updated successfully for {file}")

print("All operations completed.")