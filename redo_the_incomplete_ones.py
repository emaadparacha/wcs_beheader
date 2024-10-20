#!/usr/bin/env python3

# ------------ Imports ------------

import argparse
import os
from astropy.io import fits
import shutil
from concurrent.futures import ThreadPoolExecutor
import sip_tpv

# ------------ Functions ------------

# Beheader function
def remove_astrometry_data(fits_file):

    # Check if the file exists and if it is a FITS file, if not, skip
    if not os.path.isfile(fits_file) or not fits_file.endswith(".fits"):
        print(f"Skipping {fits_file}: Not a valid FITS file")
        return
    
    # Check if the file is corrupted
    try:
        with fits.open(fits_file) as hdul:
            pass
    except:
        print(f"Skipping {fits_file}: File is corrupted")
        return

    # Open the FITS file
    with fits.open(fits_file) as hdul:
        header = hdul[0].header
        data = hdul[0].data

        # Find the index of the line 'COMMENT Original key: "END"'
        end_index = None
        for i, card in enumerate(header.cards):
            if card[0] == 'COMMENT' and 'Original key: "END"' in card[1]:
                end_index = i
                break
        
        # If found, truncate the header after this line
        if end_index is not None:
            new_header = fits.Header(cards=header.cards[:end_index + 1])
            hdul_new = fits.HDUList([fits.PrimaryHDU(data, header=new_header)])
            hdul_new.writeto(fits_file, overwrite=True)
            print(f"Beheaded {fits_file}")

    # Return the file name
    print(f"Beheaded {fits_file}")

# SIP to TPV function
def convert_sip_to_tpv(fits_filename):
    # Open the FITS file
    with fits.open(fits_filename) as hdul:

        # Check if B_ORDER is present in the header, if not, skip
        if 'B_ORDER' not in hdul[0].header:
            print(f"Skipping {fits_filename}: B_ORDER not present in header")
            return

        # Modify the header in-place
        sip_tpv.sip_to_pv(hdul[0].header)
        
        # Save the modified file, overwriting the original file
        hdul.writeto(fits_filename, overwrite=True)
        print(f"Converted {fits_filename}")




# ------------ Main Code ------------

# Original science images directory on Hen
original_dir = "/data/flights/superbit_2023/raw_images/science_images"

# Astrometry out scratch directory on Hen
astrometry_out_dir = "/data/analysis/superbit_2023/scratch_analysis/astrometry-out"

# Text file with the incomplete files on Hen
files_not_complete_text = "/home/bit/git/wcs_beheader/files_not_complete.txt"

# Set up argument parser (take in a list of files or a directory)
parser = argparse.ArgumentParser(description='Run beheader, update WCS, and convert SIP to TPV on FITS files.')

# Set it up to take either a list of files of files or a directory (one or the other is required, not both)
parser.add_argument('--num_threads', type=int, default=30, help='Number of threads to use for solve-field command')
parser.add_argument('--mode', type=str, default='all', help='Mode to run the script in. Options are: all, behead, wcs, siptotpv')

# Parse command line arguments
args = parser.parse_args()

# Check if the mode is valid
if args.mode not in ['all', 'behead', 'wcs', 'siptotpv']:
    print("Invalid mode provided. Options are: all, behead, wcs, siptotpv")
    exit()

# Get the list of files in the directory if a directory is provided
all_files = []

# Read in the text file with the list of files
with open(files_not_complete_text, 'r') as f:
    all_files = f.read().splitlines()

# Remove any empty strings from the list or any duplicates or anything not ending in .fits
all_files = [f for f in all_files if f != '' and f.endswith('.fits')]

# Check if there are any files to process
if len(all_files) == 0:
    print("No FITS files found or in the provided list.")
    exit()

# Sort the files
all_files.sort()

# Print number of files
print(f"Number of files to process: {len(all_files)}")

# ### Beheader ###

if args.mode == 'behead' or args.mode == 'all':

    status_counter = 0

    # Number of threads to use
    num_threads = args.num_threads

    # Parallel processing of the beheader function
    print("Running beheader on all files.")
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Run the beheader commands
        executor.map(remove_astrometry_data, all_files)

# ### Update WCS ###

if args.mode == 'wcs' or args.mode == 'all':

    # Set the directory paths
    old_dir = os.getcwd()
    wcs_dir = astrometry_out_dir # Moving to /data because of space issues

    # Ensure new directory exists
    if not os.path.exists(wcs_dir):
        os.makedirs(wcs_dir)

    # Create a list to store the solve-field commands
    solve_field_cmds = []

    # For each file in the list, go to the original directory and check for the same image there, open it and get the "TRG_RA" and "TRG_DEC" values to input into the solve-field command for each image
    # Temporary fix - since the original dir does not have images ending in {}_clean.fits, remove the _clean.fits from the file name to get the original image name
    for file in all_files:
        # Get the original file name
        original_file = file.replace('_clean.fits', '.fits') # Temporary fix

        # Remove everything before the last / to get the file name
        original_file = original_file.split('/')[-1] # Temporary fix

        # If the original file does not exist, skip
        if not os.path.exists(os.path.join(original_dir, original_file)):
            print(f"Original file {original_file} does not exist. Skipping solve-field for {file}.")
            continue

        # If the original file is truncated or corrupted, skip
        try:
            with fits.open(os.path.join(original_dir, original_file)) as hdul:
                print(f"Checking original file {original_file}")
                pass
        except:
            print(f"Original file {original_file} is corrupted. Skipping solve-field for {file}.")
            continue

        # Do the same for the fits file
        try:
            with fits.open(file) as hdul:
                print(f"Checking file {file}")
                pass
        except:
            print(f"File {file} is corrupted. Skipping solve-field for {file}.")
            continue

        # Open the original file
        with fits.open(os.path.join(original_dir, original_file)) as hdul:

            # Print full path name of the original file
            print(f"Original file: {os.path.join(original_dir, original_file)}")

            # Get the target RA and DEC
            target_ra = hdul[0].header['TRG_RA']
            target_dec = hdul[0].header['TRG_DEC']

        # Construct the solve-field command
        solve_field_cmd = ['solve-field', '--scale-units', 'arcsecperpix', '--scale-low', '0.135', '--scale-high', '0.145', 
                        '--ra', str(target_ra), '--dec', str(target_dec), '--radius', '10', '--downsample', '6', '--objs', '1000', 
                        '--tweak-order', '4', '--cpulimit', '75', '--overwrite', '-D', wcs_dir, file]
        solve_field_cmd = ' '.join(solve_field_cmd)

        # Append the command to the list
        solve_field_cmds.append(solve_field_cmd)

    #print("Solve-field commands created - here they are:")
    #print(solve_field_cmds)

    # Number of threads to use
    num_threads = args.num_threads

    # Run solve-field command on all files using multiple threads
    print(f"Running solve-field on all files using {num_threads} threads.")
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Run the solve-field commands
        executor.map(os.system, solve_field_cmds)

    # Status print
    print("Solve-field completed.")

    # Files not complete list
    files_not_complete = []

    # Iterate over files to update the WCS information
    for file in all_files:
        # construct the file paths
        fits_base = os.path.basename(file)
        fits_new_base = os.path.splitext(fits_base)[0] + '.new'
        new_wcs_file = os.path.join(wcs_dir, fits_new_base)

        # Check if the new file exists before proceeding
        if not os.path.exists(new_wcs_file):

            # Add the file to the list of files not complete
            files_not_complete.append(file)

            print(f"New FITS file {new_wcs_file} does not exist. Skipping header update for {file}.")
            continue
        else:
            # open the old fits file
            with fits.open(file) as hdul_old:
                # open the new fits file
                with fits.open(new_wcs_file) as hdul_new:
                    # replace the header of the old file with that of the new file
                    hdul_old[0].header = hdul_new[0].header
                    # write the changes to a temporary fits file
                    hdul_old.writeto(file + ".temp", overwrite=True)

            # remove the old fits file
            os.remove(file)
            # rename the temporary fits file to the old fits file name
            shutil.move(file + ".temp", file)

        print(f"Header updated successfully for {file}")

    # Status print
    print("Header update completed with downsampling of 6. Number of files not complete:", len(files_not_complete))


    # Redo the solve-field command for the files that were not complete but with downsample 8
    print("Redoing solve-field for files that were not complete with downsample 8.")

    # Create a list to store the solve-field commands
    solve_field_cmds = []

    # For each file in the list, go to the original directory and check for the same image there, open it and get the "TRG_RA" and "TRG_DEC" values to input into the solve-field command for each image
    # Temporary fix - since the original dir does not have images ending in {}_clean.fits, remove the _clean.fits from the file name to get the original image name
    for file in files_not_complete:
        # Get the original file name
        original_file = file.replace('_clean.fits', '.fits') # Temporary fix

        # Remove everything before the last / to get the file name
        original_file = original_file.split('/')[-1] # Temporary fix

        # Open the original file
        with fits.open(os.path.join(original_dir, original_file)) as hdul:

            # Print full path name of the original file
            print(f"Original file: {os.path.join(original_dir, original_file)}")

            # Get the target RA and DEC
            target_ra = hdul[0].header['TRG_RA']
            target_dec = hdul[0].header['TRG_DEC']

        # Construct the solve-field command
        solve_field_cmd = ['solve-field', '--scale-units', 'arcsecperpix', '--scale-low', '0.135', '--scale-high', '0.145', 
                        '--ra', str(target_ra), '--dec', str(target_dec), '--radius', '10', '--downsample', '8', '--objs', '1000', 
                        '--tweak-order', '4', '--cpulimit', '75', '--overwrite', '-D', wcs_dir, file]
        solve_field_cmd = ' '.join(solve_field_cmd)

        # Append the command to the list
        solve_field_cmds.append(solve_field_cmd)

    #print("Solve-field commands created - here they are:")
    #print(solve_field_cmds)

    # Number of threads to use
    num_threads = args.num_threads

    # Run solve-field command on all files using multiple threads
    print(f"Running solve-field on all files using {num_threads} threads.")

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Run the solve-field commands
        executor.map(os.system, solve_field_cmds)

    # Status print
    print("Solve-field completed.")

    # Really not complete list
    files_not_complete_round_2 = []

    # Iterate over files to update the WCS information
    for file in files_not_complete:
        # construct the file paths
        fits_base = os.path.basename(file)
        fits_new_base = os.path.splitext(fits_base)[0] + '.new'
        new_wcs_file = os.path.join(wcs_dir, fits_new_base)

        # Check if the new file exists before proceeding
        if not os.path.exists(new_wcs_file):

            # Add the file to the list of files not complete
            files_not_complete_round_2.append(file)

            print(f"New FITS file {new_wcs_file} does not exist. Skipping header update for {file}.")
            continue
        else:
            # open the old fits file
            with fits.open(file) as hdul_old:
                # open the new fits file
                with fits.open(new_wcs_file) as hdul_new:
                    # replace the header of the old file with that of the new file
                    hdul_old[0].header = hdul_new[0].header
                    # write the changes to a temporary fits file
                    hdul_old.writeto(file + ".temp", overwrite=True)

            # remove the old fits file
            os.remove(file)
            # rename the temporary fits file to the old fits file name
            shutil.move(file + ".temp", file)

        print(f"Header updated successfully for {file}")

    # Status print
    print("Header update completed, for the second run. Number of files not complete:", len(files_not_complete_round_2))

    # Write the files not complete to a text file (overwriting the previous file)
    with open('files_not_complete_downsample_8.txt', 'w') as f:
        for file in files_not_complete_round_2:
            f.write(file + '\n')

    print("All operations completed.")

# ### SIP to TPV ###

if args.mode == 'siptotpv' or args.mode == 'all':

    # Number of threads to use
    num_threads = args.num_threads

    # List to store the files not complete
    files_not_complete_round_2 = []
    status_counter = 0

    # Read the files_not_complete_round_2 file if it exists and store the files in a list
    if os.path.exists('files_not_complete_downsample_8.txt'):
        with open('files_not_complete_downsample_8.txt', 'r') as f:
            files_not_complete_round_2 = f.read().splitlines()

    # Error check and print the number of files not complete
    if len(files_not_complete_round_2) == 0:
        print("No files not complete found.")
    else:
        print(f"Number of files not complete: {len(files_not_complete_round_2)}")

    # Remove the files_not_complete_round_2 from the all_files list
    all_files = [file for file in all_files if file not in files_not_complete_round_2]

    # Convert SIP to TPV
    print("Converting SIP to TPV.")
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Run the SIP to TPV commands
        executor.map(convert_sip_to_tpv, all_files)

    print("Files converted")

print("All operations completed.")

