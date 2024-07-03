import os
import argparse
from astropy.io import fits

def remove_astrometry_data(fits_file):
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

def main():
    parser = argparse.ArgumentParser(description='Remove astrometry data from FITS file headers.')
    parser.add_argument('fits_files', metavar='FITS_FILES', type=str, nargs='+', help='Paths to the FITS files to process.')
    args = parser.parse_args()

    for fits_file in args.fits_files:
        if os.path.isfile(fits_file) and fits_file.endswith(".fits"):
            remove_astrometry_data(fits_file)
        else:
            print(f"Skipping {fits_file}: Not a valid FITS file")

if __name__ == "__main__":
    main()