# This script will retrieve the Landolt (1992) stars and corresponding APASS
# data from Vizier (or some other online service).

import numpy as np
import warnings
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.table import Table

#  Define the name of the output Landolt/APASS catalog
outputFile = 'landoltStars.csv'

# Reset ROW_LIMIT property to retrieve FULL catalog
Vizier.ROW_LIMIT = -1

###########
# It should be possible to significantly increase the sample by using the more
# recent data at +50 deg. declination found in "2013AJ....146..131L"
###########
landoltStars = (Vizier.get_catalogs('J/AJ/137/4186'))[0]

outputNames  = ['_RAJ2000',
                '_DEJ2000',
                'Name',
                'Vmag_landolt',
                'e_Vmag_landolt',
                'B-V_landolt',
                'e_B-V_landolt',
                'U-B_landolt',
                'e_U-B_landolt',
                'V-R_landolt',
                'e_V-R_landolt',
                'R-I_landolt',
                'e_R-I_landolt',
                'V-I_landolt',
                'e_V-I_landolt',
                'B-V_apass',
                'e_B-V_apass',
                'Bmag_apass',
                'e_Bmag_apass',
                'Vmag_apass',
                'e_Vmag_apass',
                'g_mag_apass',
                'e_g_mag_apass',
                'r_mag_apass',
                'e_r_mag_apass',
                'i_mag_apass',
                'e_i_mag_apass']

outputDtypes = ['<f8',
                '<f8',
                'S11',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4',
                '<f4']

# Generate the initial blank table
outputTable = Table(masked = True, names = outputNames, dtype = outputDtypes)

# Loop through each Landolt star and retrieve the APASS data for that object.
for iStar, star in enumerate(landoltStars):
    # Parse the star's name in SimbadName
    simbadName = star['SimbadName'].decode('utf-8')

    # Parse this star's pointing
    RA, Dec  = star['_RAJ2000'], star['_DEJ2000']

    # # Skip over Landolt stars more than 10 degrees from the equator.
    # if np.abs(Dec) > 20.0:
    #     print('Star {0} is far from the equator... skipping'.format(simbadName))
    #     continue

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Query the APASS data for this object
        APASS = Vizier.query_object(simbadName,
            radius = 2.5*u.arcsec, catalog='APASS')

    if len(APASS) > 0:
        # Some match was found, so grab that 'catalog'
        APASS = APASS[0]
    else:
        print('Star {0} not matched in APASS... skipping'.format(simbadName))
        continue

    # If there is more than one object returned,
    if len(APASS) > 1:
        # Then find the object closest to the query point and just use that.
        matchInd = np.where(APASS['_r'].data == APASS['_r'].data.min())
        APASS    = APASS[matchInd]

    # Test if we have some of the expected data
    BVbool = (not APASS['B-V'].mask[0])   and (not APASS['e_B-V'].mask[0])
    Bbool  = (not APASS['Bmag'].mask[0])  and (not APASS['e_Bmag'].mask[0])
    Vbool  = (not APASS['Vmag'].mask[0])  and (not APASS['e_Vmag'].mask[0])
    g_bool = (not APASS['g_mag'].mask[0]) and (not APASS['e_g_mag'].mask[0])
    r_bool = (not APASS['r_mag'].mask[0]) and (not APASS['e_r_mag'].mask[0])
    i_bool = (not APASS['i_mag'].mask[0]) and (not APASS['e_i_mag'].mask[0])

    # Test if more than one band was provided
    sourceBool = np.sum([
        BVbool,
        Bbool,
        Vbool,
        g_bool,
        r_bool,
        i_bool
    ]) >= 2

    if sourceBool == True:
        print('Parsing star {0}'.format(simbadName), end=' ')

        # If there is at least some source information, then add a new row
        outputTable.add_row([0]*len(outputNames), mask = [True]*len(outputNames))

        # Determine the row index of the added null row
        iRow = len(outputTable) - 1

        # Copy the star data into the outputTable
        for col in star.columns:
            if col in outputNames:
                outputTable[col][iRow] = star[col]
            elif col + '_landolt' in outputNames:
                outputTable[col + '_landolt'][iRow] = star[col]

        # Extract the sources for the B1/2 and R1/2 magnitudes
        BVmag   = APASS['B-V'].data.data[0]
        e_BVmag = APASS['e_B-V'].data.data[0]
        Vmag    = APASS['Vmag'].data.data[0]
        e_Vmag  = APASS['e_Vmag'].data.data[0]
        Bmag    = APASS['Bmag'].data.data[0]
        e_Bmag  = APASS['e_Bmag'].data.data[0]
        g_mag   = APASS['g_mag'].data.data[0]
        e_g_mag = APASS['e_g_mag'].data.data[0]
        r_mag   = APASS['r_mag'].data.data[0]
        e_r_mag = APASS['e_r_mag'].data.data[0]
        i_mag   = APASS['i_mag'].data.data[0]
        e_i_mag = APASS['e_i_mag'].data.data[0]

        # Treat each magnitude separately
        ############################ BV magnitudes ############################
        if BVbool == True:
            outputTable['B-V_apass'][iRow]   = BVmag
            outputTable['e_B-V_apass'][iRow] = e_BVmag

        ############################ V magnitudes #############################
        if Vbool == True:
            outputTable['Vmag_apass'][iRow]   = Vmag
            outputTable['e_Vmag_apass'][iRow] = e_Vmag

        ############################ V magnitudes #############################
        if Bbool == True:
            outputTable['Bmag_apass'][iRow]   = Bmag
            outputTable['e_Bmag_apass'][iRow] = e_Bmag

        ############################ V magnitudes #############################
        if g_bool == True:
            outputTable['g_mag_apass'][iRow]   = g_mag
            outputTable['e_g_mag_apass'][iRow] = e_g_mag

        ############################ V magnitudes #############################
        if r_bool == True:
            outputTable['r_mag_apass'][iRow]   = r_mag
            outputTable['e_r_mag_apass'][iRow] = e_r_mag

        ############################ V magnitudes #############################
        if i_bool == True:
            outputTable['i_mag_apass'][iRow]   = i_mag
            outputTable['e_i_mag_apass'][iRow] = e_i_mag

# # Now that the entire Landolt catalog has been parsed, save it to disk
# percentageKept = float(len(outputTable))/float(len(landoltStars))
# print('\n{0:4.4g} percent of Landolt Stars parsed'.format(percentageKept))

# Perform the save
outputTable.write(outputFile, format='ascii.csv', overwrite=True)

print('Done!')
