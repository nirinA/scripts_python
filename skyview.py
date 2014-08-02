'''SkyView multi-wavelength FITS image retriever.

available survey are:'''
## available survey on SkyView
all_survey = ['1420MHZ', '408MHZ', '4850MHZ', 'CO2D', 'FIRST', '35MHZ', 'NH', 'NVSS', 'SUMSS', 'VLSS', 'WENSS', '2MASSH', '2MASSJ', '2MASSK', 'COBEAAM', 'COBEZSMA', 'DSS2 IR', 'IRAS 100', 'IRAS 12', 'IRAS 25', 'IRAS 60', 'IRIS 100', 'IRIS 12', 'IRIS 25', 'IRIS 60', 'SFD100M', 'SFDDUST', 'UKIDSSH', 'UKIDSSJ', 'UKIDSSK', 'UKIDSSY', 'WISE3.4', 'WISE4.6', 'WISE12', 'WISE22', 'WMAPK', 'WMAPKA', 'WMAPQ', 'WMAPV', 'WMAPW', 'WMAPILC', 'DSS', 'DSS1 BLUE', 'DSS1 RED', 'DSS2 BLUE', 'DSS2 RED', 'HALPHA', 'MELLINGER-B', 'MELLINGER-G', 'MELLINGER-R', 'NEAT', 'SDSSG', 'SDSSI', 'SDSSR', 'SDSSU', 'SDSSZ', 'SDSSDR7G', 'SDSSDR7I', 'SDSSDR7R', 'SDSSDR7U', 'SDSSDR7Z', 'SHASSA-C', 'SHASSA-CC', 'SHASSA-H', 'SHASSA-SM', 'EUVE 171 A', 'EUVE 405 A', 'EUVE 555 A', 'EUVE 83 A', 'GALEX FAR UV', 'GALEX NEAR UV', 'ROSAT WFC F1', 'ROSAT WFC F2', 'BAT FLUX 14-20', 'BAT SNR 14-20', 'BAT FLUX 20-24', 'BAT SNR 20-24', 'BAT FLUX 24-35', 'BAT SNR 24-35', 'BAT FLUX 35-50', 'BAT SNR 35-50', 'BAT FLUX 50-75', 'BAT SNR 50-75', 'BAT FLUX 75-100', 'BAT SNR 75-100', 'BAT FLUX 100-150', 'BAT SNR 100-150', 'BAT FLUX 150-195', 'BAT SNR 150-195', 'BAT SNR 14-195', 'GRANAT_SIGMA_FLUX', 'GRANAT_SIGMA_SIG', 'HEAO1A', 'HRIINT', 'INTEGRALSPI_GC', 'INTGAL1735E', 'INTGAL1735F', 'INTGAL1735S', 'INTGAL1760E', 'INTGAL1760F', 'INTGAL1760S', 'INTGAL3580E', 'INTGAL3580F', 'INTGAL3580S', 'PSPC1CNT', 'PSPC1EXP', 'PSPC1INT', 'PSPC2CNT', 'PSPC2EXP', 'PSPC2INT', 'PSPC0.6CNT', 'PSPC0.6EXP', 'PSPC0.6INT', 'RASS-CNT HARD', 'RASS-CNT BROAD', 'RASS-CNT SOFT', 'RASS-INT HARD', 'RASS-INT BROAD', 'RASS-INT SOFT', 'RASS BACKGROUND 1', 'RASS BACKGROUND 2', 'RASS BACKGROUND 3', 'RASS BACKGROUND 4', 'RASS BACKGROUND 5', 'RASS BACKGROUND 6', 'RASS BACKGROUND 7', 'RXTE ALLSKY 3-20KEV FLUX', 'RXTE ALLSKY 3-8KEV FLUX', 'RXTE ALLSKY 8-20KEV FLUX', 'RXTE ALLSKY 3-20KEV SIG', 'RXTE ALLSKY 3-8KEV SIG', 'RXTE ALLSKY 8-20KEV SIG', 'COMPTEL', 'EGRETHARD', 'EGRETSOFT,', 'EGRET3D', 'FERMI1', 'FERMI2', 'FERMI3', 'FERMI4', 'FERMI5', 'PLANCK 217', 'PLANCK 545', 'PLANCK 044', 'PLANCK 143', 'PLANCK 070', 'PLANCK 030', 'PLANCK 100', 'PLANCK 857', 'PLANCK 353']

import os
import sys
import argparse
import urllib.request
import urllib.parse

SkyView = 'http://skyview.gsfc.nasa.gov/cgi-bin/pskcall?%s'

def get_skyview(ra, dec, survey=None, output=None,
                pixels=None,
                size=None):
    ra = str(ra)
    dec = str(dec)
    Position = '%s, %s'%(ra, dec)
    if survey is None:
        survey = 'DSS'
    survey = survey.upper()
    if survey not in all_survey:
        raise NameError('%s not in SkyView survey'%survey)
    if output is None:
        output = 'skyview_%s_%s.fits'%(ra, dec)
    if os.path.isfile(output):
        raise FileExistsError('%s already exists!'%output)
    if pixels is None:
        pixels = 300
    if size is None:
        size = 0.2798933387264
    params = {'Position':Position,
              'Survey':survey,
              'Size':size,
              }
    p_encode = urllib.parse.urlencode(params)
    print('retrieving SkyView: Survey: %s, at Position: %s ...'%(survey, Position), end= '')
    with open(output, 'wb') as s:
        f = urllib.request.urlopen(SkyView%p_encode)
        s.write(f.read())
    print('done')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "%s \n%s"%(__doc__, '\n'.join(all_survey)),
        epilog = 'For more detail on SkyView, see: http://skyview.gsfc.nasa.gov')

    parser.add_argument(
        '-r', '--ra',
        required=True,
        help = 'rigth ascension in degrees'
        )

    parser.add_argument(
        '-d', '--dec',
        required=True,
        help = 'declination in degrees'
        )

    parser.add_argument(
        '-S', '--survey',
        required=False,
        help = 'name of the survey'
        )

    parser.add_argument(
        '-p', '--pixels',
        required=False,
        help = 'number of pixels'
        )

    parser.add_argument(
        '-s', '--size',
        required=False,
        help = 'imge size'
        )

    parser.add_argument(
        '-o', '--output',
        required=False,
        help = 'filename output'
        )

    res = parser.parse_args(sys.argv[1:])
    get_skyview(res.ra, res.dec,
                survey=res.survey,
                output=res.output,
                size=res.size,
                pixels=res.pixels)
