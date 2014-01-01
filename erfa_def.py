'''This module provides ERFA,\\n\\
the Essential Routine for Fundamental Astronomy,\\n\\
interface to Python\\n\\
'''

from build_extension_module import *

module_name = '_erfa'
module_filename = '_erfamodule.c'
module_doc = __doc__
include_files = ['"structseq.h"',
                 '<math.h>',
                 '<time.h>',
                 '"erfa.h"',
                 '"erfam.h"',
                 ]

module_localprototype = '' #'''/* local prototype */'''

module_localfunction = '''
/* local function */
static PyObject *
_to_py_vector(double v[3])
{
    return Py_BuildValue("(ddd)",
                           v[0], v[1], v[2]
                           );
}

static PyObject *
_to_py_matrix(double m[3][3])
{
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                           m[0][0], m[0][1], m[0][2],
                           m[1][0], m[1][1], m[1][2],
                           m[2][0], m[2][1], m[2][2]
                           );
}

static PyObject *
_to_py_astrom(eraASTROM *a)
{
    PyObject *v = PyStructSequence_New(&AstromType);
    if (v == NULL)
        return NULL;
#define SET(i,val) PyStructSequence_SET_ITEM(v, i, PyFloat_FromDouble((double) val))
    SET(0, a->pmt);
    PyStructSequence_SET_ITEM(v, 1, _to_py_vector(a->eb));
    PyStructSequence_SET_ITEM(v, 2, _to_py_vector(a->eh));
    SET(3, a->em);
    PyStructSequence_SET_ITEM(v, 4, _to_py_vector(a->v));
    SET(5, a->bm1);
    PyStructSequence_SET_ITEM(v, 6, _to_py_matrix(a->bpn));
    SET(7, a->along);
    SET(8, a->phi);
    SET(9, a->xpl);
    SET(10, a->ypl);
    SET(11, a->sphi);
    SET(12, a->cphi);
    SET(13, a->diurab);
    SET(14, a->eral);
    SET(15, a->refa);
    SET(16, a->refb);
    if (PyErr_Occurred()) {
        Py_XDECREF(v);
        return NULL;
    }
#undef SET
    return v;
}
'''
module_struct_seq = [
        ## list of 3-tuple
        ## (name, (struct field, field_doc), description))
        ('ASTROM', (
            ("pmt", "PM time interval (SSB, Julian years)"),
            ("eb", "SSB to observer (vector, au)"),
            ("eh", "Sun to observer (unit vector)"),
            ("em", "distance from Sun to observer (au)"),
            ("v", "barycentric observer velocity (vector, c"),
            ("bm1", "sqrt(1-|v|^2): reciprocal of Lorenz factor"),
            ("bpn", "bias-precession-nutation matrix"),
            ("along", "longitude + s' + dERA(DUT) (radians)"),
            ("phi", "geodetic latitude (radians)"),
            ("xpl", "polar motion xp wrt local meridian (radians)"),
            ("ypl", "polar motion yp wrt local meridian (radians)"),
            ("sphi", "sine of geodetic latitude"),
            ("cphi", "cosine of geodetic latitude"),
            ("diurab", "magnitude of diurnal aberration vector"),
            ("eral", "local Earth rotation angle (radians)"),
            ("refa", "refraction constant A (radians)"),
            ("refb", "refraction constant B (radians)")
             ),
        '''"Star-independent astrometry parameters\\n"
"(Vectors eb, eh, em and v are all with respect to BCRS axes.)"'''
        ),
        ('LDBODY', (
            ("bm", "mass of the body (solar masses)"),
            ("dl", "deflection limiter (radians^2/2)"),
            ("pv", "barycentric PV of the body (au, au/day)")
            ),
        '''"Body parameters for light deflection"'''
        )
    ]

module_add = [
        ## additional module constants 
        ('DPI', 'PyFloat_FromDouble(ERFA_DPI)'),
        ('D2PI', 'PyFloat_FromDouble(ERFA_D2PI)'),
        ('DR2D', 'PyFloat_FromDouble(ERFA_DR2D)'),
        ('DD2R', 'PyFloat_FromDouble(ERFA_DD2R)'),
        ('DR2AS', 'PyFloat_FromDouble(ERFA_DR2AS)'),
        ('DAS2R', 'PyFloat_FromDouble(ERFA_DAS2R)'),
        ('DS2R', 'PyFloat_FromDouble(ERFA_DS2R)'),
        ('TURNAS', 'PyFloat_FromDouble(ERFA_TURNAS)'),
        ('DMAS2R', 'PyFloat_FromDouble(ERFA_DMAS2R)'),
        ('DTY', 'PyFloat_FromDouble(ERFA_DTY)'),
        ('DAYSEC', 'PyFloat_FromDouble(ERFA_DAYSEC)'),
        ('DJY', 'PyFloat_FromDouble(ERFA_DJY)'),
        ('DJC', 'PyFloat_FromDouble(ERFA_DJC)'),
        ('DJM', 'PyFloat_FromDouble(ERFA_DJM)'),
        ('DJ00', 'PyFloat_FromDouble(ERFA_DJ00)'),
        ('DJM0', 'PyFloat_FromDouble(ERFA_DJM0)'),
        ('DJM00', 'PyFloat_FromDouble(ERFA_DJM00)'),
        ('DJM77', 'PyFloat_FromDouble(ERFA_DJM77)'),
        ('TTMTAI', 'PyFloat_FromDouble(ERFA_TTMTAI)'),
        ('DAU', 'PyFloat_FromDouble(ERFA_DAU)'),
        ('CMPS', 'PyFloat_FromDouble(ERFA_CMPS)'),
        ('AULT', 'PyFloat_FromDouble(ERFA_AULT)'),
        ('DC', 'PyFloat_FromDouble(ERFA_DC)'),
        ('ELG', 'PyFloat_FromDouble(ERFA_ELG)'),
        ('ELB', 'PyFloat_FromDouble(ERFA_ELB)'),
        ('TDB0', 'PyFloat_FromDouble(ERFA_TDB0)'),
        ('SRS', 'PyFloat_FromDouble(ERFA_SRS)'),
        ('WGS84', 'PyLong_FromLong(ERFA_WGS84)'),
        ('GRS80', 'PyLong_FromLong(ERFA_GRS80)'),
        ('WGS72', 'PyLong_FromLong(ERFA_WGS72)'),
         ]

module_methods = [
        ## list of 4-tuple (name, arguments, doc string, CPython function body)
        ## describing each method
        ## e.g. :
        ##('name', 'METH_VARARGS',
        ##     '''"\\nname() -> \\n"
        ##"the doc string \\n"
        ##"Given:\\n"
        ##"   x          x\\n"
        ##"Returned:\\n"
        ##"   x       \\n"
        ##"   x        "''',
        ##     '''{
        ##    return;
        ##}
        ##'''),
# astrometry
    ('ab', 'METH_VARARGS',
     '''"\\nab(pnat[3], v[3], s, bm1) -> ppr\\n"
"Apply aberration to transform natural direction into proper direction.\\n"
"Given:\\n"
"    pnat       natural direction to the source (unit vector)\\n"
"    v          observer barycentric velocity in units of c\\n"
"    s          distance between the Sun and the observer (au)\\n"
"    bm1        sqrt(1-|v|^2): reciprocal of Lorenz factor\\n"
"Returned:\\n"
"    ppr        proper direction to source (unit vector)"''',
     '''{
    double pnat[3], v[3], s, bm1, ppr[3];
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)dd",
        &pnat[0],&pnat[1],&pnat[2],&v[0],&v[1],&v[2],&s,&bm1))
        return NULL;
    eraAb(pnat, v, s, bm1, ppr);
    return Py_BuildValue("(ddd)",ppr[0], ppr[1], ppr[2]);
}
'''),
    ('apcg', 'METH_VARARGS',
     '''"\\napcg(date1, date2, ebpv[2][3], ehp[3]) -> astrom\\n"
"For a geocentric observer, prepare star-independent astrometry\\n"
"parameters for transformations between ICRS and GCRS coordinates.\\n"
"The Earth ephemeris is supplied by the caller.\\n"
"\\n"
"The parameters produced by this function are required in the\\n"
"parallax, light deflection and aberration parts of the astrometric\\n"
"transformation chain.\\n"
"Given:\\n"
"    date1      TDB as a 2-part...\\n"
"    date2      ...Julian Date\\n"
"    ebpv       Earth barycentric pos/vel (au, au/day)\\n"
"    ehp        Earth heliocentric position (au)\\n"
"Returned:\\n"
"    astrom     star-independent astrometry parameters"''',
     '''{
    double date1, date2, ebpv[2][3], ehp[3];
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "dd((ddd)(ddd))(ddd)",
                          &date1, &date2,
                          &ebpv[0][0],&ebpv[0][1],&ebpv[0][2],
                          &ebpv[1][0],&ebpv[1][1],&ebpv[1][2],
                          &ehp[0],&ehp[1],&ehp[2]))      
        return NULL;
    eraApcg(date1, date2, ebpv, ehp, &astrom);
    return _to_py_astrom(&astrom);
}
'''),
    ('apcg13', 'METH_VARARGS',
     '''"\\napcg13(date1, date2) -> astrom\\n"
"For a geocentric observer, prepare star-independent astrometry\\n"
"parameters for transformations between ICRS and GCRS coordinates.\\n"
"The caller supplies the date, and ERFA models are used to predict\\n"
"the Earth ephemeris.\\n"
"\\n"
"The parameters produced by this function are required in the\\n"
"parallax, light deflection and aberration parts of the astrometric\\n"
"transformation chain.\\n"
"Given:\\n"
"    date1      TDB as a 2-part...\\n"
"    date2      ...Julian Date\\n"
"Returned:\\n"
"    astrom     star-independent astrometry parameters"''',
     '''{
    double date1, date2;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "dd", &date1, &date2))      
        return NULL;
    eraApcg13(date1, date2, &astrom);
    return _to_py_astrom(&astrom);
}
'''),
    ('apci', 'METH_VARARGS',
     '''"\\napci(date1, date2, ebpv[2][3], ehp[3], x, y, s) -> astrom\\n"
"For a terrestrial observer, prepare star-independent astrometry\\n"
"parameters for transformations between ICRS and geocentric CIRS\\n"
"coordinates.  The Earth ephemeris and CIP/CIO are supplied by the caller.\\n"
"\\n"
"The parameters produced by this function are required in the\\n"
"parallax, light deflection, aberration, and bias-precession-nutation\\n"
"parts of the astrometric transformation chain.\\n"
"Given:\\n"
"    date1      TDB as a 2-part...\\n"
"    date2      ...Julian Date\\n"
"    ebpv       Earth barycentric pos/vel (au, au/day)\\n"
"    ehp        Earth heliocentric position (au)\\n"
"    x,y        CIP X,Y (components of unit vector)\\n"
"    s          the CIO locator s (radians)\\n"
"Returned:\\n"
"    astrom     star-independent astrometry parameters"''',
     '''{
    double date1, date2, ebpv[2][3], ehp[3], x, y, s;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "dd((ddd)(ddd))(ddd)ddd",
                          &date1, &date2,
                          &ebpv[0][0],&ebpv[0][1],&ebpv[0][2],
                          &ebpv[1][0],&ebpv[1][1],&ebpv[1][2],
                          &ehp[0],&ehp[1],&ehp[2],
                          &x, &y, &s))      
        return NULL;
    eraApci(date1, date2, ebpv, ehp, x, y, s, &astrom);
    return _to_py_astrom(&astrom);
}
'''),
    ('apci13', 'METH_VARARGS',
     '''"\\napci13(date1, date2) -> astrom, eo\\n"
"For a terrestrial observer, prepare star-independent astrometry\\n"
"parameters for transformations between ICRS and geocentric CIRS\\n"
"coordinates.  The caller supplies the date, and ERFA models are used\\n"
"to predict the Earth ephemeris and CIP/CIO.\\n"
"\\n"
"The parameters produced by this function are required in the\\n"
"parallax, light deflection and aberration parts of the astrometric\\n"
"transformation chain.\\n"
"Given:\\n"
"    date1      TDB as a 2-part...\\n"
"    date2      ...Julian Date\\n"
"Returned:\\n"
"    astrom     star-independent astrometry parameters\\n"
"    eo         equation of the origins (ERA-GST)"''',
     '''{
    double date1, date2, eo;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "dd", &date1, &date2))      
        return NULL;
    eraApci13(date1, date2, &astrom, &eo);
    return Py_BuildValue("Od", _to_py_astrom(&astrom), eo);
}
'''),
    ('apco', 'METH_VARARGS',
     '''"\\napco(date1, date2, ebpv[2][3], ehp[3], x, y, s,theta,,elong, phi, hm,xp, yp, refa, refb) -> astrom\\n"
"For a terrestrial observer, prepare star-independent astrometry\\n"
"parameters for transformations between ICRS and observed\\n"
"coordinates. The caller supplies the Earth ephemeris, the Earth\\n"
"rotation information and the refraction constants as well as the\\n"
"site coordinates.\\n"
"Given:\\n"
"    date1      TDB as a 2-part...\\n"
"    date2      ...Julian Date\\n"
"    ebpv       Earth barycentric pos/vel (au, au/day)\\n"
"    ehp        Earth heliocentric position (au)\\n"
"    x,y        CIP X,Y (components of unit vector)\\n"
"    s          the CIO locator s (radians)\\n"
"    theta      Earth rotation angle (radians)\\n"
"    elong      longitude (radians, east +ve)\\n"
"    phi        latitude (geodetic, radians)\\n"
"    hm         height above ellipsoid (m, geodetic3)\\n"
"    xp,yp      polar motion coordinates (radians)\\n"
"    sp         the TIO locator s' (radians)\\n"
"    refa       refraction constant A (radians)\\n"
"    refb       refraction constant B (radians)\\n"
"Returned:\\n"
"    astrom     star-independent astrometry parameters"''',
     '''{
    double date1, date2, ebpv[2][3], ehp[3], x, y, s;
    double theta, elong, phi, hm, xp, yp;
    double sp, refa, refb;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "dd((ddd)(ddd))(ddd)dddddddddddd",
                          &date1, &date2,
                          &ebpv[0][0],&ebpv[0][1],&ebpv[0][2],
                          &ebpv[1][0],&ebpv[1][1],&ebpv[1][2],
                          &ehp[0],&ehp[1],&ehp[2],
                          &x, &y, &s,
                          &theta, &elong, &phi, &hm, &xp, &yp,
                          &sp, &refa, &refb))      
        return NULL;
    eraApco(date1, date2, ebpv, ehp, x, y, s,
            theta, elong, phi, hm, xp, yp,
            sp, refa, refb,
            &astrom);
    return _to_py_astrom(&astrom);
}
'''),
    ('apco13', 'METH_VARARGS',
     '''"\\napco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl) -> astrom, eo\\n"
"For a terrestrial observer, prepare star-independent astrometry\\n"
"parameters for transformations between ICRS and geocentric CIRS\\n"
"coordinates.  The caller supplies the date, and ERFA models are used\\n"
"to predict the Earth ephemeris and CIP/CIO.\\n"
"\\n"
"The parameters produced by this function are required in the\\n"
"parallax, light deflection, aberration, and bias-precession-nutation\\n"
"parts of the ICRS/CIRS transformations.\\n"
"Given:\\n"
"    utc1   double     UTC as a 2-part...\\n"
"    utc2   double     ...quasi Julian Date\\n"
"    dut1   double     UT1-UTC (seconds)\\n"
"    theta  double     Earth rotation angle (radians)\\n"
"    elong  double     longitude (radians, east +ve)\\n"
"    phi    double     latitude (geodetic, radians)\\n"
"    hm     double     height above ellipsoid (m, geodetic)\\n"
"    xp,yp  double     polar motion coordinates (radians)\\n"
"    phpa   double     pressure at the observer (hPa = mB)\\n"
"    tc     double     ambient temperature at the observer (deg C)\\n"
"    rh     double     relative humidity at the observer (range 0-1)\\n"
"    wl     double     wavelength (micrometers)\\n"
"Returned:\\n"
"    astrom     star-independent astrometry parameters\\n"
"    eo         equation of the origins (ERA-GST)"''',
     '''{
    int j;
    double utc1, utc2, dut1;
    double elong, phi, hm, xp, yp;
    double phpa, tc, rh, wl, eo;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "dddddddddddd",
                                 &utc1, &utc2, &dut1,
                                 &elong, &phi, &hm, &xp, &yp,
                                 &phpa, &tc, &rh, &wl))      
        return NULL;
    j = eraApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &astrom, &eo);
    if (j == +1) {
        PyErr_SetString(_erfaError, "doubious year");
        return NULL;
    }
    else if (j == -1) {
        PyErr_SetString(_erfaError, "unacceptable date");
        return NULL;
    }
    return Py_BuildValue("Od", _to_py_astrom(&astrom), eo);
}
'''),
    ('apcs', 'METH_VARARGS',
     '''"\\napcs(date1, date2, pv[2][3], ebpv[2][3], ehp[3]) -> astrom\\n"
"For an observer whose geocentric position and velocity are known,\\n"
"prepare star-independent astrometry parameters for transformations\\n"
"between ICRS and GCRS.  The Earth ephemeris is supplied by the caller.\\n"
"\\n"
"The parameters produced by this function are required in the\\n"
"parallax, light deflection and aberration parts of the astrometric\\n"
"transformation chain.\\n"
"Given:\\n"
"    date1      TDB as a 2-part...\\n"
"    date2      ...Julian Date\\n"
"    pv         observer's geocentric pos/vel (m, m/s)\\n"
"    ebpv       Earth barycentric pos/vel (au, au/day)\\n"
"    ehp        Earth heliocentric position (au)\\n"
"Returned:\\n"
"    astrom     star-independent astrometry parameters"''',
     '''{
    double date1, date2, pv[2][3], ebpv[2][3], ehp[3];
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "dd((ddd)(ddd))((ddd)(ddd))(ddd)",
                          &date1, &date2,
                          &pv[0][0],&pv[0][1],&pv[0][2],
                          &pv[1][0],&pv[1][1],&pv[1][2],
                          &ebpv[0][0],&ebpv[0][1],&ebpv[0][2],
                          &ebpv[1][0],&ebpv[1][1],&ebpv[1][2],
                          &ehp[0],&ehp[1],&ehp[2]))      
        return NULL;
    eraApcs(date1, date2, pv, ebpv, ehp, &astrom);
    return _to_py_astrom(&astrom);
}
'''),
    ('apcs13', 'METH_VARARGS',
     '''"\\napcs13(date1, date2, pv[2][3]) -> astrom\\n"
"For an observer whose geocentric position and velocity are known,\\n"
"prepare star-independent astrometry parameters for transformations\\n"
"between ICRS and GCRS.  The Earth ephemeris is is from ERFA models.\\n"
"\\n"
"The parameters produced by this function are required in the space\\n"
"motion, parallax, light deflection and aberration parts of the astrometric\\n"
"transformation chain.\\n"
"Given:\\n"
"    date1      TDB as a 2-part...\\n"
"    date2      ...Julian Date\\n"
"    pv         observer's geocentric pos/vel (m, m/s)\\n"
"Returned:\\n"
"    astrom     star-independent astrometry parameters"''',
     '''{
    double date1, date2, pv[2][3];
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "dd((ddd)(ddd))",
                          &date1, &date2,
                          &pv[0][0],&pv[0][1],&pv[0][2],
                          &pv[1][0],&pv[1][1],&pv[1][2]))      
        return NULL;
    eraApcs13(date1, date2, pv, &astrom);
    return _to_py_astrom(&astrom);
}
'''),
##
## aper[13] defined is erfa.py
##
    ('apio', 'METH_VARARGS',
     '''"\\napio(sp,theta,elong,phi,hm,xp,yp,refa,refb) -> astrom\\n"
"For a terrestrial observer, prepare star-independent astrometry\\n"
"parameters for transformations between CIRS and observed\\n"
"coordinates.  The caller supplies the Earth orientation information\\n"
"and the refraction constants as well as the site coordinates.\\n"
"Given:\\n"
"    sp         the TIO locator s'\\n"
"    theta      Earth rotation angle (radians)\\n"
"    elong      longitude (radians, east +ve)\\n"
"    phi        latitude (geodetic, radians)\\n"
"    hm         height above ellipsoid (m, geodetic3)\\n"
"    xp,yp      polar motion coordinates (radians)\\n"
"    refa       refraction constant A (radians)\\n"
"    refb       refraction constant B (radians)\\n"
"Returned:\\n""   x       \\n"
"    astrom     star-independent astrometry parameters"''',
     '''{
    double sp, theta, elong, phi, hm, xp, yp, refa, refb;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "ddddddddd",
                          &sp, &theta, &elong, &phi, &hm,
                          &xp, &yp, &refa, &refb))      
        return NULL;
    eraApio(sp,theta,elong,phi,hm,xp,yp,refa,refb, &astrom);
    return _to_py_astrom(&astrom);
}
'''),
    ('apio13', 'METH_VARARGS',
     '''"\\napio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl) -> astrom, eo\\n"
"For a terrestrial observer, prepare star-independent astrometry\\n"
"parameters for transformations between CIRS and observed\\n"
"coordinates.  The caller supplies UTC, site coordinates, ambient air\\n"
"conditions and observing wavelength.\\n"
"Given:\\n"
"    utc1   double     UTC as a 2-part...\\n"
"    utc2   double     ...quasi Julian Date\\n"
"    dut1   double     UT1-UTC (seconds)\\n"
"    elong  double     longitude (radians, east +ve)\\n"
"    phi    double     latitude (geodetic, radians)\\n"
"    hm     double     height above ellipsoid (m, geodetic)\\n"
"    xp,yp  double     polar motion coordinates (radians)\\n"
"    phpa   double     pressure at the observer (hPa = mB)\\n"
"    tc     double     ambient temperature at the observer (deg C)\\n"
"    rh     double     relative humidity at the observer (range 0-1)\\n"
"    wl     double     wavelength (micrometers)\\n"
"Returned:\\n"
"    astrom     star-independent astrometry parameters"''',
     '''{
    int j;
    double utc1, utc2, dut1;
    double elong, phi, hm, xp, yp;
    double phpa, tc, rh, wl;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "dddddddddddd",
                                 &utc1, &utc2, &dut1,
                                 &elong, &phi, &hm, &xp, &yp,
                                 &phpa, &tc, &rh, &wl))      
        return NULL;
    j = eraApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, &astrom);
    if (j == +1) {
        PyErr_SetString(_erfaError, "doubious year");
        return NULL;
    }
    else if (j == -1) {
        PyErr_SetString(_erfaError, "unacceptable date");
        return NULL;
    }
    return _to_py_astrom(&astrom);
}
'''),
    ('atci13', 'METH_VARARGS',
     '''"\\natci13(rc, dc, pr, pd, px, rv, date1, date2) -> ri, di, eo\\n"
"Transform ICRS star data, epoch J2000.0, to CIRS.\\n"
"Given:\\n"
"    rc     ICRS right ascension at J2000.0 (radians)\\n"
"    dc     ICRS declination at J2000.0 (radians)\\n"
"    pr     RA proper motion (radians/year)\\n"
"    pd     Dec proper motion (radians/year)\\n"
"    px     parallax (arcsec)\\n"
"    rv     radial velocity (km/s, +ve if receding)\\n"
"    date1  TDB as a 2-part...\\n"
"    date2  ...Julian Date\\n"
"Returned:\\n"
"    ri,di  CIRS geocentric RA,Dec (radians)\\n"
"    eo     double*  equation of the origins (ERA-GST)"''',
     '''{
     double rc, dc, pr, pd, px, rv, date1, date2;
     double ri, di, eo;
    if (!PyArg_ParseTuple(args, "dddddddd",
                                 &rc, &dc, &pr, &pd, &px, &rv, &date1, &date2))      
        return NULL;
    eraAtci13(rc, dc, pr, pd, px, rv, date1, date2, &ri, &di, &eo);
    return Py_BuildValue("ddd", ri, di, eo);
}
'''),
    ('atciq', 'METH_VARARGS',
     '''"\\natciq( rc, dc, pr, pd, px, rv, astrom) -> ri,di\\n"
"Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed\\n"
"star-independent astrometry parameters.\\n"
"\\n"
"Use of this function is appropriate when efficiency is important and\\n"
"where many star positions are to be transformed for one date.  The\\n"
"star-independent parameters can be obtained by calling one of the\\n"
"functions apci[13], apcg[13], apco[13] or apcs[13].\\n"
"\\n"
"If the parallax and proper motions are zero the eraAtciqz function\\n"
"can be used instead.\\n"
"Given:\\n"
"    rc,dc      ICRS RA,Dec at J2000.0 (radians)\\n"
"    pr         RA proper motion (radians/year)\\n"
"    pd         Dec proper motion (radians/year)\\n"
"    px         parallax (arcsec)\\n"
"    rv         radial velocity (km/s, +ve if receding)\\n"
"    astrom     star-independent astrometry parameters\\n"
"Returned:\\n"
"    ri,di      CIRS RA,Dec (radians)\\n"''',
     '''{
    double rc, dc, pr, pd, px, rv;
    double ri, di;
    double pmt, eb[3], eh[3], em, v[3], bm1, bpn[3][3];
    double along, phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "ddddddd(ddd)(ddd)d(ddd)d((ddd)(ddd)(ddd))dddddddddd",
                          &rc, &dc, &pr, &pd, &px, &rv,
                          &pmt, &eb[0], &eb[1], &eb[2],
                          &eh[0], &eh[1], &eh[2], &em,
                          &v[0], &v[1], &v[2], &bm1,
                          &bpn[0][0],&bpn[0][1],&bpn[0][2],
                          &bpn[1][0],&bpn[1][1],&bpn[1][2],      
                          &bpn[2][0],&bpn[2][1],&bpn[2][2],
                          &along, &phi, &xpl, &ypl, &sphi, &cphi, &diurab,
                          &eral, &refa, &refb))      
        return NULL;
    astrom.pmt = pmt;
    astrom.eb[0] = eb[0];
    astrom.eb[1] = eb[1];
    astrom.eb[2] = eb[2];
    astrom.eh[0] = eh[0];
    astrom.eh[1] = eh[1];
    astrom.eh[2] = eh[2];
    astrom.em = em;
    astrom.v[0] = v[0];
    astrom.v[1] = v[1];
    astrom.v[2] = v[2];
    astrom.bm1 = bm1;
    astrom.bpn[0][0] = bpn[0][0];
    astrom.bpn[0][1] = bpn[0][1];
    astrom.bpn[0][2] = bpn[0][2];
    astrom.bpn[1][0] = bpn[1][0];
    astrom.bpn[1][1] = bpn[1][1];
    astrom.bpn[1][2] = bpn[1][2];
    astrom.bpn[2][0] = bpn[2][0];
    astrom.bpn[2][1] = bpn[2][1];
    astrom.bpn[2][2] = bpn[2][2];
    astrom.along = along;
    astrom.phi = phi;
    astrom.xpl = xpl;
    astrom.ypl = ypl;
    astrom.sphi = sphi;
    astrom.cphi = cphi;
    astrom.diurab = diurab;
    astrom.eral = eral;
    astrom.refa = refa;
    astrom.refb = refb;
    eraAtciq(rc, dc, pr, pd, px, rv, &astrom, &ri, &di);
    return Py_BuildValue("dd", ri, di);
}
'''),
##
## atciqn see erfa.py
##
    ('atciqz', 'METH_VARARGS',
     '''"\\natciqz( rc, dc, astrom) -> ri,di\\n"
"Quick ICRS to CIRS transformation, given precomputed star-\\n"
"independent astrometry parameters, and assuming zero parallax and proper motion.\\n"
"\\n"
"Use of this function is appropriate when efficiency is important and\\n"
"where many star positions are to be transformed for one date.  The\\n"
"star-independent parameters can be obtained by calling one of the\\n"
"functions apci[13], apcg[13], apco[13] or apcs[13].\\n"
"\\n"
"The corresponding function for the case of non-zero parallax and\\n"
"proper motion is atciq.\\n"
"Given:\\n"
"    rc,dc      ICRS RA,Dec at J2000.0 (radians)\\n"
"    astrom     star-independent astrometry parameters\\n"
"Returned:\\n"
"    ri,di      CIRS RA,Dec (radians)\\n"''',
     '''{
    double rc, dc;
    double ri, di;
    double pmt, eb[3], eh[3], em, v[3], bm1, bpn[3][3];
    double along, phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "ddd(ddd)(ddd)d(ddd)d((ddd)(ddd)(ddd))dddddddddd",
                          &rc, &dc,
                          &pmt, &eb[0], &eb[1], &eb[2],
                          &eh[0], &eh[1], &eh[2], &em,
                          &v[0], &v[1], &v[2], &bm1,
                          &bpn[0][0],&bpn[0][1],&bpn[0][2],
                          &bpn[1][0],&bpn[1][1],&bpn[1][2],      
                          &bpn[2][0],&bpn[2][1],&bpn[2][2],
                          &along, &phi, &xpl, &ypl, &sphi, &cphi, &diurab,
                          &eral, &refa, &refb))      
        return NULL;
    astrom.pmt = pmt;
    astrom.eb[0] = eb[0];
    astrom.eb[1] = eb[1];
    astrom.eb[2] = eb[2];
    astrom.eh[0] = eh[0];
    astrom.eh[1] = eh[1];
    astrom.eh[2] = eh[2];
    astrom.em = em;
    astrom.v[0] = v[0];
    astrom.v[1] = v[1];
    astrom.v[2] = v[2];
    astrom.bm1 = bm1;
    astrom.bpn[0][0] = bpn[0][0];
    astrom.bpn[0][1] = bpn[0][1];
    astrom.bpn[0][2] = bpn[0][2];
    astrom.bpn[1][0] = bpn[1][0];
    astrom.bpn[1][1] = bpn[1][1];
    astrom.bpn[1][2] = bpn[1][2];
    astrom.bpn[2][0] = bpn[2][0];
    astrom.bpn[2][1] = bpn[2][1];
    astrom.bpn[2][2] = bpn[2][2];
    astrom.along = along;
    astrom.phi = phi;
    astrom.xpl = xpl;
    astrom.ypl = ypl;
    astrom.sphi = sphi;
    astrom.cphi = cphi;
    astrom.diurab = diurab;
    astrom.eral = eral;
    astrom.refa = refa;
    astrom.refb = refb;
    eraAtciqz(rc, dc, &astrom, &ri, &di);
    return Py_BuildValue("dd", ri, di);
}
'''),
    ('atco13', 'METH_VARARGS',
     '''"\\natco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl) -> aob, zob, hob, dob, rob, eo\\n"
"ICRS RA,Dec to observed place.  The caller supplies UTC, site\\n"
"coordinates, ambient air conditions and observing wavelength.\\n"
"\\n"
"ERFA models are used for the Earth ephemeris, bias-precession-\\n"
"nutation, Earth orientation and refraction.\\n"
"Given:\\n"
"    rc,dc  ICRS right ascension at J2000.0 (radians)\\n"
"    pr     RA proper motion (radians/year)\\n"
"    pd     Dec proper motion (radians/year)\\n"
"    px     parallax (arcsec)\\n"
"    rv     radial velocity (km/s, +ve if receding)\\n"
"    utc1   UTC as a 2-part...\\n"
"    utc2    ...quasi Julian Date\\n"
"    dut1   UT1-UTC (seconds)\\n"
"    elong  longitude (radians, east +ve)\\n"
"    phi    latitude (geodetic, radians)\\n"
"    hm     height above ellipsoid (m, geodetic)\\n"
"    xp,yp  polar motion coordinates (radians)\\n"
"    phpa   pressure at the observer (hPa = mB)\\n"
"    tc     ambient temperature at the observer (deg C)\\n"
"    rh     relative humidity at the observer (range 0-1)\\n"
"    wl     wavelength (micrometers)\\n"
"Returned:\\n"
"    aob    observed azimuth (radians: N=0,E=90)\\n"
"    zob    observed zenith distance (radians)\\n"
"    hob    observed hour angle (radians)\\n"
"    dob    observed declination (radians)\\n"
"    rob    observed right ascension (CIO-based, radians)\\n"
"    eo     equation of the origins (ERA-GST)"''',
     '''{
    int j;
    double rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl;
    double aob, zob, hob, dob, rob, eo;
    if (!PyArg_ParseTuple(args, "dddddddddddddddddd",
                          &rc, &dc,
                          &pr, &pd, &px, &rv, &utc1, &utc2, &dut1,
                          &elong, &phi, &hm, &xp, &yp,
                          &phpa, &tc, &rh, &wl))
        return NULL;
    j = eraAtco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                  &aob, &zob, &hob, &dob, &rob, &eo);
    if (j == +1) {
        PyErr_SetString(_erfaError, "doubious year");
        return NULL;
    }
    else if (j == -1) {
        PyErr_SetString(_erfaError, "unacceptable date");
        return NULL;
    }
    return Py_BuildValue("dddddd", aob, zob, hob, dob, rob, eo);
}
'''),
    ('atic13', 'METH_VARARGS',
     '''"\\natic13(ri, di, date1, date2) -> rc, dc, eo\\n"
"Transform star RA,Dec from geocentric CIRS to ICRS astrometric.\\n"
"Given:\\n"
"    ri,di  CIRS geocentric RA,Dec (radians)\\n"
"    date1  TDB as a 2-part...\\n"
"    date2  ...Julian Date\\n"
"Returned:\\n"
"    rc,dc  ICRS astrometric RA,Dec (radians)\\n"
"    eo     equation of the origins (ERA-GST)"''',
     '''{
    double ri, di, date1, date2;
    double rc, dc, eo;
    if (!PyArg_ParseTuple(args, "dddd", &ri, &di, &date1, &date2))
        return NULL;
    eraAtic13(ri, di, date1, date2, &rc, &dc, &eo);
    return Py_BuildValue("ddd", rc, dc, eo);
}
'''),
    ('aticq', 'METH_VARARGS',
     '''"\\naticq(ri,di, astrom) -> rc, dc\\n"
"Quick CIRS RA,Dec to ICRS astrometric place, given the star-\\n"
"independent astrometry parameters.\\n"
"\\n"
"Use of this function is appropriate when efficiency is important and\\n"
"where many star positions are all to be transformed for one date.  The\\n"
"star-independent parameters can be obtained by calling one of the\\n"
"functions apci[13], apcg[13], apco[13] or apcs[13].\\n"
"Given:\\n"
"    ri,di      CIRS RA,Dec (radians)\\n"
"    astrom     star-independent astrometry parameters\\n"
"Returned:\\n"
"    rc,dc      ICRS astrometric RA,Dec (radians)\\n"''',
     '''{
    double rc, dc, ri, di;
    double pmt, eb[3], eh[3], em, v[3], bm1, bpn[3][3];
    double along, phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "ddd(ddd)(ddd)d(ddd)d((ddd)(ddd)(ddd))dddddddddd",
                          &ri, &di,
                          &pmt, &eb[0], &eb[1], &eb[2],
                          &eh[0], &eh[1], &eh[2], &em,
                          &v[0], &v[1], &v[2], &bm1,
                          &bpn[0][0],&bpn[0][1],&bpn[0][2],
                          &bpn[1][0],&bpn[1][1],&bpn[1][2],      
                          &bpn[2][0],&bpn[2][1],&bpn[2][2],
                          &along, &phi, &xpl, &ypl, &sphi, &cphi, &diurab,
                          &eral, &refa, &refb))      
        return NULL;
    astrom.pmt = pmt;
    astrom.eb[0] = eb[0];
    astrom.eb[1] = eb[1];
    astrom.eb[2] = eb[2];
    astrom.eh[0] = eh[0];
    astrom.eh[1] = eh[1];
    astrom.eh[2] = eh[2];
    astrom.em = em;
    astrom.v[0] = v[0];
    astrom.v[1] = v[1];
    astrom.v[2] = v[2];
    astrom.bm1 = bm1;
    astrom.bpn[0][0] = bpn[0][0];
    astrom.bpn[0][1] = bpn[0][1];
    astrom.bpn[0][2] = bpn[0][2];
    astrom.bpn[1][0] = bpn[1][0];
    astrom.bpn[1][1] = bpn[1][1];
    astrom.bpn[1][2] = bpn[1][2];
    astrom.bpn[2][0] = bpn[2][0];
    astrom.bpn[2][1] = bpn[2][1];
    astrom.bpn[2][2] = bpn[2][2];
    astrom.along = along;
    astrom.phi = phi;
    astrom.xpl = xpl;
    astrom.ypl = ypl;
    astrom.sphi = sphi;
    astrom.cphi = cphi;
    astrom.diurab = diurab;
    astrom.eral = eral;
    astrom.refa = refa;
    astrom.refb = refb;
    eraAticq(ri, di, &astrom, &rc, &dc);
    return Py_BuildValue("dd", rc, dc);
}
'''),
##
## aticqn see erfa.py
##
    ('atio13', 'METH_VARARGS',
     '''"\\natio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl) -> aob, zob, hob, dob, rob\\n"
"CIRS RA,Dec to observed place.  The caller supplies UTC, site\\n"
"coordinates, ambient air conditions and observing wavelength.\\n"
"Given:\\n"
"    ri     CIRS right ascension (CIO-based, radians)\\n"
"    di     CIRS declination (radians)\\n"
"    utc1   UTC as a 2-part...\\n"
"    utc2    ...quasi Julian Date\\n"
"    dut1   UT1-UTC (seconds)\\n"
"    elong  longitude (radians, east +ve)\\n"
"    phi    geodetic latitude (radians)\\n"
"    hm     height above ellipsoid (m, geodetic)\\n"
"    xp,yp  polar motion coordinates (radians)\\n"
"    phpa   pressure at the observer (hPa = mB)\\n"
"    tc     ambient temperature at the observer (deg C)\\n"
"    rh     relative humidity at the observer (range 0-1)\\n"
"    wl     wavelength (micrometers)\\n"
"Returned:\\n"
"    aob    observed azimuth (radians: N=0,E=90)\\n"
"    zob    observed zenith distance (radians)\\n"
"    hob    observed hour angle (radians)\\n"
"    dob    observed declination (radians)\\n"
"    rob    observed right ascension (CIO-based, radians)"''',
     '''{
    int j;
    double ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl;
    double aob, zob, hob, dob, rob;
    if (!PyArg_ParseTuple(args, "dddddddddddddd",
                          &ri, &di, &utc1, &utc2, &dut1,
                          &elong, &phi, &hm, &xp, &yp,
                          &phpa, &tc, &rh, &wl))
        return NULL;
    j = eraAtio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                  &aob, &zob, &hob, &dob, &rob);
    if (j == +1) {
        PyErr_SetString(_erfaError, "doubious year");
        return NULL;
    }
    else if (j == -1) {
        PyErr_SetString(_erfaError, "unacceptable date");
        return NULL;
    }
    return Py_BuildValue("ddddd", aob, zob, hob, dob, rob);
}
'''),
    ('atioq', 'METH_VARARGS',
     '''"\\natioq(ri,di, astrom) -> aob, zob, hob, dob, rob\\n"
"Quick CIRS to observed place transformation.\\n"
"\\n"
"Use of this function is appropriate when efficiency is important and\\n"
"where many star positions are all to be transformed for one date.  The\\n"
"star-independent parameters can be obtained by calling\\n"
"apio[13], or apco[13].\\n"
"Given:\\n"
"    ri,di      CIRS RA,Dec (radians)\\n"
"    astrom     star-independent astrometry parameters\\n"
"Returned:\\n"
"    aob    observed azimuth (radians: N=0,E=90)\\n"
"    zob    observed zenith distance (radians)\\n"
"    hob    observed hour angle (radians)\\n"
"    dob    observed declination (radians)\\n"
"    rob    observed right ascension (CIO-based, radians)"''',
     '''{
    double ri, di;
    double pmt, eb[3], eh[3], em, v[3], bm1, bpn[3][3];
    double along, phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb;
    double aob, zob, hob, dob, rob;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "ddd(ddd)(ddd)d(ddd)d((ddd)(ddd)(ddd))dddddddddd",
                          &ri, &di,
                          &pmt, &eb[0], &eb[1], &eb[2],
                          &eh[0], &eh[1], &eh[2], &em,
                          &v[0], &v[1], &v[2], &bm1,
                          &bpn[0][0],&bpn[0][1],&bpn[0][2],
                          &bpn[1][0],&bpn[1][1],&bpn[1][2],      
                          &bpn[2][0],&bpn[2][1],&bpn[2][2],
                          &along, &phi, &xpl, &ypl, &sphi, &cphi, &diurab,
                          &eral, &refa, &refb))      
        return NULL;
    astrom.pmt = pmt;
    astrom.eb[0] = eb[0];
    astrom.eb[1] = eb[1];
    astrom.eb[2] = eb[2];
    astrom.eh[0] = eh[0];
    astrom.eh[1] = eh[1];
    astrom.eh[2] = eh[2];
    astrom.em = em;
    astrom.v[0] = v[0];
    astrom.v[1] = v[1];
    astrom.v[2] = v[2];
    astrom.bm1 = bm1;
    astrom.bpn[0][0] = bpn[0][0];
    astrom.bpn[0][1] = bpn[0][1];
    astrom.bpn[0][2] = bpn[0][2];
    astrom.bpn[1][0] = bpn[1][0];
    astrom.bpn[1][1] = bpn[1][1];
    astrom.bpn[1][2] = bpn[1][2];
    astrom.bpn[2][0] = bpn[2][0];
    astrom.bpn[2][1] = bpn[2][1];
    astrom.bpn[2][2] = bpn[2][2];
    astrom.along = along;
    astrom.phi = phi;
    astrom.xpl = xpl;
    astrom.ypl = ypl;
    astrom.sphi = sphi;
    astrom.cphi = cphi;
    astrom.diurab = diurab;
    astrom.eral = eral;
    astrom.refa = refa;
    astrom.refb = refb;
    eraAtioq(ri, di, &astrom, &aob, &zob, &hob, &dob, &rob);
    return Py_BuildValue("ddddd", aob, zob, hob, dob, rob);
}
'''),
    ('atoc13', 'METH_VARARGS',
     '''"\\natoc13(type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl) -> rc, dc\\n"
"Observed place at a groundbased site to to ICRS astrometric RA,Dec.\\n"
"The caller supplies UTC, site coordinates, ambient air conditions\\n"
"and observing wavelength.\\n"
"Given:\\n"
"    type   type of coordinates - ''R'', ''H'' or ''A''\\n"
"    ob1    observed Az, HA or RA (radians; Az is N=0,E=90)\\n"
"    ob2    observed ZD or Dec (radians)\\n"
"    utc1   UTC as a 2-part...\\n"
"    utc2   ...quasi Julian Date\\n"
"    dut1   UT1-UTC (seconds\\n"
"    elong  longitude (radians, east +ve)\\n"
"    phi    geodetic latitude (radians)\\n"
"    hm     height above ellipsoid (m, geodetic)\\n"
"    xp,yp  polar motion coordinates (radians)\\n"
"    phpa   pressure at the observer (hPa = mB)\\n"
"    tc     ambient temperature at the observer (deg C)\\n"
"    rh     relative humidity at the observer (range 0-1)\\n"
"    wl     wavelength (micrometers)\\n"
"Returned:\\n"
"    rc,dc  ICRS astrometric RA,Dec (radians)"''',
     '''{
    int j;
    const char *type;
    double ob1, ob2, utc1, utc2, dut1;
    double elong, phi, hm, xp, yp, phpa, tc, rh, wl;
    double rc, dc;
    if (!PyArg_ParseTuple(args, "sdddddddddddddd",
                                 &type, &ob1, &ob2, &utc1, &utc2, &dut1,
                                 &elong, &phi, &hm, &xp, &yp, &phpa, &tc, &rh, &wl))      
        return NULL;
    if (strcmp("R", type) == 0 || strcmp("H", type) == 0 || strcmp("A", type) == 0) {
        j = eraAtoc13(type, ob1, ob2, utc1, utc2, dut1,
                      elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                      &rc, &dc);
    }
    else {
        PyErr_SetString(_erfaError, "unknown type of coordinates");
        return NULL;
    }
    if (j == +1) {
        PyErr_SetString(_erfaError, "doubious year");
        return NULL;
    }
    else if (j == -1) {
        PyErr_SetString(_erfaError, "unacceptable date");
        return NULL;
    }
    return Py_BuildValue("dd", rc, dc);
}
'''),
    ('atoi13', 'METH_VARARGS',
     '''"\\natoi13(type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl) -> ri, di\\n"
"Observed place at a groundbased site to to ICRS astrometric RA,Dec.\\n"
"The caller supplies UTC, site coordinates, ambient air conditions\\n"
"and observing wavelength.\\n"
"Given:\\n"
"    type   type of coordinates - ''R'', ''H'' or ''A''\\n"
"    ob1    observed Az, HA or RA (radians; Az is N=0,E=90)\\n"
"    ob2    observed ZD or Dec (radians)\\n"
"    utc1   UTC as a 2-part...\\n"
"    utc2   ...quasi Julian Date\\n"
"    dut1   UT1-UTC (seconds\\n"
"    elong  longitude (radians, east +ve)\\n"
"    phi    geodetic latitude (radians)\\n"
"    hm     height above ellipsoid (m, geodetic)\\n"
"    xp,yp  polar motion coordinates (radians)\\n"
"    phpa   pressure at the observer (hPa = mB)\\n"
"    tc     ambient temperature at the observer (deg C)\\n"
"    rh     relative humidity at the observer (range 0-1)\\n"
"    wl     wavelength (micrometers)\\n"
"Returned:\\n"
"    ri     CIRS right ascension (CIO-based, radians)\\n"
"    di     CIRS declination (radians)"''',
     '''{
    int j;
    const char *type;
    double ob1, ob2, utc1, utc2, dut1;
    double elong, phi, hm, xp, yp, phpa, tc, rh, wl;
    double ri, di;
    if (!PyArg_ParseTuple(args, "sdddddddddddddd",
                                 &type, &ob1, &ob2, &utc1, &utc2, &dut1,
                                 &elong, &phi, &hm, &xp, &yp, &phpa, &tc, &rh, &wl))      
        return NULL;
    if (strcmp("R", type) == 0 || strcmp("H", type) == 0 || strcmp("A", type) == 0) {
        j = eraAtoi13(type, ob1, ob2, utc1, utc2, dut1,
                      elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                      &ri, &di);
    }
    else {
        PyErr_SetString(_erfaError, "unknown type of coordinates");
        return NULL;
    }
    if (j == +1) {
        PyErr_SetString(_erfaError, "doubious year");
        return NULL;
    }
    else if (j == -1) {
        PyErr_SetString(_erfaError, "unacceptable date");
        return NULL;
    }
    return Py_BuildValue("dd", ri, di);
}
'''),
    ('atoiq', 'METH_VARARGS',
     '''"\\natoiq(type, ob1, ob2, astrom) -> ri, di\\n"
"Observed place at a groundbased site to to ICRS astrometric RA,Dec.\\n"
"The caller supplies UTC, site coordinates, ambient air conditions\\n"
"and observing wavelength.\\n"
"Given:\\n"
"    type   type of coordinates - ''R'', ''H'' or ''A''\\n"
"    ob1    observed Az, HA or RA (radians; Az is N=0,E=90)\\n"
"    ob2    observed ZD or Dec (radians)\\n"
"    astrom star-independent astrometry parameters\\n"
"Returned:\\n"
"    ri     CIRS right ascension (CIO-based, radians)\\n"
"    di     CIRS declination (radians)"''',
     '''{
    const char *type;
    double ob1, ob2;
    double ri, di;
    double pmt, eb[3], eh[3], em, v[3], bm1, bpn[3][3];
    double along, phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "sddd(ddd)(ddd)d(ddd)d((ddd)(ddd)(ddd))dddddddddd",
                          &type, &ob1, &ob2,
                          &pmt, &eb[0], &eb[1], &eb[2],
                          &eh[0], &eh[1], &eh[2], &em,
                          &v[0], &v[1], &v[2], &bm1,
                          &bpn[0][0],&bpn[0][1],&bpn[0][2],
                          &bpn[1][0],&bpn[1][1],&bpn[1][2],      
                          &bpn[2][0],&bpn[2][1],&bpn[2][2],
                          &along, &phi, &xpl, &ypl, &sphi, &cphi, &diurab,
                          &eral, &refa, &refb))      
        return NULL;
    astrom.pmt = pmt;
    astrom.eb[0] = eb[0];
    astrom.eb[1] = eb[1];
    astrom.eb[2] = eb[2];
    astrom.eh[0] = eh[0];
    astrom.eh[1] = eh[1];
    astrom.eh[2] = eh[2];
    astrom.em = em;
    astrom.v[0] = v[0];
    astrom.v[1] = v[1];
    astrom.v[2] = v[2];
    astrom.bm1 = bm1;
    astrom.bpn[0][0] = bpn[0][0];
    astrom.bpn[0][1] = bpn[0][1];
    astrom.bpn[0][2] = bpn[0][2];
    astrom.bpn[1][0] = bpn[1][0];
    astrom.bpn[1][1] = bpn[1][1];
    astrom.bpn[1][2] = bpn[1][2];
    astrom.bpn[2][0] = bpn[2][0];
    astrom.bpn[2][1] = bpn[2][1];
    astrom.bpn[2][2] = bpn[2][2];
    astrom.along = along;
    astrom.phi = phi;
    astrom.xpl = xpl;
    astrom.ypl = ypl;
    astrom.sphi = sphi;
    astrom.cphi = cphi;
    astrom.diurab = diurab;
    astrom.eral = eral;
    astrom.refa = refa;
    astrom.refb = refb;
    if (strcmp("R", type) == 0 || strcmp("H", type) == 0 || strcmp("A", type) == 0) {
        eraAtoiq(type, ob1, ob2, &astrom, &ri, &di);
    }
    else {
        PyErr_SetString(_erfaError, "unknown type of coordinates");
        return NULL;
    }
    return Py_BuildValue("dd", ri, di);
}
'''),
    ('ld', 'METH_VARARGS',
     '''"\\nld(bm, p[3], q[3], e[3], em, dlim) -> p1[3]\\n"
"Apply light deflection by a solar-system body, as part of\\n"
"transforming coordinate direction into natural direction.\\n"
"Given:\\n"
"    bm     mass of the gravitating body (solar masses)\\n"
"    p      direction from observer to source (unit vector)\\n"
"    q      direction from body to source (unit vector)\\n"
"    e      direction from body to observer (unit vector)\\n"
"    em     distance from body to observer (au)\\n"
"    dlim   deflection limiter\\n"
"Returned:\\n"
"    p1     observer to deflected source (unit vector)"''',
     '''{
    double bm, p[3], q[3], e[3], em, dlim, p1[3];
    if (!PyArg_ParseTuple(args, "d(ddd)(ddd)(ddd)dd",
                                 &bm,
                                 &p[0], &p[1], &p[2],
                                 &q[0], &q[1], &q[2],
                                 &e[0], &e[1], &e[2],
                                 &em, &dlim))
        return NULL;
    eraLd(bm, p, q, e, em, dlim, p1);
    return Py_BuildValue("(ddd)", p1[0], p1[1], p1[2]);
}
'''),
##
## ldn see erfa.py
##
    ('ldsun', 'METH_VARARGS',
     '''"\\nldsun(p[3], e[3], em) -> p1[3]\\n"
"Light deflection by the Sun.\\n"
"Given:\\n"
"    p      direction from observer to source (unit vector)\\n"
"    e      direction from Sun to observer (unit vector)\\n"
"    em     distance from Sun to observer (au)\\n"
"Returned:\\n"
"    p1     observer to deflected source (unit vector)"''',
     '''{
    double p[3], e[3], em, p1[3];
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)d",
                                 &p[0], &p[1], &p[2],
                                 &e[0], &e[1], &e[2],
                                 &em))
        return NULL;
    eraLdsun(p, e, em, p1);
    return Py_BuildValue("(ddd)", p1[0], p1[1], p1[2]);
}
'''),
    ('pmpx', 'METH_VARARGS',
     '''"\\npmpx(rc, dc, pr, pd, px, rv, pmt, pob[3],) -> pco[3]\\n"
"Proper motion and parallax.\\n"
"Given:\\n"
"    rc,dc  ICRS RA,Dec at catalog epoch (radians)\\n"
"    pr     RA proper motion (radians/year; Note 1)\\n"
"    pd     Dec proper motion (radians/year)\\n"
"    px     parallax (arcsec)\\n"
"    rv     radial velocity (km/s, +ve if receding)\\n"
"    pmt    proper motion time interval (SSB, Julian years)\\n"
"    pob    SSB to observer vector (au)\\n"
"Returned:\\n"
"    pco    coordinate direction (BCRS unit vector)"''',
     '''{
    double rc, dc, pr, pd, px, rv, pmt, pob[3], pco[3];
    if (!PyArg_ParseTuple(args, "ddddddd(ddd)",
                                 &rc, &dc, &pr, &pd, &px, &rv, &pmt,
                                 &pob[0], &pob[1], &pob[2]))
        return NULL;
    eraPmpx(rc, dc, pr, pd, px, rv, pmt, pob, pco);
    return Py_BuildValue("(ddd)", pco[0], pco[1], pco[2]);
}
'''),
    ('pmsafe', 'METH_VARARGS',
     '''"\\npmsafe(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b, -> ra2, dec2, pmr2, pmd2, px2, rv2)\\n"
"Star proper motion:  update star catalog data for space motion, with\\n"
"special handling to handle the zero parallax case.\\n"
"Given:\\n"
"    ra1    right ascension (radians), before\\n"
"    dec1   declination (radians), before\\n"
"    pmr1   RA proper motion (radians/year), before\\n"
"    pmd1   Dec proper motion (radians/year), before\\n"
"    px1    parallax (arcseconds), before\\n"
"    rv1    radial velocity (km/s, +ve = receding), before\\n"
"    ep1a   ''before'' epoch, part A\\n"
"    ep1b   ''before'' epoch, part B\\n"
"    ep2a   ''after'' epoch, part A\\n"
"    ep2b   ''after'' epoch, part B\\n"
"Returned:\\n"
"    ra2    double      right ascension (radians), after\\n"
"    dec2   double      declination (radians), after\\n"
"    pmr2   double      RA proper motion (radians/year), after\\n"
"    pmd2   double      Dec proper motion (radians/year), after\\n"
"    px2    double      parallax (arcseconds), after\\n"
"    rv2    double      radial velocity (km/s, +ve = receding), after"''',
     '''{
    double ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b;
    double ra2, dec2, pmr2, pmd2, px2, rv2;
    int j;
    if (!PyArg_ParseTuple(args, "dddddddddd",
                                 &ra1, &dec1, &pmr1, &pmd1, &px1,
                                 &rv1, &ep1a, &ep1b, &ep2a, &ep2b))
        return NULL;
    j = eraPmsafe(ra1, dec1, pmr1, pmd1, px1, rv1,
                  ep1a, ep1b, ep2a, ep2b,
                  &ra2, &dec2, &pmr2, &pmd2, &px2, &rv2);
    if (j == -1) {
        PyErr_SetString(_erfaError, "system error (should not occur)");
        return NULL;
    }
    else if (j == 1) {
        PyErr_SetString(_erfaError, "distance overridden");
        return NULL;
    }
    else if (j == 2) {
        PyErr_SetString(_erfaError, "excessive velocity");
        return NULL;
    }
    else if (j == 4) {
        PyErr_SetString(_erfaError, "solution didn't converge");
        return NULL;
    }
    return Py_BuildValue("dddddd", ra2, dec2, pmr2, pmd2, px2, rv2);
}
'''),
    ('pvtob', 'METH_VARARGS',
     '''"\\npvtob(elong, phi, hm, xp, yp, sp, theta) -> pv[2][3]\\n"
"Position and velocity of a terrestrial observing station.\\n"
"Given:\\n"
"    elong   double       longitude (radians, east +ve)\\n"
"    phi     double       latitude (geodetic, radians))\\n"
"    hm      double       height above ref. ellipsoid (geodetic, m)\\n"
"    xp,yp   double       coordinates of the pole (radians))\\n"
"    sp      double       the TIO locator s' (radians))\\n"
"    theta   double       Earth rotation angle (radians)\\n"
"Returned:\\n"
"    pv      position/velocity vector (m, m/s, CIRS)"''',
     '''{
    double elong, phi, hm, xp, yp, sp, theta, pv[2][3];
    if (!PyArg_ParseTuple(args, "ddddddd",
                                 &elong, &phi, &hm, &xp, &yp, &sp, &theta))
        return NULL;
    eraPvtob(elong, phi, hm, xp, yp, sp, theta, pv);
    return Py_BuildValue("((ddd)(ddd))",
                          pv[0][0], pv[0][1], pv[0][2],
                          pv[1][0], pv[1][1], pv[1][2]);
}
'''),
    ('refco', 'METH_VARARGS',
     '''"\\nrefco(phpa, tc, rh, wl) -> refa, refb\\n"
"Determine the constants A and B in the atmospheric refraction model\\n"
"dZ = A tan Z + B tan^3 Z.\\n"
"\\n"
"Z is the ''observed'' zenith distance (i.e. affected by refraction)\\n"
"and dZ is what to add to Z to give the ''topocentric'' (i.e. in vacuo)\\n"
"zenith distance.\\n"
"Given:\\n"
"   phpa   pressure at the observer (hPa = millibar)\\n"
"   tc     ambient temperature at the observer (deg C)\\n"
"   rh     relative humidity at the observer (range 0-1)\\n"
"   wl     wavelength (micrometers)\\n"
"Returned:\\n"
"   refa   tan Z coefficient (radians)\\n"
"   refb   tan^3 Z coefficient (radians)"''',
     '''{
    double phpa, tc, rh, wl;
    double refa, refb;
    if (!PyArg_ParseTuple(args, "dddd", &phpa, &tc, &rh, &wl))
        return NULL;
    eraRefco(phpa, tc, rh, wl, &refa, &refb);
    return Py_BuildValue("dd", refa, refb);
}
'''),

# astronomy library
    ('bi00', 'METH_NOARGS',
     '''"\\nbi00() -> dpsibi,depsbi,dra\\n"
"Frame bias components of IAU 2000 precession-nutation models (part of MHB2000 with additions).\\n"
"Returned:\\n"
"    dpsibi,depsbi    obliquity and correction\\n"
"    dra              the ICRS RA of the J2000.0 mean equinox"''',
     '''{
    double dpsibi, depsbi, dra;
    eraBi00(&dpsibi, &depsbi, &dra);
    if (!PyErr_Occurred()) {
        return Py_BuildValue("ddd", dpsibi, depsbi, dra);
    }
    else {
        return NULL;
    }
}
'''),
    ('bp00', 'METH_VARARGS',
     '''"\\nbp00(d1, d2) -> rb, rp, rbp\\n\\n"
"Frame bias and precession, IAU 2000.\\n"
"Given:\\n"
"    d1, d2     TT as a 2-part Julian Date\\n"
"Returned:\\n"
"    rb         frame bias matrix\\n"
"    rp         precession matrix\\n"
"    rbp        bias-precession matrix"''',
     '''{
    double d1, d2, rb[3][3], rp[3][3], rbp[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        eraBp00(d1, d2, rb, rp, rbp);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd)),((ddd)(ddd)(ddd)),((ddd)(ddd)(ddd))",
        rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
        rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
        rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2]
        );
    }
    else {
        return NULL;
    }
}
'''),
    ('bp06', 'METH_VARARGS',
     '''"\\nbp06(d1, d2) -> rb, rp, rbp\\n\\n"
"Frame bias and precession, IAU 2006.\\n"
"Given:\\n"
"    d1, d2     TT as 2-part Julian Date\\n"
"Returned:\\n"
"    rb         frame bias matrix\\n"
"    p          precession matrix)\\n"
"    rbp        bias-precession matrix"''',
     '''{
    double d1, d2, rb[3][3], rp[3][3], rbp[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        eraBp06(d1, d2, rb, rp, rbp);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd)),((ddd)(ddd)(ddd)),((ddd)(ddd)(ddd))",
        rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
        rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
        rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2]
        );
    }
    else {
        return NULL;
    }
}
'''),
    ('bpn2xy', 'METH_VARARGS',
     '''"\\nbpn2xy(rbpn) -> x, y\\n\\n"
"Extract from the bias-precession-nutation matrix\\n"
"the X,Y coordinates of the Celestial Intermediate Pole.\\n"
"Given:\\n"
"    rbpn       celestial-to-true matrix\\n"
"Returned:\\n"
"    x,y        celestial Intermediate Pole"''',
     '''{
    double x, y, rbpn[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))",
                                   &rbpn[0][0],&rbpn[0][1],&rbpn[0][2],
                                   &rbpn[1][0],&rbpn[1][1],&rbpn[1][2],
                                   &rbpn[2][0],&rbpn[2][1],&rbpn[2][2]);
    if (ok) {
        eraBpn2xy(rbpn, &x, &y);
        return Py_BuildValue("dd", x, y);
    }
    else {
        return NULL;
    }
}
'''),
    ('c2i00a', 'METH_VARARGS',
     '''"\\nc2i00a(d1, d2) -> rc2i\\n\\n"
"Form the celestial-to-intermediate matrix\\n"
"for a given date using the IAU 2000A precession-nutation model.\\n"
"Given:\\n"
"    d1, d2     TT as 2-part Julian date\\n"
"Returned:\\n"
"    rc2i       celestial-to-intermediate matrix\\n"''',
     '''{
    double d1, d2, rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        eraC2i00a(d1, d2, rc2i);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2i[0][0],rc2i[0][1],rc2i[0][2],
        rc2i[1][0],rc2i[1][1],rc2i[1][2],
        rc2i[2][0],rc2i[2][1],rc2i[2][2]
        );
    }
    else {
        return NULL;
    }
}
'''),
    ('c2i00b', 'METH_VARARGS',
     '''"\\nc2i00b(d1, d2) -> rc2i\\n\\n"
"Form the celestial-to-intermediate matrix for a given\\n"
"date using the IAU 2000B precession-nutation model.\\n"
"Given:\\n"
"    d1,d2      TT as 2-part Julian date\\n"
"Returned:\\n"
"    rc2i       celestial-to-intermediate matrix"''',
     '''{
    double d1, d2, rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        eraC2i00b(d1, d2, rc2i);
        return Py_BuildValue("((ddd)(ddd)(ddd))",
                                rc2i[0][0],rc2i[0][1],rc2i[0][2],
                                rc2i[1][0],rc2i[1][1],rc2i[1][2],
                                rc2i[2][0],rc2i[2][1],rc2i[2][2]);
    }
    else {
        return NULL;
    }
}
'''),
    ('c2i06a', 'METH_VARARGS',
     '''"\\nc2i06a(d1, d2) -> rc2i\\n\\n"
"Form the celestial-to-intermediate matrix for a given date \\n"
"using the IAU 2006 precession and IAU 200A nutation model.\\n"
"Given:\\n"
"    d1,d2      TT as 2-part Julian date\\n"
"Returned:\\n"
"    rc2i       celestial-to-intermediate matrix"''',
     '''{
    double d1, d2, rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        eraC2i06a(d1, d2, rc2i);
        return Py_BuildValue("((ddd)(ddd)(ddd))",
                                rc2i[0][0],rc2i[0][1],rc2i[0][2],
                                rc2i[1][0],rc2i[1][1],rc2i[1][2],
                                rc2i[2][0],rc2i[2][1],rc2i[2][2]);
    }
    else {
        return NULL;
    }
}
'''),
    ('c2ibpn', 'METH_VARARGS',
     '''"\\nc2ibpn(d1, d2, rbpn) -> rc2i\\n\\n"
"Form the celestial-to-intermediate matrix for a given date\\n"
"and given the bias-precession-nutation matrix. IAU 2000.\\n"
"Given:\\n"
"    d1,d2      TT as 2-part Julian date\\n"
"    rbpn       celestial-to-true matrix\\n"
"Returned:\\n"
"    rc2i       celestial-to-intermediate matrix"''',
     '''{
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    double d1, d2, rbpn[3][3], rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd((ddd)(ddd)(ddd))", &d1, &d2,
    &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22);
    if (ok) {
        rbpn[0][0] = r00;
        rbpn[0][1] = r01;
        rbpn[0][2] = r02;
        rbpn[1][0] = r10;
        rbpn[1][1] = r11;
        rbpn[1][2] = r12;
        rbpn[2][0] = r20;
        rbpn[2][1] = r21;
        rbpn[2][2] = r22;
        eraC2ibpn(d1, d2, rbpn, rc2i);
        return Py_BuildValue("((ddd)(ddd)(ddd))",
                                rc2i[0][0],rc2i[0][1],rc2i[0][2],
                                rc2i[1][0],rc2i[1][1],rc2i[1][2],
                                rc2i[2][0],rc2i[2][1],rc2i[2][2]
        );
    }
    else {
        return NULL;
    }
}
'''),
    ('c2ixy', 'METH_VARARGS',
     '''"\\nc2ixy(d1, d2, x, y) -> rc2i\\n\\n"
"Form the celestial to intermediate-frame-of-date matrix for a given\\n"
"date when the CIP X,Y coordinates are known. IAU 2000.\\n"
"Given:\\n"
"    d1,d2      TT as 2-part Julian Date\\n"
"    x, y       Celestial Intermediate Pole\\n"
"Returned:\\n"
"    rc2i       celestial-to-intermediate matrix"''',
     '''{
    double d1, d2, x, y, rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddd", &d1, &d2, &x, &y);
    if (ok) {
        eraC2ixy(d1,d2,x,y,rc2i);
        return Py_BuildValue("((ddd)(ddd)(ddd))",
                                rc2i[0][0],rc2i[0][1],rc2i[0][2],
                                rc2i[1][0],rc2i[1][1],rc2i[1][2],
                                rc2i[2][0],rc2i[2][1],rc2i[2][2]
        );    
    }
    else {
        return NULL;
    }
}
'''),
    ('c2ixys', 'METH_VARARGS',
     '''"\\nc2ixys(x, y, s) -> rc2i\\n\\n"
"Form the celestial to intermediate-frame-of-date matrix\\n"
" given the CIP X,Y and the CIO locator s.\\n"
"Given:\\n"
"    x, y       Celestial Intermediate Pole\\n"
"    s          CIO locator \\n"
"Returned:\\n"
"   rc2i        celestial-to-intermediate matrix"''',
     '''{
    double x, y, s, rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "ddd", &x, &y, &s);
    if (ok) {
        eraC2ixys(x,y,s,rc2i);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2i[0][0],rc2i[0][1],rc2i[0][2],
        rc2i[1][0],rc2i[1][1],rc2i[1][2],
        rc2i[2][0],rc2i[2][1],rc2i[2][2]
        );    
    }
    else {
        return NULL;
    }
}
'''),
    ('c2t00a', 'METH_VARARGS',
     '''"\\nc2t00a(tta,ttb,uta,utb,xp,yp) -> rc2t\\n\\n"
"Form the celestial to terrestrial matrix given the date, the UT1 and\\n"
"the polar motion, using the IAU 2000A nutation model.\\n"
"Given:\\n"
"    tta,ttb    TT as 2-part Julian Date\\n"
"    uta,utb    UT1 as 2-part Julian Date\\n"
"    xp, yp     coordinates of the pole (radians)\\n"
"Returned:\\n"
"    rc2t       celestial-to-terrestrial matrix"''',
     '''{
    double tta,ttb,uta,utb,xp,yp,rc2t[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddddd", &tta,&ttb,&uta,&utb,&xp,&yp);
    if (ok) {
        eraC2t00a(tta,ttb,uta,utb,xp,yp,rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        ); 
    }
    else {
        return NULL;
    }
}
'''),
    ('c2t00b', 'METH_VARARGS',
     '''"\\nc2t00b(tta,ttb,uta,utb,xp,yp) -> rc2t\\n\\n"
"Form the celestial to terrestrial matrix given the date, the UT1 and\\n"
"the polar motion, using the IAU 2000B nutation model.\\n"
"Given:\\n"
"    tta,ttb    TT as 2-part Julian Date\\n"
"    uta,utb    UT1 as 2-part Julian Date\\n"
"    xp, yp     coordinates of the pole (radians)\\n"
"Returned:\\n"
"    rc2t       celestial-to-terrestrial matrix"''',
     '''{
    double tta,ttb,uta,utb,xp,yp,rc2t[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddddd", &tta,&ttb,&uta,&utb,&xp,&yp);
    if (ok) {
        eraC2t00b(tta,ttb,uta,utb,xp,yp,rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        ); 
    }
    else {
        return NULL;
    }
}
'''),
    ('c2t06a', 'METH_VARARGS',
     '''"\\nc2t06a(tta,ttb,uta,utb,xp,yp) -> rc2t\\n\\n"
"Form the celestial to terrestrial matrix given the date, the UT1 and\\n"
"the polar motion, using the IAU 2006 precession and IAU 2000A nutation models.\\n"
"Given:\\n"
"    tta,ttb    TT as 2-part Julian Date\\n"
"    uta,utb    UT1 as 2-part Julian Date\\n"
"    xp, yp     coordinates of the pole (radians)\\n"
"Returned:\\n"
"    rc2t       celestial-to-terrestrial matrix"
''',
     '''{
    double tta,ttb,uta,utb,xp,yp,rc2t[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddddd", &tta,&ttb,&uta,&utb,&xp,&yp);
    if (ok) {
        eraC2t06a(tta,ttb,uta,utb,xp,yp,rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        ); 
    }
    else {
        return NULL;
    }
}
'''),
    ('c2tcio', 'METH_VARARGS',
     '''"\\nc2tcio(rc2i, era, rpom) -> rc2t\\n\\n"
"Assemble the celestial to terrestrial matrix from CIO-based components\\n"
"(the celestial-to-intermediate matrix, the Earth Rotation Angle and the polar motion matrix)\\n"
"Given:\\n"
"    rc2        celestial-to-intermediate matrix\\n"
"    era        Earth rotation angle\\n"
"    rpom       polar-motion matrix\\n"
"Returned:t\\n"
"    rc2t       celestial-to-terrestrial matrix"''',
     '''{
    double c00,c01,c02,c10,c11,c12,c20,c21,c22,rc2i[3][3];
    double era, rc2t[3][3];
    double p00,p01,p02,p10,p11,p12,p20,p21,p22,rpom[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))d((ddd)(ddd)(ddd))",
    &c00,&c01,&c02,&c10,&c11,&c12,&c20,&c21,&c22,
    &era,
    &p00,&p01,&p02,&p10,&p11,&p12,&p20,&p21,&p22);
    if (ok) {
        rc2i[0][0] = c00;
        rc2i[0][1] = c01;
        rc2i[0][2] = c02;
        rc2i[1][0] = c10;
        rc2i[1][1] = c11;
        rc2i[1][2] = c12;
        rc2i[2][0] = c20;
        rc2i[2][1] = c21;
        rc2i[2][2] = c22;
        rpom[0][0] = p00;
        rpom[0][1] = p01;
        rpom[0][2] = p02;
        rpom[1][0] = p10;
        rpom[1][1] = p11;
        rpom[1][2] = p12;
        rpom[2][0] = p20;
        rpom[2][1] = p21;
        rpom[2][2] = p22;
        eraC2tcio(rc2i, era, rpom, rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        );
    }
    else {
        return NULL;
    }
}
'''),
    ('c2teqx', 'METH_VARARGS',
     '''"\\nc2teqx(rbpn, gst, rpom -> rc2t\\n\\n"
"Assemble the celestial to terrestrial matrix from equinox-based\\n"
"components (the celestial-to-true matrix, the Greenwich Apparent\\n"
"Sidereal Time and the polar motion matrix).\\n"
"Given:\\n"
"    rbpn       celestial-to-true matrix\\n"
"    gst        Greenwich (apparent) Sidereal Time\\n"
"    rpom       polar-motion matrix\\n"
"Returned:\\n"
"    rc2t       celestial-to-terrestrial matrix"''',
     '''{
    double c00,c01,c02,c10,c11,c12,c20,c21,c22,rc2i[3][3];
    double gst, rc2t[3][3];
    double p00,p01,p02,p10,p11,p12,p20,p21,p22,rpom[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))d((ddd)(ddd)(ddd))",
    &c00,&c01,&c02,&c10,&c11,&c12,&c20,&c21,&c22,
    &gst,
    &p00,&p01,&p02,&p10,&p11,&p12,&p20,&p21,&p22);
    if (ok) {
        rc2i[0][0] = c00;
        rc2i[0][1] = c01;
        rc2i[0][2] = c02;
        rc2i[1][0] = c10;
        rc2i[1][1] = c11;
        rc2i[1][2] = c12;
        rc2i[2][0] = c20;
        rc2i[2][1] = c21;
        rc2i[2][2] = c22;
        rpom[0][0] = p00;
        rpom[0][1] = p01;
        rpom[0][2] = p02;
        rpom[1][0] = p10;
        rpom[1][1] = p11;
        rpom[1][2] = p12;
        rpom[2][0] = p20;
        rpom[2][1] = p21;
        rpom[2][2] = p22;
        eraC2teqx(rc2i, gst, rpom, rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        );
    }
    else {
        return NULL;
    }
}
'''),
    ('c2tpe', 'METH_VARARGS',
     '''"\\nc2tpe(tta,ttb,uta,utb,dpsi,deps,xp,yp) -> rc2t\\n\\n"
"Form the celestial to terrestrial matrix given the date, the UT1,\\n"
"the nutation and the polar motion.  IAU 2000.\\n"
"Given:\\n"
"    tta,ttb    TT as 2-part Julian Date\\n"
"    uta,utb    UT1 as 2-part Julian Date\\n"
"    dpsi,deps  nutation\\n"
"    xp,yp      coordinates of the pole\\n"
"Returned:\\n"
"    rc2t       celestial-to-terrestrial matrix"''',
     '''{
    double tta,ttb,uta,utb,dpsi,deps,xp,yp,rc2t[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddddddd", &tta,&ttb,&uta,&utb,&dpsi,&deps,&xp,&yp);
    if (ok) {
        eraC2tpe(tta,ttb,uta,utb,dpsi,deps,xp,yp,rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        ); 
    }
    else {
        return NULL;
    }
}
'''),
    ('c2txy', 'METH_VARARGS',
     '''"\\nc2txy(tta,ttb,uta,utb,x,y,xp,yp) -> rc2t\\n\\n"
"Form the celestial to terrestrial matrix given the date, the UT1,\\n"
"the CIP coordinates and the polar motion.  IAU 2000."
"Given:\\n"
"    tta,ttb    TT as 2-part Julian Date\\n"
"    uta,utb    UT1 as 2-part Julian Date\\n"
"    x,y        Celestial Intermediate Pole\\n"
"    xp,yp      coordinates of the pole\\n"
"Returned:\\n"
"    rc2t       celestial-to-terrestrial matrix"''',
     '''{
    double tta,ttb,uta,utb,x,y,xp,yp,rc2t[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddddddd", &tta,&ttb,&uta,&utb,&x,&y,&xp,&yp);
    if (ok) {
        eraC2txy(tta,ttb,uta,utb,x,y,xp,yp,rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        ); 
    }
    else {
        return NULL;
    }
}
'''),
    ('cal2jd', 'METH_VARARGS',
     '''"\\ncal2jd(year, month, day) -> 2400000.5,djm\\n\\n"
"Gregorian Calendar to Julian Date.\\n"
"Given:\\n"
"    year, month      day in Gregorian calendar\\n"
"Returned:\\n"
"    2400000.5,djm    MJD zero-point and Modified Julian Date for 0 hrs"''',
     '''{
    int iy, im, id;
    double dmj0, dmj;
    int ok, status;
    ok = PyArg_ParseTuple(args, "iii", &iy, &im, &id);
    if (ok) {
        status = eraCal2jd(iy, im, id, &dmj0, &dmj);
    }
    else {
        return NULL;
    }
    if (status < 0){
        if (status == -1){
            PyErr_SetString(_erfaError, "bad year");
            return NULL;
        }
        else if (status == -2){
            PyErr_SetString(_erfaError, "bad month");
            return NULL;
        }
        else if (status == -3){
            PyErr_SetString(_erfaError, "bad day");
            return NULL;
        }
    }
    return Py_BuildValue("dd", dmj0, dmj);
}
'''),
    ('d2dtf', 'METH_VARARGS',
     '''"\\nd2dtf(n, d1, d2 | scale) -> year, month, day, hours, minutes, seconds, fraction\\n\\n"
"Format for output a 2-part Julian Date (or in the case of UTC a\\n"
"quasi-JD form that includes special provision for leap seconds).\\n"
"Given:\\n"
"    ndp        resolution\\n"
"    d1,d2      time as a 2-part Julian Date\\n"
"    scale      optional time scale ID, default to ''UTC''\\n"
"Returned:\\n"
"   iy,im,id    year, month, day in Gregorian calendar\\n"
"   ihmsf       hours, minutes, seconds, fraction\\n"

''',
     '''{
    double d1, d2;
    int n, ok, status;
    int iy, im, id, ihmsf[4];
    char *scale = "UTC";
    ok = PyArg_ParseTuple(args, "idd|s", &n, &d1, &d2, &scale);
    if (ok) {
        status = eraD2dtf(scale, n, d1, d2, &iy, &im, &id, ihmsf);
    }
    else {
        return NULL;
    }
    if (status > 0){
        PyErr_SetString(_erfaError, "doubious year: date predates UTC or too far in the future");
        return NULL;
    }
    if (status < 0){
        PyErr_SetString(_erfaError, "unaceptable date");
        return NULL;
    }
    return Py_BuildValue("iiiiiii", iy,im,id,ihmsf[0],ihmsf[1],ihmsf[2],ihmsf[3]);
}
'''),
    ('dat', 'METH_VARARGS',
     '''"\\ndat(y,m,d,fd) -> delta t\\n\\n"
"For a given UTC date, calculate delta(AT) = TAI-UTC.\\n"
"Given:\\n"
"    y          year\\n"
"    m          month\\n"
"    d          day\\n"
"    fd         fraction of day\\n"
"Returned:\\n"
"    deltat     TAI minus UTC, seconds"''',
     '''{
    int y, m, d, status;
    double fd, deltat;
    if (! PyArg_ParseTuple(args, "iiid", &y, &m, &d, &fd))
        return NULL;
    status = eraDat(y,m,d,fd,&deltat);
    if (status > 0){
        PyErr_SetString(_erfaError, "doubious year: date before UTC:1960 January 1.0.");
        return NULL;
    }
    else if (status < 0){
        if (status == -1){
            PyErr_SetString(_erfaError, "unaceptable date, bad year");
            return NULL;
        }
        else if (status == -2){
            PyErr_SetString(_erfaError, "unaceptable date, bad month");
            return NULL;
        }
        else if (status == -3){
            PyErr_SetString(_erfaError, "unaceptable date, bad day");
            return NULL;
        }      
        else if (status == -4){
            PyErr_SetString(_erfaError, "bad fraction day, should be < 1.");
            return NULL;
        }      
    }
    return PyFloat_FromDouble(deltat);
}
'''),
    ('dtdb', 'METH_VARARGS',
     '''"\\ndtdb(d1, d2, ut, elon, u, v) -> TDB-TT\\n\\n"
"An approximation to TDB-TT, the difference between barycentric\\n"
"dynamical time and terrestrial time, for an observer on the Earth.\\n"
"Given:\\n"
"    d1,d2      TDB date as 2-part Julian Date\\n"
"    ut         UT1 universal time as fraction of one day\\n"
"    elong      longitude (east positive, radians)\\n"
"    u          distance from Earth spin axis (km)\\n"
"    v          distance north of equatorial plane (km)"''',
     '''{
    double tdbtt, jd1, jd2, ut1, elon, u, v;
    if (! PyArg_ParseTuple(args, "dddddd", &jd1, &jd2, &ut1, &elon, &u, &v))
        return NULL;
    tdbtt = eraDtdb(jd1, jd2, ut1, elon, u, v);
    return PyFloat_FromDouble(tdbtt);
}
'''),
    ('dtf2d', 'METH_VARARGS',
     '''"\\ndtf2d(y,m,d,hr,mn,sec |scale) -> d1,d2\\n\\n"
"Encode date and time fields into 2-part Julian Date (or in the case\\n"
"of UTC a quasi-JD form that includes special provision for leap seconds).\\n"
"Given:\\n"
"    y,m,d   year, month, day in Gregorian calendar\\n"
"    hr,mn   hour, minute\\n"
"    sec     seconds\\n"
"    scale   optional time scale ID, default ''UTC''\\n"
"Returned:\\n"
"   d1,d2     2-part Julian Date"''',
     '''{
    int iy, im, id, ihr, imn, ok, status;
    double sec, dj1, dj2;
    char *scale = "UTC";
    ok = PyArg_ParseTuple(args, "iiiiid|s", &iy,&im,&id,&ihr,&imn,&sec,&scale);
    if (ok) {
        status = eraDtf2d(scale, iy, im, id, ihr, imn, sec, &dj1, &dj2);
        if (status == 3) {
            PyErr_SetString(_erfaError, "both of next two");
            return NULL;
        }
        else if (status == 2) {
            PyErr_SetString(_erfaError, "time is after end of day");
            return NULL;
        }
        else if (status == 1) {
            PyErr_SetString(_erfaError, "dubious year");
            return NULL;
        }
        else if (status == -1) {
            PyErr_SetString(_erfaError, "bad year");
            return NULL;
        }
        else if (status == -2) {
            PyErr_SetString(_erfaError, "bad month");
            return NULL;
        }
        else if (status == -3) {
            PyErr_SetString(_erfaError, "bad day");
            return NULL;
        }
        else if (status == -4) {
            PyErr_SetString(_erfaError, "bad hour");
            return NULL;
        }
        else if (status == -5) {
            PyErr_SetString(_erfaError, "bad minute");
            return NULL;
        }
        else if (status == -6) {
            PyErr_SetString(_erfaError, "bad second (<0)");
            return NULL;
        }
        else {
            return  Py_BuildValue("dd",dj1, dj2);
        }
    }
    else {
        return NULL;
    }
}
'''),
    ('ee00', 'METH_VARARGS',
     '''"\\nee00(d1,d2,epsa,dpsi) -> ee\\n\\n"
"The equation of the equinoxes compatible with IAU 2000 resolutions,\\n"
"given the nutation in longitude and the mean obliquity.\\n"
"Given:\\n"
"    d1,d2      TT as 2-part Julian Date\\n"
"    epsa       mean obliquity\\n"
"    dpsi       nutation in longitude\\n"
"Returned:\\n"
"    ee         equation of the equinoxes"''',
     '''{
    double d1,d2,epsa,dpsi,ee;
    if (!PyArg_ParseTuple(args, "dddd", &d1,&d2,&epsa,&dpsi)) {
        return NULL;
    }
    ee = eraEe00(d1,d2,epsa,dpsi);
    return PyFloat_FromDouble(ee);
}
'''),
    ('ee00a', 'METH_VARARGS',
     '''"\\nee00a(d1,d2) -> ee\\n\\n"
"equation of the equinoxes, compatible with IAU 2000 resolutions.\\n"
"Given:\\n"
"    d1,d2      TT as 2-part Julian Date\\n"
"Returned:\\n"
"    ee         equation of the equinoxes"''',
     '''{
    double d1,d2,ee;
    if (!PyArg_ParseTuple(args, "dd", &d1,&d2)) {
        return NULL;
    }
    ee = eraEe00a(d1,d2);
    return Py_BuildValue("d",ee); //PyFloat_FromDouble(ee);
}
'''),
    ('ee00b', 'METH_VARARGS',
     '''"\\nee00b(d1,d2) -> ee\\n\\n"
"Equation of the equinoxes, compatible with IAU 2000 resolutions\\n"
"but using the truncated nutation model IAU 2000B.\\n"
"Given:\\n"
"    d1,d2      TT as 2-part Julian Date\\n"
"Returned:\\n"
"    ee         equation of the equinoxes"''',
     '''{
    double d1,d2,ee;
    if (!PyArg_ParseTuple(args, "dd", &d1,&d2)) {
        return NULL;
    }
    ee = eraEe00b(d1,d2);
    return Py_BuildValue("d",ee); //PyFloat_FromDouble(ee);
}
'''),
    ('ee06a', 'METH_VARARGS',
     '''"\\nee06a(d1,d2) -> ee\\n\\n"
"Equation of the equinoxes, compatible with IAU 2000 resolutions and \\n"
"IAU 2006/2000A precession-nutation.\\n"
"Given:\\n"
"    d1,d2      TT as 2-part Julian Date\\n"
"Returned:\\n"
"    ee         equation of the equinoxes"''',
     '''{
    double d1,d2,ee;
    if (!PyArg_ParseTuple(args, "dd", &d1,&d2)) {
        return NULL;
    }
    ee = eraEe06a(d1,d2);
    return Py_BuildValue("d",ee); //PyFloat_FromDouble(ee);
}
'''),
    ('eect00', 'METH_VARARGS',
     '''"\\neect00(d1,d2) -> ct\\n\\n"
"Equation of the equinoxes complementary terms, consistent with IAU 2000 resolutions.\\n"
"Given:\\n"
"    d1,d2      TT as 2-part Julian Date\\n"
"Returned:\\n"
"    ct        complementary terms"''',
     '''{
    double d1,d2,ct;
    if (!PyArg_ParseTuple(args, "dd", &d1,&d2)) {
        return NULL;
    }
    ct = eraEect00(d1,d2);
    return Py_BuildValue("d",ct);
}
'''),
    ('eform', 'METH_VARARGS',
     '''"\\neform(n) -> a, f\\n\\n"
"Earth reference ellipsoids.\\n"
"Given:\\n"
"    n          ellipsoid identifier\\n\\n"
"    n=1 WGS84\\n"
"    n=2 GRS80\\n"
"    n=3 WGS72\\n"
"Returned:\\n"
"    a          equatorial radius (meters)\\n"
"    f          flattening"''',
     '''{
    int n, status;
    double a, f;
    if (!PyArg_ParseTuple(args, "i", &n)) {
        return NULL;
    }
    status = eraEform(n, &a, &f);
    if (status) {
        PyErr_SetString(_erfaError, "illegal identifier; n should be 1,2 or 3");
        return NULL;
    }
    return Py_BuildValue("dd",a,f);
}
'''),
    ('eo06a', 'METH_VARARGS',
     '''"\\neo06a(d1, d2) -> eo\\n"
"equation of the origins, IAU 2006 precession and IAU 2000A nutation.\\n"
"Given:\\n"
"    d1,d2      TT as 2-part Julian Date\\n"
"Returned:\\n"
"    eo         equation of the origins (radians)"''',
     '''{
    double d1, d2, eo;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eo = eraEo06a(d1, d2);
    return Py_BuildValue("d",eo);
}
'''),
    ('eors', 'METH_VARARGS',
     '''"\\neors(rnpb, s) -> eo\\n\\n"
"Equation of the origins, given the classical NPB matrix and the quantity s\\n"
"Given:\\n"
"    rnpb       classical nutation x precession x bias matrix\\n"
"    s          CIO locator\\n"
"Returned:\\n"
"    eo         equation of the origins in radians"''',
     '''{
    double rnpb[3][3], s, eo;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))d",
                                 &rnpb[0][0], &rnpb[0][1], &rnpb[0][2],
                                 &rnpb[1][0], &rnpb[1][1], &rnpb[1][2],
                                 &rnpb[2][0], &rnpb[2][1], &rnpb[2][2], &s)) {
        return NULL;
    }        
    eo = eraEors(rnpb, s);
    return Py_BuildValue("d",eo);
}
'''),
    ('epb', 'METH_VARARGS',
     '''"\\nepb(d1, d2) -> b\\n\\n"
"Julian Date to Besselian Epoch\\n"
"Given:\\n"
"    d1,d2      2-part Julian Date\\n"
"Returned:\\n"
"    b          Besselian Epoch."''',
     '''{
    double d1, d2, b;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    b = eraEpb(d1,d2);
    return Py_BuildValue("d",b);
}
'''),
    ('epb2jd', 'METH_VARARGS',
     '''"\\nepb2jd(epb) -> 2400000.5 djm\\n\\n"
"Given:\\n"
"    epb        Besselian Epoch,\\n"
"Returned:\\n"
"    djm        Modified Julian Date"''',
     '''{
    double epb, jd0, jd1;
    if (!PyArg_ParseTuple(args, "d", &epb)) {
        return NULL;
    }
    eraEpb2jd(epb, &jd0, &jd1);
    return Py_BuildValue("dd",jd0,jd1);
}
'''),
    ('epj', 'METH_VARARGS',
     '''"\\nepj(d1, d2) -> b\\n\\n"
"Julian Date to Julian Epoch.\\n"
"Given:\\n"
"    d1,d2      2-part Julian Date\\n"
"Returned:\\n"
"    b          Julian Epoch"''',
     '''{
    double d1, d2, j;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    j = eraEpj(d1,d2);
    return Py_BuildValue("d",j);
}
'''),
    ('epj2jd', 'METH_VARARGS',
     '''"\\nepj2jd(epj) -> 2400000.5 djm\\n\\n"
"Julian Epoch to Julian Date\\n"
"Given:\\n"
"    epj        Julian Epoch\\n"
"Returned:\\n"
"    djm        Modified Julian Date"''',
     '''{
    double epj, jd0, jd1;
    if (!PyArg_ParseTuple(args, "d", &epj)) {
        return NULL;
    }
    eraEpj2jd(epj, &jd0, &jd1);
    return Py_BuildValue("dd",jd0,jd1);
}
'''),
    ('epv00', 'METH_VARARGS|METH_KEYWORDS',
     '''"\\nepv00(d1,d2 |prec=0) -> pvh, pvb\\n\\n"
"Earth position and velocity, heliocentric and barycentric,\\n"
"with respect to the Barycentric Celestial Reference System.\\n"
"Given:\\n"
"    d1,d2      TDB as 2-part Julian Date\\n"
"    prec!=0    optional, to compute outside the range 1900-2100 AD\\n"
"Returned:\\n"
"    pvh        heliocentric Earth position/velocity\\n"
"    pvb        barycentric Earth position/velocity"''',
     '''{
    static char *kwlist[] = {"dj1", "dj2", "prec", NULL};
    double dj1, dj2;
    int prec=0, status;
    double pvh[2][3] = {{[0]=0.},{[0]=0.}};
    double pvb[2][3] = {{[0]=0.},{[0]=0.}};
    if (! PyArg_ParseTupleAndKeywords(args, kwds, "dd|i", kwlist, &dj1, &dj2, &prec))
        return NULL;
    status = eraEpv00(dj1, dj2, pvh, pvb);
    if (status == 1) {
        if (prec) {
           return Py_BuildValue("((ddd)(ddd))((ddd)(ddd))",
           pvh[0][0], pvh[0][1], pvh[0][2], pvh[1][0], pvh[1][1], pvh[1][2],
           pvb[0][0], pvb[0][1], pvb[0][2], pvb[1][0], pvb[1][1], pvb[1][2]);
        }
        else {
            PyErr_SetString(_erfaError, "date outside the range 1900-2100 AD, set prec!=0 to force computation");
            return NULL;
        }
    }
    else {
       return Py_BuildValue("((ddd)(ddd))((ddd)(ddd))",
       pvh[0][0], pvh[0][1], pvh[0][2], pvh[1][0], pvh[1][1], pvh[1][2],
       pvb[0][0], pvb[0][1], pvb[0][2], pvb[1][0], pvb[1][1], pvb[1][2]);
    }
}
'''),
    ('eqeq94', 'METH_VARARGS',
     '''"\\neqeq94(d1,d2) -> ee\\n\\n"
"Equation of the equinoxes, IAU 1994 model.\\n"
"Given:\\n"
"    d1,d2      TDB as 2-part Julian Date\\n"
"Returned:\\n"
"    ee         equation of the equinoxes"''',
     '''{
    double d1, d2, ee;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    ee = eraEqeq94(d1,d2);
    return Py_BuildValue("d",ee);
}
'''),
    ('era00', 'METH_VARARGS',
     '''"\\nera00(d1,d2) -> era\\n\\n"
"Earth rotation angle (IAU 2000 model).\\n"
"Given:\\n"
"    d1,d2      UT1 as 2-part Julian Date (d1,d2)\\n"
"Returned:\\n"
"    era         Earth rotation angle (radians), range 0-2pi"''',
     '''{
    double d1, d2, era;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }   
    era = eraEra00(d1,d2);
    return Py_BuildValue("d",era);
}
'''),
    ('fad03', 'METH_VARARGS',
     '''"\\nfad03(t) -> d\\n\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean elongation of the Moon from the Sun.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    d          mean elongation of the Moon from the Sun, radians."''',
     '''{
    double t,d;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    d = eraFad03(t);
    return Py_BuildValue("d",d);
}
'''),
    ('fae03', 'METH_VARARGS',
     '''"\\nfae03(t) -> e\\n\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean longitude of Earth.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    e          mean longitude of Earth, radians."''',
     '''{
    double t,e;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    e = eraFae03(t);
    return Py_BuildValue("d",e);
}
'''),
    ('faf03', 'METH_VARARGS',
     '''"\\nfaf03(t) -> f\\n\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean longitude of the Moon minus mean longitude of the ascending node.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    f          mean longitude of the Moon minus mean longitude of the ascending node, radians."''',
     '''{
    double t,f;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    f = eraFaf03(t);
    return Py_BuildValue("d",f);
}
'''),
    ('faju03', 'METH_VARARGS',
     '''"\\nfaju03(t) -> l\\n\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean longitude of Jupiter.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    l          mean longitude of Jupiter., in radians."''',
     '''{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFaju03(t);
    return Py_BuildValue("d",l);
}
'''),
    ('fal03', 'METH_VARARGS',
     '''"\\nfal03(t) -> l\\n\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean anomaly of the Moon.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    l          mean anomaly of the Moon, in radians."''',
     '''{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFal03(t);
    return Py_BuildValue("d",l);
}
'''),
    ('falp03', 'METH_VARARGS',
     '''"\\nfal03(t) -> lp\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean anomaly of the Sun.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned\\n"
"    lp         mean anomaly of the Sun, in radians."''',
     '''{
    double t,lp;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    lp = eraFalp03(t);
    return Py_BuildValue("d",lp);
}
'''),
    ('fama03', 'METH_VARARGS',
     '''"\\nfama03(t) -> l\\n\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean longitude of Mars.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    l          mean longitude of Mars, in radians."''',
     '''{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFama03(t);
    return Py_BuildValue("d",l);
}
'''),
    ('fame03', 'METH_VARARGS',
     '''"\\nfame03(t) -> l\\n\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean longitude of Mercury.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    l          mean longitude of Mercury, in radians."''',
     '''{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFame03(t);
    return Py_BuildValue("d",l);
}
'''),
    ('fane03', 'METH_VARARGS',
     '''"\\nfane03(t) -> l\\n\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean longitude of Neptune.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    l          mean longitude of Neptune, in radians."''',
     '''{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFane03(t);
    return Py_BuildValue("d",l);
}
'''),
    ('faom03', 'METH_VARARGS',
     '''"\\nfaom03(t) -> l\\n\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean longitude of Moon's ascending node.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    l          mean longitude of Moon's ascending node, in radians."''',
     '''{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFaom03(t);
    return Py_BuildValue("d",l);
}
'''),
    ('fapa03', 'METH_VARARGS',
     '''"\\nfapa03(t) -> l\\n\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"general accumulated precession in longitude.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    l          general accumulated precession in longitude, in radians."''',
     '''{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFapa03(t);
    return Py_BuildValue("d",l);
}
'''),
    ('fasa03', 'METH_VARARGS',
     '''"\\nfasa03(t) -> l\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean longitude of Saturn.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    l          mean longitude of Saturn, in radians."''',
     '''{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFasa03(t);
    return Py_BuildValue("d",l);
}
'''),
    ('faur03', 'METH_VARARGS',
     '''"\\nfaur03(t) -> l\\n\\n"
"Fundamental argument, IERS Conventions (2003):\\n"
"mean longitude of Uranus.\\n"
"Given:\\n"
"    t          TDB as Julian centuries since J2000.,\\n"
"Returned:\\n"
"    l          mean longitude of Uranus, in radians."''',
     '''{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFaur03(t);
    return Py_BuildValue("d",l);
}
'''),
    ('fave03', 'METH_VARARGS',
     '''"\\nfaver03(t) -> l\\n\\n"
"Fundamental argument, IERS Conventions (2003)\\n"
"mean longitude of Venus."
"Given:\\n"
"    t          TDB as Julian centuries since J2000.0\\n"
"Returned:\\n"
"    l          mean longitude of Venus, in radians."''',
     '''{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFave03(t);
    return Py_BuildValue("d",l);
}
'''),
    ('fk52h', 'METH_VARARGS',
     '''"\\nfk52h(r5, d5, dr5, dd5, px5,rv5) -> rh, dh, drh, ddh, pxh, rvh)\\n\\n"
"Transform FK5 (J2000.0) star data into the Hipparcos system.\\n"
"Given (all FK5, equinox J2000.0, epoch J2000.0):\\n"
"    r5         RA (radians)\\n"
"    d5         Dec (radians)\\n"
"    dr5        proper motion in RA (dRA/dt, rad/Jyear)\\n"
"    dd5        proper motion in Dec (dDec/dt, rad/Jyear)\\n"
"    px5        parallax (arcsec)\\n"
"    rv5        radial velocity (km/s, positive = receding)\\n"
"Returned (all Hipparcos, epoch J2000.0):\\n"
"    rh         RA (radians)\\n"
"    dh         Dec (radians)\\n"
"    drh        proper motion in RA (dRA/dt, rad/Jyear)\\n"
"    ddh        proper motion in Dec (dDec/dt, rad/Jyear)\\n"
"    pxh        parallax (arcsec)\\n"
"    rvh        radial velocity (km/s, positive = receding)"''',
     '''{
    double r5, d5, dr5, dd5, px5, rv5;
    double rh, dh, drh, ddh, pxh, rvh;
    if (!PyArg_ParseTuple(args, "dddddd", &r5, &d5, &dr5, &dd5, &px5, &rv5)) {
        return NULL;
    }
    eraFk52h(r5, d5, dr5, dd5, px5, rv5,
             &rh, &dh, &drh, &ddh, &pxh, &rvh);
    return Py_BuildValue("dddddd", rh, dh, drh, ddh, pxh, rvh);
}
'''),
    ('fk5hip', 'METH_NOARGS',
     '''"\\nfk5hip() -> r5h, s5h\\n\\n"
"FK5 to Hipparcos rotation and spin.\\n"
"Returned:\\n"
"    r5h        r-matrix: FK5 rotation wrt Hipparcos \\n"
"    s5         r-vector: FK5 spin wrt Hipparcos."''',
     '''{
    double r5h[3][3], s5h[3];
    eraFk5hip(r5h, s5h);
    return Py_BuildValue("((ddd)(ddd)(ddd))(ddd)",
        r5h[0][0],r5h[0][1],r5h[0][2],
        r5h[1][0],r5h[1][1],r5h[1][2],
        r5h[2][0],r5h[2][1],r5h[2][2],
        s5h[0],s5h[1],s5h[2]);
}
'''),
    ('fk5hz', 'METH_VARARGS',
     '''"\\nfk5hz(r5, d5, d1, d2) -> rh, dh\\n\\n"
"Transform an FK5 (J2000.0) star position into the system of the\\n"
"Hipparcos catalogue, assuming zero Hipparcos proper motion.\\n"
"Given:\\n"
"    r5         FK5 RA (radians), equinox J2000.0, at date d1,d2\\n"
"    d5         FK5 Dec (radians), equinox J2000.0, at date d1,d2\\n"
"    d1,d2      TDB date \\n"
"Returned:\\n"
"    rh         Hipparcos RA (radians)\\n"
"    dh         Hipparcos Dec (radians)"''',
     '''{
    double r5,d5,d1,d2,rh,dh;
    if (!PyArg_ParseTuple(args, "dddd", &r5, &d5, &d1, &d2)) {
        return NULL;
    }
    eraFk5hz(r5,d5,d1,d2,&rh,&dh);
    return Py_BuildValue("dd", rh, dh);
}
'''),
    ('fw2m', 'METH_VARARGS',
     '''"\\nfw2m(gamb, phib, psi, eps) -> r\\n\\n"
"Form rotation matrix given the Fukushima-Williams angles.\\n"
"Given:\\n"
"    gamb       F-W angle gamma_bar (radians)\\n"
"    phib       F-W angle phi_bar (radians)\\n"
"    si         F-W angle psi (radians)\\n"
"    eps        F-W angle epsilon (radians)\\n"
"Returned:\\n"
"    r          rotation matrix"''',
     '''{
    double gamb, phib, psi, eps, r[3][3];
    if (!PyArg_ParseTuple(args, "dddd", &gamb, &phib, &psi, &eps)) {
        return NULL;
    }
    eraFw2m(gamb, phib, psi, eps, r);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
        r[0][0],r[0][1],r[0][2],
        r[1][0],r[1][1],r[1][2],
        r[2][0],r[2][1],r[2][2]);
}
'''),
    ('fw2xy', 'METH_VARARGS',
     '''"\\nfw2xy(gamb, phib, psi, eps) -> x, y\\n\\n"
"Form CIP X,Y given Fukushima-Williams bias-precession-nutation angles.\\n"
"Given:\\n"
"    gamb       F-W angle gamma_bar (radians)\\n"
"    phib       F-W angle phi_bar (radians)\\n"
"    psi        F-W angle psi (radians)\\n"
"    eps        F-W angle epsilon (radians)\\n"
"Returned:\\n"
"    x,y        CIP X,Y (radians)"''',
     '''{
    double gamb, phib, psi, eps, x, y;
    if (!PyArg_ParseTuple(args, "dddd", &gamb, &phib, &psi, &eps)) {
        return NULL;
    }
    eraFw2xy(gamb, phib, psi, eps, &x, &y);
    return Py_BuildValue("dd", x, y);
}
'''),
    ('gc2gd', 'METH_VARARGS',
     '''"\\ngc2gd(n, xyz[3]) -> elong, phi, height\\n\\n"
"Transform geocentric coordinates to geodetic\\n"
"using the specified reference ellipsoid.\\n"
"Given:\\n"
"    n          ellipsoid identifier\\n"
"    xyz        geocentric vector\\n"
"Returned:\\n"
"    elong      longitude (radians, east +ve)\\n"
"    phi        latitude (geodetic, radians)\\n"
"    height     height above ellipsoid (geodetic)"''',
    '''{
    double xyz[3], elong, phi, height;
    int n, status;
    if (!PyArg_ParseTuple(args, "n(ddd)", &n, &xyz[0], &xyz[1], &xyz[2])) {
        return NULL;
    }
    status = eraGc2gd(n, xyz, &elong, &phi, &height);
    if (status == -1) {
        PyErr_SetString(_erfaError, "illegal identifier; n should be 1,2 or 3");
        return NULL;
    }
    else if (status == -2) {
        PyErr_SetString(_erfaError, "internal error");
        return NULL;
    }
    else {
        return Py_BuildValue("ddd", elong, phi, height);
    }
}
'''),
    ('gc2gde', 'METH_VARARGS',
     '''"\\ngc2gde(a, f, xyz) -> elong, phi, height\\n\\n"
"Transform geocentric coordinates to geodetic\\n"
" for a reference ellipsoid of specified form.\\n"
"Given:\\n"
"    a          equatorial radius\\n"
"    f          flattening\\n"
"    xyz        geocentric vector\\n"
"Returned:\\n"
"    elong      longitude (radians, east +ve)\\n"
"    phi        latitude (geodetic, radians)\\n"
"    height     height above ellipsoid (geodetic)"''',
     '''{
    double xyz[3], a, f, elong, phi, height;
    int status;
    if (!PyArg_ParseTuple(args, "dd(ddd)", &a, &f, &xyz[0], &xyz[1], &xyz[2])) {
        return NULL;
    }    
    status = eraGc2gde(a, f, xyz, &elong, &phi, &height);
    if (status == -1) {
        PyErr_SetString(_erfaError, "illegal f");
        return NULL;
    }
    else if (status == -2) {
        PyErr_SetString(_erfaError, "illegal a");
        return NULL;
    }
    else {
        return Py_BuildValue("ddd", elong, phi, height);
    }
}
'''),
    ('gd2gc', 'METH_VARARGS',
     '''"\\ngd2gc(n, elong, phi, height) -> xyz\\n\\n"
"Transform geodetic coordinates to geocentric\\n"
" using the specified reference ellipsoid.\\n"
"Given:\\n"
"    n          ellipsoid identifier\\n"
"    elong      longitude (radians, east +ve)\\n"
"    phi        latitude (geodetic, radians)\\n"
"    height     height above ellipsoid (geodetic in meters)\\n"
"  Returned:\\n"
"    xyz        geocentric vector in meters"''',
     '''{
    double elong, phi, height, xyz[3];
    int n, status;
    if (!PyArg_ParseTuple(args, "iddd", &n, &elong, &phi, &height)) {
        return NULL;
    }
    status = eraGd2gc(n, elong, phi, height, xyz);
    if (status == -1) {
        PyErr_SetString(_erfaError, "illegal identifier; n should be 1,2 or 3");
        return NULL;
    }
    else if (status == -2) {
        PyErr_SetString(_erfaError, "illegal case");
        return NULL;
    }
    else {
        return Py_BuildValue("ddd", xyz[0], xyz[1], xyz[2]);
    }
}
'''),    
    ('gd2gce', 'METH_VARARGS',
     '''"\\ngd2gce(a, f, elong, phi, height) -> xyz\\n\\n"
"Transform geodetic coordinates to geocentric for a reference\\n"
" for a reference ellipsoid of specified form\\n"
"Given:\\n"
"    a          equatorial radius in meters\\n"
"    f          flattening\\n"
"    elong      longitude (radians, east +ve)\\n"
"    phi        latitude (geodetic, radians)\\n"
"    height     height above ellipsoid (geodetic in meters)\\n"
"  Returned:\\n"
"    xyz        geocentric vector in meters"''',
     '''{
    double a, f, elong, phi, height, xyz[3];
    int status;
    if (!PyArg_ParseTuple(args, "ddddd", &a, &f, &elong, &phi, &height)) {
        return NULL;
    }
    status = eraGd2gce(a, f, elong, phi, height, xyz);
    if (status == -1) {
        PyErr_SetString(_erfaError, "illegal case");
        return NULL;
    }
    else {
        return Py_BuildValue("ddd", xyz[0], xyz[1], xyz[2]);
    }
}
'''),
    ('gmst00', 'METH_VARARGS',
     '''"\\ngmst00(uta, utb, tta, ttb) -> gmst\\n\\n"
"Greenwich mean sidereal time\\n"
"(model consistent with IAU 2000 resolutions).\\n"
"Given:\\n"
"    uta,utb    UT1 as a 2-part Julian Date\\n"
"    tta,ttb    TT as a 2-part Julian Date\\n"
"Returned:\\n"
"    g          Greenwich mean sidereal time (radians)"''',
     '''{
    double uta, utb, tta, ttb, g;
    if (!PyArg_ParseTuple(args, "dddd", &uta, &utb, &tta, &ttb)) {
        return NULL;
    }
    g = eraGmst00(uta, utb, tta, ttb);
    return Py_BuildValue("d", g);
}
'''),
    ('gmst06', 'METH_VARARGS',
     '''"\\ngmst06(uta, utb, tta, ttb) -> gmst\\n\\n"
"Greenwich mean sidereal time\\n"
"(model consistent with IAU 2006resolutions).\\n"
"Given:\\n"
"    uta,utb    UT1 as a 2-part Julian Date\\n"
"    tta,ttb    TT as a 2-part Julian Date\\n"
"Returned:\\n"
"    g          Greenwich mean sidereal time (radians)"''',
     '''{
    double uta, utb, tta, ttb, g;
    if (!PyArg_ParseTuple(args, "dddd", &uta, &utb, &tta, &ttb)) {
        return NULL;
    }
    g = eraGmst06(uta, utb, tta, ttb);
    return Py_BuildValue("d", g);
}
'''),
    ('gmst82', 'METH_VARARGS',
     '''"\\ngmst82(d1, d2) -> gmst\\n\\n"
"Universal Time to Greenwich mean sidereal time (IAU 1982 model)\\n"
"Given:\\n"
"    d1,d2      UT1 as a 2-part Julian Date\\n"
"Returned:\\n"
"    g          Greenwich mean sidereal time (radians)"''',
     '''{
    double dj1, dj2, g;
    if (!PyArg_ParseTuple(args, "dd", &dj1, &dj2)) {
        return NULL;
    }
    g = eraGmst82(dj1, dj2);
    return Py_BuildValue("d", g);
}
'''),
    ('gst00a', 'METH_VARARGS',
     '''"\\ngst00a(uta, utb, tta, ttb) -> gast\\n\\n"
"Greenwich apparent sidereal time\\n"
"(model consistent with IAU 2000resolutions).\\n"
"Given:\\n"
"    uta,utb    UT1 as a 2-part Julian Date\\n"
"    tta,ttb    TT as a 2-part Julian Date\\n"
"Returned:\\n"
"    g          Greenwich apparent sidereal time (radians)"''',
     '''{
    double uta, utb, tta, ttb, g;
    if (!PyArg_ParseTuple(args, "dddd", &uta, &utb, &tta, &ttb)) {
        return NULL;
    }
    g = eraGst00a(uta, utb, tta, ttb);
    return Py_BuildValue("d", g);
}
'''),
    ('gst00b', 'METH_VARARGS',
     '''"\\ngst00b(uta, utb) -> gast\\n\\n"
"Greenwich apparent sidereal time (model consistent with IAU 2000\\n"
"resolutions but using the truncated nutation model IAU 2000B).\\n"
"Given:\\n"
"    uta,utb    UT1 as a 2-part Julian Date\\n"
"Returned:\\n"
"    g          Greenwich apparent sidereal time (radians)"''',
     '''{
    double uta, utb, g;
    if (!PyArg_ParseTuple(args, "dd", &uta, &utb)) {
        return NULL;
    }
    g = eraGst00b(uta, utb);
    return Py_BuildValue("d", g);
}
'''),
    ('gst06', 'METH_VARARGS',
     '''"\\ngst06(uta, utb, tta, ttb, rnpb[3][3]) -> gast\\n\\n"
"Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.\\n"
"Given:\\n"
"    uta,utb    UT1 as a 2-part Julian Date\\n"
"    tta,ttb    TT as a 2-part Julian Date\\n"
"    rnpb       nutation x precession x bias matrix\\n"
"Returned:\\n"
"    g          Greenwich apparent sidereal time"''',
     '''{
    double uta, utb, tta, ttb, rnpb[3][3], g;
    if (!PyArg_ParseTuple(args, "dddd((ddd)(ddd)(ddd))",
                                 &uta, &utb, &tta, &ttb,
                                 &rnpb[0][0], &rnpb[0][1], &rnpb[0][2],
                                 &rnpb[1][0], &rnpb[1][1], &rnpb[1][2],
                                 &rnpb[2][0], &rnpb[2][1], &rnpb[2][2])) {
        return NULL;
    }        
    g = eraGst06(uta, utb, tta, ttb, rnpb);
    return Py_BuildValue("d",g);
}
'''),
    ('gst06a', 'METH_VARARGS',
     '''"\\ngst06a(uta, utb, tta, ttb) -> gast\\n\\n"
"Greenwich apparent sidereal time\\n"
"(model consistent with IAU 2000and 2006 resolutions).\\n"
"Given:\\n"
"    uta,utb    UT1 as a 2-part Julian Date\\n"
"    tta,ttb    TT as a 2-part Julian Date\\n"
"Returned:\\n"
"    g          Greenwich apparent sidereal time (radians)"''',
     '''{
    double uta, utb, tta, ttb, g;
    if (!PyArg_ParseTuple(args, "dddd", &uta, &utb, &tta, &ttb)) {
        return NULL;
    }
    g = eraGst06a(uta, utb, tta, ttb);
    return Py_BuildValue("d", g);
}
'''),
    ('gst94', 'METH_VARARGS',
     '''"\\ngst94(uta, utb) -> gast\\n\\n"
"Greenwich apparent sidereal time\\n"
"(consistent with IAU 1982/94 resolutions)\\n"
"Given:\\n"
"    uta,utb    UT1 as a 2-part Julian Date\\n"
"Returned:\\n"
"    g          Greenwich mean sidereal time (radians)"''',
     '''{
    double uta, utb, g;
    if (!PyArg_ParseTuple(args, "dd", &uta, &utb)) {
        return NULL;
    }
    g = eraGst94(uta, utb);
    return Py_BuildValue("d", g);
}
'''),
    ('h2fk5', 'METH_VARARGS',
     '''"\\nh2fk5(rh, dh, drh, ddh, pxh, rvh) -> r5, d5, dr5, dd5, px5, rv5\\n\\n"
"Transform Hipparcos star data into the FK5 (J2000.0) system.\\n"
"Given (all Hipparcos, epoch J2000.0):\\n"
"     rh    RA (radians)\\n"
"     dh    Dec (radians)\\n"
"     drh   proper motion in RA (dRA/dt, rad/Jyear)\\n"
"     ddh   proper motion in Dec (dDec/dt, rad/Jyear)\\n"
"     pxh   parallax (arcsec)\\n"
"     rvh   radial velocity (km/s, positive = receding)\\n"
"Returned (all FK5, equinox J2000.0, epoch J2000.0):\\n"
"     r5    RA (radians)\\n"
"     d5    Dec (radians)\\n"
"     dr5   proper motion in RA (dRA/dt, rad/Jyear)\\n"
"     dd5   proper motion in Dec (dDec/dt, rad/Jyear)\\n"
"     px5   parallax (arcsec)\\n"
"     rv5   radial velocity (km/s, positive = receding)"''',
     '''{
    double rh, dh, drh, ddh, pxh, rvh;
    double r5, d5, dr5, dd5, px5, rv5;
    if (!PyArg_ParseTuple(args, "dddddd", &rh, &dh, &drh, &ddh, &pxh, &rvh)) {
        return NULL;
    }
    eraH2fk5(rh, dh, drh, ddh, pxh, rvh, &r5, &d5, &dr5, &dd5, &px5, &rv5);
    return Py_BuildValue("dddddd", r5, d5, dr5, dd5, px5, rv5);
}
'''),
    ('hfk5z', 'METH_VARARGS',
     '''"\\nhf5kz(rh, dh, d1, d2) -> r5, d5, dr5, dd5\\n\\n"
"Transform a Hipparcos star position into FK5 J2000.0, assuming\\n"
"zero Hipparcos proper motion.\\n"
"Given:\\n"
"     rh    Hipparcos RA (radians)\\n"
"     dh    Hipparcos Dec (radians)\\n"
"     d1,d2 TDB date\\n"
"Returned (all FK5, equinox J2000.0, date date1+date2):\\n"
"     r5    RA (radians)\\n"
"     d5    Dec (radians)\\n"
"     dr5   FK5 RA proper motion (rad/year)\\n"
"     dd5   Dec proper motion (rad/year)"''',
     '''{
    double rh, dh, d1, d2;
    double r5, d5, dr5, dd5;
    if (!PyArg_ParseTuple(args, "dddd", &rh, &dh, &d1, &d2)) {
        return NULL;
    }
    eraHfk5z(rh, dh, d1, d2, &r5, &d5, &dr5, &dd5);
    return Py_BuildValue("dddd", r5, d5, dr5, dd5);    
}
'''),
    ('jd2cal', 'METH_VARARGS',
     '''"\\njd2cal(dj1, dj2) -> year, month, day, fraction of day\\n\\n"
"Julian Date to Gregorian year, month, day, and fraction of a day.\\n"
"Given:\\n"
"     dj1,dj2   2-part Julian Date\\n"
"Returned:\\n"
"     iy    year\\n"
"     im    month\\n"
"     id    day\\n"
"     fd    fraction of day"''',
    '''{
    double dj1, dj2, fd;
    int iy, im, id, status;
    if (!PyArg_ParseTuple(args, "dd", &dj1, &dj2)) {
        return NULL;
    }
    status = eraJd2cal(dj1, dj2, &iy, &im, &id, &fd);
    if (status) {
        PyErr_SetString(_erfaError, "unacceptable date");
        return NULL;
    }
    return Py_BuildValue("iiid", iy, im, id, fd);    
}
'''),
    ('jdcalf', 'METH_VARARGS',
     '''"\\njdcalf(n, d1, d2) -> y, m, d, fd\\n\\n"
"Julian Date to Gregorian Calendar, expressed in a form convenient\\n"
"for formatting messages:  rounded to a specified precision.\\n"
"Given:\\n"
"     n         number of decimal places of days in fraction\\n"
"     d1,d2     2-part Julian Date\\n"
"Returned:\\n"
"     y         year\\n"
"     m         month\\n"
"     d         day\\n"
"     fd        fraction of day"''',
    '''{
    double dj1, dj2;
    int ndp, status, iymdf[4];
    if (!PyArg_ParseTuple(args, "idd", &ndp, &dj1, &dj2)) {
        return NULL;
    }
    status = eraJdcalf(ndp, dj1, dj2, iymdf);
    if (status == -1) {
        PyErr_SetString(_erfaError, "date out of range");
        return NULL;
    }
    if (status == 1) {
        PyErr_SetString(_erfaError, "n not in 0-9 (interpreted as 0)");
        return NULL;
    }
    else {
        return Py_BuildValue("iiii", iymdf[0],iymdf[1],iymdf[2],iymdf[3]);
    }
}
'''),
    ('num00a', 'METH_VARARGS',
     '''"\\nnum00a(d1, d2) -> rmatn\\n\\n"
"Form the matrix of nutation for a given date, IAU 2000A model.\\n"
"Given:\\n"
"     d1,d2     TT as a 2-part Julian Date\\n"
"Returned:\\n"
"     rmatn     nutation matrix"''',
     '''{
    double d1, d2, rmatn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNum00a(d1, d2, rmatn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            rmatn[0][0],rmatn[0][1],rmatn[0][2],
                            rmatn[1][0],rmatn[1][1],rmatn[1][2],
                            rmatn[2][0],rmatn[2][1],rmatn[2][2]);    
}
'''),
    ('num00b', 'METH_VARARGS',
     '''"\\nnum00b(d1, d2) -> rmatn\\n\\n"
"Form the matrix of nutation for a given date, IAU 2000B model.\\n"
"Given:\\n"
"     d1,d2     TT as a 2-part Julian Date\\n"
"Returned:\\n"
"     rmatn     nutation matrix"''',
     '''{
    double d1, d2, rmatn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNum00b(d1, d2, rmatn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatn[0][0],rmatn[0][1],rmatn[0][2],
    rmatn[1][0],rmatn[1][1],rmatn[1][2],
    rmatn[2][0],rmatn[2][1],rmatn[2][2]);    
}
'''),
    ('num06a', 'METH_VARARGS',
     '''"\\nnum06a(d1, d2) -> rmatn\\n\\n"
"Form the matrix of nutation for a given date, IAU 2006/2000A model.\\n"
"Given:\\n"
"     d1,d2     TT as a 2-part Julian Date\\n"
"Returned:\\n"
"     rmatn     nutation matrix"''',
     '''{
    double d1, d2, rmatn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNum06a(d1, d2, rmatn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatn[0][0],rmatn[0][1],rmatn[0][2],
    rmatn[1][0],rmatn[1][1],rmatn[1][2],
    rmatn[2][0],rmatn[2][1],rmatn[2][2]);    
}
'''),
    ('numat', 'METH_VARARGS',
     '''"\\nnumat(epsa, dpsi, deps) -> rmatn\\n\\n"
"Form the matrix of nutation.\\n"
"Given:\\n"
"     epsa          mean obliquity of date\\n"
"     dpsi,deps     nutation\\n"
"Returned:\\n"
"     rmatn         nutation matrix"''',
     '''{
    double epsa, dpsi, deps, rmatn[3][3];
    if (!PyArg_ParseTuple(args, "ddd", &epsa, &dpsi, &deps)) {
        return NULL;
    }
    eraNumat(epsa, dpsi, deps, rmatn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatn[0][0],rmatn[0][1],rmatn[0][2],
    rmatn[1][0],rmatn[1][1],rmatn[1][2],
    rmatn[2][0],rmatn[2][1],rmatn[2][2]);    
}
'''),
    ('nut00a', 'METH_VARARGS',
     '''"\\nnut00a(d1, d2) -> dpsi, deps\\n\\n"
"Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation\\n"
"with free core nutation omitted).\\n"
"Given:\\n"
"    d1,d2      TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   dpsi,deps   nutation, luni-solar + planetary\\n"''',
     '''{
    double d1, d2, dpsi, deps;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNut00a(d1, d2, &dpsi, &deps);
    return Py_BuildValue("dd", dpsi, deps);
}
'''),
    ('nut00b', 'METH_VARARGS',
     '''"\\nnut00b(d1, d2) -> dpsi, deps\\n\\n"
"Nutation, IAU 2000B.\\n"
"Given:\\n"
"    d1,d2       TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   dpsi,deps   nutation, luni-solar + planetary\\n"''',
     '''{
    double d1, d2, dpsi, deps;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNut00b(d1, d2, &dpsi, &deps);
    return Py_BuildValue("dd", dpsi, deps);
}
'''),
    ('nut06a', 'METH_VARARGS',
     '''"\\nnut06a(d1, d2) -> dpsi, deps\\n\\n"
"IAU 2000A nutation with adjustments to match the IAU 2006 precession.\\n"
"Given:\\n"
"    d1,d2      TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   dpsi,deps   nutation, luni-solar + planetary\\n"''',
     '''{
    double d1, d2, dpsi, deps;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNut06a(d1, d2, &dpsi, &deps);
    return Py_BuildValue("dd", dpsi, deps);
}
'''),
    ('nut80', 'METH_VARARGS',
     '''"\\nnut80(d1, d2) -> dpsi, deps\\n\\n"
"Nutation, IAU 1980 model.\\n"
"Given:\\n"
"    d1,d2      TT as a 2-part Julian Date\\n"
"Returned:\\n"
"    dpsi        nutation in longitude (radians)\\n"
"    deps        nutation in obliquity (radians)\\n"''',
     '''{
    double d1, d2, dpsi, deps;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNut80(d1, d2, &dpsi, &deps);
    return Py_BuildValue("dd", dpsi, deps);
}
'''),
    ('nutm80', 'METH_VARARGS',
     '''"\\nnutm80(d1, d2) -> rmatn\\n\\n"
"Form the matrix of nutation for a given date, IAU 1980 model.\\n"
"Given:\\n"
"    d1,d2      TDB as a 2-part Julian Date\\n"
"Returned:\\n"
"    rmatn[     nutation matrix"''',
     '''{
    double d1, d2, rmatn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNutm80(d1, d2, rmatn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatn[0][0],rmatn[0][1],rmatn[0][2],
    rmatn[1][0],rmatn[1][1],rmatn[1][2],
    rmatn[2][0],rmatn[2][1],rmatn[2][2]);    
}
'''),
    ('obl06', 'METH_VARARGS',
     '''"\\nobl06(d1, d2) -> obl\\n\\n"
"Mean obliquity of the ecliptic, IAU 2006 precession model.\\n"
"Given:\\n"
"    d1,d2      TT as a 2-part Julian Date\\n"
"Returned:\\n"
"    obl        obliquity of the ecliptic (radians)"''',
     '''{
    double d1, d2, obl;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    obl = eraObl06(d1, d2);
        return Py_BuildValue("d", obl);
}
'''),
    ('obl80', 'METH_VARARGS',
     '''"\\nobl80(d1, d2) -> obl\\n\\n"
"Mean obliquity of the ecliptic, IAU 1980 model.\\n"
"Given:\\n"
"    d1,d2      TT as a 2-part Julian Date\\n"
"Returned:\\n"
"    obl        obliquity of the ecliptic (radians)"''',
     '''{
    double d1, d2, obl;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    obl = eraObl80(d1, d2);
    return Py_BuildValue("d", obl);
}
'''),
    ('p06e', 'METH_VARARGS',
     '''"\\np06e(d1, d2) -> eps0,psia,oma,bpa,bqa,pia,bpia,epsa,chia,za,zetaa,thetaa,pa,gam,phi,psi\\n\\n"
"Precession angles, IAU 2006, equinox based.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   eps0   epsilon_0   obliquity at J2000.0\\n"
"   psia   psi_A       luni-solar precession\\n"
"   oma    omega_A     inclination of equator wrt J2000.0 ecliptic\\n"
"   bpa    P_A         ecliptic pole x, J2000.0 ecliptic triad\\n"
"   bqa    Q_A         ecliptic pole -y, J2000.0 ecliptic triad\\n"
"   pia    pi_A        angle between moving and J2000.0 ecliptics\\n"
"   bpia   Pi_A        longitude of ascending node of the ecliptic\\n"
"   epsa   epsilon_A   obliquity of the ecliptic\\n"
"   chia   chi_A       planetary precession\\n"
"   za     z_A         equatorial precession: -3rd 323 Euler angle\\n"
"   zetaa  zeta_A      equatorial precession: -1st 323 Euler angle\\n"
"   thetaa theta_A     equatorial precession: 2nd 323 Euler angle\\n"
"   pa     p_A         general precession\\n"
"   gam    gamma_J2000 J2000.0 RA difference of ecliptic poles\\n"
"   phi    phi_J2000   J2000.0 codeclination of ecliptic pole\\n"
"   psi    psi_J2000   longitude difference of equator poles, J2000.0"''',
     '''{
    double d1,d2, eps0,psia,oma,bpa,bqa,pia,bpia,epsa,chia,za,zetaa,thetaa,pa,gam,phi,psi;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraP06e(d1,d2,
    &eps0,&psia,&oma,&bpa,&bqa,&pia,&bpia,&epsa,
    &chia,&za,&zetaa,&thetaa,&pa,&gam,&phi,&psi);
    return Py_BuildValue("dddddddddddddddd",
    eps0,psia,oma,bpa,bqa,pia,bpia,epsa,chia,za,zetaa,thetaa,pa,gam,phi,psi);
}
'''),
    ('pb06', 'METH_VARARGS',
     '''"\\npb06(d1, d2) -> bzeta, bz, btheta\\n\\n"
"Forms three Euler angles which implement general\\n"
"precession from epoch J2000.0, using the IAU 2006 model.\\n"
"Framebias (the offset between ICRS and mean J2000.0) is included.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   bzeta   1st rotation: radians cw around z\\n"
"   bz      3rd rotation: radians cw around z\\n"
"   btheta  2nd rotation: radians ccw around y"''',
     '''{
    double d1, d2, bzeta, bz, btheta;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPb06(d1, d2, &bzeta, &bz, &btheta);
    return Py_BuildValue("ddd", bzeta, bz, btheta);
}
'''),
    ('pfw06', 'METH_VARARGS',
     '''"\\npfw06(d1, d2) -> gamb, phib, psib, epsa\\n\\n"
"Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   gamb    F-W angle gamma_bar (radians)\\n"
"   phib    F-W angle phi_bar (radians)\\n"
"   psib    F-W angle psi_bar (radians)\\n"
"   epsa    F-W angle epsilon_A (radians)"''',
     '''{
    double d1, d2, gamb, phib, psib, epsa;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPfw06(d1,d2, &gamb, &phib, &psib, &epsa);
    return Py_BuildValue("dddd", gamb, phib, psib, epsa);
}
'''),
    ('plan94', 'METH_VARARGS',
     '''"\\nplan94(d1, d2, np) -> pv\\n\\n"
"Approximate heliocentric position and velocity of a nominated major\\n"
"planet:  Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or\\n"
"Neptune (but not the Earth itself).\\n"
"Given:\\n"
"    d1         TDB date part A\\n"
"    d2         TDB date part B\\n"
"    np         planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars,\\n"
"                       5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)\\n"
"Returned:\\n"
"    pv         planet p,v (heliocentric, J2000.0, AU,AU/d)"''',
     '''{
    double d1, d2, pv[2][3];
    int np, status;
    if (!PyArg_ParseTuple(args, "ddi", &d1, &d2, &np)) {
        return NULL;
    }
    status = eraPlan94(d1, d2, np, pv);
    if (status == -1){
        PyErr_SetString(_erfaError, "illegal np,  not in range(1,8) for planet");
        return NULL;
    }
/*    else if (status == 1){
        PyErr_SetString(_erfaError, "year outside range(1000:3000)");
        return NULL;
    }*/
    else if (status == 2){
        PyErr_SetString(_erfaError, "computation failed to converge");
        return NULL;
    }
    else {
        return Py_BuildValue("(ddd)(ddd)",
        pv[0][0], pv[0][1], pv[0][2],
        pv[1][0], pv[1][1], pv[1][2]);
    }
}
'''
     ),
    ('pmat00', 'METH_VARARGS',
     '''"\\npmat00(d1, d2) -> rbp\\n\\n"
"Precession matrix (including frame bias) from GCRS to a specified\\n"
"date, IAU 2000 model.\\n"
"Given:\\n"
"   d1,d2       TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   rbp         bias-precession matrix"''',
     '''{
    double d1, d2, rbp[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPmat00(d1, d2, rbp);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rbp[0][0],rbp[0][1],rbp[0][2],
    rbp[1][0],rbp[1][1],rbp[1][2],
    rbp[2][0],rbp[2][1],rbp[2][2]);    
}
'''),
    ('pmat06', 'METH_VARARGS',
     '''"\\npmat06(d1, d2) -> rbp\\n\\n"
"Precession matrix (including frame bias) from GCRS to a specified\\n"
"date, IAU 2006 model.\\n"
"Given:\\n"
"   d1,d2       TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   rbp         bias-precession matrix"''',
     '''{
    double d1, d2, rbp[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPmat06(d1, d2, rbp);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rbp[0][0],rbp[0][1],rbp[0][2],
    rbp[1][0],rbp[1][1],rbp[1][2],
    rbp[2][0],rbp[2][1],rbp[2][2]);    
}
'''),
    ('pmat76', 'METH_VARARGS',
     '''"\\npmat76(d1, d2) -> rmatp\\n\\n"
"Precession matrix from J2000.0 to a specified date, IAU 1976 model.\\n"
"Given:\\n"
"   d1,d2       TT ending date as a 2-part Julian Date\\n"
"Returned:\\n"
"   rmatp       precession matrix, J2000.0 -> d1+d2"''',
     '''{
    double d1, d2, rmatp[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPmat76(d1, d2, rmatp);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatp[0][0],rmatp[0][1],rmatp[0][2],
    rmatp[1][0],rmatp[1][1],rmatp[1][2],
    rmatp[2][0],rmatp[2][1],rmatp[2][2]);    
}
'''),
    ('pn00', 'METH_VARARGS',
     '''"\\npn00(d1,d2,dpsi,deps) -> epsa,rb,rp,rbp,rn,rbpn\\n\\n"
"Precession-nutation, IAU 2000 model:  a multi-purpose function,\\n"
"supporting classical (equinox-based) use directly and CIO-based\\n"
"use indirectly.\\n"
"Given:\\n"
"   d1,d2       TT as a 2-part Julian Date\\n"
"   dpsi,deps   nutation\\n"
"Returned:\\n"
"   epsa        mean obliquity\\n"
"   rb          frame bias matrix\\n"
"   rp          precession matrix\\n"
"   rbp         bias-precession matrix\\n"
"   rn          nutation matrix\\n"
"   rbpn        GCRS-to-true matrix"''',
     '''{
    double d1,d2,dpsi,deps,epsa;
    double rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dddd", &d1, &d2, &dpsi, &deps)) {
        return NULL;
    }
    eraPn00(d1,d2,dpsi,deps,&epsa,rb,rp,rbp,rn,rbpn);
    return Py_BuildValue("d((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
    epsa,
    rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
    rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
    rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2],
    rn[0][0],rn[0][1],rn[0][2],rn[1][0],rn[1][1],rn[1][2],rn[2][0],rn[2][1],rn[2][2],
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}
'''),
    ('pn00a', 'METH_VARARGS',
     '''"\\npn00a(d1,d2) -> dpsi,deps,epsa,rb,rp,rbp,rn,rbpn\\n\\n"
"Precession-nutation, IAU 2000A model:  a multi-purpose function,\\n"
"supporting classical (equinox-based) use directly and CIO-based\\n"
"use indirectly.\\n"
"Given:\\n"
"   d1,d2       TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   dpsi,deps   nutation\\n"
"   epsa        mean obliquity\\n"
"   rb          frame bias matrix\\n"
"   rp          precession matrix\\n"
"   rbp         bias-precession matrix\\n"
"   rn          nutation matrix\\n"
"   rbpn        GCRS-to-true matrix"''',
     '''{
    double d1,d2,dpsi,deps,epsa;
    double rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPn00a(d1,d2,&dpsi,&deps,&epsa,rb,rp,rbp,rn,rbpn);
    return Py_BuildValue("ddd((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
    dpsi,deps,epsa,
    rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
    rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
    rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2],
    rn[0][0],rn[0][1],rn[0][2],rn[1][0],rn[1][1],rn[1][2],rn[2][0],rn[2][1],rn[2][2],
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}
'''),
    ('pn00b', 'METH_VARARGS',
     '''"\\npn00b(d1,d2) -> dpsi,deps,epsa,rb,rp,rbp,rn,rbpn\\n\\n"
"Precession-nutation, IAU 2000B model:  a multi-purpose function,\\n"
"supporting classical (equinox-based) use directly and CIO-based\\n"
"use indirectly.\\n"
"Given:\\n"
"   d1,d2       TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   dpsi,deps   nutation\\n"
"   epsa        mean obliquity\\n"
"   rb          frame bias matrix\\n"
"   rp          precession matrix\\n"
"   rbp         bias-precession matrix\\n"
"   rn[         nutation matrix\\n"
"   rbpn        GCRS-to-true matrix"''',
     '''{
    double d1,d2,dpsi,deps,epsa;
    double rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPn00b(d1,d2,&dpsi,&deps,&epsa,rb,rp,rbp,rn,rbpn);
    return Py_BuildValue("ddd((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
    dpsi,deps,epsa,
    rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
    rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
    rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2],
    rn[0][0],rn[0][1],rn[0][2],rn[1][0],rn[1][1],rn[1][2],rn[2][0],rn[2][1],rn[2][2],
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}
'''),
    ('pn06', 'METH_VARARGS',
     '''"\\npn06(d1,d2,dpsi,deps) -> epsa,rb,rp,rbp,rn,rbpn\\n\\n"
"Precession-nutation, IAU 2006 model:  a multi-purpose function,\\n"
"supporting classical (equinox-based) use directly and CIO-based\\n"
"use indirectly.\\n"
"Given:\\n"
"   d1,d2       TT as a 2-part Julian Date\\n"
"   dpsi,deps   nutation\\n"
"Returned:\\n"
"   epsa        mean obliquity\\n"
"   rb          frame bias matrix\\n"
"   rp          precession matrix\\n"
"   rbp         bias-precession matrix\\n"
"   rn          nutation matrix\\n"
"   rbpn        GCRS-to-true matrix"''',
     '''{
    double d1,d2,dpsi,deps,epsa;
    double rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dddd", &d1, &d2, &dpsi, &deps)) {
        return NULL;
    }
    eraPn06(d1,d2,dpsi,deps,&epsa,rb,rp,rbp,rn,rbpn);
    return Py_BuildValue("d((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
    epsa,
    rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
    rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
    rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2],
    rn[0][0],rn[0][1],rn[0][2],rn[1][0],rn[1][1],rn[1][2],rn[2][0],rn[2][1],rn[2][2],
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}
'''),
    ('pn06a', 'METH_VARARGS',
     '''"\\npn06a(d1,d2) -> dpsi,deps,epsa,rb,rp,rbp,rn,rbpn\\n\\n"
"Precession-nutation, IAU 2006/2000A model:  a multi-purpose function,\\n"
"supporting classical (equinox-based) use directly and CIO-based\\n"
"use indirectly.\\n"
"Given:\\n"
"   d1,d2       TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   dpsi,deps   nutation\\n"
"   epsa        mean obliquity\\n"
"   rb          frame bias matrix\\n"
"   rp          precession matrix\\n"
"   rbp         bias-precession matrix\\n"
"   rn          nutation matrix\\n"
"   rbpn        GCRS-to-true matrix"''',
     '''{
    double d1,d2,dpsi,deps,epsa;
    double rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPn06a(d1,d2,&dpsi,&deps,&epsa,rb,rp,rbp,rn,rbpn);
    return Py_BuildValue("ddd((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
    dpsi,deps,epsa,
    rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
    rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
    rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2],
    rn[0][0],rn[0][1],rn[0][2],rn[1][0],rn[1][1],rn[1][2],rn[2][0],rn[2][1],rn[2][2],
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}
'''),
    ('pnm00a', 'METH_VARARGS',
     '''"\\npnm00a(d1, d2) -> rbpn\\n\\n"
"Form the matrix of precession-nutation for a given date (including\\n"
"frame bias), equinox-based, IAU 2000A model.\\n"
"Given:\\n"
"   d1,d2       TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   rbpn        classical NPB matrix"''',
     '''{
    double d1,d2,rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPnm00a(d1, d2, rbpn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}
'''),
    ('pnm00b', 'METH_VARARGS',
     '''"\\npnm00b(d1, d2) -> rbpn\\n\\n"
"Form the matrix of precession-nutation for a given date (including\\n"
"frame bias), equinox-based, IAU 2000B model.\\n"
"Given:\\n"
"   d1,d2       TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   rbpn        classical NPB matrix"''',
     '''{
    double d1,d2,rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPnm00b(d1, d2, rbpn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}
'''),
    ('pnm06a', 'METH_VARARGS',
     '''"\\npnm06a(d1, d2) -> rbpn\\n\\n"
"Form the matrix of precession-nutation for a given date (including\\n"
"frame bias), IAU 2006 precession and IAU 2000A nutation models.\\n"
"Given:\\n"
"   d1,d2       TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   rbpn        classical NPB matrix"''',
     '''{
    double d1,d2,rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPnm06a(d1, d2, rbpn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}
'''),
    ('pnm80', 'METH_VARARGS',
     '''"\\npnm80(d1, d2) -> rmatp\\n\\n"
"Form the matrix of precession/nutation for a given date, IAU 1976\\n"
"precession model, IAU 1980 nutation model.\\n"
"Given:\\n"
"   d1,d2           TDB date as a 2-part Julian Date\\n"
"Returned:\\n"
"   rmatp           combined precession/nutation matrix"''',
     '''{
    double d1, d2, rmatp[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPnm80(d1, d2, rmatp);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatp[0][0],rmatp[0][1],rmatp[0][2],
    rmatp[1][0],rmatp[1][1],rmatp[1][2],
    rmatp[2][0],rmatp[2][1],rmatp[2][2]);    
}
'''),
    ('pom00', 'METH_VARARGS',
     '''"\\npom00(xp, yp, sp) -> rpom\\n\\n"
"Form the matrix of polar motion for a given date, IAU 2000.\\n"
"Given:\\n"
"   xp,yp       coordinates of the pole (radians)\\n"
"   sp          the TIO locator s' (radians)\\n"
"Returned:\\n"
"   rpom        polar-motion matrix"''',
     '''{    
    double xp, yp, sp, rpom[3][3];
    if (!PyArg_ParseTuple(args, "ddd", &xp, &yp, &sp)) {
        return NULL;
    }
    eraPom00(xp, yp, sp, rpom);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rpom[0][0],rpom[0][1],rpom[0][2],
    rpom[1][0],rpom[1][1],rpom[1][2],
    rpom[2][0],rpom[2][1],rpom[2][2]);    
}
'''),
    ('pr00', 'METH_VARARGS',
     '''"\\npr00(d1,d2) -> dpsipr,depspr\\n\\n"
"Precession-rate part of the IAU 2000 precession-nutation models\\n"
"(part of MHB2000).\\n"
"Given:\\n"
"   d1,d2           TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   dpsipr,depspr   precession corrections"''',
     '''{
    double d1, d2, dpsipr, depspr;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPr00(d1, d2, &dpsipr, &depspr);
    return Py_BuildValue("dd", dpsipr, depspr);
}
'''),
    ('prec76', 'METH_VARARGS',
     '''"\\nprec76(ep01, ep02, ep11, ep12) -> zeta, z, theta\\n\\n"
"IAU 1976 precession model.\\n"
"Given:\\n"
"   ep01,ep02   TDB starting epoch as 2-part Julian Date\\n"
"   ep11,ep12   TDB ending epoch as 2-part Julian Date\\n"
"Returned:\\n"
"   zeta        1st rotation: radians cw around z\\n"
"   z           3rd rotation: radians cw around z\\n"
"   theta       2nd rotation: radians ccw around y"''',
     '''{
    double ep01, ep02, ep11, ep12, zeta, z, theta;
    if (!PyArg_ParseTuple(args, "dddd", &ep01, &ep02, &ep11, &ep12)) {
        return NULL;
    }
    eraPrec76(ep01, ep02, ep11, ep12, &zeta, &z, &theta);
    return Py_BuildValue("ddd", zeta, z, theta);
}
'''),
    ('pvstar', 'METH_VARARGS',
     '''"\\npvstar(pv[2][3]) -> ra, dec, pmr, pmd, px, rv\\n\\n"
"Convert star position+velocity vector to catalog coordinates.\\n"
"Given:\\n"
"   pv[2][3]    pv-vector (AU, AU/day)\\n"
"Returned:\\n"
"   ra          right ascension (radians)\\n"
"   dec         declination (radians)\\n"
"   pmr         RA proper motion (radians/year)\\n"
"   pmd         Dec proper motion (radians/year)\\n"
"   px          parallax (arcsec)\\n"
"   rv          radial velocity (km/s, positive = receding)"''',
     '''{
    double pv[2][3], ra, dec, pmr, pmd, px, rv;
    int status;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))",
                                   &pv[0][0], &pv[0][1], &pv[0][2],
                                   &pv[1][0], &pv[1][1], &pv[1][2])) {
        return NULL;
    }
    status = eraPvstar(pv, &ra, &dec, &pmr, &pmd, &px, &rv);
    if (status == -1) {
        PyErr_SetString(_erfaError, "superluminal speed");
        return NULL;
    }
    else if (status == -2) {
        PyErr_SetString(_erfaError, "null position vector");
        return NULL;
    }
    else {
    return Py_BuildValue("dddddd", ra, dec, pmr, pmd, px, rv);
    }
}'''),
    ('s00', 'METH_VARARGS',
     '''"\\ns00(d1, d2, x, y) -> s\\n\\n"
"The CIO locator s, positioning the Celestial Intermediate Origin on\\n"
"the equator of the Celestial Intermediate Pole, given the CIP's X,Y\\n"
"coordinates.  Compatible with IAU 2000A precession-nutation.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"   x,y     CIP coordinates\\n"
"Returned:\\n"
"   s       the CIO locator s in radians"''',
     '''{
    double d1, d2, x, y, s;
    if (!PyArg_ParseTuple(args, "dddd", &d1, &d2, &x, &y)) {
        return NULL;
    }
    s = eraS00(d1, d2, x, y);
    return Py_BuildValue("d", s);
}
'''),
    ('s00a', 'METH_VARARGS',
     '''"\\ns00a(d1, d2) -> s\\n\\n"
"The CIO locator s, positioning the Celestial Intermediate Origin on\\n"
"the equator of the Celestial Intermediate Pole, using the IAU 2000A\\n"
"precession-nutation model.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   s       the CIO locator s in radians"''',
     '''{
    double d1, d2,s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    s = eraS00a(d1, d2);
    return Py_BuildValue("d", s);
}
'''),
    ('s00b', 'METH_VARARGS',
     '''"\\ns00b(d1, d2) -> s\\n\\n"
"The CIO locator s, positioning the Celestial Intermediate Origin on\\n"
"the equator of the Celestial Intermediate Pole, using the IAU 2000B\\n"
"precession-nutation model.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   s       the CIO locator s in radians"''',
     '''{
    double d1, d2,s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    s = eraS00b(d1, d2);
    return Py_BuildValue("d", s);
}
'''),
    ('s06', 'METH_VARARGS',
     '''"\\ns06(d1, d2, x, y) -> s\\n\\n"
"The CIO locator s, positioning the Celestial Intermediate Origin on\\n"
"the equator of the Celestial Intermediate Pole, given the CIP's X,Y\\n"
"coordinates.  Compatible with IAU 2006/2000A precession-nutation.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"   x,y     CIP coordinates\\n"
"Returned:\\n"
"   s       the CIO locator s in radians"''',
     '''{
    double d1, d2, x, y, s;
    if (!PyArg_ParseTuple(args, "dddd", &d1, &d2, &x, &y)) {
        return NULL;
    }
    s = eraS06(d1, d2, x, y);
    return Py_BuildValue("d", s);
}
'''),
    ('s06a', 'METH_VARARGS',
     '''"\\ns06a(d1, d2) -> s\\n\\n"
"The CIO locator s, positioning the Celestial Intermediate Origin on\\n"
"the equator of the Celestial Intermediate Pole, using the IAU 2006\\n"
"precession and IAU 2000A nutation model.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   s       the CIO locator s in radians"''',
     '''{
    double d1, d2,s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    s = eraS06a(d1, d2);
    return Py_BuildValue("d", s);
}
'''),
    ('sp00', 'METH_VARARGS',
     '''"\\nsp00(d1, d2) -> s\\n\\n"
"The TIO locator s', positioning the Terrestrial Intermediate Origin\\n"
"on the equator of the Celestial Intermediate Pole.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   s       the TIO locator s' in radians"''',
     '''{
    double d1, d2, s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    s = eraSp00(d1, d2);  
    return Py_BuildValue("d", s);
}
'''),
    ('starpm', 'METH_VARARGS',
     '''"\\nstarpm(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b) -> ra2, dec2, pmr2, pmd2, px2, rv2\\n\\n"
"Star proper motion:  update star catalog data for space motion.\\n"
"Given:\\n"
"   ra1     right ascension (radians), before\\n"
"   dec1    declination (radians), before\\n"
"   pmr1    RA proper motion (radians/year), before\\n"
"   pmd1    Dec proper motion (radians/year), before\\n"
"   px1     parallax (arcseconds), before\\n"
"   rv1     radial velocity (km/s, +ve = receding), before\\n"
"   ep1a    ''before'' epoch, part A \\n"
"   ep1b    ''before'' epoch, part B\\n"
"   ep2a    ''after'' epoch, part A\\n"
"   ep2b    ''after'' epoch, part B\\n"
"Returned:\\n"
"   ra2     right ascension (radians), after\\n"
"   dec2    declination (radians), after\\n"
"   pmr2    RA proper motion (radians/year), after\\n"
"   pmd2    Dec proper motion (radians/year), after\\n"
"   px2     parallax (arcseconds), after\\n"
"   rv2     radial velocity (km/s, +ve = receding), after"''',
     '''{
   double ra1, dec1, pmr1, pmd1, px1, rv1;
   double ep1a, ep1b, ep2a, ep2b;
   double ra2, dec2, pmr2, pmd2, px2, rv2;
   int status;
   if (!PyArg_ParseTuple(args, "dddddddddd",
        &ra1, &dec1, &pmr1, &pmd1, &px1, &rv1,
        &ep1a, &ep1b, &ep2a, &ep2b)) {
        return NULL;
    }
    status = eraStarpm(ra1, dec1, pmr1, pmd1, px1, rv1,
                 ep1a, ep1b, ep2a, ep2b,
                 &ra2, &dec2, &pmr2, &pmd2, &px2, &rv2);
    if (status) {
        if (status == -1) {
            PyErr_SetString(_erfaError, "system error, normally should NOT occur");
            return NULL;
        }
        else if (status == 1) {
            PyErr_SetString(_erfaError, "distance overriden, extremely small (or zero or negative) parallax");
            return NULL;
        }
        else if (status == 2) {
            PyErr_SetString(_erfaError, "excessive velocity");
            return NULL;
        }
        else if (status == 4) {
            PyErr_SetString(_erfaError, "solution didn't converge");
            return NULL;
        }
        else {
            PyErr_SetString(_erfaError, "binary logical OR of other error above");
            return NULL;
        }
    }
    else {
        return Py_BuildValue("dddddd", ra2, dec2, pmr2, pmd2, px2, rv2);
    }
}
'''),
    ('starpv', 'METH_VARARGS',
     '''"\\nstarpv(ra1, dec1, pmr1, pmd1, px1, rv1) -> pv\\n\\n"
"Convert star catalog coordinates to position+velocity vector.\\n"
"Given:\\n"
"   ra      right ascension (radians)\\n"
"   dec     declination (radians)\\n"
"   pmr     RA proper motion (radians/year)\\n"
"   pmd     Dec proper motion (radians/year)\\n"
"   px      parallax (arcseconds)\\n"
"   rv  radial velocity (km/s, positive = receding)\\n"
"Returned:\\n"
"   pv      pv-vector (AU, AU/day)"''',
     '''{
   double ra, dec, pmr, pmd, px, rv;
   double pv[2][3];
   int status;
   if (!PyArg_ParseTuple(args, "dddddd",
                                &ra, &dec, &pmr, &pmd, &px, &rv)) {
        return NULL;
    }
    status = eraStarpv(ra, dec, pmr, pmd, px, rv, pv);
    if (status) {
        if (status == 1) {
            PyErr_SetString(_erfaError, "distance overriden, extremely small (or zero or negative) parallax");
            return NULL;
        }
        else if (status == 2) {
            PyErr_SetString(_erfaError, "excessive velocity");
            return NULL;
        }
        else if (status == 4) {
            PyErr_SetString(_erfaError, "solution didn't converge");
            return NULL;
        }
        else {
            PyErr_SetString(_erfaError, "binary logical OR of other error above");
            return NULL;
        }
    }
    else {
        return Py_BuildValue("(ddd)(ddd)",
            pv[0][0],pv[0][1],pv[0][2],
            pv[1][0],pv[1][1],pv[1][2]);
    }
}
'''),
    ('taitt', 'METH_VARARGS',
     '''"\\ntaitt(tai1, tai2) -> tt1, tt2\\n\\n"
"Time scale transformation:  International Atomic Time, TAI, to\\n"
"Terrestrial Time, TT.\\n"
"Given:\\n"
"   tai1,tai2   TAI as a 2-part Julian Date\\n"
"Returned:\\n"
"   tt1,tt2     TT as a 2-part Julian Date"''',
     '''{
    double tai1, tai2, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tai1, &tai2)) {
        return NULL;
    }
    status = eraTaitt(tai1, tai2, &tt1, &tt2);
    if (status) {
        PyErr_SetString(_erfaError, "internal error...");
        return NULL;
    }
    return Py_BuildValue("dd", tt1, tt2);
}
'''),
    ('taiut1', 'METH_VARARGS',
     '''"\\ntaiut1(tai1, tai2, dta) -> ut11, ut12\\n\\n"
"Time scale transformation:  International Atomic Time, TAI, to\\n"
"Universal Time, UT1..\\n"
"The argument dta, i.e. UT1-TAI, is an observed quantity, and is\\n"
"available from IERS tabulations.\\n"
"Given:\\n"
"   tai1,tai2   TAI as a 2-part Julian Date\\n"
"   dta         UT1-TAI in seconds\\n"
"Returned:\\n"
"   ut11,ut12     TT as a 2-part Julian Date"''',
     '''{
    double tai1, tai2, dta, ut11, ut12;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &tai1, &tai2, &dta)) {
        return NULL;
    }
    status = eraTaiut1(tai1, tai2, dta, &ut11, &ut12);
    if (status) {
        PyErr_SetString(_erfaError, "internal error...");
        return NULL;
    }
    return Py_BuildValue("dd", ut11, ut12);
}
'''),
    ('taiutc', 'METH_VARARGS',
     '''"\\ntaiutc(tai1, tai2) -> utc1, utc2\\n\\n"
"Time scale transformation:  International Atomic Time, TAI, to\\n"
"Coordinated Universal Time, UTC.\\n"
"Given:\\n"
"   tai1,tai2   TAI as a 2-part Julian Date\\n"
"Returned:\\n"
"   utc1,utc2   TT as a 2-part Julian Date"''',
     '''{
    double tai1, tai2, utc1, utc2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tai1, &tai2)) {
        return NULL;
    }
    status = eraTaiutc(tai1, tai2, &utc1, &utc2);
    if (status) {
        if (status == 1) {
            PyErr_SetString(_erfaError, "dubious year");
            return NULL;
        }
        else if (status == -1) {
            PyErr_SetString(_erfaError, "unacceptable date");
            return NULL;
        }
    }
    return Py_BuildValue("dd", utc1, utc2);
}
'''),
    ('tcbtdb', 'METH_VARARGS',
     '''"\\ntcbtdb(tcb1, tcb2) -> tdb1, tdb2\\n\\n"
"Time scale transformation:  Barycentric Coordinate Time, TCB, to\\n"
"Barycentric Dynamical Time, TDB.\\n"
"Given:\\n"
"   tcb1,tcb2   TCB as a 2-part Julian Date\\n"
"Returned:\\n"
"   tdb1,tdb2   TDB as a 2-part Julian Date"''',
     '''{
    double tcb1, tcb2, tdb1, tdb2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tcb1, &tcb2)) {
        return NULL;
    }
    status = eraTcbtdb(tcb1, tcb2, &tdb1, &tdb2);
    if (status) {
        PyErr_SetString(_erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tdb1, tdb2);
}
'''),
    ('tcgtt', 'METH_VARARGS',
     '''"\\ntcgtt(tcg1, tcg2) -> tt1, tt2\\n\\n"
"Time scale transformation:  Geocentric Coordinate Time, TCG, to\\n"
"Terrestrial Time, TT.\\n"
"Given:\\n"
"   tcg1,tcg2   TCG as a 2-part Julian Date\\n"
"Returned:\\n"
"   tt1,tt2   TT as a 2-part Julian Date"''',
     '''{
    double tcg1, tcg2, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tcg1, &tcg2)) {
        return NULL;
    }
    status = eraTcgtt(tcg1, tcg2, &tt1, &tt2);
    if (status) {
        PyErr_SetString(_erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tt1, tt2);
}
'''),
    ('tdbtcb', 'METH_VARARGS',
     '''"\\ntdbtcb(tdb1, tdb2) -> tcb1, tcb2\\n\\n"
"Time scale transformation:  Barycentric Dynamical Time, TDB, to\\n"
"Barycentric Coordinate Time, TCB.\\n"
"Given:\\n"
"   tdb1,tdb2   TDB as a 2-part Julian Date\\n"
"Returned:\\n"
"   tcb1,tcb2   TCB as a 2-part Julian Date"''',
     '''{
    double tdb1, tdb2, tcb1, tcb2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tdb1, &tdb2)) {
        return NULL;
    }
    status = eraTdbtcb(tdb1, tdb2, &tcb1, &tcb2);
    if (status) {
        PyErr_SetString(_erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tcb1, tcb2);
}
'''),
    ('tdbtt', 'METH_VARARGS',
     '''"\\ntdbtt(tdb1, tdb2, dtr) -> tt1, tt2\\n\\n"
"Time scale transformation: Barycentric Dynamical Time, TDB, to\\n"
"Terrestrial Time, TT.\\n"
"Given:\\n"
"   tdb1,tdb2   TDB as a 2-part Julian Date\\n"
"   dtr         TDB-TT in seconds\\n"
"Returned:\\n"
"   tt1,tt2   TT as a 2-part Julian Date"''',
     '''{
    double tdb1, tdb2, dtr, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &tdb1, &tdb2, &dtr)) {
        return NULL;
    }
    status = eraTdbtt(tdb1, tdb2, dtr, &tt1, &tt2);
    if (status) {
        PyErr_SetString(_erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tt1, tt2);
}
'''),
    ('tttai', 'METH_VARARGS',
     '''"\\ntttai(tt1, tt2) -> tai1, tai2\\n\\n"
"Time scale transformation: Terrestrial Time, TT, to\\n"
"International Atomic Time, TAI.\\n"
"Given:\\n"
"   tt1,tt2     TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   tai1,tai2   TAI as a 2-part Julian Date"''',
     '''{
    double tai1, tai2, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tt1, &tt2)) {
        return NULL;
    }
    status = eraTttai(tt1, tt2, &tai1, &tai2);
    if (status) {
        PyErr_SetString(_erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tai1, tai2);
}
'''),
    ('tttcg', 'METH_VARARGS',
     '''"\\ntttcg(tt1, tt2) -> tcg1, tcg2\\n\\n"
"Time scale transformation: Terrestrial Time, TT, to Geocentric\\n"
"Coordinate Time, TCG.\\n"
"Given:\\n"
"   tt1,tt2     TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   tcg1,tcg2   TCG as a 2-part Julian Date"''',
     '''{
    double tcg1, tcg2, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tt1, &tt2)) {
        return NULL;
    }
    status = eraTttcg(tt1, tt2, &tcg1, &tcg2);
    if (status) {
        PyErr_SetString(_erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tcg1, tcg2);
}
'''),
    ('tttdb', 'METH_VARARGS',
     '''"\\ntttdb(tt1, tt2, dtr) -> tdb1, tdb2\\n\\n"
"Time scale transformation: Terrestrial Time, TT, to\\n"
"Barycentric Dynamical Time, TDB.\\n"
"Given:\\n"
"   tt1,tt2     TT as a 2-part Julian Date\\n"
"   dtr         TDB-TT in seconds\\n"
"Returned:\\n"
"   tdb1,tdb2   TDB as a 2-part Julian Date"''',
     '''{
    double tdb1, tdb2, dtr, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &tt1, &tt2, &dtr)) {
        return NULL;
    }
    status = eraTttdb(tt1, tt2, dtr, &tdb1, &tdb2);
    if (status) {
        PyErr_SetString(_erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tdb1, tdb2);
}
'''),
    ('ttut1', 'METH_VARARGS',
     '''"\\nttut1(tt1, tt2, dt) -> ut11, ut12\\n\\n"
"Time scale transformation: Terrestrial Time, TT, to\\n"
"Universal Time UT1.\\n"
"Given:\\n"
"   tt1,tt2     TT as a 2-part Julian Date\\n"
"   dt          TT-UT1 in seconds\\n"
"Returned:\\n"
"   ut11,ut12   UT1 as a 2-part Julian Date"''',
     '''{
    double ut11, ut12, dt, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &tt1, &tt2, &dt)) {
        return NULL;
    }
    status = eraTtut1(tt1, tt2, dt, &ut11, &ut12);
    if (status) {
        PyErr_SetString(_erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", ut11, ut12);
}
'''),
    ('ut1tai', 'METH_VARARGS',
     '''"\\nut1tai(ut11, ut12, dta) -> tai1, tai2\\n\\n"
"Time scale transformation: Universal Time, UT1, to\\n"
"International Atomic Time, TAI.\\n"
"The argument dta, i.e. UT1-TAI, is an observed quantity, and is\\n"
"available from IERS tabulations.\\n"
"Given:\\n"
"   ut11,ut12     TT as a 2-part Julian Date\\n"
"   dta         UT1-TAI in seconds\\n"
"Returned:\\n"
"   tai1,tai2   TAI as a 2-part Julian Date"''',
     '''{
    double tai1, tai2, dta, ut11, ut12;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &ut11, &ut12, &dta)) {
        return NULL;
    }
    status = eraUt1tai(ut11, ut12, dta, &tai1, &tai2);
    if (status) {
        PyErr_SetString(_erfaError, "internal error...");
        return NULL;
    }
    return Py_BuildValue("dd", tai1, tai2);
}
'''),
    ('ut1tt', 'METH_VARARGS',
     '''"\\nut1tt(ut11, ut12, dt) -> tt1, tt2\\n\\n"
"Time scale transformation: Universal Time, UT1, to\\n"
"Terrestrial Time, TT.\\n"
"Given:\\n"
"   ut11,ut12   UT1 as a 2-part Julian Date\\n"
"   dt          TT-UT1 in seconds, dt is classical Delta T\\n"
"Returned:\\n"
"   tt1,tt2     TT as a 2-part Julian Date"''',
     '''{
    double tt1, tt2, dt, ut11, ut12;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &ut11, &ut12, &dt)) {
        return NULL;
    }
    status = eraUt1tt(ut11, ut12, dt, &tt1, &tt2);
    if (status) {
        PyErr_SetString(_erfaError, "internal error...");
        return NULL;
    }
    return Py_BuildValue("dd", tt1, tt2);
}
'''),
    ('ut1utc', 'METH_VARARGS',
     '''"\\nut1utc(ut11, ut12, dut1) -> utc1, utc2\\n\\n"
"Time scale transformation: Universal Time, UT1, to\\n"
"Coordinated Universal Time, UTC.\\n"
"Given:\\n"
"   ut11,ut12   UT1 as a 2-part Julian Date\\n"
"   dut1        UT1-UTC in seconds, Delta UT1\\n"
"Returned:\\n"
"   utc1,utc2   UTC as a 2-part Julian Date"''',
     '''{
    double utc1, utc2, dut1, ut11, ut12;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &ut11, &ut12, &dut1)) {
        return NULL;
    }
    status = eraUt1utc(ut11, ut12, dut1, &utc1, &utc2);
    if (status) {
        if (status == 1) {
            PyErr_SetString(_erfaError, "dubious year");
            return NULL;
        }
        else if (status == -1) {
            PyErr_SetString(_erfaError, "unacceptable date");
            return NULL;
        }
    }
    return Py_BuildValue("dd", utc1, utc2);
}
'''),
    ('utctai', 'METH_VARARGS',
     '''"\\nutctai(utc1, utc2) -> tai1, tai2\\n\\n"
"Time scale transformation: Coordinated Universal Time, UTC, to\\n"
"International Atomic Time, TAI.\\n"
"Given:\\n"
"   utc1,uc12   UTC as a 2-part Julian Date\\n"
"Returned:\\n"
"   tai1,tai2   TAI as a 2-part Julian Date"''',
     '''{
    double utc1, utc2, tai1, tai2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &utc1, &utc2)) {
        return NULL;
    }
    status = eraUtctai(utc1, utc2, &tai1, &tai2);
    if (status) {
        if (status == 1) {
            PyErr_SetString(_erfaError, "dubious year");
            return NULL;
        }
        else if (status == -1) {
            PyErr_SetString(_erfaError, "unacceptable date");
            return NULL;
        }
    }
    return Py_BuildValue("dd", tai1, tai2);
}
'''),
    ('utcut1', 'METH_VARARGS',
     '''"\\nutcut1(utc1, utc2, dut1) -> ut11, ut12\\n\\n"
"Time scale transformation: Coordinated Universal Time, UTC, to\\n"
"Universal Time, UT1.\\n"
"Given:\\n"
"   utc1,utc2   UTC as a 2-part Julian Date\\n"
"   dut1        UT1-UTC in seconds, Delta UT1\\n"
"Returned:\\n"
"   ut11,ut12   UT1 as a 2-part Julian Date"''',
     '''{
    double utc1, utc2, dut1, ut11, ut12;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &utc1, &utc2, &dut1)) {
        return NULL;
    }
    status = eraUtcut1(utc1, utc2, dut1, &ut11, &ut12);
    if (status) {
        if (status == 1) {
            PyErr_SetString(_erfaError, "dubious year");
            return NULL;
        }
        else if (status == -1) {
            PyErr_SetString(_erfaError, "unacceptable date");
            return NULL;
        }
    }
    return Py_BuildValue("dd", ut11, ut12);
}
'''),
    ('xy06', 'METH_VARARGS',
     '''"\\nxy06(d1, d2) -> x, y\\n\\n"
"X,Y coordinates of celestial intermediate pole from series based\\n"
"on IAU 2006 precession and IAU 2000A nutation.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   x,y     CIP X,Y coordinates"''',
     '''{
    double d1, d2, x, y;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraXy06(d1, d2, &x, &y);
    return Py_BuildValue("dd", x, y);
}
'''),
    ('xys00a', 'METH_VARARGS',
     '''"\\nxys00a(d1, d2) -> x, y, s\\n\\n"
"For a given TT date, compute the X,Y coordinates of the Celestial\\n"
"Intermediate Pole and the CIO locator s, using the IAU 2000A\\n"
"precession-nutation model.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   x,y     Celestial Intermediate Pole\\n"
"   s       the CIO locator s"''',
     '''{
    double d1, d2, x, y,  s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraXys00a(d1, d2, &x, &y, &s);
    return Py_BuildValue("ddd", x, y, s);
}
'''),
    ('xys00b', 'METH_VARARGS',
     '''"\\nxys00b(d1, d2) -> x, y, s\\n\\n"
"For a given TT date, compute the X,Y coordinates of the Celestial\\n"
"Intermediate Pole and the CIO locator s, using the IAU 2000B\\n"
"precession-nutation model.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   x,y     Celestial Intermediate Pole\\n"
"   s       the CIO locator s"''',
     '''{
    double d1, d2, x, y,  s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraXys00b(d1, d2, &x, &y, &s);
    return Py_BuildValue("ddd", x, y, s);
}
'''),
    ('xys06a', 'METH_VARARGS',
     '''"\\nxys06a(d1, d2) -> x, y, s\\n\\n"
"For a given TT date, compute the X,Y coordinates of the Celestial\\n"
"Intermediate Pole and the CIO locator s, using the IAU 2006\\n"
"precession and IAU 2000A nutation model.\\n"
"Given:\\n"
"   d1,d2   TT as a 2-part Julian Date\\n"
"Returned:\\n"
"   x,y     Celestial Intermediate Pole\\n"
"   s       the CIO locator s"''',
     '''{
    double d1, d2, x, y,  s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraXys06a(d1, d2, &x, &y, &s);
    return Py_BuildValue("ddd", x, y, s);
}
'''),

# vector-matrix library
    ('a2af', 'METH_VARARGS',
     '''"\\na2af(n, a) -> +/-, d, m, s, f\\n\\n"
"Decompose radians into degrees, arcminutes, arcseconds, fraction.\\n"
"Given:\\n"
"   n           resolution\\n"
"   a           angle in radians\\n"
"Returned:\\n"
"   sign        '+' or '-'\\n"
"   d           degrees\\n"
"   m           arcminutes\\n"
"   s           arcseconds\\n"
"   f           fraction"''',
     '''{
    int ndp, idmsf[4];
    char sign;
    double a;
    if (!PyArg_ParseTuple(args, "id", &ndp, &a)) {
        return NULL;
    }
    eraA2af(ndp, a, &sign, idmsf);
#if PY_VERSION_HEX >= 0x03000000
    return Py_BuildValue("Ciiii",sign,idmsf[0],idmsf[1],idmsf[2],idmsf[3]);
#else
    return Py_BuildValue("ciiii",sign,idmsf[0],idmsf[1],idmsf[2],idmsf[3]);
#endif
}
'''),
    ('a2tf', 'METH_VARARGS',
     '''"\\na2tf(n, a) -> +/-, h, m, s, f\\n\\n"
"Decompose radians into hours, minutes, seconds, fraction.\\n"
"Given:\\n"
"   n           resolution\\n"
"   a           angle in radians\\n"
"Returned:\\n"
"   sign        '+' or '-'\\n"
"   h           hours\\n"
"   m           minutes\\n"
"   s           seconds\\n"
"   f           fraction"''',
     '''{
    int ndp, ihmsf[4];
    char sign;
    double a;
    if (!PyArg_ParseTuple(args, "id", &ndp, &a)) {
        return NULL;
    }
    eraA2tf(ndp, a, &sign, ihmsf);
#if PY_VERSION_HEX >= 0x03000000
    return Py_BuildValue("Ciiii",sign,ihmsf[0],ihmsf[1],ihmsf[2],ihmsf[3]);
#else
    return Py_BuildValue("ciiii",sign,ihmsf[0],ihmsf[1],ihmsf[2],ihmsf[3]);
#endif
}
'''),
    ('af2a', 'METH_VARARGS',
     '''"\\naf2a(d, m, s) -> r\\n"
"Convert degrees, arcminutes, arcseconds to radians.\\n"
"Given:\\n"
/*"    sign       '-' = negative, otherwise positive\\n"*/
"    d          degrees\\n"
"    m          arcminutes\\n"
"    s          arcseconds\\n"
"Returned:\\n"
"    r          angle in radians"''',
     '''{
    int ideg, iamin;
    double asec, rad;
    int sign = '+';
    if (!PyArg_ParseTuple(args, "iid", &ideg, &iamin, &asec)) {
        return NULL;
    }
    if (ideg < 0) {
        sign = '-';
    }
    eraAf2a(sign, ideg, iamin, asec, &rad);
    return Py_BuildValue("d", rad);
}
'''),
    ('anp', 'METH_O',
     '''"\\nanp(a) -> 0 <= a < 2pi\\n\\n"
"Normalize angle into the range 0 <= a < 2pi.\\n"
"Given:\\n"
"    a          angle (radians)\\n"
"Returned:\\n"
"    a          angle in range 0-2pi"''',
     '''{
    double a = PyFloat_AsDouble(arg);
    if ((a == -1.0) && PyErr_Occurred()) {
        PyErr_SetString(_erfaError, "cannot convert angle to float!");
        return NULL;
    }
    return Py_BuildValue("d", eraAnp(a));
}
'''),
    ('anpm', 'METH_O',
     '''"\\nanpm(a) -> -pi <= a < +pi\\n\\n"
"Normalize angle into the range -pi <= a < +pi.\\n"
"Given:\\n"
"    a          angle (radians)\\n"
"Returned:\\n"
"    a          angle in range 0-2pi"''',
     '''{
    double a = PyFloat_AsDouble(arg);
    if ((a == -1.0) && PyErr_Occurred()) {
        PyErr_SetString(_erfaError, "cannot convert angle to float!");
        return NULL;
    }
    return Py_BuildValue("d", eraAnpm(a));
}
'''),
    ('c2s', 'METH_VARARGS',
     '''"\\nc2s(p) -> theta, phi\\n\\n"
"P-vector to spherical coordinates.\\n"
"Given:\\n"
"    p          p-vector\\n"
"Returned:\\n"
"    theta      longitude angle (radians)\\n"
"    phi        latitude angle (radians)"''',
     '''{
    double theta, phi, p[3];
    if (!PyArg_ParseTuple(args, "(ddd)", &p[0], &p[1], &p[2])) {
        return NULL;
    }
    eraC2s(p, &theta, &phi);
    return Py_BuildValue("dd", theta, phi);
}
'''),
    ('cp', 'METH_VARARGS',
     '''"\\ncp(p) -> c\\n\\n"
"Copy a p-vector.\\n"
"Given:\\n"
"   p           p-vector to be copied\\n"
"  Returned:\\n"
"   c           copy"''',
     '''{
    double p[3], c[3];
    if (!PyArg_ParseTuple(args, "(ddd)", &p[0], &p[1], &p[2])) {
        return NULL;
    }
    eraCp(p, c);
    return Py_BuildValue("ddd", c[0], c[1], c[2]);
}
'''),
    ('cpv', 'METH_VARARGS',
     '''"\\ncp(pv) -> c\\n\\n"
"Copy a position/velocity vector.\\n"
"Given:\\n"
"   pv          position/velocity vector to be copied\\n"
"  Returned:\\n"
"   c           copy"''',
     '''{
    double pv[2][3], c[2][3];
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
                                  &pv[0][0], &pv[0][1], &pv[0][2],
                                  &pv[1][0], &pv[1][1], &pv[1][2])) {
        return NULL;
    }
    eraCpv(pv, c);
    return Py_BuildValue("(ddd)(ddd)",
                           c[0][0],c[0][1],c[0][2],
                           c[1][0],c[1][1],c[1][2]);    
}
'''),
    ('cr', 'METH_VARARGS',
     '''"\\ncr(r) -> c\\n\\n"
"Copy an r-matrix.\\n"
"Given:\\n"
"   r           r-matrix to be copied\\n"
"  Returned:\\n"
"   c           copy"''',
     '''{
    double r[3][3], c[3][3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))",
                                   &r[0][0], &r[0][1], &r[0][2],
                                   &r[1][0], &r[1][1], &r[1][2],
                                   &r[2][0], &r[2][1], &r[2][2])) {
        return NULL;
    }
    eraCr(r, c);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            c[0][0],c[0][1],c[0][2],
                            c[1][0],c[1][1],c[1][2],
                            c[2][0],c[2][1],c[2][2]);    
}
'''),
    ('d2tf', 'METH_VARARGS',
     '''"\\nd2tf(n, d) -> +/-,h, m, s, f\\n\\n"
"Decompose days to hours, minutes, seconds, fraction.\\n"
"Given:\\n"
"    n          resolution\\n"
"    d          interval in days\\n"
"Returned:\\n"
"    sign       '+' or '-'\\n"
"    h          hours\\n"
"    m          minutes\\n"
"    s          seconds\\n"
"    f          fraction"''',
     '''{
    double d;
    char sign = '+';
    int ndp, df[4];
    if (! PyArg_ParseTuple(args, "id",&ndp, &d)) {
        return NULL;
    }
    eraD2tf(ndp, d, &sign, df);
#if PY_VERSION_HEX >= 0x03000000
    return Py_BuildValue("Ciiii", sign,df[0],df[1],df[2],df[3]);
#else
    return Py_BuildValue("ciiii", sign,df[0],df[1],df[2],df[3]);
#endif
}
'''),
##
## Ir see erfa.py
##
    ('p2pv', 'METH_VARARGS',
     '''"\\np2pv(p) -> pv\\n\\n"
"Extend a p-vector to a pv-vector by appending a zero velocity.\\n"
"Given:\\n"
"    p          p-vector\\n"
"Returned:\\n"
"    pv         pv-vector"''',
     '''{
    double p[3], pv[2][3];
    if (!PyArg_ParseTuple(args, "(ddd)", &p[0], &p[1], &p[2])) {
        return NULL;
    }
    eraP2pv(p, pv);
    return Py_BuildValue("(ddd)(ddd)",
    pv[0][0],pv[0][1],pv[0][2],pv[1][0],pv[1][1],pv[1][2]);
}
'''),
    ('p2s', 'METH_VARARGS',
     '''"\\np2s(p) -> theta, phi, r\\n\\n"
"P-vector to spherical polar coordinates.\\n"
"Given:\\n"
"    p          p-vector\\n"
"Returned:\\n"
"    theta      longitude angle (radians)\\n"
"    phi        latitude angle (radians)\\n"
"    r          radial distance"''',
     '''{
    double p[3], theta, phi, r;
    if (!PyArg_ParseTuple(args, "(ddd)", &p[0], &p[1], &p[2])) {
        return NULL;
    }
    eraP2s(p, &theta, &phi, &r);
    return Py_BuildValue("ddd", theta, phi, r);
}
'''),
    ('pap', 'METH_VARARGS',
     '''"\\npap(a, b) -> theta\\n\\n"
"Position-angle from two p-vectors.\\n"
"Given:\\n"
"   a       direction of reference point\\n"
"   b       direction of point whose PA is required\\n"
"Returned:\\n"
"   theta   position angle of b with respect to a (radians)"''',
     '''{
    double a[3], b[3], theta;
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
                                  &a[0], &a[1], &a[2],
                                  &b[0], &b[1], &b[2])) {
        return NULL;
    }
    theta = eraPap(a, b);
    return Py_BuildValue("d", theta);
}
'''),
    ('pas', 'METH_VARARGS',
     '''"\\npas(al, ap, bl, bp) -> p\\n\\n"
"Position-angle from spherical coordinates.\\n"
"Given:\\n"
"   al  longitude of point A (e.g. RA) in radians\\n"
"   ap  latitude of point A (e.g. Dec) in radians\\n"
"   bl  longitude of point B\\n"
"   bp  latitude of point B\\n"
"Returned:\\n"
"   p   position angle of B with respect to A"''',
    '''{
    double al, ap, bl, bp, p;
    if (!PyArg_ParseTuple(args, "dddd",&al, &ap, &bl, &bp)) {
        return NULL;
    }
    p = eraPas(al, ap, bl, bp);
    return Py_BuildValue("d", p);
}
'''),
    ('pdp', 'METH_VARARGS',
     '''"\\npdp(a, b -> a.b\\n\\n"
"p-vector inner (=scalar=dot) product.\\n"
"Given:\\n"
"   a       first p-vector\\n"
"   b       second p-vector\\n"
"Returned:\\n"
"   ab      a . b"''',
    '''{
    double a[3], b[3], ab;
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
                                  &a[0], &a[1], &a[2],
                                  &b[0], &b[1], &b[2])) {
        return NULL;
    }
    ab = eraPdp(a,b);
    return Py_BuildValue("d", ab);
}
'''),
    ('pm', 'METH_VARARGS',
     '''"\\npm(p) -> modulus\\n\\n"
"Modulus of p-vector.\\n"
"Given:\\n"
"   p       p-vector\\n"
"Returned:\\n"
"   m       modulus"''',
     '''{
    double p[3], m;
    if (!PyArg_ParseTuple(args, "(ddd)",&p[0], &p[1], &p[2])) {
        return NULL;
    }
    m = eraPm(p);
    return Py_BuildValue("d", m);
}
'''),
    ('pmp', 'METH_VARARGS',
     '''"\\npmp(a, b) -> amb = a-b\\n\\n"
"P-vector subtraction.\\n"
"Given:\\n"
"   a           first p-vector\\n"
"   b           second p-vector\\n"
"Returned:\\n"
"   amb         a - b"''',
     '''{
    double a[3], b[3], amb[3];
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
                                  &a[0], &a[1], &a[2],
                                  &b[0], &b[1], &b[2])) {
        return NULL;
    }
    eraPmp(a, b, amb);
    return Py_BuildValue("ddd", amb[0], amb[1], amb[2]);
}
'''),
    ('pn', 'METH_VARARGS',
     '''"\\npn(p) -> r,u\\n\\n"
"Convert a p-vector into modulus and unit vector.\\n"
"Given:\\n"
"   p           p-vector\\n"
"Returned:\\n"
"   r           modulus\\n"
"   u           unit vector"''',
     '''{
    double p[3], r, u[3];
    if (!PyArg_ParseTuple(args, "(ddd)",&p[0], &p[1], &p[2])) {
        return NULL;
    }
    eraPn(p, &r, u);
    return Py_BuildValue("d(ddd)", r, u[0], u[1], u[2]);    
}
'''),
    ('ppp', 'METH_VARARGS',
     '''"\\nppp(a, b) -> apb = a+b\\n\\n"
"P-vector addition.\\n"
"Given:\\n"
"   a           first p-vector\\n"
"   b           second p-vector\\n"
"Returned:\\n"
"   apb         a + b"''',
     '''{
    double a[3], b[3], apb[3];
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
                                  &a[0], &a[1], &a[2],
                                  &b[0], &b[1], &b[2])) {
        return NULL;
    }
    eraPpp(a, b, apb);
    return Py_BuildValue("ddd", apb[0], apb[1], apb[2]);
}
'''),
    ('ppsp', 'METH_VARARGS',
     '''"\\nppsp(a, s, b) -> apsb = a + s*b\\n\\n"
"P-vector plus scaled p-vector.\\n"
"Given:\\n"
"   a           first p-vector\\n"
"   s           scalar (multiplier for b)\\n"
"   b           second p-vector\\n"
"Returned:\\n"
"   apsb        a + s*b"''',
     '''{
    double s, a[3], b[3], apsb[3];
    if (!PyArg_ParseTuple(args, "(ddd)d(ddd)",
                                  &a[0], &a[1], &a[2],
                                  &s,
                                  &b[0], &b[1], &b[2])) {
        return NULL;
    }
    eraPpsp(a, s, b, apsb);
    return Py_BuildValue("ddd", apsb[0], apsb[1], apsb[2]);
}
'''),
    ('pv2p', 'METH_VARARGS',
     '''"\\npv2p(pv) -> p\\n\\n"
"Discard velocity component of a pv-vector.\\n"
"Given:\\n"
"   pv          pv-vector\\n"
"Returned:\\n"
"   p           p-vector"''',
     '''{
    double pv[2][3], p[3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))",
                                   &pv[0][0], &pv[0][1], &pv[0][2],
                                   &pv[1][0], &pv[1][1], &pv[1][2])) {
        return NULL;
    }
    eraPv2p(pv, p);
    return Py_BuildValue("ddd", p[0],p[1],p[2]);    
}
'''),
    ('pv2s', 'METH_VARARGS',
     '''"\\npv2s(pv) -> theta, phi, r, td, pd, rd\\n\\n"
"Convert position/velocity from Cartesian to spherical coordinates.\\n"
"Given:\\n"
"   pv          pv-vector\\n"
"Returned:\\n"
"   theta       longitude angle (radians)\\n"
"   phi         latitude angle (radians)\\n"
"   r           radial distance\\n"
"   td          rate of change of theta\\n"
"   pd          rate of change of phi\\n"
"   rd          rate of change of r"''',
     '''{
    double pv[2][3], theta, phi, r, td, pd, rd;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))",
                                   &pv[0][0], &pv[0][1], &pv[0][2],
                                   &pv[1][0], &pv[1][1], &pv[1][2])) {
        return NULL;
    }
    eraPv2s(pv, &theta, &phi, &r, &td, &pd, &rd);
    return Py_BuildValue("dddddd", theta, phi, r, td, pd, rd);
}
'''),
    ('pvdpv', 'METH_VARARGS',
     '''"\\npvdpv(a, b) -> adb = a.b\\n\\n"
"Inner (=scalar=dot) product of two pv-vectors.\\n"
"Given:\\n"
"   a           first pv-vector\\n"
"   b           second pv-vector\\n"
"Returned:\\n"
"   adb         a . b"''',
     '''{
    double pva[2][3], pvb[2][3], adb[2];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))((ddd)(ddd))",
                                   &pva[0][0], &pva[0][1], &pva[0][2],
                                   &pva[1][0], &pva[1][1], &pva[1][2],
                                   &pvb[0][0], &pvb[0][1], &pvb[0][2],
                                   &pvb[1][0], &pvb[1][1], &pvb[1][2])) {
        return NULL;
    }
    eraPvdpv(pva, pvb, adb);
    return Py_BuildValue("dd", adb[0], adb[1]);
}
'''),
    ('pvm', 'METH_VARARGS',
     '''"\\npvm(pv) -> r,s\\n\\n"
"Modulus of pv-vector.\\n"
"Given:\\n"
"   pv          pv-vector\\n"
"Returned:\\n"
"   r           modulus of position component\\n"
"   s           modulus of velocity component"''',
     '''{
    double pv[2][3], r, s;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))",
                                   &pv[0][0], &pv[0][1], &pv[0][2],
                                   &pv[1][0], &pv[1][1], &pv[1][2])) {
        return NULL;
    }
    eraPvm(pv, &r, &s);
    return Py_BuildValue("dd", r, s);
}
'''),
    ('pvmpv', 'METH_VARARGS',
     '''"\\npvmpv(a, b) -> amb = a-b\\n\\n"
"Subtract one pv-vector from another.\\n"
"Given:\\n"
"   a           first pv-vector\\n"
"   b           second pv-vector\\n"
"Returned:\\n"
"   amb         a - b"''',
     '''{
    double pva[2][3], pvb[2][3], amb[2][3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))((ddd)(ddd))",
                                   &pva[0][0], &pva[0][1], &pva[0][2],
                                   &pva[1][0], &pva[1][1], &pva[1][2],
                                   &pvb[0][0], &pvb[0][1], &pvb[0][2],
                                   &pvb[1][0], &pvb[1][1], &pvb[1][2])) {
        return NULL;
    }
    eraPvmpv(pva, pvb, amb);
    return Py_BuildValue("(ddd)(ddd)",
                           amb[0][0], amb[0][1], amb[0][2],
                           amb[1][0], amb[1][1], amb[1][2]);
}
'''),
    ('pvppv', 'METH_VARARGS',
     '''"\\npvppv(a, b) -> apb = a+b\\n\\n"
"Add one pv-vector to another.\\n"
"Given:\\n"
"   a           first pv-vector\\n"
"   b           second pv-vector\\n"
"Returned:\\n"
"   apb         a + b"''',
     '''{
    double pva[2][3], pvb[2][3], apb[2][3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))((ddd)(ddd))",
                                   &pva[0][0], &pva[0][1], &pva[0][2],
                                   &pva[1][0], &pva[1][1], &pva[1][2],
                                   &pvb[0][0], &pvb[0][1], &pvb[0][2],
                                   &pvb[1][0], &pvb[1][1], &pvb[1][2])) {
        return NULL;
    }
    eraPvppv(pva, pvb, apb);
    return Py_BuildValue("(ddd)(ddd)",
                           apb[0][0], apb[0][1], apb[0][2],
                           apb[1][0], apb[1][1], apb[1][2]);
}
'''),
    ('pvu', 'METH_VARARGS',
     '''"\\npvu(dt, pv) -> upv\\n\\n"
"Update a pv-vector.\\n"
"Given:\\n"
"   dt          time interval\\n"
"   pv          pv-vector\\n"
"Returned:\\n"
"   upv         p updated, v unchanged"''',
     '''{
    double pv[2][3], dt, upv[2][3];
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd))", &dt,
                                   &pv[0][0], &pv[0][1], &pv[0][2],
                                   &pv[1][0], &pv[1][1], &pv[1][2])) {
        return NULL;
    }
    eraPvu(dt, pv, upv);
    return Py_BuildValue("(ddd)(ddd)",
                           upv[0][0], upv[0][1], upv[0][2],
                           upv[1][0], upv[1][1], upv[1][2]);
}
'''),
    ('pvup', 'METH_VARARGS',
     '''"\\npvup(dt, pv) -> p\\n\\n"
"Update a pv-vector, discarding the velocity component.\\n"
"Given:\\n"
"   dt          time interval\\n"
"   pv          pv-vector\\n"
"Returned:\\n"
"   p           p-vector"''',
     '''{
    double pv[2][3], dt, p[3];
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd))", &dt,
                                   &pv[0][0], &pv[0][1], &pv[0][2],
                                   &pv[1][0], &pv[1][1], &pv[1][2])) {
        return NULL;
    }
    eraPvup(dt, pv, p);
    return Py_BuildValue("ddd", p[0], p[1], p[2]);
}
'''),
    ('pvxpv', 'METH_VARARGS',
     '''"\\npvxpv(a, b) -> axb = a x b\\n\\n"
"Outer (=vector=cross) product of two pv-vectors.\\n"
"Given:\\n"
"   a           first pv-vector\\n"
"   b           second pv-vector\\n"
"Returned:\\n"
"   axb         a x b"''',
     '''{
    double pva[2][3], pvb[2][3], axb[2][3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))((ddd)(ddd))",
                                   &pva[0][0], &pva[0][1], &pva[0][2],
                                   &pva[1][0], &pva[1][1], &pva[1][2],
                                   &pvb[0][0], &pvb[0][1], &pvb[0][2],
                                   &pvb[1][0], &pvb[1][1], &pvb[1][2])) {
        return NULL;
    }
    eraPvxpv(pva, pvb, axb);
    return Py_BuildValue("(ddd)(ddd)",
                           axb[0][0], axb[0][1], axb[0][2],
                           axb[1][0], axb[1][1], axb[1][2]);
}
'''),
    ('pxp', 'METH_VARARGS',
     '''"\\npxp(a, b) -> axb = a x b\\n\\n"
"p-vector outer (=vector=cross) product.\\n"
"Given:\\n"
"   a           first p-vector\\n"
"   b           second p-vector\\n"
"Returned:\\n"
"   axb         a x b"''',
     '''{
    double a[3], b[3], axb[3];
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
                                  &a[0], &a[1], &a[2],
                                  &b[0], &b[1], &b[2])) {
        return NULL;
    }
    eraPxp(a, b, axb);
    return Py_BuildValue("ddd", axb[0], axb[1], axb[2]);
}
'''),
    ('rm2v', 'METH_VARARGS',
     '''"\\nrm2v(r) -> w\\n\\n"
"Express an r-matrix as an r-vector.\\n"
"Given:\\n"
"   r          rotation matrix\\n"
"Returned:\\n"
"   w          rotation vector"''',
     '''{
    double r[3][3], w[3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))",
                                   &r[0][0], &r[0][1], &r[0][2],
                                   &r[1][0], &r[1][1], &r[1][2],
                                   &r[2][0], &r[2][1], &r[2][2])) {
        return NULL;
    }
    eraRm2v(r, w);
    return Py_BuildValue("ddd", w[0], w[1], w[2]);    
}
'''),
    ('rv2m', 'METH_VARARGS',
     '''"\\nrv2m(w) -> r\\n\\n"
"Form the r-matrix corresponding to a given r-vector.\\n"
"Given:\\n"
"   w           rotation vector\\n"
"Returned:\\n"
"   r           rotation matrix"''',
     '''{    
    double r[3][3], w[3];
    if (!PyArg_ParseTuple(args, "(ddd)", &w[0], &w[1], &w[2])) {
        return NULL;
    }
    eraRv2m(w, r);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            r[0][0],r[0][1],r[0][2],
                            r[1][0],r[1][1],r[1][2],
                            r[2][0],r[2][1],r[2][2]);        
}
'''),
    ('rx', 'METH_VARARGS',
     '''"\\nrx(phi, r) -> r\\n\\n"
"Rotate an r-matrix about the x-axis.\\n"
"Given:\\n"
"   phi         angle (radians)\\n"
"Given and returned:\\n"
"   r           r-matrix, rotated"''',
     '''{
    double r[3][3], phi;
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd)(ddd))", &phi,
                                   &r[0][0], &r[0][1], &r[0][2],
                                   &r[1][0], &r[1][1], &r[1][2],
                                   &r[2][0], &r[2][1], &r[2][2])) {
        return NULL;
    }
    eraRx(phi, r);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            r[0][0],r[0][1],r[0][2],
                            r[1][0],r[1][1],r[1][2],
                            r[2][0],r[2][1],r[2][2]);    
}
'''),
    ('rxp', 'METH_VARARGS',
     '''"\\nrxp(r, p) -> rp\\n\\n"
"Multiply a p-vector by an r-matrix.\\n"
"Given:\\n"
"   r           r-matrix\\n"
"   p           p-vector\\n"
"Returned:\\n"
"   rp          r * p"''',
     '''{
    double r[3][3], p[3], rp[3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))(ddd)", 
                                   &r[0][0], &r[0][1], &r[0][2],
                                   &r[1][0], &r[1][1], &r[1][2],
                                   &r[2][0], &r[2][1], &r[2][2],
                                   &p[0], &p[1], &p[2])) {
        return NULL;
    }
    eraRxp(r, p, rp);
    return Py_BuildValue("ddd", rp[0], rp[1], rp[2]);
}
'''),
    ('rxpv', 'METH_VARARGS',
     '''"\\nrxpv(r, pv) -> rpv\\n\\n"
"Multiply a pv-vector by an r-matrix.\\n"
"Given:\\n"
"   r           r-matrix\\n"
"   pv          pv-vector\\n"
"Returned:\\n"
"   rpv         r * pv"''',
     '''{
    double r[3][3], pv[2][3], rpv[2][3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))((ddd)(ddd))", 
                                   &r[0][0], &r[0][1], &r[0][2],
                                   &r[1][0], &r[1][1], &r[1][2],
                                   &r[2][0], &r[2][1], &r[2][2],
                                   &pv[0][0], &pv[0][1], &pv[0][2],
                                   &pv[1][0], &pv[1][1], &pv[1][2])) {
        return NULL;
    }
    eraRxpv(r, pv, rpv);
    return Py_BuildValue("(ddd)(ddd)",
                          rpv[0][0], rpv[0][1], rpv[0][2],
                          rpv[1][0], rpv[1][1], rpv[1][2]);
}
'''),
    ('rxr', 'METH_VARARGS',
     '''"\\nrxr(a, b -> atb\\n\\n"
"Multiply two r-matrices.\\n"
"Given:\\n"
"   a           first r-matrix\\n"
"   b           second r-matrix\\n"
"Returned:\\n"
"   atb         a * b"''',
     '''{
    double a[3][3], b[3][3], atb[3][3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
                                   &a[0][0], &a[0][1], &a[0][2],
                                   &a[1][0], &a[1][1], &a[1][2],
                                   &a[2][0], &a[2][1], &a[2][2],
                                   &b[0][0], &b[0][1], &b[0][2],
                                   &b[1][0], &b[1][1], &b[1][2],
                                   &b[2][0], &b[2][1], &b[2][2])) {
        return NULL;
    }
    eraRxr(a, b, atb);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                           atb[0][0],atb[0][1],atb[0][2],
                           atb[1][0],atb[1][1],atb[1][2],
                           atb[2][0],atb[2][1],atb[2][2]);    
}
'''),
    ('ry', 'METH_VARARGS',
     '''"\\nry(theta, r) -> r\\n\\n"
"Rotate an r-matrix about the y-axis.\\n"
"Given:\\n"
"   theta       angle (radians)\\n"
"Given and returned:\\n"
"   r           r-matrix, rotated"''',
     '''{
    double r[3][3], theta;
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd)(ddd))", &theta,
                                 &r[0][0], &r[0][1], &r[0][2],
                                 &r[1][0], &r[1][1], &r[1][2],
                                 &r[2][0], &r[2][1], &r[2][2])) {
        return NULL;
    }
    eraRy(theta, r);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            r[0][0],r[0][1],r[0][2],
                            r[1][0],r[1][1],r[1][2],
                            r[2][0],r[2][1],r[2][2]);    
}
'''),
    ('rz', 'METH_VARARGS',
     '''"\\nrz(psi, r) -> r\\n\\n"
"Rotate an r-matrix about the z-axis.\\n"
"Given:\\n"
"   psi         angle (radians)\\n"
"Given and returned:\\n"
"   r           r-matrix, rotated"''',
     '''{
    double r[3][3], psi;
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd)(ddd))", &psi,
                                 &r[0][0], &r[0][1], &r[0][2],
                                 &r[1][0], &r[1][1], &r[1][2],
                                 &r[2][0], &r[2][1], &r[2][2])) {
        return NULL;
    }
    eraRz(psi, r);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            r[0][0],r[0][1],r[0][2],
                            r[1][0],r[1][1],r[1][2],
                            r[2][0],r[2][1],r[2][2]);    
}
'''),
    ('s2c', 'METH_VARARGS',
     '''"\\ns2c(theta, phi) -> c\\n\\n"
"Convert spherical coordinates to Cartesian.\\n"
"Given:\\n"
"    theta   longitude angle (radians)\\n"
"    phi     latitude angle (radians)\\n"
"Returned:\\n"
"    c       direction cosines"''',
     '''{
    double theta, phi, c[3];
    if (!PyArg_ParseTuple(args, "dd", &theta, &phi)) {
        return NULL;
    }
    eraS2c(theta, phi, c);
    return Py_BuildValue("ddd", c[0], c[1], c[2]);
}
'''),
    ('s2p', 'METH_VARARGS',
     '''"\\ns2p(theta, phi, r) -> p\\n\\n"
"Convert spherical polar coordinates to p-vector.\\n"
"Given:\\n"
"   theta   longitude angle (radians)\\n"
"   phi     latitude angle (radians)\\n"
"   r       radial distance\\n"
"Returned:\\n"
"   p       direction cosines"''',
     '''{
    double theta, phi, r, p[3];
    if (!PyArg_ParseTuple(args, "ddd", &theta, &phi, &r)) {
        return NULL;
    }
    eraS2p(theta, phi, r, p);
    return Py_BuildValue("ddd", p[0], p[1], p[2]);
}
'''),
    ('s2pv', 'METH_VARARGS',
     '''"\\ns2pv(theta, phi, r, td, pd, rd) -> pv\\n\\n"
"Convert position/velocity from spherical to Cartesian coordinates.\\n"
"Given:\\n"
"   theta   longitude angle (radians)\\n"
"   phi     latitude angle (radians)\\n"
"   r       radial distance\\n"
"   td      rate of change of theta\\n"
"   pd      rate of change of phi\\n"
"   rd      rate of change of r\\n"
"Returned:\\n"
"   pv      pv-vector"''',
     '''{
    double theta, phi, r, td, pd, rd, pv[2][3];
    if (!PyArg_ParseTuple(args, "dddddd",
        &theta, &phi, &r, &td, &pd, &rd)) {
        return NULL;
    }
    eraS2pv(theta, phi, r, td, pd, rd, pv);
    return Py_BuildValue("(ddd)(ddd)",
                        pv[0][0], pv[0][1], pv[0][2],
                        pv[1][0], pv[1][1], pv[1][2]);
}
'''),
    ('s2xpv', 'METH_VARARGS',
     '''"\\ns2xpv(s1, s2, pv) -> spv\\n\\n"
"Multiply a pv-vector by two scalars.\\n"
"Given:\\n"
"   s1          scalar to multiply position component by\\n"
"   s2          scalar to multiply velocity component by\\n"
"   pv          pv-vector\\n"
"Returned:\\n"
"   spv         pv-vector: p scaled by s1, v scaled by s2"''',
     '''{
    double s1, s2, pv[2][3], spv[2][3];
    if (!PyArg_ParseTuple(args, "dd((ddd)(ddd))", &s1, &s2,
                                 &pv[0][0], &pv[0][1], &pv[0][2],
                                 &pv[1][0], &pv[1][1], &pv[1][2])) {
        return NULL;
    }
    eraS2xpv(s1, s2, pv, spv);
    return Py_BuildValue("(ddd)(ddd)",
                        spv[0][0], spv[0][1], spv[0][2],
                        spv[1][0], spv[1][1], spv[1][2]);
}
'''),
    ('sepp', 'METH_VARARGS',
     '''"\\nsepp(a, b) -> s\\n\\n"
"Angular separation between two p-vectors.\\n"
"Given:\\n"
"   a       first p-vector (not necessarily unit length)\\n"
"   b       second p-vector (not necessarily unit length)\\n"
"Returned:\\n"
"   s       angular separation (radians, always positive)"''',
     '''{
    double a[3], b[3], s;
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
                                  &a[0], &a[1], &a[2],
                                  &b[0], &b[1], &b[2])) {
        return NULL;
    }
    s = eraSepp(a, b);
    return Py_BuildValue("d", s);
}
'''),
    ('seps', 'METH_VARARGS',
    '''"seps(al, ap, bl, bp) -> s\\n\\n"
"Angular separation between two sets of spherical coordinates.\\n"
"Given:\\n"
"   al  first longitude (radians)\\n"
"   ap  first latitude (radians)\\n"
"   bl  second longitude (radians)\\n"
"   bp  second latitude (radians)\\n"
"Returned:\\n"
"   s   angular separation (radians)"''',
    '''{
    double al, ap, bl, bp, s;
    if (!PyArg_ParseTuple(args, "dddd",
        &al, &ap, &bl, &bp)) {
        return NULL;
    }
    s = eraSeps(al, ap, bl, bp);
    return Py_BuildValue("d", s);
}
'''),
    ('sxp', 'METH_VARARGS',
     '''"\\nsxp(s, p) -> sp\\n\\n"
"Multiply a p-vector by a scalar.\\n"
"Given:\\n"
"   s           scalar\\n"
"   p           p-vector\\n"
"Returned:\\n"
"   sp          s * p"''',
     '''{
    double s, p[3], sp[3];
    if (!PyArg_ParseTuple(args, "d(ddd)", &s,
                                 &p[0], &p[1], &p[2])) {
        return NULL;
    }
    eraSxp(s, p, sp);
    return Py_BuildValue("ddd", sp[0], sp[1], sp[2]);
}
'''),
    ('sxpv', 'METH_VARARGS',
     '''"\\nsxpv(s, pv) -> spv\\n\\n"
"Multiply a pv-vector by a scalar.\\n"
"Given:\\n"
"   s           scalar\\n"
"   pv          pv-vector\\n"
"Returned:\\n"
"   spv         s * pv"''',
     '''{
    double s, pv[2][3], spv[2][3];
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd))", &s,
                                 &pv[0][0], &pv[0][1], &pv[0][2],
                                 &pv[1][0], &pv[1][1], &pv[1][2])) {
        return NULL;
    }
    eraSxpv(s, pv, spv);
    return Py_BuildValue("(ddd)(ddd)",
                        spv[0][0], spv[0][1], spv[0][2],
                        spv[1][0], spv[1][1], spv[1][2]);
}
'''),
    ('tf2a', 'METH_VARARGS',
     '''"\\ntf2a(hour, min, sec) -> rad\\n\\n"
"Convert hours, minutes, seconds to radians.\\n"
"Given:\\n"
"   hour    hours\\n"
"   min     minutes\\n"
"   sec     seconds\\n"
"Returned:\\n"
"   rad     angle in radians"''',
     '''{
    int ihour, imin, status;
    double sec, rad;
    char s = '+';
    if (!PyArg_ParseTuple(args, "iid", &ihour, &imin, &sec)) {
        return NULL;
    }
    if (ihour < 0) {
        s = '-';
        /*ihour = -1 * ihour;*/
    }
    status = eraTf2a(s, ihour, imin, sec, &rad);
    if (status == 1) {
        PyErr_WarnEx(PyExc_Warning, "hour outside range 0-23", 1);
    }
    else if (status == 2) {
        PyErr_WarnEx(PyExc_Warning, "min outside range 0-59", 1);
    }
    else if (status == 3) {
        PyErr_WarnEx(PyExc_Warning, "sec outside range 0-59", 1);
    }
    return Py_BuildValue("d", rad);
}
'''),
    ('tf2d', 'METH_VARARGS',
     '''"\\ntf2d(hour, min, sec) -> days\\n\\n"
"Convert hours, minutes, seconds to days.\\n"
"Given:\\n"
"   hour    hours\\n"
"   min     minutes\\n"
"   sec     seconds\\n"
"Returned:\\n"
"   days    interval in days"''',
     '''{
    int ihour, imin, status;
    double sec, days;
    char s = '+';
    if (!PyArg_ParseTuple(args, "iid", &ihour, &imin, &sec)) {
        return NULL;
    }
    if (ihour < 0)  {
        s = '-';
        ihour = -1 * ihour;
    }
    status = eraTf2d(s, ihour, imin, sec, &days);
    if (status == 1) {
        PyErr_WarnEx(PyExc_Warning, "hour outside range 0-23", 1);
    }
    else if (status == 2) {
        PyErr_WarnEx(PyExc_Warning, "min outside range 0-59", 1);
    }
    else if (status == 3) {
        PyErr_WarnEx(PyExc_Warning, "sec outside range 0-59", 1);
    }
    return Py_BuildValue("d", days);
}
'''),
    ('tr', 'METH_VARARGS',
     '''"\\ntr(r) -> rt\\n\\n"
"Transpose an r-matrix.\\n"
"Given:\\n"
"   r           r-matrix\\n"
"Given and returned:\\n"
"   rt          transpose"''',
     '''{
    double r[3][3], rt[3][3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))",
                                   &r[0][0], &r[0][1], &r[0][2],
                                   &r[1][0], &r[1][1], &r[1][2],
                                   &r[2][0], &r[2][1], &r[2][2])) {
        return NULL;
    }
    eraTr(r, rt);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            rt[0][0],rt[0][1],rt[0][2],
                            rt[1][0],rt[1][1],rt[1][2],
                            rt[2][0],rt[2][1],rt[2][2]);    
}
'''),
    ('trxp', 'METH_VARARGS',
     '''"\\ntrxp(r, p) -> trp\\n\\n"
"Multiply a p-vector by the transpose of an r-matrix.\\n"
"Given:\\n"
"   r           r-matrix\\n"
"   p           p-vector\\n"
"Given and returned:\\n"
"   trp         r * p"''',
     '''{
    double r[3][3], p[3], trp[3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))(ddd)", 
                                   &r[0][0], &r[0][1], &r[0][2],
                                   &r[1][0], &r[1][1], &r[1][2],
                                   &r[2][0], &r[2][1], &r[2][2],
                                   &p[0], &p[1], &p[2])) {
        return NULL;
    }
    eraTrxp(r, p, trp);
    return Py_BuildValue("ddd", trp[0], trp[1], trp[2]);    
}
'''),
    ('trxpv', 'METH_VARARGS',
     '''"\\ntrxpv(r, pv) -> trpv\\n\\n"
"Multiply a pv-vector by the transpose of an r-matrix.\\n"
"Given:\\n"
"   r           r-matrix\\n"
"   pv          pv-vector\\n"
"Given and returned:\\n"
"   trpv        r * pv"''',
     '''{
    double r[3][3], pv[2][3], trpv[2][3];
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))((ddd)(ddd))", 
                                   &r[0][0], &r[0][1], &r[0][2],
                                   &r[1][0], &r[1][1], &r[1][2],
                                   &r[2][0], &r[2][1], &r[2][2],
                                   &pv[0][0], &pv[0][1], &pv[0][2],
                                   &pv[1][0], &pv[1][1], &pv[1][2])) {
        return NULL;
    }
    eraTrxpv(r, pv, trpv);
    return Py_BuildValue("(ddd)(ddd)",
                        trpv[0][0], trpv[0][1], trpv[0][2],
                        trpv[1][0], trpv[1][1], trpv[1][2]);
}
'''),
##
## Zp, Zpv and Zr see erfa.py
##
    ]


module = Builder(module_name, include_files,
                 module_localprototype, module_localfunction,
                 module_add, module_methods,
                 module_doc, module_filename,
                 module_struct_seq)
module.create()
