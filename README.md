# MCML

Last update: 3/1/2000

MCML is a Monte Carlo simulation program for Multi-layered Turbid
Media with an infinitely narrow photon beam as the light source. The
simulation is specified by an input text file called, for example,
"sample.mci", which can be altered by any simple text editor.  The
output is another text file called, for example, "sample.mco". (The
names are arbitrary.)

CONV is a convolution program which uses the MCML output file to
convolve for photon beams of variable size or shape (Gaussian or flat
field).  CONV can provide a variety of output formats (reflectance,
transmission, iso-fluence contours, etc.) which are compatible with
standard graphics applications.


HISTORY
======================================================================
MCML 1.2 corrected a bug that sometimes led to a sqrt() domain error
due to the finite machine precision.  However, this error has not
affected previously computed results because if this error happened,
the program would have stopped execution and yielded no output.

MCML 1.2.1 corrected a bug that sometimes led to memory problems when
the grid size was too small.

MCML 1.2.2 corrected a bug that was introduced in 1.2.1.

CONV 1.1 corrected a bug of the qtrap() function in the convnr.c file.
The integration by the original qtrap() sometimes converged incorrectly 
to zero when the Monte Carlo data was noisy.

======================================================================

Lihong Wang, Ph.D.<br />
Bren Professor of Medical Engineering and Electrical Engineering<br />
Affiliated Faculty of Physics and Applied Physics<br />
Andrew and Peggy Cherng Medical Engineering Leadership Chair<br />
Executive Officer (aka Department Chair) of Medical Engineering<br />
Andrew and Peggy Cherng Department of Medical Engineering<br />
Department of Electrical Engineering<br />
Division of Engineering and Applied Science<br />
California Institute of Technology<br />
1200 E. California Blvd., MC 138-78<br />
Pasadena, CA 91125<br />
Office: 205 Keck Labs (via 207)<br />
Work: 626-395-1959<br />
Fax: 626-395-1347<br />
Email: LVW@Caltech.edu<br />
Web: http://COILab.Caltech.edu<br />


Steven L. Jacques, Ph.D.<br />
Professor<br />
Oregon Medical Laser Center<br />
Providence/St. Vincent Hospital<br />
9205 SW Barnes Rd.<br />
Portland, OR 97225<br />
Tel:	503-216-4092<br />
Fax:	503-291-2422<br />
Email:	sjacques@ece.ogi.edu<br />
URL:	http://omlc.ogi.edu/staff/jacques.html<br />
