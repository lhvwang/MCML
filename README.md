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
Lihong Wang, Ph.D.
Associate Professor
Biomedical Engineering Program
234C Zachry Engineering Center
Texas A&M University
College Station, TX 77843-3120
Tel:	979-847-9040
Fax:	979-845-4450
Email:	LWang@tamu.edu
URL:	http://oilab.tamu.edu

Steven L. Jacques, Ph.D.
Professor
Oregon Medical Laser Center
Providence/St. Vincent Hospital
9205 SW Barnes Rd.
Portland, OR 97225
Tel:	503-216-4092
Fax:	503-291-2422
Email:	sjacques@ece.ogi.edu
URL:	http://omlc.ogi.edu/staff/jacques.html
