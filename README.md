# stable-stats
Software system to convert phase data to web page with frequency stability information.

John Ackermann   N8UR
Copyright 2009 - 2016
Licensed under GPLV3 -- All other rights reserved

stable-stats consists of a perl script supported by a couple of bash scripts and "template" files.  It relies on the "grace" plotting program (http://plasma-gate.weizmann.ac.il/Grace/) which ought to be more popular than it is.

It expects to read data files with two fields -- the first an MJD datestamp, and the second a phase value.  These files are created by external programs.

TO DO:

1.  Move html and text into external text files and do variable substitution; this will simplify the .pl code dramatically
2.  Improve statistics code; add additional ADEV types and validate performance
3.  Better documentation
