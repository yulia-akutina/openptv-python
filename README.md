openptv-python
==============

PyPTV is the Python version of [OpenPTV](http://www.openptv.net) which is 
Traits-based GUI and core algorithms in C (see <https://github.com/yosefm/openptv/tree/fb_submission>) 

## Getting started with OpenPTV-python on Mac OS X

This document is only for the Mac OS X users that can use pre-compiled ```.so``` files


1. Install EPD-free 32bit version from http://epd-free.enthought.com/ (includes Traits and wxPython)
2. Download the tarball from ...
3. Open the archive, get to the unzipped folder, open the terminal and run:

     python pyptv_gui/pyptv_gui.py test/
     
If you are not on Mac OS X and you need to compile the software, make sure you get the Cython and run:

    python pyptv_gui/setup.py build_ext --inplace
    
you can write us to openptv at gmail dot com and we'll provide pre-compiled version for your platform.

4. Follow the instructions in our screencast tutorials:
  *  Tutorial 1: http://youtu.be/S2fY5WFsFwo
  *  Tutorial 2: http://www.youtube.com/watch?v=_JxFxwVDSt0
  *  Tutorial 3: http://www.youtube.com/watch?v=z1eqFL5JIJc


## If you want to contribute, see http://www.openptv.net

