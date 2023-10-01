# LEGO plate generator

This is a Python script that generates arbitrary shaped LEGO plates.
The input is a SVG with a single path. It outputs a file for OpenSCAD.

Run as:

>python3 shape.svg shape.scad

This will create the file shape.scad .

When running openscad the file shape_parts needs to be in the search path
