===============================
Running a geometry optimization
===============================
In the following we will assume to have a valid user input file for the water
molecule called ``h2o.inp``, e.g. like this

.. include:: h2o_geopt.inp
    :literal:

A geometry optimization can be run by adding GeometryOptimizer section to any normal .inp file and setting the run keyword to true:

.. include:: geopt_section.inp
    :literal:
