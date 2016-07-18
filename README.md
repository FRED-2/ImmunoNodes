ImmunoNodes - Commandline tools and KINME extension of FRED 2
=============================================================

What is ImmunoNodes:
--------------------

ImmunoNodes is a collection of immunoinformatics command line tools and community KNIME nodes (http://www.knime.org), written in Python with FRED 2 (http://fred-2.github.io/).
It offers unified interfaces to many popular immunoinformatics methods and covers a wide variety of standard applications, 
such as: (neo-)epitope, cleavage site, and TAP prediction, as well as HLA genotyping, and epitope-based vaccine design.
Building onto of KNIME, ImmunoNodes is a powerful and intuitive toolbox to develop complex workflows, even without any programming experience. 




Licensing:
---------

The ImmunoNodes source code is published under a 3-clause BSD license. The licenses of the individual applications 
integrated into ImmunoNodes can be found in LICENSE.

* DTU(CBS) Prediction tools (NetMHC etc): Academic use only (please contact the administration of CBS for non-academic use http://www.cbs.dtu.dk/)
* LKH Solver: Academic use only (please contact the author of LKH for non-academic use http://www.akira.ruc.dk/~keld/research/LKH/)
* CBC Solver: Eclipse Public License v1.0
* Docker: Apache License v2.0
* Generic KNIME Nodes: Apache License v2.0


Prerequisite:
------------
ImmunNodes is platform independent by using Docker images (https://www.docker.com/), the CPU must support VT-x or AMD-v.

For KNIME integration you need:
* KNIME >= 3.1: http://www.knime.org
* Docker >= 1.9: https://www.docker.com/
* Generic KNIME Node (with Docker support): (https://github.com/genericworkflownodes/GenericKnimeNodes)

For commandline tool usage you can either use the provided Docker image and call the command lines within the Docker container (see our tutorial of how to do this). 
Or you have to install FRED 2 and all its dependencies (see http://fred-2.github.io/getting-started/)


Documentation Resources
-----------------------


Contact
=======
Benjamin Schubert
schubert@informatik.uni-tuebingen.de
University of Tübingen, Applied Bioinformatics,
Center for Bioinformatics, Quantitative Biology Center,
and Dept. of Computer Science,
Sand 14, 72076 Tübingen, Germany