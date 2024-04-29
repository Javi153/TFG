
/**
 * \mainpage CPSP - Constraint-based Protein Structure Prediction
 *
 * This document provides reference information about the 
 * <A HREF="http://www.bioinf.uni-freiburg.de/SW/CPSP/">Constraint-based
 * Protein Structure Prediction</A> Approach of Rolf Backofen and Sebastian Will.
 *
 * <CENTER>Version 2.4.6</CENTER>
 * 
 * \section SecGoals The goal
 * 
 * This C++ programming library provides fast and easy extensible implementation
 * of the constraint-based approach of Rolf Backofen and Sebastian Will to 
 * predict proven optimal structures of simplified 3D-lattice proteins in the
 * HP-model.
 *
 * Furthermore it will collect additional tools related to this.
 * 
 * \section SecFeatures CPSP Library Features
 * 
 * - A set of tools for direct 3D-lattice HP-model concerning tasks
 * - Full library implementation of the CPSP approach of R. Backofen and S. Will
 * - Interfacing classes to access H-core databases
 * - Complete API for all classes etc.
 * 
 * \section SecTools Included Tools in 2.4.6
 * 
 * - \ref SecHPstruct
 * - \ref SecHPstructSC
 * - \ref SecHPdeg
 * - \ref SecHPoptdeg
 * - \ref SecHPoptdegSC
 * - \ref SecHPdesign
 * - \ref SecHPnnet
 * - \ref SecHPrand
 * - \ref SecHPcompress
 * - \ref SecHPconvert
 * - \ref SecHPview
 * - \ref SecHPviewSC
 * 
 * \section SecDependencies Dependencies
 * 
 * The CPSP depends on the following libraries:
 * 
 * - the <A HREF="http://www.bioinf.uni-freiburg.de/SW/BIU/">Bioinformatic Utility</A>
 *   library (BIU) (version 1.3.0)
 * - the <A HREF="http://www.gecode.org">Gecode</A> Constraint Programming
 *   Library (version 1.3.0)
 * 
 * Additionally, the documentation also features the following parts:
 *  - \ref PageInstallation
 *
 * The following lists and indices are available
 *  - <a class="el" href="annotated.html">List of all classes including brief documentation</a>
 *  - <a class="el" href="namespaces.html">List of all namespaces including brief documentation</a>
 *  - <a class="el" href="files.html">List of all files</a>
 *  - <a class="el" href="hierarchy.html">Class hierarchy</a>
 *  - <a class="el" href="classes.html">Alphabetical class index</a>
 *  - <a class="el" href="namespacemembers.html">Namespace members</a>
 *  - <a class="el" href="functions.html">Class members</a>
 *  - <a class="el" href="globals.html">File members</a>
 *
 * Contact : http://www.bioinf.uni-freiburg.de/
 */



/**
 * \page PageInstallation Installation
 * 
 * Follow this steps to create a local version of the CPSP:
 * 
 * - download and unzip the <A HREF="http://www.bioinf.uni-freiburg.de/SW/CPSP/">CPSP source package</A> distribution
 * 
 * - run "./configure --help" to get an overview of further possible parameters
 * - run "./configure" with your desired parameters to generate the Makefiles
 * - "make" to generate the library
 * 
 * - "make test" to test if everything works fine
 * 
 * - if so run "make install" to install the headers, libraries and tools
 * 
 * - "make doc" to generate the html documentation (if selected during the
 *   configure call)
 * - "make doxygen-doc" to generate the whole documentation
 *
 * See the file INSTALL in the package root path for further informations.
 * 
 */
 
/**
 * \page PageTests Test cases and usage examples
 * 
 * To check if the CPSP does function as expected a test case was implemented
 * for each non-abstracted class. 
 * 
 * The source code can be used as an example for the usage of the CPSP-classes.
 * 
 * All files can be found in a subfolder named "tests" in the CPSP root
 * directory. They are named in the format "testCLASSNAME.cc".
 * 
 */
 
/**
 * \page PageTools Tools of the CPSP package
 * 
 * \section SecHPstruct HPstruct - 3D-lattice structure prediction
 * 
 * HPstruct is a tool to enumerate optimal 3D-lattice structures in the common
 * HP-model. The structures are given in absolute move string representation
 * to describe the structure relatively. Each internal neighboring vector  
 * between two monomers along the sequence (an absolute move) is represented by 
 * letters. We are using the following encoding:
 * 
 * - F/B = +- x
 * - L/R = +- y
 * - U/D = +- z
 * 
 * To calculate the structures a precalculated H-core database is neccessary and
 * used. Currently only a file-based database is available.
 * 
 * The approach uses the generic <A HREF="http://www.gecode.org">Gecode</A>
 * constraint programming framework to implement the CPSP method by R. Backofen
 * and S. Will.
 * 
 * To face the highly degenerated structure space of lattice proteins one can
 * exclude symmetric solutions directly. For that, the symmetry breaking approach
 * introduced by S. Will and R. Backofen is used and embedded into the Gecode
 * framework.
 * 
 * In addition, it is possible to constrain the resulting structures. 
 * A restriction of the equal absolute move string positions between all 
 * structures to a maximal value can be defined or it can be forced that all 
 * structures derived by one core differ in at least k monomer positions. Both
 * can be used to get a small set of different samples of the optimal structures 
 * if their number is high.   
 *   
 * For a complete list of the available parameters run the tool using '-help'.
 * 
 * 
 * \section SecHPstructSC HPstructSC - structure prediction in side chain models
 * 
 * HPstructSC extends CPSP approach to be utilized within HP side chain models 
 * in 3D lattices similar to HPstruct. Within side chain models each amino acid
 * is represented by two monomers, one for the backbone atoms and one 
 * representing the side chain atoms. 
 * The structure encoding follows the 
 * an extended absolute move string encoding:
 * 
 * e.g. (U)F(R)D(R)
 *
 * The absolute positioning of the side chain monomer according to the backbone
 * monomer is given in brackets. Inbetween the absolute positioning of the 
 * backbone monomers to each other is given as in the standard HP model.
 * 
 * 
 * \section SecHPdeg HPdeg - 3D-lattice degeneracy prediction
 * 
 * The degeneracy of a sequence S is the number of optimal structures that S can
 * adopt in a specific lattice. This can be calculated using the CPSP approach 
 * of R. Backofen and S. Will. 
 * 
 * The resulting tool HPdeg calculates the degeneracy of a given sequence in a
 * specified lattices using the CPSP library.
 * 
 * To handle high degenerated sequences as well and to allow testing for a 
 * maximal degeneracy this can be constrained to an upper bound.
 *   
 * For a complete list of the available parameters run the tool using '-help'.
 * 
 * 
 * 
 * \section SecHPoptdeg HPoptdeg - Search for low degenerated HP-sequences
 * 
 * The degeneracy of HP-sequences forms funnel-like structures in the sequence
 * space. Local search algorithms are therefore a possibility to find local 
 * minima.
 * 
 * HPoptdeg performs a Monte-Carlo search in the sequence space and finds low
 * degenerated HP-sequences.
 *   
 * For a complete list of the available parameters run the tool using '-help'.
 * 
 * 
 * \section SecHPoptdegSC HPoptdegSC - Search for low degenerated HP-sequences
 * in side chain models
 *
 * The extension of HPoptdeg to side chain HP models in 3D-lattices (see
 * \ref SecHPstructSC)
 * 
 * 
 * \section SecHPdesign HPdesign - HP-sequence design for a given structure
 * 
 * The problem HPdesign is facing is about the design of HP-sequences that fold
 * optimal into a given structure and have a degeneracy below a given upper 
 * bound.
 * 
 * The approach first uses a precalculated database of H-cores to detect 
 * sequences that can adopt the structure as an optimal one. Afterwards the 
 * degeneracy of the sequences is checked using the CPSP approach of R. 
 * Backofen and S. Will.
 * 
 * The level of suboptimal H-cores taken into account can be restricted to speed
 * up the search. If no sequence is found you should increase this level to 
 * take more sequences for tests into account.
 * 
 * Additionally, the H-content of the sequence can be constrained in order to
 * restrict the enumerated sequences. 
 *   
 * For a complete list of the available parameters run the tool using '-help'.
 * 
 * 
 * 
 * \section SecHPnnet HPnnet - Neutral nets of HP-sequences
 * 
 * A neutral net for a given sequence S and its only optimal structure X
 * includes all sequences S' that can adopt X as their only optimal structure 
 * too. Additionally, all sequences in S' have to be direct or indirect 
 * neighbors of S. Two sequences are neighbored if they differ only in one 
 * sequence position.
 * 
 * HPnnet uses for its calculation the CPSP approach of R. Backofen and S. Will
 * in order to check the degeneracy of a sequence neighbor and to compare its
 * optimal structure to X if degeneracy is 1. Per default symmetric structures
 * are excluded but can be included on demand.
 * 
 * To weaken the degeneracy criteria one can increase the maximal value allowed.
 * 
 * For a complete list of the available parameters run the tool using '-help'.
 * 
 * 
 * 
 * \section SecHPrand HPrand - Random HP-sequence generation
 * 
 * HPrand generates random HP-sequences of a given length.
 * 
 * The number and the H-monomer content of the output sequences can be 
 * constrained.
 *   
 * For a complete list of the available parameters run the tool using '-help'.
 * 
 * 
 *
 * \section SecHPcompress HPcompress - HP-sequence (de-)compression
 *
 * HPcompress allows the conversion of HP-sequences between normal/expanded
 * representation and a compressed one.
 *
 * e.g. HHHHPPPPPH <--> 4H5PH
 * 
 * For a complete list of the available parameters run the tool using '-help'.
 * 
 * 
 *
 *
 * \section SecHPconvert HPconvert - Lattice structure representation conversion
 *
 * HPconvert converts lattice structures between different formats.
 *
 * Currently supported representation formats are:
 * - Absolute move string
 * - Relative move string
 * - Absolute monomer positions given in XYZ-file format
 *
 * The move string representation follows the encoding:
 * - F/B = +- x
 * - L/R = +- y
 * - U/D = +- z
 *
 * The given structure is not validated (check if connected and selfavoiding).
 * For invalid structures a normal tool execution cant be guarantied.
 *
 * The XYZ-file format looks like that:
 *
 * <PRE>
 * # Beginning with '#' marks a comment line
 * # The lattice positions x,y and z of each point are given 
 * # in integer coding e.g.
 * 0 1 0
 * 1 1 0
 * 1 1 -1 
 * # EOF #</PRE>
 *
 * 
 * For a complete list of the available parameters run the tool using '-help'.
 * 
 * \section SecHPview HPview - HP lattice protein visualization
 *
 * HPview creates an output in CML-file format of a sequence/structure that can be
 * viewed with molecule viewers like Chime or Jmol.
 * 
 * The structure is NOT validated (if connected and selfavoiding). If it is invalid
 * normal execution cant be guarantied.
 * 
 * The moves are encoded using:
 *  F/B : +-x
 *  L/R : +-y
 *  U/D : +-z
 * 
 * Currently supported viewers are:
 *  - <A HREF="http://jmol.sourceforge.net"> Jmol </A> (http://jmol.sourceforge.net)
 * 
 * For a complete list of the available parameters run the tool using '-help'.
 * 
 * 
 * \section SecHPviewSC HPviewSC - HP side chain lattice protein visualization
 *
 * The extension of HPview to side chain HP models in 3D-lattices (see
 * \ref SecHPstructSC)
 * 
 * 
 * 
 */
 
