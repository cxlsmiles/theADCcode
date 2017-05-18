=============================
Condensed Instructions
=============================

in .../rand/rand_750x750-asymm
	compile
in .../matvec
	compile
in .../Lanczos
        Lanczos -> v#.#.#
in .../Lanczos/v#.#.#
	lapack_
	compile
in .../Tests/test_util
        -llapack
in .../Tests/test_util/lib
        compile
in .../Tests/test_util/inspect
        compile
in .../Tests/test#
        compile
        test [write]

=============================
=============================



=============================
To-Do List
=============================

REMINDER: the M*v operation should have += character on the target vector(s)!

1) Technical changes
X	a) Replace Resizable with an STL container
X	b) Get rid of static variables (interfere between different instances of the class)
	c) Merge in adcman changes (see ADCman_merge/readme.txt)
X	d) double==double (in lanczos.h)
	e) //prints in lanczos.h
	f) use error throwing instead of error messages
	g) how to unify output mechanism between streams (overload function protocol to take either?)

2) Design changes
	a) "wrap" the demands that I make on the matrix and vector classes using traits that associate the proper
	   vector and scalar types with a given matrix type
	   it is at that level that functions for taking the real part and conjugating must be defined (probably via other traits of those types)
	b) move diag-converged one level up and make time-step and other interfaces (engine makes the matrix, and Lanczos manages it)
	   Move DiagConverged up one level.  To this end, GetEigenvals should automatically diagonalize the matrix, unless it hasn't increased
	   since last time (no need to call diagonalized directly) and residuals can be computed by * operator directly between projection
	   and vector?
	c) similarly allow the user to insert their own band diagonalizer

3) New Features
	a) Make it so that Reset(NULL,0) puts us in the fresh from the constructor state (deallocate engine, since it can't handle 0 block size)
	   why did i want this?
	b) projective lanczos for eigenstate finding (projects out converged states)
	   Make methods to streamline the following process:  As one of the eigenvectors converges, it introduces noise the procedure.  I would
	   like to restart, using the current vectors as the guess vectors, except the one that is nearly converged.  This involves restarting
	   with a smaller block size, for which facilities largely already exist.  Now the major question is what to do with the converged vector.
	   We want to keep the entire procedure for the other vector strictly orthogonal to it (should this be done in the Lanczos, Lanczos Engine
	   or LancMat level).

4) Things that will probably never get done
	a) implement the object-manager mechanism to talk to Michael and Evgeny's context
	b) Rewrite insides on top of spaces and expressions libraries

5) Documentation changes
	a) in 2. How to use the 'Lanczos' one can mention that the class can be used to approximate the spectrum of the matrix rather than obtaining
	   the exact eigenvalues and eigenvectors.
	b) for LancVec add:
		- LinCombo must now act properly even if the two input arguments are *this, for all values of the coefficients
		- an operator=(double) method needs to work at least for zero 
	c) for LancMat add:
	   the method 'void operator()(LancVec&, const LancVec&)' now changes to 'int operator()(LancVec**, LancVec**, int)'
	   the last input argument tells the matrix class how many vectors to act on in the second argument and *add to*
	   the locations in the first argument.  The matrix class may do fewer and return the number of vectors acted
	   on, but may not do more.

=============================
=============================



=============================
History
=============================

The patch slot "0.0.1" etc is reserved for bug fixes, no new features.

===========
version 0.0
===========

Version 0.0 is the first largely stable version whose results are verifiably compatible with the tested 
beta code from the q-chem application.  This code is based on the modified version from Yasen's ADC-2-hole
code, and has had substantial changes to the front and back end.

All future versions will derive from this code.  As far as can be told, this code should be functional
in both applications that use it so far (q-chem and Yasen's), with the minor exception that the q-chem
application needs to have the "M*v" operation changed so that it has "+=" character on the target vectors.

Since the "sandox" (Yasen's) code on which this was based was largely free of bugs, once minor modifications
were made so that all routines compiled for both real and complex numbers, this version is the same as
the "sandbox" code, but it is named, so that it is clear that it is the stable starting point, and all
future bug fixes will proceed from this point, never touching/running the sandbox code again.

From now on, this code will be used to generate the reference data, even as potential as-of-yet
undiscovered bugs are fixed here.  Of course, we will continue to check the q-chem (ADC) beta code against 
this data for consistency, possibly fixing bugs there, for some time to come.

!!! For all intents and purposes this revision is now static, save bug fixes, and the Sandbox code is dead. !!!

===========
version 0.1
===========

Version 0.1 is now the working code where new features (well, new to this branch, some are already in 
Yasen's code and mine) will be merged in and issues of coding style will be cleaned up.

All future versions will derive from this code.  As far as can be told, this code should be functional
in both applications that use it so far (q-chem and Yasen's), with the minor exception that the q-chem
application needs to have the "M*v" operation changed so that it has "+=" character on the target vectors.

Yasen's changes to the front-end lanczos.h file (found in toMerge-from-Yasen/myLanc) have been merged in
have been merged now ... we will assume that we got everything that is important, and skip the suggested
step of scouring emails for suggestions and bug reports.

In spite of still not having merged some stuff from the ADC code, it is time to close this revision
so that Yasen can play with some other improvements now.  The dangling issues will be transplanted to
the v0.2 directory.

===========
version 0.2
===========

Version 0.2 is now the working code.  As of writing this sentence, it is the same as v0.1, except for
the structure of this documentation.  Yasen plans to clean up a few things, mostly under point 1) on the
above To-Do list.  Either in this revision or in 0.3, the ADCman changes will be merged in.

===========
version 1.0
===========

Version 1.0 will start the major overhaul of the backend.


=============================
=============================
