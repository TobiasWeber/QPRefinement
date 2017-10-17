# QPRefinement

This soplex version uses iterative refinement to solve quadratic optimization problems (QPs) to arbitrary precision.
Therefore more and more refined QPs are solved iteratively to increase the solution accuracy until a specified tolerance is reached.
The QPs are solved by the QP-solver qpOASES which is interfaced by the software.

## Installation

To use this sofware you need qpOASES. Follow the instructions on their download page:

[qpOASES download](https://projects.coin-or.org/qpOASES "download qpOASES here and come back after installation (compiliation)").

It is important that you check the `make_linux.mk` file in the main qpOASES directory and ensure that the `-D__USE_LONG_INTEGERS__` flag is not part of the `CPPFLAGS` variable.
Otherwise remove the flag and run make again in your main qpOASES directory to compile qpOASES.

Now you will need to download or clone this repo to your computer. Inside the repo run:

    make QPOASES=dense

You will be asked to enter the full path to your qpOASES install directory (e.g. /home/user/some/path/to/qpOASES/).
Afterwards you can run:

    make test SETTINGS=qpir_reliable TEST=marosmeszarossmall

to check the installation and solve some small QP examples. To test the programm further you will need more QP test problems in `.qps` file format. You can download Maros and Meszaros QP library for that:

[Maros and Meszaros QP library](http://www.cuter.rl.ac.uk/Problems/marmes.shtml "The Maros and Meszaros Convex Quadratic Programming Test Problem Set").

## Usage

To solve any QP problem given as `.qps` file you can run:

    bin/soplex -Q --loadset=settings/qpir_reliable.set /some/path/to/marosmeszaros/QPTEST.SIF

Here we solve the trivial `QPTEST.SIF` (this is qps format!) problem vom the Maros and Meszaros QP library. To indicate qprefinement we use the -Q flag and we load the default qprefinement setting file `qpir_reliable.set`.
