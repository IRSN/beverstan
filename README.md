The package
===========

The **beverstan** package for R is a private package performing Bayesian
inference for some Extreme-Value models. For now, these are simply
"time-varying" models with GEV margins for block maxima, usually annual
maxima. The function `TVGEVBayes` has an interface which is similar to
that of `NSGEV::TVGEV`. It further allows to use historical data via its
argument `timeMAXdata`, see the function help page `?TVGEVBayes` and the
examples therein.

The Bayesian inference relies on MCMC provided by Stan used through the
**rstan** package.

**beverstan** depends on the private packages **bever** and **NSGEV** as
well as on public packages available on the CRAN, see the `DESCRIPION`
file.

INSTALLATION
============

Cloning the repository
----------------------

If you do not have yet a local `beverstan` repository, use `git clone`
to clone the kergp repository

`git clone https://github.com/yvesdeville/beverstan/`

This will create a `beverstan` subdirectory of the current directory,
i.e. the directory from which the git command was issued. Installation
on Unix and MacOs systems

Install on Linux and MacOS systems
----------------------------------

With these sytems you can install a package from its source. Move to the
parent directory of your cloned repository and use the following command
from a terminal to create a tarball source file

`R CMD build beverstan`

This will produce a source tarball `beverstan_x.y.z` where `x`, `y` and
`z` stand for the major, minor and patch version numbers. Then you can
install from a command line

`R CMD INSTALL beverstan_x.y.z`

Note that you must formerly have installed all the packages required by
`beverstan` installed.

If you are using the **RStudio** IDE, you can alternatively use menus.
`Install` and `pre-compile for Windows`.

Install and pre-compile for Windows
-----------------------------------

In order to install the package from its source, you must have a
suitable Windows plateform with the **Rtools** installed. Then you can
proceed as Unix or MacOS users, with a build step from command line.

If you can not (or do not want to) install the **Rtools** you may get a
trusted binary from a friend or collegue next to you.

If you have the Rtools installed, you can create a binary. Using a
terminal, move if necessary by using cd to the directory containing the
source tarball and R command, and then type

`R CMD INSTALL --build beverstan_x.y.z`

This will create a `.zip` file that can be used on a Windows plateform
which may not be equipped with the **Rtools**. For instance, with
**RStudio** you can use the menu `Tools/Install Packages` and select
`Install from:`.
