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

Installation
============

Cloning the repository
----------------------

If you do not have yet a local `beverstan` repository, use `git clone`
to clone the kergp repository

`git clone https://github.com/yvesdeville/beverstan/`

This will create a `beverstan` subdirectory of the current directory,
i.e. the directory from which the git command was issued.

If you already have cloned the repos, the use the following command from
within the `beverstan` directory

`git pull`

This will update your local copy.

Install on Linux and MacOS systems
----------------------------------

With these sytems you can install a package from its source. Move to the
parent directory of your cloned repository and use the following command
from a terminal to create a tarball source file

`R CMD build beverstan`

This will produce a source tarball `beverstan_x.y.z.tar.gz` where `x`,
`y` and `z` stand for the major, minor and patch version numbers. Then
you can install from a command line

`R CMD INSTALL beverstan_x.y.z.tar.gz`

Note that you must formerly have installed all the packages required by
`beverstan` installed.

Install and pre-compile for Windows
-----------------------------------

### Rtools

In order to install the package from its source, you must have a
suitable Windows plateform with the
[**Rtools**](https://cran.r-project.org/bin/windows/Rtools/) installed.
Then you can proceed as Unix or MacOS users, with a build step from
command line.

If you can not (or do not want to) install the **Rtools** you may get a
trusted binary from a friend or collegue next to you.

### Compiling with the Rtools

If you have the Rtools installed, you can create a binary from the
source tarball `beverstan_x.y.z.tar.gz`. Using a terminal, move if
necessary by using `cd` to the directory containing the source tarball,
and then type

`R CMD INSTALL --build beverstan_x.y.z.tar.gz`

This will create a `.zip` file that can be used on a Windows plateform
which may not be equipped with the **Rtools**. For instance, with
**RStudio** you can use the menu `Tools/Install Packages` and select
`Install from:`.

### RStudio

If you are using the **RStudio** IDE, you can alternatively use menus:
`Install` and `pre-compile for Windows` to install the source tarball
`beverstan_x.y.z.tar.gz`.

If you are you can also create a project using the Menu **New
Project...** / **Existing Directory** and then use the package
developpement tools **Build**.

All systems
-----------

If you are on a Linux, MacOS system or if you are on WIndow and have the
**Rtools** installed, you can also install the package **devtools**
package.

    devtools::install_github("yvesdeville/beverstan", auth_token = myToken)

where the `myToken` is a character object containing your authorisation
token. Of course you must keep your token secret. So do not write it in
R programs

Managing R versions (all systems)
---------------------------------

Remind that when the package is compiled with `R-a.b.c` the resulting
binary is intended to be used with R versions having the same major `a`
and the same minor `b` but a possibly different patch number `c`. So if
you want to use the package on several versions of R with different
major or minor, you will have to re-compile it to avoid warnings or
errors.
