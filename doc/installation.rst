.. _installation:

============
Installation
============

Installation of WhatsHap is easiest if you use Conda.


Installation with Conda
-----------------------

If you do not have Conda, follow the `the Bioconda installation
instructions <https://bioconda.github.io/user/install.html#getting-started>`_.

If you already have Conda (for example, miniconda or Anaconda) installed,
ensure you have enabled both the “bioconda” and “conda-forge” channels::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge


Then WhatsHap can be installed with this command::

    conda install whatshap

For a smaller download and installation size, add ``nomkl`` to this command.
This will avoid installing the (quite large) “Intel Math Kernel Library”,
which WhatsHap does not use anyway::

    conda install whatshap nomkl


Installation with pip
---------------------

Before you can ``pip install``, you need to install dependencies that pip cannot
install for you. WhatsHap is implemented in C++ and Python. You need to have a
C++ compiler, Python 3.6 or later and the corresponding Python header files.
In Ubuntu, installing the packages ``build-essential`` and ``python3-dev`` will
take care of all required dependencies.

WhatsHap can then be installed with pip::

    pip3 install --user whatshap

This installs WhatsHap into ``$HOME/.local/bin``.  Then add
``$HOME/.local/bin`` to your ``$PATH`` and run the tool::

    export PATH=$HOME/.local/bin:$PATH
    whatshap --help

Alternatively, you can also install WhatsHap into a virtual environment if you
are familiar with that.


Installing an unreleased development version
--------------------------------------------

If you for want to use the most recent development version of
WhatsHap, you can install it in the following way. These instructions will
create a virtual environment in the directory ``whatshap-env`` that contains
WhatsHap. Simply delete that directory to uninstall the software. Other WhatsHap
versions you may have installed in other locations remain unaffected. ::

    python3 -m venv whatshap-env
    whatshap-env/bin/pip install --upgrade pip
    whatshap-env/bin/pip install git+https://github.com/whatshap/whatshap

You can then run WhatsHap like this::

    whatshap-env/bin/whatshap --version

You should see a version number like ``0.18.dev119+g5ba23de``, which means that
this is going to become version 0.18, with 119 commits ahead of the previous
version (0.17).
