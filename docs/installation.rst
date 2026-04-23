.. highlight:: shell

============
Installation
============


Stable release
--------------

To install Cigar Math, run this command in your terminal:

.. code-block:: console

    $ pip install cigarmath

This is the preferred method to install Cigar Math, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


Optional dependencies
---------------------

To use SAM/BAM file I/O via ``cigarmath.io``, install with the ``io`` extra:

.. code-block:: console

    $ pip install cigarmath[io]

This installs `pysam`_, which requires a Unix-like system (Linux or macOS).

.. _pysam: https://pysam.readthedocs.io


Development install
-------------------

To work on Cigar Math locally (tests, linting, docs tooling), install the
``dev`` extra. This is the supported way to get a consistent developer
environment (CI, tox, and Read the Docs use the same dependency groups).

.. code-block:: console

    $ pip install -e ".[dev,io]"


From sources
------------

The sources for Cigar Math can be downloaded from the `Github repo`_.

.. code-block:: console

    $ git clone https://github.com/DamLabResources/cigarmath.git

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ pip install .

Or with optional dependencies:

.. code-block:: console

    $ pip install ".[dev,io]"

.. _Github repo: https://github.com/DamLabResources/cigarmath
