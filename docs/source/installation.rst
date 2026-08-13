.. _installation:

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Installation for Triple-Axis Data Visualization (TAVI) development
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

* Create and activate a virtual environment with [Pixi](https://pixi.sh/). Prerequisites: Pixi installation e.g. for Linux:

   .. code-block:: bash

      $ curl -fsSL https://pixi.sh/install.sh | sh

* Setup/Update the environment of Tavi (make sure you are in Tavi project folder where pyproject.toml is)

   .. code-block:: bash

      $ pixi install

* Activate pixi environment

   .. code-block:: bash

      $ pixi shell

* Start the tool

   .. code-block:: bash

      $ tavi

  The same entry point is also available as a pixi task, which does not require an
  activated shell:

   .. code-block:: bash

      $ pixi run tavi

Useful development tasks
------------------------

All defined in the ``[tool.pixi.tasks]`` table of ``pyproject.toml``:

.. code-block:: bash

   $ pixi run test               # unit tests
   $ pixi run integration-test   # integration tests only
   $ pixi run build-docs         # build this documentation into docs/_build/html
   $ pixi run clean-all          # remove caches, build and doc artifacts
