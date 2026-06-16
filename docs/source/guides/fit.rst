.. _fit_package:

+++++++++++
Fit Package
+++++++++++

The fit package wraps a fitting backend (currently ``lmfit``) behind a small,
declarative API. Callers describe the models they want with ``ModelName`` enums
and a per-model options dict, and receive a ``FitResult`` that exposes fitted
components by prefix.


Overview
========

Source module:

.. code-block:: text

   tavi.library.fit.fit

Selecting a backend:

.. autoclass:: tavi.library.fit.fit.FitPackage
   :members:
   :undoc-members:

Available models:

.. autoclass:: tavi.library.fit.fit.ModelName
   :members:
   :undoc-members:

Primary class:

.. autoclass:: tavi.library.fit.fit.Fit
   :members:
   :undoc-members:


Result Objects
==============

A fit returns a ``FitResult`` aggregating one ``ComponentResult`` per model
component, keyed by its prefix.

.. autoclass:: tavi.library.fit.fit.FitResult
   :members:

.. autoclass:: tavi.library.fit.fit.ComponentResult
   :members:


Model Specification
===================

Models are passed as a list of ``(ModelName, options)`` tuples. Each entry adds
one component to the composite model; ``options`` is forwarded to the backend
(for example ``dict(guess=True)`` to auto-initialize parameters).

.. code-block:: python

   from tavi.library.fit.fit import FitPackage, ModelName

   model_dict = [
       (ModelName.Gaussian, dict(guess=True)),
       (ModelName.Linear, dict()),
   ]


Minimal Example
===============

.. code-block:: python

   import numpy as np
   from tavi.library.fit.fit import Fit, FitPackage, ModelName

   x = np.linspace(-5, 5, 100)
   y = np.exp(-(x**2)) + 0.01 * x

   fit = Fit(package=FitPackage.lmfit)
   result = fit.fit(x, y, [(ModelName.Gaussian, dict(guess=True))])

   # Single-component fits expose scalar parameters directly.
   print(result.center, result.center_err)
   print(result.redchi)

   # With multiple components, read each by its prefix instead:
   #   result["g1_"].values["center"]
