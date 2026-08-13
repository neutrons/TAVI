.. _cfsignatures:

Application Programming Interface
###################################

Generated function signatures and class documentation for the modules that make
up TAVI. The narrative descriptions live in :doc:`taviclasses` and the
:ref:`maintenance tutorials <maintainence_tutorials>`; this page is the reference
listing.

.. note::

   These directives import the modules they document. The Qt-backed frontend
   modules require a working ``qtpy``/PySide6 installation in the docs
   environment, which ``pixi run build-docs`` provides.

Backend
=======

Models
------

.. autoclass:: tavi.backend.model.tavi_project_model.TaviProjectModel
   :members:

.. autoclass:: tavi.backend.model.plot_model.PlotModel
   :members:

.. autoclass:: tavi.backend.model.application_model.ApplicationModel
   :members:

.. automodule:: tavi.backend.model.plot_resolver
   :members:

Model Interfaces
----------------

.. automodule:: tavi.backend.model.interface.model_interface
   :members:

.. automodule:: tavi.backend.model.interface.tavi_project_interface
   :members:

.. automodule:: tavi.backend.model.interface.plot_model_interface
   :members:

.. automodule:: tavi.backend.model.interface.application_model_interface
   :members:

Classification
--------------

.. autoclass:: tavi.backend.classification.raw_scan_classifier.RawScanClassifier
   :members:

.. autoclass:: tavi.backend.classification.rule_based_classifier.RuleBasedClassifier
   :members:

.. autoclass:: tavi.backend.classification.rule_set.rule_set.RuleSet
   :members:

.. autoclass:: tavi.backend.classification.rule_set.ornl_spice_rule_set.ORNLSpiceRuleSet
   :members:

.. autoclass:: tavi.backend.classification.rule.interface.rule_interface.RuleInterface
   :members:

Frontend
========

Presenters
----------

.. autoclass:: tavi.frontend.presenter.abstract_presenter.AbstractPresenter
   :members:

.. autoclass:: tavi.frontend.presenter.main_presenter.MainPresenter
   :members:

.. autoclass:: tavi.frontend.presenter.file_menu_presenter.FileMenuPresenter
   :members:

.. autoclass:: tavi.frontend.presenter.load_raw_scan_presenter.LoadRawScanPresenter
   :members:

.. autoclass:: tavi.frontend.presenter.plotter_presenter.PlotterPresenter
   :members:

.. autoclass:: tavi.frontend.presenter.data_file_presenter.DataFilePresenter
   :members:

.. autoclass:: tavi.frontend.presenter.error_presenter.ErrorPresenter
   :members:

Views
-----

.. autoclass:: tavi.frontend.view.main_view.TaviView
   :members:

.. autoclass:: tavi.frontend.view.project_view.ProjectView
   :members:

.. autoclass:: tavi.frontend.view.project_view.TreeViewWidget
   :members:

.. autoclass:: tavi.frontend.view.plotter_view.Plot1DView
   :members:

.. autoclass:: tavi.frontend.view.data_file_view.DataFileView
   :members:

.. autoclass:: tavi.frontend.view.filter_view.FilterView
   :members:

.. autoclass:: tavi.frontend.view.menubar_view.MainMenuBar
   :members:

.. autoclass:: tavi.frontend.view.file_menu_view.FileMenu
   :members:

.. autoclass:: tavi.frontend.view.error_view.ErrorView
   :members:

Widgets
-------

.. autoclass:: tavi.frontend.widget.tavi_message_box.TaviMessageBox
   :members:

Library
=======

Data Model
----------

.. automodule:: tavi.library.data.scan
   :members:

.. automodule:: tavi.library.data.plot
   :members:

.. automodule:: tavi.library.data.tavi_data
   :members:

.. automodule:: tavi.library.data.model_response
   :members:

.. automodule:: tavi.library.data.enum.raw_scan_type
   :members:

.. automodule:: tavi.library.data.enum.preset_type
   :members:

.. automodule:: tavi.library.data.enum.rebin_mode
   :members:

Storage and Loading
-------------------

.. autoclass:: tavi.library.storage.interface.file_store_interface.FileStoreInterface
   :members:

.. autoclass:: tavi.library.storage.local_file_store.LocalFileStore
   :members:

.. autoclass:: tavi.library.storage.controller.raw_scan_load_controller.RawScanLoadController
   :members:

.. autoclass:: tavi.library.storage.loader.loader_registry.LoaderRegistry
   :members:

.. autoclass:: tavi.library.storage.loader.interface.loader_interface.LoaderInterface
   :members:

.. autoclass:: tavi.library.storage.loader.interface.base.AbstractLoader
   :members:

.. autoclass:: tavi.library.storage.loader.ornl_spice_loader.ORNLSpiceLoader
   :members:
   :no-index:

.. autoclass:: tavi.library.storage.loader.default_loader.DefaultLoader
   :members:

Instrument and Components
-------------------------

.. automodule:: tavi.library.Instrument.instrument
   :members:

.. autoclass:: tavi.library.component.crystal.Crystal
   :members:

.. autoclass:: tavi.library.component.collimators.Collimators
   :members:

.. autoclass:: tavi.library.component.goniometer.Goniometer
   :members:

Geometry
--------

.. autoclass:: tavi.library.geometry.oriented_lattice.OrientedLattice
   :members:

.. autoclass:: tavi.library.geometry.sample.Sample
   :members:

Experiment
----------

.. autoclass:: tavi.library.experiment.experiment.Experiment
   :members:
   :no-index:

.. automodule:: tavi.library.experiment.enum
   :members:

.. automodule:: tavi.library.experiment.peak
   :members:

.. automodule:: tavi.library.experiment.utilities
   :members:

Resolution
----------

.. automodule:: tavi.library.resolution.resolution
   :members:
   :no-index:

.. autoclass:: tavi.library.resolution.ellipsoid.ResolutionEllipsoid
   :members:
   :no-index:

.. autoclass:: tavi.library.resolution.cooper_nathans.CooperNathans
   :members:

Fitting
-------

.. automodule:: tavi.library.fit.fit
   :members:
   :no-index:

Plotting
--------

.. autofunction:: tavi.library.plot.browser.browse_scans
   :no-index:

.. automodule:: tavi.library.plot.plot_ellipse
   :members:
   :no-index:

UB Algorithm
------------

.. autoclass:: tavi.library.ubalgorithm.ub.UBAlgorithm
   :members:
   :no-index:

TAS
---

.. automodule:: tavi.library.tas.triple_axis
   :members:
   :no-index:

Meta
====

Events
------

.. autoclass:: tavi.meta.event.event_broker.EventBroker
   :members:

.. automodule:: tavi.meta.event.event_interface
   :members:

.. automodule:: tavi.meta.event.type.model_event
   :members:

.. automodule:: tavi.meta.event.type.presenter_event
   :members:

.. automodule:: tavi.meta.event.type.exception_event
   :members:

Exceptions and Recovery
-----------------------

.. autoclass:: tavi.meta.exception.tavi_exception.TaviError
   :members:

.. autoclass:: tavi.meta.exception.recoverable.base.RecoverableError
   :members:

.. autoclass:: tavi.meta.exception.nonrecoverable.base.NonRecoverableError
   :members:

.. autoclass:: tavi.meta.exception.recovery_service.RecoveryService
   :members:

Threading
---------

.. automodule:: tavi.meta.multithreading.proxy
   :members:

.. automodule:: tavi.meta.multithreading.worker_pool
   :members:

.. automodule:: tavi.meta.multithreading.signal
   :members:

Logging and Configuration
-------------------------

.. automodule:: tavi.meta.logging
   :members:

.. automodule:: tavi.configuration
   :members:
