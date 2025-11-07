.. _cfsignatures:

Application Programming Interface
###################################

It includes the function signatures and classes of the various modules.

Hyspecppt
----------

.. automodule:: hyspecppt
   :members:

Hppt
------
It contains the majority of the software's functionality organized in a Model-View-Presenter architectural pattern.
Tenchnically hppt is included as a widget in the base layout of the MainWindow.

------
Model
------

.. autoclass:: hyspecppt.hppt.hppt_model.HyspecPPTModel
   :members:

.. autoclass:: hyspecppt.hppt.hppt_model.CrosshairParameters
   :members:

.. autoclass:: hyspecppt.hppt.hppt_model.SingleCrystalParameters
   :members:

----------
Presenter
----------

.. automodule:: hyspecppt.hppt.hppt_presenter
   :members:

-----
View
-----

.. autoclass:: hyspecppt.hppt.hppt_view.HyspecPPTView
   :members:

.. autoclass:: hyspecppt.hppt.hppt_view.ExperimentWidget
   :members:

.. autoclass:: hyspecppt.hppt.hppt_view.SelectorWidget
   :members:

.. autoclass:: hyspecppt.hppt.hppt_view.CrosshairWidget
   :members:

.. autoclass:: hyspecppt.hppt.hppt_view.SingleCrystalWidget
   :members:

.. autoclass:: hyspecppt.hppt.hppt_view.PlotWidget
   :members:


View Field Validators
----------------------

.. automodule:: hyspecppt.hppt.hppt_view_validators
   :members:

Help
-----

The help functionality, includes a button that directs the user to the documentation.
The help button is added in the MainWindow.

.. automodule:: hyspecppt.help.help_model
   :members:


Configuration
--------------

The configuration functionality includes settings that are retrieved though a configuration template.
The mechanism is created in a way to be updated automatically, when a new setting is added,
while allowing the users to update the default ones in their environment.

.. automodule:: hyspecppt.configuration
   :members:
