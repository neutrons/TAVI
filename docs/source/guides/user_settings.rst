.. _user_settings:

++++++++++++++
User Settings
++++++++++++++

The ``settings.yaml`` file stores user-specific application preferences and state in the user's application home directory. This YAML configuration file allows TAVI to persist user settings across sessions.

Structure
=========

The ``settings.yaml`` file follows this structure:

.. code-block:: yaml

    TAVI:
      recent:
        projects:
          - project1.h5
          - project2.h5

Key Components
==============

**TAVI**
    Top-level namespace for all TAVI-specific settings.

**recent.projects**
    A list of recently opened project file paths. This list is used to populate the recent projects menu in the GUI, allowing users to quickly access previously opened files.

Usage
=====

The ``settings.yaml`` file is:

- Created automatically when TAVI starts if it doesn't exist
- Read by ``tavi_project_model.py`` to initialize the application state
- Managed by ``LocalFileStore`` for file I/O operations
- Used to pass the ``recent_projects`` list to the UI presenter

Users typically don't need to manually edit this file, as TAVI manages it automatically during operation.
