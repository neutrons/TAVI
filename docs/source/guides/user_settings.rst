.. _user_settings:

++++++++++++++
User Settings
++++++++++++++

TAVI keeps two separate files under the user's application home directory
(``Config["user.application.home"]``, ``~/.TAVI/`` by default). They are managed
by different mechanisms and should not be confused.

``settings.yaml`` — application state
=====================================

Stores user-specific state that persists across sessions. Read and written
through ``LocalFileStore``'s user-data helpers.

Structure
---------

.. code-block:: yaml

    TAVI:
      recent:
        projects:
          - project1.h5
          - project2.h5

Key Components
--------------

**TAVI**
    Top-level namespace for all TAVI-specific settings.

**recent.projects**
    A list of recently opened project file paths. Populates the *Recent Projects*
    submenu of the File menu.

How it is loaded
----------------

``TaviProjectModel._get_recent_projects`` performs the read:

.. code-block:: python

    self.filestore.write_user_data_file("settings.yaml", Resource.read("default_settings.yml"))
    raw_settings_yml = self.filestore.read_user_data_file("settings.yaml")
    settings = YAML().load(raw_settings_yml)
    return dict(settings)["TAVI"]["recent"]["projects"]

.. warning::

   The first line **overwrites** ``~/.TAVI/settings.yaml`` with the packaged
   defaults from ``src/tavi/resources/default_settings.yml`` on every read. This
   is flagged in the source as demo scaffolding:

   .. code-block:: python

      # TODO: Demo purposes only. Remove this line when settings.yaml is actually used.

   So today the file is not really persistent state — the recent-projects list
   always comes back as the placeholder entries shipped in
   ``default_settings.yml`` (``untitled.h5``, ``adfasdf.h5``, ``1234.h5``), and
   nothing ever appends a genuinely opened project to it. Editing the file by
   hand has no lasting effect until that write is removed.

When it is read
---------------

Reads are event-driven rather than happening at import time:

#. ``MainPresenter`` finishes wiring the UI and publishes ``DownstreamReadyEvent``.
#. ``TaviProjectModel.sync_on_ready`` calls ``emit_sync_recent_projects()``.
#. That publishes ``SyncRecentProjects(recent_projects=[...])``.
#. ``FileMenuPresenter.sync_recent_projects`` forwards the list to
   ``FileMenu.init_recent_projects``, which adds one ``QAction`` per entry.

The indirection exists because the model is constructed before any presenter has
subscribed; ``DownstreamReadyEvent`` is the signal that consumers are ready for
startup state to be pushed.

``configuration.ini`` — configuration template
==============================================

A separate, unrelated mechanism handled by
:class:`tavi.configuration.Configuration`. On startup, ``tavi.__main__.execute``
constructs it and exits with status ``-1`` if ``is_valid()`` returns ``False``.

Structure
---------

.. code-block:: ini

    [global.other]
    help_url = https://github.com/neutrons/python_project_template/blob/main/README.md

Behaviour
---------

- Lives at ``~/.TAVI/configuration.ini`` (``CONFIG_PATH_FILE``).
- Created by copying ``src/tavi/configuration_template.ini`` if it does not exist,
  creating the parent directory as needed.
- On every startup, ``validate()`` walks the template and writes back any section
  or field the user's file is missing. Existing user values are preserved, so new
  settings appear automatically on upgrade without discarding local edits.
- ``get_data(section, name)`` reads a value, casting ``"True"``/``"False"`` to
  ``bool`` and ``"None"`` to ``None``.

.. note::

   The shipped template still points ``help_url`` at the
   ``python_project_template`` README rather than at TAVI's documentation.

Application defaults
====================

Neither file holds the runtime defaults used by the library (file size limits,
instrument names, classification flags). Those live in
``src/tavi/resources/application.yml`` and are read through
``neutrons_standard.config.Config`` with dotted keys, e.g.
``Config["library.filestore.raw.size-limit"]``. They are package resources, not
user-editable settings; tests override them with ``Config_override``
(see :doc:`../design/filestore`).
