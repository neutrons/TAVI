.. _integration_testing:

+++++++++++++++++++++++++++++++++++
Integration Testing Guide
+++++++++++++++++++++++++++++++++++

This guide outlines the methodology for writing and maintaining integration tests in TAVI. Integration tests verify how multiple components work together as a system, complementing unit tests that validate individual components.

Overview
========

Integration tests in TAVI are located in ``tests/integration/`` and serve to:

- Verify complete workflows end-to-end
- Test interaction between multiple components (backend, frontend, storage, etc.)
- Ensure configuration and initialization work correctly in a realistic environment
- Catch issues that unit tests cannot detect

Integration tests carry the ``integration`` pytest marker, which they inherit by
subclassing ``IntegrationTest`` — an explicit ``@pytest.mark.integration``
decorator is not needed. ``pyproject.toml`` sets ``addopts = "-m 'not
(integration)'"``, so a plain ``pytest`` run skips them. Run them with:

.. code-block:: bash

    pixi run integration-test   # pytest -m integration

Running All Tests (Including Integration)
==========================================

To run both unit and integration tests:

.. code-block:: bash

    pytest -m ""

Test Structure
==============

Integration tests inherit from ``IntegrationTest`` base class provided by the ``neutrons_standard`` library and use pytest with Qt support via ``pytest-qt``.

Basic Test Template
-------------------

.. code-block:: python

    import pytest
    from neutrons_standard.test.integration.test_base import IntegrationTest
    from neutrons_standard.test.integration.test_summary import TestSummary

    class TestMyWorkflow(IntegrationTest):
        """Subclassing IntegrationTest applies the `integration` marker."""

        @pytest.fixture(scope="function", autouse=True)  # noqa: PT003
        def setup(self):
            """Initialize state for each test function."""
            pass

        def test_workflow_name(self, qtbot, qapp, monkeypatch):
            """Test a complete workflow with multiple steps."""
            self.test_summary = (
                TestSummary.builder()
                .step("Initialize components")
                .step("Perform action")
                .step("Verify result")
                .build()
            )

            # Test implementation here

            self.test_summary.SUCCESS()

Key Principles
==============

1. **Reusable Subtest Methods**

   Break down test workflows into reusable subtest methods that can be composed in different ways. Each subtest method encapsulates a specific action or verification step and can execute independently given the appropriate component state.

   **Good Example:**

   .. code-block:: python

       class TestLoadRawScans(IntegrationTest):

           def _subtest_load_folder(self, state, qtbot):
               """Subtest: Load a folder via file dialog."""
               def auto_accept():
                   for w in QApplication.topLevelWidgets():
                       if isinstance(w, QFileDialog):
                           w.selectFile(Resource.getPath("/inputs/integration/load/datafiles"))
                           w.accept()
                           break

               QTimer.singleShot(0, auto_accept)
               qtbot.mouseClick(state.load_button)
               qtbot.waitUntil(lambda: len(state.project_model.tavi_data.raw_scans) > 0)

           def _subtest_verify_scans_visible(self, state):
               """Subtest: Verify loaded scans appear in the project view."""
               assert len(state.project_model.tavi_data.raw_scans) > 0
               assert len(state.tree_widget.uuid_map) > 0

           def test_load_and_display_scans(self, qtbot, qapp):
               """Main test: Load scans and verify they're displayed."""
               state = ApplicationState(qtbot, qapp)
               state.gui.show()

               self._subtest_load_folder(state, qtbot)
               self._subtest_verify_scans_visible(state)

           def test_load_multiple_folders(self, qtbot, qapp):
               """Main test: Load multiple folders sequentially."""
               state = ApplicationState(qtbot, qapp)
               state.gui.show()

               self._subtest_load_folder(state, qtbot)
               self._subtest_load_folder(state, qtbot)  # Reuse the same subtest
               assert len(state.project_model.tavi_data.raw_scans) > 1

   **Benefits:**

   - Subtests are composable and reusable across multiple test methods
   - Main test methods orchestrate which subtests run based on the scenario
   - State is explicitly managed by the main test method
   - Tests become more modular and easier to maintain
   - Reduces code duplication for common operations

2. **Test Isolation**

   Each integration test should be independent:

   - Use fixtures with appropriate scope to initialize state fresh for each test
   - Clean up GUI widgets and temporary resources after each test
   - Monkeypatch external interactions (file dialogs, message boxes) to avoid GUI popups

   .. code-block:: python

       def test_workflow(self, qtbot, qapp, monkeypatch):
           def fail_on_exec(*args, **kwargs):
               raise AssertionError("Unexpected message box exec() called")

           monkeypatch.setattr(TaviMessageBox, "exec", fail_on_exec)

           # Test implementation

           file_menu_view.close()  # Always clean up

3. **Simulate User Interactions**

   Use pytest-qt features to simulate realistic user behavior:

   - ``qtbot.mouseClick()`` for button clicks
   - ``qtbot.waitUntil()`` for waiting on asynchronous events
   - ``QTimer.singleShot()`` to schedule operations
   - File dialog auto-acceptance via monkeypatching

   .. code-block:: python

       def auto_accept_file_dialog():
           for w in QApplication.topLevelWidgets():
               if isinstance(w, QFileDialog):
                   w.selectFile(Resource.getPath("/inputs/integration/load/datafiles"))
                   w.accept()
                   break

       QTimer.singleShot(0, auto_accept_file_dialog)
       qtbot.mouseClick(button)

4. **Use TestSummary for Documentation**

   Document test steps using ``TestSummary`` to make test intent clear:

   .. code-block:: python

       self.test_summary = (
           TestSummary.builder()
           .step("Open Folder Browser")
           .step("Load Raw Scans into ProjectView")
           .step("Verify data loaded correctly")
           .build()
       )

5. **Wait for Asynchronous Operations**

   Integration tests often need to wait for backend processing:

   .. code-block:: python

       # Wait for a condition to be true (with timeout)
       qtbot.waitUntil(lambda: file_menu_view.isVisible())

       # Wait for a fixed duration (use sparingly)
       qtbot.wait(50)  # 50 milliseconds

Test Data Management
====================

Two separate trees hold test data, and they are not interchangeable:

``tests/resources/`` — **what integration tests actually read.**
    Resolved through ``neutrons_standard.config.Resource.getPath()``, whose paths
    are relative to this directory. The load test's input lives at
    ``tests/resources/inputs/integration/load/`` with ``datafiles/`` and
    ``UBConf/`` subdirectories, referenced as
    ``Resource.getPath("/inputs/integration/load/datafiles")``. This tree also
    carries the test-environment ``application.yml``, ``default_settings.yml`` and
    ``logging_config.json``, which override the packaged ones because
    ``tests/conftest.py`` sets ``env=test`` before ``neutrons_standard.init``.

``test_data/`` — **a large corpus of real experiment data**, used ad hoc by
notebooks, scripts, and some unit tests. It holds full SPICE experiment folders
(``IPTS32912_HB1A_exp1031/``, ``exp815/``, ...) and HDF5 files such as
``tavi_test_exp1031.h5``. It is not on the ``Resource`` search path.

When adding new integration tests that require test data:

1. Store minimal, representative fixtures under ``tests/resources/inputs/``
2. Document the data format and contents
3. Reference them with ``Resource.getPath()`` — never a hardcoded path
4. Keep file sizes small; ``library.filestore.raw.size-limit`` is 1 MB, so any
   larger file is silently skipped by ``fetch_files_at`` and the test will look
   like it loaded nothing

Git-LFS Consideration
---------------------

As test data grows, consider using Git LFS (Git Large File Storage):

- Beneficial when: binary test files exceed 50-100 MB collectively
- Advantages: keeps Git repository lightweight, maintains full history
- Implementation: add patterns to ``.gitattributes`` for large files
- Example: ``*.h5 filter=lfs diff=lfs merge=lfs -text``

Future plans may include:

- Migration of large HDF5 test files to LFS
- Automated test data generation from smaller base files
- Docker-based test environment with pre-built test data

Common Patterns
===============

Loading Application Components
-------------------------------

.. code-block:: python

    def test_application_startup(self, qtbot, qapp):
        # Mirror tavi.__main__.execute: filestore first, then the models.
        filestore = LocalFileStore()
        tavi_project_model = TaviProjectModel(filestore)
        plot_model = PlotModel(
            tavi_project_model.get_plots_handle(),
            tavi_project_model.get_raw_scans_handle(),
        )
        application_model = ApplicationModel(filestore)

        # Wire up models through a presenter
        dict_of_model = {
            "TaviProjectProxy": TaviProjectProxy(tavi_project_model),
            ApplicationModelInterface.__name__: ApplicationModelProxy(application_model),
            PlotModelInterface.__name__: PlotModelProxy(plot_model),
        }

        main_presenter = MainPresenter(dict_of_model)
        main_presenter.safe_exit = False
        gui = main_presenter._view
        gui.show()
        qtbot.addWidget(gui)

        # Now test interactions with the initialized GUI

``MainPresenter`` requires all three keys — it indexes ``model_dict`` for
``"TaviProjectProxy"``, ``ApplicationModelInterface.__name__`` and
``PlotModelInterface.__name__`` — so omitting one raises ``KeyError`` during
construction. ``PlotModel`` must be built from the project model's live handles,
which is why the ordering matters.

Testing File Operations
-----------------------

.. code-block:: python

    def test_file_loading(self, qtbot, qapp, monkeypatch):
        # Set up the state
        state = ApplicationState(qtbot, qapp)

        # Mock file dialog to automatically accept and select a file
        def auto_accept():
            for w in QApplication.topLevelWidgets():
                if isinstance(w, QFileDialog):
                    w.selectFile(Resource.getPath("/inputs/integration/load/datafiles"))
                    w.accept()
                    break

        QTimer.singleShot(0, auto_accept)

        # Trigger the file open action by clicking its rect in the popped-up menu
        file_menu_view = state.presenter.file_menu_presenter._view
        action = file_menu_view.load_folder_action
        file_menu_view.popup(state.gui.mapToGlobal(QPoint(10, 10)))
        qtbot.waitUntil(lambda: file_menu_view.isVisible())
        qtbot.wait(50)  # allow layout to complete before reading geometry

        rect = file_menu_view.actionGeometry(action)
        assert rect.isValid()
        qtbot.mouseClick(file_menu_view, Qt.LeftButton, pos=rect.center())

        # Verify against the tree view's uuid_map. The uuid is the md5 of the
        # file's text content, so it is stable for a fixed fixture file.
        tree_widget = state.presenter.project_view.tree_widget
        qtbot.wait(50)
        assert UUID(value="3d682ef8c6633c0dd9bdad1f35439a7c") in tree_widget.uuid_map

        file_menu_view.close()  # always clean up

Best Practices
==============

✓ **Do:**

- Create reusable state-based components for common setup
- Document test steps using TestSummary
- Keep individual tests focused on a single workflow
- Use monkeypatching to control external dependencies
- Wait for conditions rather than using fixed delays
- Clean up GUI resources after each test
- Add comments explaining non-obvious GUI interactions

✗ **Don't:**

- Create monolithic tests that test multiple unrelated workflows
- Use fixed delays (``time.sleep()``) for waiting on conditions
- Leave GUI windows open after tests complete
- Import test utilities directly from integration modules (use fixtures)
- Write tests that depend on the order of execution
- Hardcode file paths; use ``Resource.getPath()`` instead

Debugging Integration Tests
============================

Running with Verbose Output
----------------------------

.. code-block:: bash

    pytest tests/integration/test_load_raw_scans.py -v -s

Running a Single Test
---------------------

.. code-block:: bash

    pytest tests/integration/test_load_raw_scans.py::TestLoadRawScans::test_load_ORNL_Spice -v -s

Capturing Screenshots
---------------------

For debugging GUI issues, pytest-qt can capture screenshots:

.. code-block:: python

    def test_with_screenshot(self, qtbot, qapp):
        # ... test code ...
        path = qtbot.screenshot(gui)  # takes the widget; returns the file path
        print(path)                   # written under pytest's tmp dir

Troubleshooting
---------------

**Test hangs waiting for dialogs:**
- Ensure file dialogs are being mocked with ``monkeypatch.setattr()``
- Check that ``QTimer.singleShot(0, ...)`` is scheduled before triggering the action
- Use ``qtbot.waitUntil()`` with timeout instead of infinite waits

**GUI elements not found:**
- Add ``qtbot.wait(50)`` to allow layout to complete after showing widgets
- Use ``qtbot.waitUntil()`` to wait for elements to become visible
- Check that widgets are being properly initialized in the fixture

**Singleton state pollution:**

- ``tests/conftest.py`` calls ``reset_Singletons()`` before each test, class and
  module — but only for tests that do *not* carry the ``integration`` keyword.
  Integration tests deliberately keep singleton state across the run, since they
  exercise the real application wiring.
- That means ``LoaderRegistry``, ``RawScanLoadController``, ``EventBroker``,
  ``RecoveryService``, ``TaviProjectModel`` and ``WorkerPool`` persist between
  integration tests. Constructor arguments passed by a later test are ignored —
  the first construction wins.
- Ensure proper cleanup in test fixtures if singletons are used
- Use separate test classes for tests with conflicting singleton state

Integration Test Examples
=========================

See ``tests/integration/test_load_raw_scans.py`` for a complete example of:

- Setting up application state with multiple models
- Mocking file dialogs and message boxes
- Simulating user interactions with the GUI
- Using TestSummary to document test steps
- Verifying application state changes

Contributing Integration Tests
==============================

When adding a new integration test:

1. Create a new test file in ``tests/integration/`` named ``test_<feature_name>.py``
2. Inherit from ``IntegrationTest`` base class
3. It will inherit the ``@pytest.mark.integration`` decorator this way.
4. Extract reusable state setup into separate classes
5. Document test steps with TestSummary
6. Add comments for non-obvious GUI manipulations
7. Ensure tests clean up after themselves
8. Run tests locally: ``pixi run integration-test``
9. Verify unit tests still pass: ``pixi run test``

Future Enhancements
===================

**Planned Improvements:**

- **Git-LFS Integration:** Large binary test files will be managed with Git LFS for faster clones and reduced repository size
- **Test Data Factories:** Utilities to programmatically generate test data instead of storing large files
- **Snapshot Testing:** Visual regression testing for GUI components
- **Performance Baselines:** Track integration test execution time to catch performance regressions
- **CI/CD Optimization:** Parallel execution of independent integration tests
- **Test Video Recording:** Automated recording of test execution for debugging failed tests

See Also
========

- :ref:`maintainence_tutorials`
- `pytest-qt Documentation <https://pytest-qt.readthedocs.io/>`_
- `pytest Documentation <https://docs.pytest.org/>`_
- ``tests/integration/`` - Integration test examples
- ``tests/conftest.py`` - Test fixtures and configuration
