.. _Taviclasses:

Core Classes
#############################

TAVI is organized in a Model-View-Presenter pattern. Models own state and never
touch Qt; views own widgets and never touch domain objects; presenters wire the
two together and subscribe to events. Cross-layer communication goes through the
``EventBroker`` singleton rather than direct references — see
:doc:`../guides/event_broker`.

Models
++++++

``TaviProjectModel`` owns the project's data (``TaviData``) and drives loading.
``PlotModel`` turns a focused scan into a ``Plot`` composition. ``ApplicationModel``
handles application-level concerns such as writing error logs.

.. mermaid::

 classDiagram
    TaviProjectModel "1" --> "1" TaviData
    TaviProjectModel "1" --> "1" RawScanLoadController
    TaviData "1" --> "*" RawScan
    TaviData "1" --> "*" Plot
    PlotModel ..> TaviData : live handles

    class TaviProjectModel{
        -Filestore filestore
        -TaviData tavi_data
        -EventBroker _event_broker
        -RawScanLoadController raw_scan_load_controller
        +get_plots_handle() dict
        +get_raw_scans_handle() dict
        +load_raw_scan_from_folder(folder) ModelResponse
        +emit_sync_recent_projects()
        -_handle_focus_event(FocusEvent)
        -_handle_save_plot_event(SavePlotEvent)
    }

    class TaviData{
        +dict~UUID,RawScan~ raw_scans
        +dict~UUID,Plot~ plots
        +fetch_by_uuid(uuid) RawScan|Plot
    }

    class PlotModel{
        -dict~UUID,Plot~ _plots
        -dict~UUID,RawScan~ _raw_scans
        -Plot _last_plot
        +update_fields(PlotFields) ModelResponse
        -_handle_raw_scan_focus_event(RawScanFocusEvent)
    }

    class ApplicationModel{
        -FileStoreInterface filestore
        +write_error_log(message) ModelResponse
    }

    class RawScan{
        +UUID uuid
        +ScanData data
        +ScanMetadata metadata
        +TaviMetadata tavimeta
        +Provenance prov
    }

    class Plot{
        +UUID uuid
        +list~PlotSeries~ series
    }

``TaviProjectModel``, ``ApplicationModel`` and ``RawScanLoadController`` are
singletons (``@Singleton`` from ``neutrons_standard.decorators.singleton``).
``PlotModel`` is not — it is constructed in ``tavi.__main__`` with live handles
into ``TaviProjectModel``'s storage, obtained via ``get_plots_handle()`` and
``get_raw_scans_handle()``.

.. note::

   ``PlotModel.__init__`` annotates its first parameter as ``list[Plot]``, but
   ``get_plots_handle()`` returns ``TaviData.plots``, which is a
   ``dict[UUID, Plot]``. The annotation is wrong, not the call — the attribute
   holds a dict. It is currently harmless because ``self._plots`` is stored and
   never read; ``PlotModel`` works off ``_last_plot`` and ``_raw_scans`` instead.

Every model method returns a ``ModelResponse`` (``code``, optional ``message``),
because models are invoked through a ``Proxy`` that runs them on a worker thread
and discards return values — see `Threading`_.

Views
+++++

``TaviView`` is the ``QMainWindow``. It does not create its sub-views; the
presenter layer owns them and passes them in through ``build_ui``.

.. mermaid::

 classDiagram
    TaviView "1" --> "1" MainMenuBar
    TaviView "1" --> "1" ProjectView
    TaviView "1" --> "1" Plot1DView
    TaviView "1" --> "1" DataFileView
    TaviView "1" --> "1" FilterView
    MainMenuBar "1" --> "1" FileMenu
    ProjectView "1" --> "1" TreeViewWidget

    class TaviView{
        +Signal exit_requested
        +install_menu_bar(menu_bar)
        +build_ui(project_view, plot_view, data_file_view, filter_view)
        +closeEvent(event)
        +force_close()
        +exit_message_box() bool
    }

    class ProjectView{
        -TreeViewWidget tree_widget
        +add_raw_scan(uuid, name, path)
        +add_plot(uuid, name, path)
        +get_selected_items() list~UUID~
        +hookup_select_signal(callback)
    }

    class TreeViewWidget{
        +dict~str,StandardItem~ path_map
        +dict~UUID,StandardItem~ uuid_map
        +Signal selected_signal
        +add_item_at_path(uuid, name, path)
        +remove_entry(index)
        +show_context_menu(position)
    }

    class Plot1DView{
        +Signal fields_focus_changed
        +Signal render_plots_signal
        +Signal plot_clicked
        +append_plot(x, y, err, ...)
        +get_plot_fields() PlotFields
        +sync_preset_fields(...)
        +reset_controls_to_defaults()
    }

    class DataFileView{
        +populate_columns(data)
        +populate_variables(names)
        +populate_metadata(metadata)
        +clear_data()
    }

``TaviView.closeEvent`` never closes the window itself. It emits
``exit_requested`` and ignores the event, leaving the decision to the presenter,
which calls ``force_close()`` once exit is approved.

Presenters
++++++++++

``MainPresenter`` constructs every other presenter, collects their views, and
hands them to ``TaviView.build_ui``. Each sub-presenter derives from
``AbstractPresenter``, which calls ``init_view()`` during construction and
exposes the result through ``view()``.

.. mermaid::

 classDiagram
    MainPresenter "1" --> "1" TaviView
    MainPresenter "1" --> "1" FileMenuPresenter
    MainPresenter "1" --> "1" LoadRawScanPresenter
    MainPresenter "1" --> "1" PlotterPresenter
    MainPresenter "1" --> "1" DataFilePresenter
    MainPresenter "1" --> "1" ErrorPresenter
    AbstractPresenter <|-- LoadRawScanPresenter
    AbstractPresenter <|-- PlotterPresenter
    AbstractPresenter <|-- DataFilePresenter
    AbstractPresenter <|-- ErrorPresenter

    class AbstractPresenter{
        -ViewInterface _view
        +init_view()
        +view() ViewInterface
    }

    class MainPresenter{
        -bool safe_exit
        -TaviView _view
        +exit() bool
    }

    class FileMenuPresenter{
        +handle_load_folder(folder)
        +sync_recent_projects(SyncRecentProjects)
        +exit()
    }

    class LoadRawScanPresenter{
        -dict~UUID,tuple~ inventory
        +update_treeview_data(RawScanAppendEvent)
        +update_plot_treeview_data(PlotAppendEvent)
        +handle_selection_event()
    }

    class PlotterPresenter{
        -Plot _current_plot
        +handle_plot_focus(PlotFocusEvent)
        +handle_raw_scan_focus(RawScanFocusEvent)
        +handle_fields_changed()
        +handle_plot_clicked()
    }

    class DataFilePresenter{
        +handle_raw_scan_focus(RawScanFocusEvent)
    }

    class ErrorPresenter{
        -RecoveryService recovery_service
        +handle_nonrecoverable_exception(NonRecoverableError)
    }

Presenters hold **no** domain state. Everything they need to render arrives
inside the event that triggered them; they never keep a live handle into a
model's storage and never call a model synchronously for data. The reasoning is
in :doc:`frontend/plot_data_model`.

Startup sequence
++++++++++++++++

``tavi.__main__.execute`` builds the models, wraps each in its ``Proxy``, and
hands the resulting dict to ``MainPresenter``.

.. mermaid::

    sequenceDiagram
        participant Main as __main__.execute
        participant Model as TaviProjectModel
        participant Presenter as MainPresenter
        participant Broker as EventBroker
        participant View as TaviView

        Main->>Main: Configuration() and validate
        Main->>Main: LocalFileStore()
        Main->>Model: TaviProjectModel(filestore)
        Main->>Main: PlotModel(plots_handle, raw_scans_handle)
        Main->>Main: ApplicationModel(filestore)
        Main->>Main: wrap each model in its Proxy
        Main->>Presenter: MainPresenter(dict_of_model)
        Presenter->>View: TaviView() and install_menu_bar()
        Presenter->>Presenter: construct sub-presenters (each creates its view)
        Presenter->>View: build_ui(project, plot, data_file, filter)
        Presenter->>Broker: publish DownstreamReadyEvent
        Broker->>Model: sync_on_ready
        Model->>Broker: publish SyncRecentProjects(recent_projects)
        Broker->>Presenter: FileMenuPresenter populates Recent Projects menu
        Main->>View: show()

``DownstreamReadyEvent`` exists because the models are constructed before the
presenters subscribe. Publishing it once the UI is fully wired lets the model
push its startup state (currently the recent-projects list) downstream without
the presenters having to poll.

Loading a folder
++++++++++++++++

.. mermaid::

    sequenceDiagram
        participant User
        participant FileMenu
        participant Presenter as FileMenuPresenter
        participant Model as TaviProjectModel
        participant Controller as RawScanLoadController
        participant Broker as EventBroker
        participant Tree as ProjectView

        User->>FileMenu: File → Load Experiment Folder
        FileMenu->>FileMenu: QFileDialog (directory mode)
        FileMenu->>Presenter: load_folder([path, ...])
        Presenter->>Model: load_raw_scan_from_folder(path)
        Note right of Model: runs on a worker thread via TaviProjectProxy
        Model->>Controller: load_folder(path)
        Controller-->>Model: list[RawScan]
        Model->>Model: store each scan in TaviData.raw_scans
        loop per scan
            Model->>Broker: publish RawScanAppendEvent
            Broker->>Tree: LoadRawScanPresenter.update_treeview_data
        end

The loading machinery itself — classification, loader selection, batching — is
documented in :doc:`raw_scan_load_controller` and :doc:`load/index`.

Threading
+++++++++

Models are never called directly from a presenter. ``tavi.__main__`` wraps each
model in a class generated by :func:`tavi.meta.multithreading.proxy.Proxy`, which
replaces every abstract method of the model's interface with a call that:

#. creates a ``Worker`` through the ``WorkerPool`` singleton (max 8 threads,
   excess work queued),
#. submits it to a thread, and
#. returns ``None`` immediately.

Models therefore cannot return values to a caller — a synchronous return would
block the GUI thread. They publish events instead, which is why every model
method's declared return type (``ModelResponse``) is consumed by the ``Worker``
rather than by the presenter.

``Worker.run`` is also the single exception-capture point for backend code: any
exception raised on the worker thread is converted into a ``NonRecoverableError``
carrying the captured stack trace and published as an ``ExceptionEvent``. See
:doc:`recovery_service`.
