.. _Taviclasses:

HPPT
#############################

The software is organized in a Model-View-Presenter pattern.
The main (Tavi) model, view and presenter classes and their components interactions are described here.


Tavi Model
+++++++++++++++

The TaviModel encapsulates the backend functionality of data calculations. The object fields are updated
every time there are new valid values received from the user (front end).

.. mermaid::

 classDiagram
    TaviProject "1" -->"1" TaviData
    TaviData "1" --> "1" Scan

    class TaviProject{
        -TaviData tavi_data
        +load_raw_scan_from_folder()
    }


    class TaviData{
        -RawDataPtr rawdataptr:dict[str, Scan]
        -CombineDataPtr combinedataptr:dict[str, np.array]
        -process_selected_data:list[str]
        -show_selected_data:list[str]
        -FitPtr fitptr
        -PlotPtr plotptr
    }

    class Scan{
        -Data raw_data
        -MetaData meta_data
        -ErrorMessage error_message
    }

Tavi View
+++++++++++++++


.. mermaid::

 classDiagram
    TaviView "1" -->"1" MainWindow

    class TaviView{
        -MainWindow:main_window
        -MainMenuBar:menu_bar
        +install_menu_bar()
        +closeEvent()
        +force_close()
        +exit_message_box()

    }

    class MainWindow{
        -ProjectView: load_view
        -Placeholder for other view widgets:Other
    }



Tavi Presenter
++++++++++++++++++++++

.. mermaid::

 classDiagram
    MainPresenter "1" -->"1" TaviProject
    MainPresenter "1" -->"1" TaviView
    TaviView "1" -->"2" MainWindow

    class MainPresenter{
        -TaviProject:model
        -TaviView:view
        +menu_bar()
        +install_menu_bar()
        +load_raw_scan_presenter(project_view, model_dict["TaviProjectProxy"])
        +exit()
    }

    class TaviProject{
        # from above
        -TaviData
        +load_raw_scan_from_folder()
    }

    class TaviView{
        # from above
        -MainWindow:main_window
        -MainMenuBar:menu_bar
        +install_menu_bar()
        +closeEvent()
        +force_close()
        +exit_message_box()

    }

    class MainWindow{
        # from above
        -ProjectView: load_view
        -Placeholder for other view widgets:Other
    }

The Hppt Model and View are unaware of one another. The Presenter is the connecting link that has a direct access and interacts with both.
The Presenter describes the main workflows that require communication and coordination between the Model and View through the Presenter. Additionally, the widgets' data initialization come from the model initialization and passed to the View.
Any value processing and/or filtering to match the requirements and logic of the View and Model side should happen on the Presenter.


#. Application Start - TaviView Initialization. All default settings (if there is any) are retrieved from the settings file.

    .. mermaid::

        sequenceDiagram
            participant View
            participant Presenter
            participant Model

            Note over View,Model:  TaviView Initialization
            Note left of View: User clicks "File"->"Load Experiment Folder"
            View->>View: Pop up a QDialog that allows user to view files and traverse system directories
            View->>Presenter: Return a single string of folder directory path
            Presenter->>Model: Pass experiment folder directory to load_raw_scan_from_folder(folder)
            Note right of Model: load_raw_scan_from_folder(folder) perform loading

            Model->>View: TO BE FILLED WITH BACKEND STORY

#. PLACEHOLDER MERMAID DIAGRAM

    * Valid Status:

        .. mermaid::

            sequenceDiagram
                participant View
                participant Presenter
                participant Model

                Note over View,Model: Updates due to switching to Single Crystal Mode
                Note left of View: User selects the Single Crystal radio button
                View->>Presenter: Trigger the update
                Presenter->>View: Update fields' visibility for single crystal case(field_visibility_in_SC)
                Note left of View: Display the SingleCrystalWidget
                Note left of View: Make crosshair_mod_q_edit field readonly
                Presenter->>Model: Set crosshair parameters with the experiment_type="SingleCrystal" (set_crosshair_data)
                Note right of Model: Store the crosshair data
                Presenter->>Model: Get crosshair parameters with the experiment_type="SingleCrystal"(get_crosshair_data)
                Note right of Model: Calculate the mod_q from the single crystal parameters and return it with the delta_e value
                Model->>Presenter: Return the crosshair data
                Presenter->>View: Return the crosshair data
                Note left of View: Display the data in the crosshair_widget
                Presenter->>Model: Get experiment parameters (get_experiment_data)
                Presenter->>View: Set experiment parameters (experiment_widget.set_values)
                Note left of View: Display experiment parameters values
                Note left of View: handle_field_values_update is triggered
                Presenter->>Model: Get single crystal parameters (get_single_crystal_data)
                Presenter->>Model: Get momentum transfer angle
                Model->>Presenter: Return momentum transfer angle
                Presenter->>View: Set single crystal parameters (singlecrystal_widget.set_parameters)
                Note left of View: Display SingleCrystal parameters values
                Note left of View: handle_field_values_update is triggered
                Presenter->>View: Display momentum transfer angle
                Note left of View: Update momentum transfer angle

    * Invalid Status:

        .. mermaid::

            sequenceDiagram
                participant View
                participant Presenter
                participant Model

                Note over View,Model: Updates due to switching to Single Crystal Mode
                Note left of View: User selects the Single Crystal radio button
                View->>Presenter: Triggers the update
                Presenter->>View: Update fields' visibility for single crystal case (field_visibility_in_SC)
                Note left of View: Display the SingleCrystalWidget
                Note left of View: Make crosshair_mod_q_edit field readonly
                Presenter->>Model: Set crosshair parameters with the experiment_type="SingleCrystal" (set_crosshair_data)
                Note right of Model: Store the crosshair data
                Presenter->>Model: Get crosshair parameters with the experiment_type="SingleCrystal"(get_crosshair_data)
                Note right of Model: Calculate the mod_q from the single crystal parameters and return it with the delta_e value
                Model->>Presenter: Return the crosshair data
                Presenter->>View: Return the crosshair data
                Presenter->>Model: Get experiment parameters (get_experiment_data)
                Presenter->>View: Set experiment parameters (experiment_widget.set_values)
                Note left of View: Display experiment parameters values
                Note left of View: handle_field_values_update is triggered
                Presenter->>Model: Get single crystal parameters (get_single_crystal_data)
                Presenter->>View: Set single crystal parameters (singlecrystal_widget.set_parameters)
                Note left of View: Display SingleCrystal parameters values
                Note left of View: handle_field_values_update is triggered


Experiment Settings
++++++++++++++++++++++
Placveholder for tavi settings, leaving HYSPECPPT's template here for reference. Update accordingly.

    * Experiment type options
        .. code-block:: bash

            class ExperimentType(Enum):
                POWDER = "Powder"
                SINGLECRYSTAL = "Single Crystal"
    * plot type options
        .. code-block:: bash

            class PlotType(Enum):
                ALPHA = "alpha_s"
                COSALPHA = "cos^2(alpha_s)"
                COSALPHAPLUS1 = "1+cos^2(alpha_s))/2"
                COS2ALPHA = (cos" + square + alpha + subscript_s + "-sin" + square + alpha + subscript_s + ")",
    * DEFAULT_MODE:dict =
        * experiment_type="single_crystal"
    * DEFAULT_CROSSHAIR: dict =
        * delta_e = 0
        * mod_q = 0
    * DEFAULT_EXPERIMENT:dict =
        * plot_type = PlotType.COS_2_ALPHA_S
        * incident_energy_e = 20
        * detector_tank_angle_s = 30
        * polarization_direction_angle_p = 0
    * DEFAULT_LATTICE:dict =
        * a = 1
        * b = 1
        * c = 1
        * alpha = 90
        * beta = 90
        * gamma = 90
        * h = 0
        * k = 0
        * l = 0
    * MAX_MODQ = 15 -- maximum momentum transfer
    * N_POINTS = 200 -- number of points in the plot
    * TANK_HALF_WIDTH = 30.0 -- tank half-width
