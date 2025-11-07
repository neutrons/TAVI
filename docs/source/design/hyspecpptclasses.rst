.. _hyspecpptclasses:

HPPT
#############################

The software is organized in a Model-View-Presenter pattern.
The main (HyspecPPT) model, view and presenter classes and their components interactions are described here.


HyspecPPT Model
+++++++++++++++

The HyspecPPTModel encapsulates the backend functionality of data calculations. The object fields are updated
every time there are new valid values received from the user (front end).

.. mermaid::

 classDiagram
    HyspecPPTModel "1" -->"1" CrosshairParameters
    CrosshairParameters "1" -->"1" SingleCrystalParameters

    class HyspecPPTModel{
        +float incident_energy_ei
        +float detector_tank_angle_s2
        +float polarization_direction_angle_p
        +enum 'PlotType' plot_type
        +CrosshairParameters crosshair_parameters
        +set_single_crystal_parameters(params: dict[str, float])
        +get_single_crystal_parameters()
        +set_crosshair(current_experiment_type: str, DeltaE: float = None, modQ: float = None)
        +get_crosshair()
        +get_ang_Q_beam()
        +set_experiment_parameters(Ei: float, S2: float, alpha_p: float, plot_type: str)
        +get_experiment_parameters()
        +check_plot_update(deltaE)
        +calculate_graph_data()
    }


    class CrosshairParameters{
        +float crosshair_delta_e
        +float crosshair_mod_q
        +experiment_type: str
        +SingleCrystalParameters single_crystal_parameters
        +set_crosshair(current_experiment_type: str, DeltaE: float = None, modQ: float = None)
        +get_crosshair()
        +get_experiment_type()
    }

    class SingleCrystalParameters{
        +float lattice_a
        +float lattice_b
        +float lattice_c
        +float lattice_alpha
        +float lattice_beta
        +float lattice_gamma
        +float q_rlu_h
        +float q_rlu_k
        +float q_rlu_l
        +set_single_crystal_parameters(params: dict[str, float])
        +get_single_crystal_parameters()
        +calculate_crosshair_mod_q()
    }


HyspecPPT View
+++++++++++++++


.. mermaid::

 classDiagram
    HyspecPPTView "1" -->"1" ExperimentWidget
    HyspecPPTView "1" -->"1" CrosshairWidget
    HyspecPPTView "1" -->"1" SingleCrystalWidget
    HyspecPPTView "1" -->"1" SelectorWidget
    HyspecPPTView "1" -->"1" PlotWidget

    class HyspecPPTView{
        +ExperimentWidget:experiment_widget
        +SingleCrystalWidget:sc_widget
        +CrosshairWidget:crosshair_widget
        +SelectorWidget:selection_widget
        +PlotWidget:plot_widget

    }

    class PlotWidget{
        +matplotlib.Figure: figure
        +matplotlib.FigureCanvas: static_canvas
        +matplotlib.NavigationToolbar: toolbar
        +matplotlib.colorbar.Colorbar: heatmap
        +float:eline_data
        +float:qline_data
        +matplotlib.lines.Line2D:eline
        +matplotlib.lines.Line2D:qline
        +update_plot_crosshair(crosshair_data: dict)
        +update_crosshair(eline: float, qline: float)
        +update_plot(q_min: list[float], q_max: list[float], energy_transfer: list[float], q2d: list[list[float]], e2d: list[list[float]], scharpf_angle: list[list[float]],plot_label: str)
        +set_axes_meta_and_draw_plot()
    }

    class SelectorWidget{
        +str:powder_label
        +QRadioButton: powder_rb
        +str:sc_label
        +QRadioButton: sc_rb
        +selector_init(selected_label: str)
        +sc_toggle()
        +get_selected_mode_label()
    }

    class ExperimentWidget{
        +QLabel:incident_energy_ei_label
        +QLineEdit:incident_energy_ei_edit
        +QLabel:detector_tank_angle_s2_label
        +QLineEdit:detector_tank_angle_s2_edit
        +QLabel:polarization_direction_angle_p_label
        +QLineEdit:polarization_direction_angle_p_edit
        +QLabel:plot_type_label
        +QComboBox:plot_type_combobox
        +initializeCombo(options: list[str])
        +validate_inputs(*_, **__)
        +validate_all_inputs()
        +set_values(values: dict[str, Union[float, str]])
    }

    class CrosshairWidget{
        +QLabel:crosshair_delta_e_label
        +QLineEdit:crosshair_delta_e_edit
        +QLabel:crosshair_mod_q_label
        +QLineEdit:crosshair_mod_q_edit
        +QLabel:angle_q_z_label
        +QLineEdit:angle_q_z_edit
        +set_mod_q_enabled(state: bool)
        +set_values(values: dict[str, float])
        +validate_inputs(*_, **__)
        +validation_status_all_inputs()
        +validate_all_inputs()
        +set_QZ_values(angle: float)
    }

    class SingleCrystalWidget{
        +QLabel:lattice_a_label
        +QLineEdit:lattice_a_edit
        +QLabel:latticeb_label
        +QLineEdit:lattice_b_edit
        +QLabel:lattice_c_label
        +QLineEdit:lattice_c_edit
        +QLabel:lattice_alpha_label
        +QLineEdit:lattice_alpha_edit
        +QLabel:lattice_beta_label
        +QLineEdit:lattice_beta_edit
        +QLabel:lattice_gamma_label
        +QLineEdit:lattice_gamma_edit
        +QLabel:q_rlu_h_label
        +QLineEdit:q_rlu_h_edit
        +QLabel:q_rlu_k_label
        +QLineEdit:q_rlu_k_edit
        +QLabel:q_rlu_l_label
        +QLineEdit:q_rlu_l_edit
        +set_values(values: dict[str, float])
        +validate_inputs(*_, **__)
        +validate_angles()
        +validate_all_inputs()
    }



HyspecPPT Presenter
++++++++++++++++++++++

.. mermaid::

 classDiagram
    HyspecPPTPresenter "1" -->"1" HyspecPPTModel
    HyspecPPTPresenter "1" -->"1" HyspecPPTView

    class HyspecPPTPresenter{
        -HyspecPPTModel:model
        -HyspecPPTView:view
        +handle_field_values_update()
        +handle_switch_to_powder()
        +handle_switch_to_sc()
        +handle_QZ_angle()
    }

    class HyspecPPTModel{
        #from above
    }

    class HyspecPPTView{
        #from above
    }

The Hppt Model and View are unaware of one another. The Presenter is the connecting link that has a direct access and interacts with both.
The Presenter describes the main workflows that require communication and coordination between the Model and View through the Presenter. Additionally, the widgets' data initialization come from the model initialization and passed to the View.
Any value processing and/or filtering to match the requirements and logic of the View and Model side should happen on the Presenter.


#. Application Start - HyspecPPTView Initialization. All default values are retrieved from the settings file.

    .. mermaid::

        sequenceDiagram
            participant View
            participant Presenter
            participant Model

            Note over View,Model:  HyspecPPTView Initialization
            Presenter->>Model: A. Get Experiment parameters
            Presenter->>View: Set Experiment parameters (experiment_widget.set_values)
            Note left of View: Display Experiment parameters values
            Note left of View: experiment_parameters_update is triggered

            Presenter->>Model: B. Get SingleCrystal parameters
            Note left of Presenter: Get the available plot types from the experiment_settings file
            Presenter->>View: Set SingleCrystal parameters (singlecrystal_widget.set_parameters)
            Note left of View: Display SingleCrystal parameters values
            Note left of View: handle_field_values_update is triggered

            Presenter->>Model: C. Get default experiment mode (Single Crystal)
            Presenter->>View: Set experiment mode (selection_widget.selector_init)
            Note left of View: Workflow continues for selecting experiment type = Single Crystal
            Note left of View: handle_field_values_update is triggered

#. This describes the sequence of events happening among M-V-P when CrosshairWidget parameters are updated in order to see a new plot : handle_field_values_update()

    * Valid Status with Replot:

        .. mermaid::

            sequenceDiagram
                participant View
                participant Presenter
                participant Model

                Note over View,Model: Plot draw due to any CrosshairWidget parameter update
                Note left of View: User updates a parameter at CrosshairWidget
                Note left of View: Check the validation status of all CrosshairWidget parameters (CrosshairWidget.validate_all_inputs)
                View->>Presenter: Emit the valid signal and pass the crosshair parameters
                Presenter->>View: Get the experiment type
                Presenter->>Model: Send crosshair_delta_e to decide on replot
                Note right of Model: Calculate replot based on new delta_e value, and previous Emin, delta_e values
                Model->>Presenter: Returns the replot to True
                Presenter->>Model: Set crosshair data (set_crosshair_data)
                Note right of Model: Store the crosshair data
                Presenter->>Model: Get momentum transfer angle
                Model->>Presenter: Return momentum transfer angle
                Presenter->>Model: Calculate plot data (calculate_graph_data)
                Note right of Model: Calculate plot dictionary data
                Model->>Presenter: Return graph data dictionary
                Presenter->>View: Return graph data (plot_widget.update_plot)
                Note left of View: Draw the (colormap) heatmap
                Presenter->>View: Return graph data (plot_widget.update_crosshair)
                Note left of View: Draw the crosshair
                Presenter->>View: Display momentum transfer angle
                Note left of View: Update momentum transfer angle


    * Valid Status without Replot:

        .. mermaid::

            sequenceDiagram
                participant View
                participant Presenter
                participant Model

                Note over View,Model: Plot draw due to any CrosshairWidget parameter update
                Note left of View: User updates a parameter at CrosshairWidget
                Note left of View: Check the validation status of all CrosshairWidget parameters (CrosshairWidget.validate_all_inputs)
                View->>Presenter: Emit the valid signal and pass the crosshair parameters
                Presenter->>Model: Send crosshair_delta_e to decide on replot
                Note right of Model: Calculate replot based on new delta_e value, and previous Emin, delta_e values
                Model->>Presenter: Returns the replot to False
                Presenter->>Model: Set crosshair data (set_crosshair_data)
                Note right of Model: Store the crosshair data
                Presenter->>Model: Get momentum transfer angle
                Model->>Presenter: Return momentum transfer angle
                Presenter->>View: Return graph data (plot_widget.update_crosshair)
                Note left of View: Draw the crosshair
                Presenter->>View: Display momentum transfer angle
                Note left of View: Update momentum transfer angle

    * Invalid Status:

    .. mermaid::

        sequenceDiagram
            participant View
            participant Presenter
            participant Model

            Note over View,Model: CrosshairWidget parameter update
            Note left of View: User updates a parameter at CrosshairWidget
            Note Left of View: Check the validation status of all CrosshairWidget parameters (CrosshairWidget.validate_all_inputs)
            Note Left of View: Red borders appear (validate_inputs) no signal is emitted


#. This describes the sequence of events happening among M-V-P when ExperimentWidget parameters are updated in order to see a new plot : handle_field_values_update()

    * Valid Status:

        .. mermaid::

            sequenceDiagram
                participant View
                participant Presenter
                participant Model

                Note over View,Model: Plot draw due to any ExperimentWidget parameter update
                Note left of View: User updates a parameter at ExperimentWidget
                Note left of View: Check the validation status of all ExperimentWidget parameters (ExperimentWidget.validate_all_inputs)
                View->>Presenter: Emit the valid signal and pass the experiment parameters
                Presenter->>Model: Set the parameters (set_experiment_data)
                Presenter->>Model: Calculate plot data (calculate_graph_data)
                Note right of Model: Calculate plot dictionary data
                Model->>Presenter: Return graph data dictionary
                Presenter->>View: Return graph data (plot_widget.update_plot)
                Presenter->>Model: Get momentum transfer angle
                Model->>Presenter: Return momentum transfer angle
                Note left of View: Draw the (colormap) heatmap
                Presenter->>View: Display momentum transfer angle
                Note left of View: Update momentum transfer angle

    * Invalid Status:

    .. mermaid::

        sequenceDiagram
            participant View
            participant Presenter
            participant Model

            Note over View,Model: ExperimentWidget parameter update
            Note left of View: User updates a parameter at ExperimentWidget
            Note left of View: Check the validation status of all ExperimentWidget parameters (ExperimentWidget.validate_all_inputs)
            Note Left of View: Red borders appear (validate_inputs) no signal is emitted


#. This describes the sequence of events happening among M-V-P when Single Crystal parameters are updated in order to see a new plot : handle_field_values_update()

    * Valid Status:

        .. mermaid::

            sequenceDiagram
                participant View
                participant Presenter
                participant Model

                Note over View,Model: Plot draw due to any SingleCrystalWidget parameter update
                Note left of View: User updates a parameter at SingleCrystalWidget
                Note left of View: Check the validation status of all SingleCrystalWidget parameters (SingleCrystalWidget.validate_all_inputs)
                View->>Presenter: Emit the valid signal and pass the single crystal parameters
                Presenter->>Model: Set the parameters (set_single_crystal_data)
                Presenter->>Model: Get the new crosshair data (get_crosshair_data)
                Presenter->>Model: Get momentum transfer angle
                Model->>Presenter: Return momentum transfer angle
                Presenter->>View: Display the crosshair data (crosshair_widget.set_values)
                Note left of Presenter: Check the validation status of all crosshair_widget parameters (CrosshairWidget.validation_status_all_inputs) is valid
                Presenter->>View: Return graph data (plot_widget.update_crosshair)
                Note left of View: Draw the crosshair
                Presenter->>View: Display momentum transfer angle
                Note left of View: Update momentum transfer angle

    * Invalid Status:

        .. mermaid::

            sequenceDiagram
                participant View
                participant Presenter
                participant Model

                Note over View,Model: Plot draw due to any SingleCrystalWidget parameter update
                Note left of View: User updates a parameter at SingleCrystalWidget
                Note left of View: Check the validation status of all SingleCrystalWidget parameters (SingleCrystalWidget.validate_all_inputs)
                View->>Presenter: Emit the valid signal and pass the single crystal parameters
                Presenter->>Model: Set the parameters (set_single_crystal_data)
                Presenter->>Model: Get the new crosshair data (get_crosshair_data)
                Presenter->>View: Display the crosshair data (crosshair_widget.set_values)
                Note left of Presenter: Check the validation status of all crosshair_widget parameters (CrosshairWidget.validation_status_all_inputs) is invalid
                Note left of Presenter: Nothing

#. This describes the sequence of events happening among M-V-P when user selects the "Powder" mode : handle_switch_to_powder()

    .. mermaid::

        sequenceDiagram
            participant View
            participant Presenter
            participant Model

            Note over View,Model: Updates due to switching to Powder Mode
            Note left of View: User selects the Powder radio button
            View->>Presenter: Trigger the update
            Presenter->>View: Update fields' visibility for powder case(field_visibility_in_Powder)
            Note left of View: Hide the SingleCrystalWidget
            Note left of View: Make crosshair_mod_q_edit field editable
            Presenter->>Model: Set crosshair parameters with the experiment_type="powder" (set_crosshair_data)
            Note right of Model: Store the crosshair data
            Presenter->>Model: Get crosshair parameters for the experiment_type="powder"(get_crosshair_data)
            Note right of Model:  Return the mod_q and the delta_e values
            Model->>Presenter: Return the crosshair data
            Presenter->>View: Return the crosshair data
            Note left of View: Display the data in the crosshair_widget
            Presenter->>View: Return graph data (plot_widget.update_crosshair)
            Note left of View: Draw the crosshair
            Presenter->>Model: Get experiment parameters (get_experiment_data)
            Presenter->>View: Set experiment parameters (experiment_widget.set_values)
            Note left of View: Display experiment parameters values
            Note left of View: handle_field_values_update is triggered


#. This describes the sequence of events happening among M-V-P when user selects the "Single Crystal" mode : handle_switch_to_sc()

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

The parameters' default values for the application are stored in a file, experiment_settings.py, next to the model file. They are imported
in the HyspecPPT Model file and used during the Experiment object's initialization and data calculations. The options for experiment and plot types are used in HyspecPPT Model and View files.
More specifically the parameters with their values are:

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
