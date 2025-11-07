.. _hyspecpptfields:

=======================
Fields and Validation
=======================

The main functionality of the tool is to display a plot, based on the filled in parameters, automatically; no plot button exists.

Overall, users can:
   * plot crosshair and heatmap from Powder experiment
   * plot crosshair and heatmap from Single Crystal experiment
   * switch between Powder and Single Crystal modes while keeping the previous valid state
   * click on "Help" button that opens up a readthedocs user documentation
   * show the angle between momentum transfer and incident beam at the crosshair position. A value of `nan` means that conservation of energy and momentum does not permit measurement at that point.

Fields
--------

Below are the fields of SingleCrystal and Powder experiment models
For the purpose of ensuring only valid fields are passed to the Model, custom qt validators are included that perform more complex checks.

.. list-table:: Common Fields
  :header-rows: 1

  * - Field
    - Type
    - Value Origin
    - Default
    - Additional validation
    - Mandatory
  * - Polarization Type
    - String
    - predefined choices:[Powder, Single Crystal]
    - Powder
    -
    - yes
  * - Ei (Incident Energy) - meV
    - Double
    -
    - 20
    - 0 < Ei < 100
    - yes
  * - S2 (HYSPEC Detector Tank Angle)
    - Double
    -
    - 30
    - -100 < S2< 100 && (S2 > 30 || S2 < -30)
    - yes
  * - Ap (Polarization Direction Angle)
    - Double
    -
    - 0
    - -180 < Ap < 180
    - yes
  * - Delta E
    - Double
    -
    - 0
    -
    - yes
  * - mod Q (\|Q\|)
    - Double
    -
    - 0
    - 0 <= \|Q\| <=10
    - yes
  * - Plot Type
    - String
    - predefined choices: :math:`[ \alpha_s, \cos^2(\alpha_s),  (1+\cos^2(\alpha_s))/2 ]`
    - :math:`\cos^2(\alpha_s)`
    -
    - yes


Below are the additional fields of SingleCrystal Model


.. list-table:: Single Crystal Model Additional Fields
  :header-rows: 1

  * - Field
    - Type
    - Default
    - Additional validation
    - Mandatory
  * - a
    - Double
    - 1
    - 1 < a < 100
    - yes
  * - b
    - Double
    - 1
    - 1 < b < 100
    - yes
  * - c
    - Double
    - 1
    - 1 < c < 100
    - yes
  * - alpha
    - Double
    - 90
    - 30 < alpha < 150
    - yes
  * - beta
    - Double
    - 90
    - 30 < beta < 150
    - yes
  * - gamma
    - Double
    - 90
    - 30 < gamma < 150
    - yes
  * - H
    - Double
    - 0
    - -100 < H < 100
    - yes
  * - K
    - Double
    - 0
    - -100 < K < 100
    - yes
  * - L
    - Double
    - 0
    - -100 < L < 100
    - yes

Validation
----------

Regarding validation, if all fields are valid, then the front end (View) triggers the backend (Model) to send the current parameters and receive the new data to plot the graphs (heatmap and/or crosshair).
If a user types an invalid value, then a red border appears on the related field(s).


Front end side validation includes:
   * required fields
   * field types
   * threshold limits: Ei, S2, Ap, modQ, and single crystal parameters


Backend side validation includes:
  * Emin recalculation due to DeltaE update


Inter-Field Validations
------------------------

Polarization Type:
  * If Polarization Type set to "Single Crystal", all parameters of Single Crystal block are required. The valid default/model-stored values are set at the appropriate fields.
  * If Polarization Type set to "Powder", all parameters of Single Crystal block are hidden and not required. The valid default/model-stored values are set at the appropriate fields.

modQ (\|Q\|):
  * If Polarization Type set to "Single Crystal", it is a read-only field. The value is returned from the backend after all Single crystal parameters are filled in.
  * If Polarization Type set to "Powder", user can fill the value in.

alpha, beta , gamma; all the below conditions need to be met for the fields to be valid:
  * All three angles' values are less than 360 degrees: (alpha + beta+ gamma) <=360
  * They can form a triangle: (alpha + beta) < gamma and (alpha + beta) < gamma and (beta + gamma) < alpha

Emin, DeltaE:
  * A change in deltaE can trigger an update to Emin; the process is hidden to the user
