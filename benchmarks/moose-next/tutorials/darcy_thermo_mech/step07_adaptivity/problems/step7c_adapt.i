[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 30
  ny = 3
  xmax = 0.304 # Length of test chamber
  ymax = 0.0257 # Test chamber radius
  uniform_refine = 3
[]

[Variables]
  [pressure]
  []
  [temperature]
    initial_condition = 300 # Start at room temperature
  []
[]

[AuxVariables]
  [velocity_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [velocity_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [velocity_z]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [darcy_pressure]
    type = DarcyPressure
    variable = pressure
  []
  [heat_conduction]
    type = ADHeatConduction
    variable = temperature
  []
  [heat_conduction_time_derivative]
    type = ADHeatConductionTimeDerivative
    variable = temperature
  []
  [heat_convection]
    type = DarcyAdvection
    variable = temperature
    pressure = pressure
  []
[]

[AuxKernels]
  [velocity_x]
    type = DarcyVelocity
    variable = velocity_x
    component = x
    execute_on = timestep_end
    pressure = pressure
  []
  [velocity_y]
    type = DarcyVelocity
    variable = velocity_y
    component = y
    execute_on = timestep_end
    pressure = pressure
  []
  [velocity_z]
    type = DarcyVelocity
    variable = velocity_z
    component = z
    execute_on = timestep_end
    pressure = pressure
  []
[]

[BCs]
  [inlet]
    type = DirichletBC
    variable = pressure
    boundary = left
    value = 4000 # (Pa) From Figure 2 from paper.  First data point for 1mm spheres.
  []
  [outlet]
    type = DirichletBC
    variable = pressure
    boundary = right
    value = 0 # (Pa) Gives the correct pressure drop from Figure 2 for 1mm spheres
  []
  [inlet_temperature]
    type = FunctionDirichletBC
    variable = temperature
    boundary = left
    function = 'if(t<0,350+50*t,350)'
  []

  [outlet_temperature]
    type = HeatConductionOutflow
    variable = temperature
    boundary = right
  []
[]

[Materials]
  [column]
    type = PackedColumn
    radius = 1
    temperature = temperature
  []
[]

[Problem]
  type = FEProblem
  coord_type = RZ
  rz_coord_axis = X
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  automatic_scaling = true

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

  end_time = 100
  dt = 0.25
  start_time = -1

  steady_state_tolerance = 1e-5
  steady_state_detection = true

  [TimeStepper]
    type = FunctionDT
    function = 'if(t<0,0.1,0.25)'
  []
[]

[Outputs]
  exodus = true
[]

[Adaptivity]
  marker = error_frac
  max_h_level = 3
  [Indicators]
    [temperature_jump]
      type = GradientJumpIndicator
      variable = temperature
      scale_by_flux_faces = true
    []
  []
  [Markers]
    [error_frac]
      type = ErrorFractionMarker
      coarsen = 0.15
      indicator = temperature_jump
      refine = 0.7
    []
  []
[]
