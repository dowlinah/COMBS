CDF      
      
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes      	   num_elem      
num_el_blk        num_node_sets         num_side_sets         num_el_in_blk1        num_nod_per_el1       num_side_ss1      num_side_ss2      num_side_ss3      num_side_ss4      num_nod_ns1       num_nod_ns2       num_nod_ns3       num_nod_ns4       num_nod_var       num_glo_var       num_info  "         api_version       @�
=   version       @�
=   floating_point_word_size            	file_size               int64_status             title         )picard_master_2subapps_out_sub_no_fail.e       maximum_name_length                 "   
time_whole                            h�   	eb_status                             	�   eb_prop1               name      ID              	�   	ns_status         	                    	�   ns_prop1      	         name      ID              	�   	ss_status         
                    	�   ss_prop1      
         name      ID              	�   coordx                      H      	�   coordy                      H      
   eb_names                       $      
`   ns_names      	                 �      
�   ss_names      
                 �         
coor_names                         D      �   node_num_map                    $      �   connect1                  	elem_type         QUAD4         @      �   elem_num_map                          4   elem_ss1                          D   side_ss1                          L   elem_ss2                          T   side_ss2                          \   elem_ss3                          d   side_ss3                          l   elem_ss4                          t   side_ss4                          |   node_ns1                          �   node_ns2                          �   node_ns3                          �   node_ns4                          �   vals_nod_var1                          H      h�   vals_nod_var2                          H      i0   name_nod_var                       D      �   name_glo_var                       $      �   vals_glo_var                             ix   info_records                      [�                                                                       ?�      ?�              ?�      ?�      ?�              ?�                      ?�      ?�              ?�      ?�      ?�      ?�                                          bottom                           left                             right                            top                              top                              right                            bottom                           left                                                                                                                            	                                             	                                                                                          	         	u                                v                                  elem_average_value                  ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                    -i picard_master.i MultiApps/active=sub sub_no_fail Outputs/file_base=picard... _master_2subapps_out sub:Outputs/file_base=picard_master_2subapps_out_sub sub... _no_fail:Outputs/file_base=picard_master_2subapps_out_sub_no_fail --error-unu... sed --error-override --no-gdb-backtrace### Input File ###                                                                                                         []                                                                                 inactive                       =                                                 element_order                  = AUTO                                            order                          = AUTO                                            side_order                     = AUTO                                            type                           = GAUSS                                         []                                                                                                                                                                [AuxVariables]                                                                                                                                                      [./u]                                                                              block                        = INVALID                                           family                       = LAGRANGE                                          inactive                     =                                                   initial_condition            = INVALID                                           order                        = FIRST                                             outputs                      = INVALID                                           initial_from_file_timestep   = LATEST                                            initial_from_file_var        = INVALID                                         [../]                                                                          []                                                                                                                                                                [BCs]                                                                                                                                                               [./left_v]                                                                         boundary                     = left                                              control_tags                 = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 inactive                     =                                                   isObjectAction               = 1                                                 matrix_tags                  = system                                            type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = v                                                 vector_tags                  = nontime                                           diag_save_in                 = INVALID                                           save_in                      = INVALID                                           seed                         = 0                                                 value                        = 1                                               [../]                                                                                                                                                             [./right_v]                                                                        boundary                     = right                                             control_tags                 = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 inactive                     =                                                   isObjectAction               = 1                                                 matrix_tags                  = system                                            type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = v                                                 vector_tags                  = nontime                                           diag_save_in                 = INVALID                                           save_in                      = INVALID                                           seed                         = 0                                                 value                        = 0                                               [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      inactive                       =                                                 isObjectAction                 = 1                                               type                           = Transient                                       abort_on_solve_fail            = 0                                               compute_initial_residual_before_preset_bcs = 0                                   contact_line_search_allowed_lambda_cuts = 2                                      contact_line_search_ltol       = INVALID                                         control_tags                   =                                                 dt                             = 0.1                                             dtmax                          = 1e+30                                           dtmin                          = 2e-14                                           enable                         = 1                                               end_time                       = 1e+30                                           l_abs_step_tol                 = -1                                              l_max_its                      = 10000                                           l_tol                          = 1e-05                                           line_search                    = default                                         line_search_package            = petsc                                           max_xfem_update                = 4294967295                                      mffd_type                      = wp                                              n_startup_steps                = 0                                               nl_abs_step_tol                = 1e-50                                           nl_abs_tol                     = 1e-10                                           nl_max_funcs                   = 10000                                           nl_max_its                     = 50                                              nl_rel_step_tol                = 1e-50                                           nl_rel_tol                     = 1e-08                                           no_fe_reinit                   = 0                                               num_steps                      = 2                                               petsc_options                  = INVALID                                         petsc_options_iname            = '-pc_type -pc_hypre_type'                       petsc_options_value            = 'hypre boomeramg'                               picard_abs_tol                 = 1e-50                                           picard_max_its                 = 1                                               picard_rel_tol                 = 1e-08                                           relaxation_factor              = 1                                               relaxed_variables              =                                                 reset_dt                       = 0                                               restart_file_base              =                                                 scheme                         = implicit-euler                                  snesmf_reuse_base              = 1                                               solve_type                     = PJFNK                                           splitting                      = INVALID                                         ss_check_tol                   = 1e-08                                           ss_tmin                        = 0                                               start_time                     = 0                                               steady_state_detection         = 0                                               steady_state_start_time        = 0                                               steady_state_tolerance         = 1e-08                                           time_period_ends               = INVALID                                         time_period_starts             = INVALID                                         time_periods                   = INVALID                                         timestep_tolerance             = 2e-14                                           trans_ss_check                 = 0                                               update_xfem_at_timestep_begin  = 0                                               use_multiapp_dt                = 0                                               verbose                        = 0                                             []                                                                                                                                                                [Kernels]                                                                                                                                                           [./diff_v]                                                                         inactive                     =                                                   isObjectAction               = 1                                                 type                         = Diffusion                                         block                        = INVALID                                           control_tags                 = Kernels                                           diag_save_in                 = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 matrix_tags                  = system                                            save_in                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 0                                                 variable                     = v                                                 vector_tags                  = nontime                                         [../]                                                                                                                                                             [./force_v]                                                                        inactive                     =                                                   isObjectAction               = 1                                                 type                         = CoupledForce                                      block                        = INVALID                                           coef                         = 1                                                 control_tags                 = Kernels                                           diag_save_in                 = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 matrix_tags                  = system                                            save_in                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 0                                                 v                            = u                                                 variable                     = v                                                 vector_tags                  = nontime                                         [../]                                                                          []                                                                                                                                                                [Mesh]                                                                             inactive                       =                                                 displacements                  = INVALID                                         block_id                       = INVALID                                         block_name                     = INVALID                                         boundary_id                    = INVALID                                         boundary_name                  = INVALID                                         construct_side_list_from_node_list = 0                                           ghosted_boundaries             = INVALID                                         ghosted_boundaries_inflation   = INVALID                                         isObjectAction                 = 1                                               second_order                   = 0                                               skip_partitioning              = 0                                               type                           = GeneratedMesh                                   uniform_refine                 = 0                                               allow_renumbering              = 1                                               bias_x                         = 1                                               bias_y                         = 1                                               bias_z                         = 1                                               centroid_partitioner_direction = INVALID                                         construct_node_list_from_side_list = 1                                           control_tags                   =                                                 dim                            = 2                                               elem_type                      = INVALID                                         enable                         = 1                                               gauss_lobatto_grid             = 0                                               ghosting_patch_size            = INVALID                                         max_leaf_size                  = 10                                              nemesis                        = 0                                               nx                             = 2                                               ny                             = 2                                               nz                             = 1                                               parallel_type                  = DEFAULT                                         partitioner                    = default                                         patch_size                     = 40                                              patch_update_strategy          = never                                           xmax                           = 1                                               xmin                           = 0                                               ymax                           = 1                                               ymin                           = 0                                               zmax                           = 1                                               zmin                           = 0                                             []                                                                                                                                                                [Mesh]                                                                           []                                                                                                                                                                [Outputs]                                                                          append_date                    = 0                                               append_date_format             = INVALID                                         checkpoint                     = 0                                               color                          = 1                                               console                        = 1                                               controls                       = 0                                               csv                            = 0                                               dofmap                         = 0                                               execute_on                     = 'INITIAL TIMESTEP_END'                          exodus                         = 1                                               file_base                      = picard_master_2subapps_out_sub_no_fail          gmv                            = 0                                               gnuplot                        = 0                                               hide                           = INVALID                                         inactive                       =                                                 interval                       = 1                                               nemesis                        = 0                                               output_if_base_contains        = INVALID                                         perf_graph                     = 0                                               print_linear_residuals         = 1                                               print_mesh_changed_info        = 0                                               print_perf_log                 = 0                                               show                           = INVALID                                         solution_history               = 0                                               sync_times                     =                                                 tecplot                        = 0                                               vtk                            = 0                                               xda                            = 0                                               xdr                            = 0                                             []                                                                                                                                                                [Postprocessors]                                                                                                                                                    [./elem_average_value]                                                             inactive                     =                                                   isObjectAction               = 1                                                 type                         = ElementAverageValue                               allow_duplicate_execution_on_initial = 0                                         block                        = INVALID                                           control_tags                 = Postprocessors                                    enable                       = 1                                                 execute_on                   = 'INITIAL TIMESTEP_END'                            implicit                     = 1                                                 outputs                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 0                                                 variable                     = v                                               [../]                                                                          []                                                                                                                                                                [Variables]                                                                                                                                                         [./v]                                                                              block                        = INVALID                                           eigen                        = 0                                                 family                       = LAGRANGE                                          inactive                     =                                                   initial_condition            = INVALID                                           order                        = FIRST                                             outputs                      = INVALID                                           scaling                      = 1                                                 initial_from_file_timestep   = LATEST                                            initial_from_file_var        = INVALID                                         [../]                                                                          []                                                                                                                                                                                                                                                 ?�������                                                                        ?�      ?�      ?�      ?�                      ?�������?�              ?�      ?�333334                                                                        ?�      ?�      ?�      ?�                      ?�������?�              ?�      