CDF      
      
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes      	   num_elem      
num_el_blk        num_node_sets         num_side_sets         num_el_in_blk1        num_nod_per_el1       num_side_ss1      num_side_ss2      num_side_ss3      num_side_ss4      num_nod_ns1       num_nod_ns2       num_nod_ns3       num_nod_ns4       num_nod_var       num_info  �         api_version       @�
=   version       @�
=   floating_point_word_size            	file_size               int64_status             title         ad_value_out.e     maximum_name_length                 "   
time_whole                            �   	eb_status                             	`   eb_prop1               name      ID              	d   	ns_status         	                    	h   ns_prop1      	         name      ID              	x   	ss_status         
                    	�   ss_prop1      
         name      ID              	�   coordx                      H      	�   coordy                      H      	�   eb_names                       $      
8   ns_names      	                 �      
\   ss_names      
                 �      
�   
coor_names                         D      d   node_num_map                    $      �   connect1                  	elem_type         QUAD4         @      �   elem_num_map                             elem_ss1                             side_ss1                          $   elem_ss2                          ,   side_ss2                          4   elem_ss3                          <   side_ss3                          D   elem_ss4                          L   side_ss4                          T   node_ns1                          \   node_ns2                          h   node_ns3                          t   node_ns4                          �   vals_nod_var1                          H      �   vals_nod_var2                          H      �\   vals_nod_var3                          H      ��   vals_nod_var4                          H      ��   name_nod_var                       �      �   info_records                      |�                                                                       ?�      ?�              ?�      ?�      ?�              ?�                      ?�      ?�              ?�      ?�      ?�      ?�                                          bottom                           left                             right                            top                              top                              right                            bottom                           left                                                                                                                            	                                             	                                                                                          	         	u                                u_jac                            v                                v_jac                            ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                    /Users/slauae/projects/moose/test/moose_test-devel -i ad_value.i --error --e... rror-unused --error-override --no-gdb-backtrace### Version Info ###              Framework Information:                                                           MOOSE Version:           git commit 3febedaa6d on 2019-05-14                     LibMesh Version:         daae6d65b81ab245f48fb51f6d581c0af2a53e29                PETSc Version:           3.9.4                                                   Current Time:            Tue May 14 11:34:37 2019                                Executable Timestamp:    Tue May 14 11:23:53 2019                                                                                                                                                                                                  ### Input File ###                                                                                                                                                []                                                                                 inactive                       = (no_default)                                    element_order                  = AUTO                                            order                          = AUTO                                            side_order                     = AUTO                                            type                           = GAUSS                                         []                                                                                                                                                                [AuxVariables]                                                                                                                                                      [./u_jac]                                                                          block                        = INVALID                                           family                       = LAGRANGE                                          inactive                     = (no_default)                                      initial_condition            = INVALID                                           order                        = FIRST                                             outputs                      = INVALID                                           initial_from_file_timestep   = LATEST                                            initial_from_file_var        = INVALID                                         [../]                                                                                                                                                             [./v_jac]                                                                          block                        = INVALID                                           family                       = LAGRANGE                                          inactive                     = (no_default)                                      initial_condition            = INVALID                                           order                        = FIRST                                             outputs                      = INVALID                                           initial_from_file_timestep   = LATEST                                            initial_from_file_var        = INVALID                                         [../]                                                                          []                                                                                                                                                                [BCs]                                                                                                                                                               [./left]                                                                           boundary                     = left                                              control_tags                 = INVALID                                           displacements                = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 inactive                     = (no_default)                                      isObjectAction               = 1                                                 matrix_tags                  = system                                            type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = u                                                 vector_tags                  = nontime                                           diag_save_in                 = INVALID                                           save_in                      = INVALID                                           seed                         = 0                                                 value                        = 0                                               [../]                                                                                                                                                             [./left_v]                                                                         boundary                     = left                                              control_tags                 = INVALID                                           displacements                = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 inactive                     = (no_default)                                      isObjectAction               = 1                                                 matrix_tags                  = system                                            type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = v                                                 vector_tags                  = nontime                                           diag_save_in                 = INVALID                                           save_in                      = INVALID                                           seed                         = 0                                                 value                        = 0                                               [../]                                                                                                                                                             [./right]                                                                          boundary                     = right                                             control_tags                 = INVALID                                           displacements                = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 inactive                     = (no_default)                                      isObjectAction               = 1                                                 matrix_tags                  = system                                            type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = u                                                 vector_tags                  = nontime                                           diag_save_in                 = INVALID                                           save_in                      = INVALID                                           seed                         = 0                                                 value                        = 1                                               [../]                                                                                                                                                             [./right_v]                                                                        boundary                     = right                                             control_tags                 = INVALID                                           displacements                = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 inactive                     = (no_default)                                      isObjectAction               = 1                                                 matrix_tags                  = system                                            type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = v                                                 vector_tags                  = nontime                                           diag_save_in                 = INVALID                                           save_in                      = INVALID                                           seed                         = 0                                                 value                        = 1                                               [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      inactive                       = (no_default)                                    isObjectAction                 = 1                                               type                           = Steady                                          accept_on_max_picard_iteration = 0                                               compute_initial_residual_before_preset_bcs = 0                                   contact_line_search_allowed_lambda_cuts = 2                                      contact_line_search_ltol       = INVALID                                         control_tags                   = (no_default)                                    disable_picard_residual_norm_check = 0                                           enable                         = 1                                               l_abs_step_tol                 = -1                                              l_max_its                      = 10000                                           l_tol                          = 1e-05                                           line_search                    = default                                         line_search_package            = petsc                                           max_xfem_update                = 4294967295                                      mffd_type                      = wp                                              nl_abs_step_tol                = 1e-50                                           nl_abs_tol                     = 1e-50                                           nl_max_funcs                   = 10000                                           nl_max_its                     = 50                                              nl_rel_step_tol                = 1e-50                                           nl_rel_tol                     = 1e-08                                           petsc_options                  = INVALID                                         petsc_options_iname            = '-pc_type -pc_hypre_type'                       petsc_options_value            = 'hypre boomeramg'                               picard_abs_tol                 = 1e-50                                           picard_force_norms             = 0                                               picard_max_its                 = 1                                               picard_rel_tol                 = 1e-08                                           relaxation_factor              = 1                                               relaxed_variables              = (no_default)                                    restart_file_base              = (no_default)                                    snesmf_reuse_base              = 1                                               solve_type                     = NEWTON                                          splitting                      = INVALID                                         time                           = 0                                               update_xfem_at_timestep_begin  = 0                                             []                                                                                                                                                                [Kernels]                                                                                                                                                           [./ad_value_test]                                                                  inactive                     = (no_default)                                      isObjectAction               = 1                                                 type                         = ADValueTest                                       block                        = INVALID                                           control_tags                 = Kernels                                           diag_save_in                 = u_jac                                             displacements                = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 matrix_tags                  = system                                            save_in                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 0                                                 variable                     = u                                                 vector_tags                  = nontime                                         [../]                                                                                                                                                             [./diff]                                                                           inactive                     = (no_default)                                      isObjectAction               = 1                                                 type                         = Diffusion                                         block                        = INVALID                                           control_tags                 = Kernels                                           diag_save_in                 = INVALID                                           displacements                = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 matrix_tags                  = system                                            save_in                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 0                                                 variable                     = u                                                 vector_tags                  = nontime                                         [../]                                                                                                                                                             [./diff_v]                                                                         inactive                     = (no_default)                                      isObjectAction               = 1                                                 type                         = Diffusion                                         block                        = INVALID                                           control_tags                 = Kernels                                           diag_save_in                 = INVALID                                           displacements                = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 matrix_tags                  = system                                            save_in                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 0                                                 variable                     = v                                                 vector_tags                  = nontime                                         [../]                                                                                                                                                             [./value_test_v]                                                                   inactive                     = (no_default)                                      isObjectAction               = 1                                                 type                         = ValueTest                                         block                        = INVALID                                           control_tags                 = Kernels                                           diag_save_in                 = v_jac                                             displacements                = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 matrix_tags                  = system                                            save_in                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 0                                                 variable                     = v                                                 vector_tags                  = nontime                                         [../]                                                                          []                                                                                                                                                                [Mesh]                                                                             displacements                  = INVALID                                         inactive                       = (no_default)                                    use_displaced_mesh             = 1                                               include_local_in_ghosting      = 0                                               output_ghosting                = 0                                               block_id                       = INVALID                                         block_name                     = INVALID                                         boundary_id                    = INVALID                                         boundary_name                  = INVALID                                         construct_side_list_from_node_list = 0                                           ghosted_boundaries             = INVALID                                         ghosted_boundaries_inflation   = INVALID                                         isObjectAction                 = 1                                               second_order                   = 0                                               skip_partitioning              = 0                                               type                           = GeneratedMesh                                   uniform_refine                 = 0                                               allow_renumbering              = 1                                               bias_x                         = 1                                               bias_y                         = 1                                               bias_z                         = 1                                               centroid_partitioner_direction = INVALID                                         construct_node_list_from_side_list = 1                                           control_tags                   = (no_default)                                    dim                            = 2                                               elem_type                      = INVALID                                         enable                         = 1                                               gauss_lobatto_grid             = 0                                               ghosting_patch_size            = INVALID                                         max_leaf_size                  = 10                                              nemesis                        = 0                                               nx                             = 2                                               ny                             = 2                                               nz                             = 1                                               parallel_type                  = DEFAULT                                         partitioner                    = default                                         patch_size                     = 40                                              patch_update_strategy          = never                                           xmax                           = 1                                               xmin                           = 0                                               ymax                           = 1                                               ymin                           = 0                                               zmax                           = 1                                               zmin                           = 0                                             []                                                                                                                                                                [Mesh]                                                                           []                                                                                                                                                                [Mesh]                                                                           []                                                                                                                                                                [Outputs]                                                                          append_date                    = 0                                               append_date_format             = INVALID                                         checkpoint                     = 0                                               color                          = 1                                               console                        = 1                                               controls                       = 0                                               csv                            = 0                                               dofmap                         = 0                                               execute_on                     = 'INITIAL TIMESTEP_END'                          exodus                         = 1                                               file_base                      = INVALID                                         gmv                            = 0                                               gnuplot                        = 0                                               hide                           = INVALID                                         inactive                       = (no_default)                                    interval                       = 1                                               nemesis                        = 0                                               output_if_base_contains        = INVALID                                         perf_graph                     = 0                                               print_linear_residuals         = 1                                               print_mesh_changed_info        = 0                                               print_perf_log                 = 0                                               show                           = INVALID                                         solution_history               = 0                                               sync_times                     = (no_default)                                    tecplot                        = 0                                               vtk                            = 0                                               xda                            = 0                                               xdr                            = 0                                             []                                                                                                                                                                [Preconditioning]                                                                                                                                                   [./smp]                                                                            inactive                     = (no_default)                                      isObjectAction               = 1                                                 type                         = SMP                                               control_tags                 = Preconditioning                                   coupled_groups               = INVALID                                           enable                       = 1                                                 full                         = 0                                                 ksp_norm                     = unpreconditioned                                  mffd_type                    = wp                                                off_diag_column              = 'v u'                                             off_diag_row                 = 'v u'                                             pc_side                      = default                                           petsc_options                = INVALID                                           petsc_options_iname          = INVALID                                           petsc_options_value          = INVALID                                           solve_type                   = INVALID                                         [../]                                                                          []                                                                                                                                                                [Variables]                                                                                                                                                         [./u]                                                                              block                        = INVALID                                           eigen                        = 0                                                 family                       = LAGRANGE                                          inactive                     = (no_default)                                      initial_condition            = INVALID                                           order                        = FIRST                                             outputs                      = INVALID                                           scaling                      = 1                                                 initial_from_file_timestep   = LATEST                                            initial_from_file_var        = INVALID                                         [../]                                                                                                                                                             [./v]                                                                              block                        = INVALID                                           eigen                        = 0                                                 family                       = LAGRANGE                                          inactive                     = (no_default)                                      initial_condition            = INVALID                                           order                        = FIRST                                             outputs                      = INVALID                                           scaling                      = 1                                                 initial_from_file_timestep   = LATEST                                            initial_from_file_var        = INVALID                                         [../]                                                                          []                                                                                                                                                                                                                                                                                                                                                                                        ?�              ?�.���.?�.���-        ?�������?�������?�.���.        ?�������                                                                                ?�.���.?�.���-        ?�������?�������?�.���.        ?���������q�q���q�q���q�q���q�q���q�q���q�q���q�q���q�q���q�q�