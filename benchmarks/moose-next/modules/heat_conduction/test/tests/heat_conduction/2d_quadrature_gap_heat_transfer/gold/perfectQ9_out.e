CDF      
      
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes      2   num_elem      
num_el_blk        num_node_sets         num_side_sets         num_el_in_blk1        num_nod_per_el1    	   num_el_in_blk2        num_nod_per_el2    	   num_side_ss1      num_side_ss2      num_side_ss3      num_side_ss4      num_nod_ns1       num_nod_ns2       num_nod_ns3       num_nod_ns4       num_nod_var       num_elem_var      num_info  	         api_version       @�
=   version       @�
=   floating_point_word_size            	file_size               int64_status             title         perfectQ9_out.e    maximum_name_length                 &   
time_whole                            f�   	eb_status                             
�   eb_prop1               name      ID              
�   	ns_status         	                    
�   ns_prop1      	         name      ID              
�   	ss_status         
                    
�   ss_prop1      
         name      ID              
�   coordx                     �      
�   coordy                     �      �   eb_names                       D         ns_names      	                 �      \   ss_names      
                 �      �   
coor_names                         D      d   node_num_map                    �      �   connect1                  	elem_type         QUAD9         �      p   connect2                  	elem_type         QUAD9         �          elem_num_map                           �   elem_ss1                          �   side_ss1                          �   elem_ss2                          �   side_ss2                          �   elem_ss3                          �   side_ss3                          �   elem_ss4                          �   side_ss4                          �   node_ns1                          �   node_ns2                             node_ns3                             node_ns4                          ,   vals_nod_var1                         �      f�   name_nod_var                       $      @   name_elem_var                          D      d   vals_elem_var1eb1                                 h,   vals_elem_var2eb1                                 hL   vals_elem_var1eb2                                 hl   vals_elem_var2eb2                                 h�   elem_var_tab                             �   info_records                      S�      �                                                            ?�      ?�                      ?�      ?�              ?�      ?�      ?�              ?�      ?�              ?�      ��      ��      ��      ��      ��      ��      ��      ��      ��      ��      @      @      @       @       @      @      @       @      @      @      @       @      @      @       @      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ��                      ��      ��              ��      ��      ��      ?�      ?�      ?�      ?�      ?�      ?�              ��              ��      ��      ��      ?�      ?�      ?�      ?�      ��                      ��      ��              ��      ��      ��      ?�      ?�      ?�      ?�      ?�      ?�              ��              ��      ��      ��      ?�      ?�      ?�      ?�      left                             right                                                                                                                                                                  leftright                        leftleft                         rightright                       rightleft                                                                                                                       	   
                                                                      !   "   #   $   %   &   '   (   )   *   +   ,   -   .   /   0   1   2                           	      
                                                                                                    !   "      #   $      %   &   '      (         )   *       +   ,   -   .      $   /   )   '   0   1   +   2                                                                                       
            #   %                  )   *   ,   /   1temp                                qpoint_penetration               paired_temp                                    ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                   -i              ��      ��                                                       perfectQ9.i                                                                      --no-gdb-backtrace                                                                                                                                                ### Version Info ###                                                             Framework Information:                                                           MOOSE version:           git commit 0dc7d2b on 2017-02-02                        PETSc Version:           3.7.4                                                   Current Time:            Thu Feb  2 09:34:02 2017                                Executable Timestamp:    Thu Feb  2 09:23:38 2017                                                                                                                                                                                                  ### Input File ###                                                                                                                                                []                                                                                 initial_from_file_timestep     = LATEST                                          initial_from_file_var          = INVALID                                         block                          = INVALID                                         coord_type                     = XYZ                                             fe_cache                       = 0                                               kernel_coverage_check          = 1                                               material_coverage_check        = 1                                               name                           = 'MOOSE Problem'                                 restart_file_base              = INVALID                                         rz_coord_axis                  = Y                                               type                           = FEProblem                                       use_legacy_uo_aux_computation  = INVALID                                         use_legacy_uo_initialization   = INVALID                                         control_tags                   = INVALID                                         enable                         = 1                                               error_on_jacobian_nonzero_reallocation = 0                                       force_restart                  = 0                                               near_null_space_dimension      = 0                                               null_space_dimension           = 0                                               petsc_inames                   =                                                 petsc_options                  = INVALID                                         petsc_values                   =                                                 solve                          = 1                                               transpose_null_space_dimension = 0                                               use_nonlinear                  = 1                                             []                                                                                                                                                                [BCs]                                                                                                                                                               [./left]                                                                           boundary                     = leftleft                                          control_tags                 = INVALID                                           enable                       = 1                                                 implicit                     = 1                                                 type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = temp                                              diag_save_in                 = INVALID                                           save_in                      = INVALID                                           seed                         = 0                                                 value                        = 300                                             [../]                                                                                                                                                             [./right]                                                                          boundary                     = rightright                                        control_tags                 = INVALID                                           enable                       = 1                                                 implicit                     = 1                                                 type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = temp                                              diag_save_in                 = INVALID                                           save_in                      = INVALID                                           seed                         = 0                                                 value                        = 400                                             [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      type                           = Steady                                          compute_initial_residual_before_preset_bcs = 0                                   control_tags                   =                                                 enable                         = 1                                               l_abs_step_tol                 = -1                                              l_max_its                      = 10000                                           l_tol                          = 1e-05                                           line_search                    = default                                         nl_abs_step_tol                = 1e-50                                           nl_abs_tol                     = 1e-50                                           nl_max_funcs                   = 10000                                           nl_max_its                     = 50                                              nl_rel_step_tol                = 1e-50                                           nl_rel_tol                     = 1e-08                                           no_fe_reinit                   = 0                                               petsc_options                  = INVALID                                         petsc_options_iname            = INVALID                                         petsc_options_value            = INVALID                                         restart_file_base              =                                                 solve_type                     = PJFNK                                           splitting                      = INVALID                                                                                                                          [./Quadrature]                                                                     element_order                = AUTO                                              order                        = THIRD                                             side_order                   = AUTO                                              type                         = GAUSS                                           [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      _fe_problem                    = 0x7fbce2040000                                  _fe_problem_base               = 0x7fbce2040000                                                                                                                   [./Quadrature]                                                                   [../]                                                                          []                                                                                                                                                                [GlobalParams]                                                                     order                          = FIRST                                         []                                                                                                                                                                [Kernels]                                                                                                                                                           [./hc]                                                                             type                         = HeatConduction                                    block                        = INVALID                                           control_tags                 = Kernels                                           diag_save_in                 = INVALID                                           diffusion_coefficient        = thermal_conductivity                              diffusion_coefficient_dT     = thermal_conductivity_dT                           eigen_kernel                 = 0                                                 enable                       = 1                                                 implicit                     = 1                                                 save_in                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 1                                                 variable                     = temp                                            [../]                                                                          []                                                                                                                                                                [Materials]                                                                                                                                                         [./hcm]                                                                            type                         = HeatConductionMaterial                            block                        = 'left right'                                      boundary                     = INVALID                                           compute                      = 1                                                 control_tags                 = Materials                                         enable                       = 1                                                 implicit                     = 1                                                 output_properties            = INVALID                                           outputs                      = none                                              seed                         = 0                                                 specific_heat                = 1                                                 specific_heat_temperature_function =                                             temp                         = INVALID                                           thermal_conductivity         = 1                                                 thermal_conductivity_temperature_function =                                      use_displaced_mesh           = 0                                               [../]                                                                          []                                                                                                                                                                [Mesh]                                                                             displacements                  = INVALID                                         block_id                       = INVALID                                         block_name                     = INVALID                                         boundary_id                    = INVALID                                         boundary_name                  = INVALID                                         construct_side_list_from_node_list = 0                                           ghosted_boundaries             = INVALID                                         ghosted_boundaries_inflation   = INVALID                                         patch_size                     = 40                                              second_order                   = 0                                               skip_partitioning              = 0                                               type                           = FileMesh                                        uniform_refine                 = 0                                               centroid_partitioner_direction = INVALID                                         construct_node_list_from_side_list = 1                                           control_tags                   =                                                 dim                            = 1                                               distribution                   = DEFAULT                                         enable                         = 1                                               file                           = perfectQ9.e                                     ghost_point_neighbors          = 0                                               nemesis                        = 0                                               num_ghosted_layers             = 1                                               parallel_type                  = DEFAULT                                         partitioner                    = default                                         patch_update_strategy          = never                                         []                                                                                                                                                                [Mesh]                                                                           []                                                                                                                                                                [Outputs]                                                                          append_date                    = 0                                               append_date_format             = INVALID                                         checkpoint                     = 0                                               color                          = 1                                               console                        = 1                                               controls                       = 0                                               csv                            = 0                                               dofmap                         = 0                                               execute_on                     = 'INITIAL TIMESTEP_END'                          exodus                         = 1                                               file_base                      = INVALID                                         gmv                            = 0                                               gnuplot                        = 0                                               hide                           = INVALID                                         interval                       = 1                                               nemesis                        = 0                                               output_if_base_contains        = INVALID                                         print_linear_residuals         = 1                                               print_mesh_changed_info        = 0                                               print_perf_log                 = 0                                               show                           = INVALID                                         solution_history               = 0                                               sync_times                     =                                                 tecplot                        = 0                                               vtk                            = 0                                               xda                            = 0                                               xdr                            = 0                                             []                                                                                                                                                                [ThermalContact]                                                                                                                                                    [./left_to_right]                                                                  conductivity_master_name     = thermal_conductivity                              conductivity_name            = thermal_conductivity                              gap_aux_type                 = GapValueAux                                       master                       = rightleft                                         normal_smoothing_distance    = INVALID                                           normal_smoothing_method      = INVALID                                           order                        = SECOND                                            quadrature                   = 1                                                 slave                        = leftright                                         tangential_tolerance         = INVALID                                           type                         = GapHeatTransfer                                   variable                     = temp                                              warnings                     = 0                                                 appended_property_name       =                                                   disp_x                       = INVALID                                           disp_y                       = INVALID                                           disp_z                       = INVALID                                           save_in                      = INVALID                                           contact_pressure             = INVALID                                           gap_conductivity             = 1                                                 gap_conductivity_function    = INVALID                                           gap_conductivity_function_variable = INVALID                                   [../]                                                                          []                                                                                                                                                                [Variables]                                                                                                                                                         [./temp]                                                                           block                        = INVALID                                           eigen                        = 0                                                 family                       = LAGRANGE                                          initial_condition            = INVALID                                           order                        = SECOND                                            outputs                      = INVALID                                           scaling                      = 1                                                 initial_from_file_timestep   = LATEST                                            initial_from_file_var        = INVALID                                         [../]                                                                          []                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ��      ��                                                                                                                      ?�      @t�UUR@&@t�UUU8�@sʪ���,@sʪ��О@t�UUUͼ@tP  =@sʪ��E�@tP   M�@tP   �3@t�UUX��@sʪ��O@t�UU[�?@tP  C�@sʪ��B2@tP  ��@r�   &?@r�   &?@sEUUW�R@r�   &?@sEUUW{�@sEUUW@r�   &?@sEUUX�M@r�   &?@sEUUW�#@x�����P@x�����P@w�UU[b@w�UUY��@x�����P@xz���$?@w�UUUp	@xz�����@xz����@x�����P@w�UU]��@x�����P@xz����@w�UUYd^@xz���h�@vꪪ���@vꪪ�9Q@wp  �@vꪪ��_@wo���Q�@wp  �T@vꪪ��@wp  �"@vꪪ�$@wp  ]ҿ�      ��                      @vꪪ�c�@vꪪ�}�                                                                @t�UUUX@t�UUZ#$