CDF      
      
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes      	   num_elem      
num_el_blk        num_node_sets         num_side_sets         num_el_in_blk1        num_nod_per_el1       num_side_ss1      
num_df_ss1        num_side_ss2      
num_df_ss2        num_nod_ns1       num_nod_ns2       num_nod_var       num_info   �         api_version       @��H   version       @��H   floating_point_word_size            	file_size               title         out.e      maximum_name_length                    
time_whole                        2H   	eb_status                         4   eb_prop1               name      ID          8   	ns_status         	                <   ns_prop1      	         name      ID          D   	ss_status         
                L   ss_prop1      
         name      ID          T   coordx                      H  \   coordy                      H  �   eb_names                       $  �   ns_names      	                 D     ss_names      
                 D  T   
coor_names                         D  �   connect1                  	elem_type         QUAD4         @  �   elem_num_map                      	   elem_ss1                      	,   side_ss1                      	4   dist_fact_ss1                          	<   elem_ss2                      	\   side_ss2                      	d   dist_fact_ss2                          	l   node_ns1                      	�   node_ns2                      	�   vals_nod_var1                          H  2P   name_nod_var                       $  	�   info_records                      (�  	�                                      ?�      ?�              ?�      ?�      ?�              ?�                      ?�      ?�              ?�      ?�      ?�      ?�                                                                                                                                                                                                                                                                                                   	                                                                                                                         	u                                   ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                   -i                                                                               2d_diffusion_test.i                                                                                                                                               ### Input File ###                                                                                                                                                                                                                                 [Mesh]                                                                             second_order                   = 0                                               uniform_refine                 = 0                                               file                           = square.e                                        nemesis                        = 0                                               patch_size                     = 40                                              skip_partitioning              = 0                                             []                                                                                                                                                                [Variables]                                                                        [./u]                                                                              family                       = LAGRANGE                                          initial_condition            = 0                                                 order                        = FIRST                                             scaling                      = 1                                                 initial_from_file_timestep   = 2                                               [../]                                                                          []                                                                                                                                                                [Kernels]                                                                          [./diff]                                                                           type                         = Diffusion                                         execute_on                   = residual                                          start_time                   = -1.79769e+308                                     stop_time                    = 1.79769e+308                                      variable                     = u                                               [../]                                                                          []                                                                                                                                                                [BCs]                                                                              [./left]                                                                           type                         = DirichletBC                                       boundary                     = 1                                                 execute_on                   = residual                                          value                        = 0                                                 variable                     = u                                               [../]                                                                                                                                                             [./right]                                                                          type                         = DirichletBC                                       boundary                     = 2                                                 execute_on                   = residual                                          value                        = 1                                                 variable                     = u                                               [../]                                                                          []                                                                                                                                                                [Materials]                                                                        [./empty]                                                                          type                         = EmptyMaterial                                     block                        = 1                                                 execute_on                   = residual                                        [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      l_abs_step_tol                 = -1                                              l_max_its                      = 10000                                           l_tol                          = 1e-05                                           nl_abs_step_tol                = 1e-50                                           nl_abs_tol                     = 1e-50                                           nl_max_funcs                   = 10000                                           nl_max_its                     = 50                                              nl_rel_step_tol                = 1e-50                                           nl_rel_tol                     = 1e-08                                           no_fe_reinit                   = 0                                               petsc_options                  = -snes_mf_operator                               scheme                         = backward-euler                                  type                           = Steady                                          execute_on                     = residual                                      []                                                                                                                                                                [Output]                                                                           exodus                         = 1                                               file_base                      = out                                             gmv                            = 0                                               gnuplot_format                 = ps                                              interval                       = 1                                               iteration_plot_start_time      = 1.79769e+308                                    max_pps_rows_screen            = 0                                               nemesis                        = 0                                               output_displaced               = 0                                               output_initial                 = 1                                               perf_log                       = 1                                               postprocessor_csv              = 0                                               postprocessor_ensight          = 0                                               postprocessor_gnuplot          = 0                                               postprocessor_screen           = 1                                               print_linear_residuals         = 1                                               screen_interval                = 1                                               show_setup_log_early           = 0                                               tecplot                        = 0                                               tecplot_binary                 = 0                                               xda                            = 0                                             []                                                                                                                                                                [setup_quadrature]                                                                 order                          = AUTO                                            type                           = GAUSS                                         []                                                                                                                                                                [setup_dampers]                                                                  []                                                                                                                                                                [no_action]                                                                      []                                                                                                                                                                [init_problem]                                                                   []                                                                                                                                                                [copy_nodal_vars]                                                                  initial_from_file_timestep     = 2                                             []                                                                                                                                                                [check_integrity]                                                                []                                                                                                                                                                [no_action]                                                                                                                                                      ?�              ?�   7��?�   7��        ?�   7��?�   7��?�   7��        ?�   7��