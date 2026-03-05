#==============================================================================#
# Config:       Default phyddle config file                                    #
# Authors:      Michael Landis and Ammon Thompson                              #
# Date:         230804                                                         #
# Description:  Simple birth-death and equal-rates CTMC model in R using ape   #
#==============================================================================#


args = {

    #-------------------------------#
    # Project organization          #
    #-------------------------------#
    'dir'     : './',

    #-------------------------------#
    # Simulate Step settings        #
    #-------------------------------#
    'sim_command'        : 'bash SIRM_sim.sh',   # Simulation command to run single job (see documentation)

    'start_idx'          : 1,                         # Start index for simulated training replicates
    #'end_idx'            : 500,                      # End index for simulated training replicates
    'end_idx'            : 500000,                      # End index for simulated training replicates
    'sim_batch_size'     : 1000,                        # Number of replicates per simulation command
    #'sim_batch_size'     : 10,                        # Number of replicates per simulation command

#    'num_proc'           : 10,                   # Number of cores for multiprocessing (-N for all but N)
    'num_proc'           : 100,                   # Number of cores for multiprocessing (-N for all but N)
    #'use_cuda'           : 'F',
    #-------------------------------#
    # Format Step settings          #
    #-------------------------------#
    'num_char'          : 5,                # number of evolutionary characters
    'num_states'        : 2,                # number of states per character
    'tree_width'        : 50,
    'tree_encode'        : 'serial',             # Encoding strategy for tree     
    'brlen_encode'      : 'height_brlen',   # how to encode phylo brlen? height_only or height_brlen
    'char_encode'       : 'integer',        # how to encode discrete states? one_hot or integer
    'asr_est'            : 'T',
    'char_format'       : 'csv',
    'min_num_taxa'       : 50,                   # Minimum number of taxa allowed when formatting
    'max_num_taxa'       : 50,                  # Maximum number of taxa allowed when formatting
    'shuffle_test'       : 'F', 

    'char_format'        : 'nexus',              # File format for character data 

    'param_data'         : {"time_of_peak_prev_0" : 'num', 
                            "time_of_peak_prev_1" : 'num', 
                            "time_of_peak_prev_2" : 'num', 
                            "time_of_peak_prev_3" : 'num', 
                            "time_of_peak_prev_4" : 'num' 
                           },

    'prop_test' : .05, 
    'prop_val' : .05, 
    'num_epochs' : 500, 
    'num_early_stop' : 25, 
    #'weight_decay' : 0.001,

    #-------------------------------#
    # Train                         #
    #-------------------------------#
    #'trn_batch_size'     : 100,                 # Training batch sizes

    #'phy_kernel_stride'  : [7, 8],               # Kernel sizes for stride convolutional layers
                                                 #     for phylogenetic state input
#    'use_parallel'       : 'F',            # Use parallelization? (recommended)
    'map_tip_states'    : {0 : [1, 0, 0, 0, 0], 
                           1 : [0, 1, 0, 0, 0], 
                           2 : [0, 0, 1, 0, 0],
                           3 : [0, 0, 0, 1, 0],
                           4 : [0, 0, 0, 0, 1]},
    'rb_nexus': 'T',
}
