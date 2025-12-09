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

    #'use_parallel'       : 'F',            # Use parallelization? (recommended)
    #-------------------------------#
    # Simulate Step settings        #
    #-------------------------------#
    'sim_command'       : 'Rscript ../sim_varyTree.R 4', # exact command string, argument is output file prefix

    'start_idx'          : 1,                         # Start index for simulated training replicates
    'end_idx'            : 150000,                      # End index for simulated training replicates
    'sim_batch_size'     : 500,                        # Number of replicates per simulation command

    #-------------------------------#
    # Format Step settings          #
    #-------------------------------#
    'num_char'          : 1,                # number of evolutionary characters
    'num_states'        : 2,                # number of states per character
    'tree_width'        : 32,
    'tree_encode'       : 'extant',         # use model with serial or extant tree
    'brlen_encode'      : 'height_brlen',   # how to encode phylo brlen? height_only or height_brlen
    'char_encode'       : 'integer',        # how to encode discrete states? one_hot or integer
    'param_est'         : {                # model parameters to predict (labels)
                            'asr_node_state'     : 'cat'

    # If you change the values in format for estimate, the true values change but not the inferred values in terms of what variables are printed out
                          },
    'param_data'        : {                      # Known model parameters to treat as aux. data
        'asr_node_label'        : 'num'
    },
    'asr_one'            : 'T',

    'char_format'       : 'csv',
    'min_num_taxa'       : 4,                  # Minimum number of taxa allowed when formatting
    'max_num_taxa'       : 4,                  # Maximum number of taxa allowed when formatting
    'shuffle_test'       : 'F', 

    #-------------------------------#
    # Train                         #
    #-------------------------------#
    #'trn_batch_size'     : 4000,                 # Training batch sizes

    #'phy_kernel_stride'  : [7, 8],              # Kernel sizes for stride convolutional layers
                                                 #     for phylogenetic state input
    #'use_parallel'       : 'F',                 # Use parallelization? (recommended)
}
