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

    #    'use_parallel'       : 'F',            # Use parallelization? (recommended)
    'num_proc'           : 50,             # Number of cores for multiprocessing
    #-------------------------------#
    # Simulate Step settings        #
    #-------------------------------#
    'sim_command'       : 'Rscript sim_bisse.R' , # exact command string, argument is output file prefix

    'start_idx'          : 1,                         # Start index for simulated training replicates
    #'start_idx'          : 2501,                         # Start index for simulated training replicates
    #'end_idx'            : 3000,                      # End index for simulated training replicates
    'end_idx'            : 50000,                      # End index for simulated training replicates
    #'sim_batch_size'     : 50,                        # Number of replicates per simulation command
    'sim_batch_size'     : 1000,                        # Number of replicates per simulation command

    #-------------------------------#
    # Format Step settings          #
    #-------------------------------#
    'num_char'          : 1,                # number of evolutionary characters
    'num_states'        : 2,                # number of states per character
    'tree_width'        : 50,
    'tree_encode'       : 'extant',         # use model with serial or extant tree
    'brlen_encode'      : 'height_brlen',   # how to encode phylo brlen? height_only or height_brlen
    'char_encode'       : 'integer',        # how to encode discrete states? one_hot or integer

    ## If you change the values in format for estimate, the true values change but not the inferred values in terms of what variables are printed out
    #                      },
    'asr_est'            : 'T',

    'char_format'       : 'csv',
    'min_num_taxa'       : 50,                  # Minimum number of taxa allowed when formatting
    'max_num_taxa'       : 50,                  # Maximum number of taxa allowed when formatting
    'shuffle_test'       : 'F', 

}
