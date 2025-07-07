# This script will simulate a batch of phylogenetic trees under a 
# density-dependent GeoSSE model. Tree label starts from start_idx to (start_idx + batch_size - 1)

# have validated script for lemma 3 (sim7)
# fixed bug where num_alive_lineage is not update when saved to a file
# add randomization on only_density, model type, and start range to generate a tree 
# remove "sp." from tree label to match with taxon label on character data
# add header on the character data
# changed categorical label using integers (i.e. for model types, dd types etc)
# save rate scalars in log scale
# Convert character data at tips (multiple binary characters -> single multi-state characters)

#####################
# Load julia packages
#####################
using Combinatorics # Package for combinatorics/permutations
using OrderedCollections # Package for ordering dictionary 
using DataStructures # Package for preserving dictionary order by its insertion order 
using Distributions # Package for using probability distributions 
using DataFrames # Package for dataframe
using CSV # Package to write to csv

import Random

####################
# Command line arguments 
####################
args          = ARGS    
num_tax       = parse(Int,args[1])   
out_path      = map(x->string(x), args)[2]
out_prefix    = map(x->string(x), args)[3]
start_idx     = parse(Int,args[4]) # start tree id
batch_size    = parse(Int,args[5]) # batch size 

#############
# File system
#############
out_dir  = out_path
################
# Initial setting
################
min_taxa       = num_tax; # Minimum number of extant taxa to simulate 
max_taxa       = num_tax; #rand(Uniform(100,150)); # Maximum number of extant taxa to simulate
max_time       = 100 #rand(Uniform(70,100)); # Maximum simulation time 
# stop_condition = "time" # If time is the stopping condition
#:stop_condition = "numtaxa" # If the final number of extant taxa is the stopping condition
stop_condition = "numtaxa" # If the final number of extant taxa and time are the stopping conditions
global resimulate = true # whether to resimulate in case of tree extinction or not
global timebin    = "100" # time bin for tracking number of species across region/range in continuous time 

#############
# Model setup
#############
# Define state space 
# num_regions = 2; # Number of regions 
num_regions = 3; # Number of regions 
num_ranges  = 2^num_regions - 1; # Number of possible ranges given num_regions

#################################################
# Create dictionary for region occupancy status (present/absent) and range labels
#################################################
# e.g. (1,0,0) means a species present only in region A, so it has range {A}.
# ranges_inv = {
#     (1,0,0):1,
#     (0,1,0):2,
#     (0,0,1):3,
#     (1,1,0):4,
#     (1,0,1):5,
#     (0,1,1):6,
#     (1,1,1):7
# }
dict_values = zeros(Int64,0); # Intialize empty vector for the values field of the dictionary
for i in 1:num_ranges
    push!(dict_values,i); # Populate the values field 
end

dict_keys = Vector{Int}[]; # Create an empty vector of vectors for the keys field of the dictionary
for j in 1:num_regions
  zero_vec = zeros(Int64,num_regions); # Create a zero vector of length num_regions
  zero_vec[1:j] = fill(1,j); # Populating the vector with 1's for indicating presence/absence of species
  range_comb = unique(permutations(zero_vec)); # different permutations of (1,...) for range size <= num_regions 
  for k in 1:length(range_comb)
    push!(dict_keys,range_comb[k]); # populating keys field of the dictionary
  end
end

ranges_inv = Dict(dict_keys .=> dict_values) # Append the keys and values to the dictionary 
ranges_inv = sort(ranges_inv; byvalue=true) # sort the dictionary by its values field
# ranges = [ [1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], [1,1,0], [1,1,1] ]


# Define event space as dictionaries (for checking if we have a valid cladogenesis/anagenesis during simulation)
################################################
# Dictionary for within-region speciation (evt_w)
################################################
# key is ancestral range, value is list of possible child ranges that differ from ancestral range
# for example:
#     (1,1,0):[(1,0,0),(0,1,0)],
# valid events:
#     AB->AB,A (budding cladogenesis)
#     AB->A,AB (budding cladogenesis)
#     AB->AB,B (budding cladogenesis)
#     AB->B,AB (budding cladogenesis)
# evt_w = {
#     (1,0,0):[(1,0,0)],
#     (0,1,0):[(0,1,0)],
#     (0,0,1):[(0,0,1)],
#     (1,1,0):[(1,0,0),(0,1,0)],
#     (1,0,1):[(1,0,0),(0,0,1)],
#     (0,1,1):[(0,1,0),(0,0,1)],
#     (1,1,1):[(1,0,0),(0,1,0),(0,0,1)]
# }

evt_w       = OrderedDict(); # Create an empty dictionary for within-region speciation 
dict_keys_w = copy(dict_keys); # keys field of the dictionary evt_w

for k in 1:length(dict_keys_w)
    if count(==(1),dict_keys_w[k])==1 # if the ancestral range is of size 1 (endemic)
        evt_w[dict_keys_w[k]] = dict_keys_w[k]; # assign possible range for children = ancestral range 
    elseif count(==(1),dict_keys_w[k])>1 # if the ancestral range is of size > 1 (widespread)
        idx = findall(dict_keys_w[k].== 1 ); # find locations of 1's in ancestral species range k
        child_range_w = Vector{Int}[]; # initialize the vector of vectors to store all possible child ranges from ancestral of range k
        for l in 1:length(idx) # loop over each 1 in ancestral species range k for within-region speciation in that location
            child_range_idx = zeros(Int64,num_regions)
            child_range_idx[idx[l]] = 1 # indicates that speciation happens in location idx[l] from the parent species with range k
            push!(child_range_w,child_range_idx) # append to possible range for child species from ancenstral species with range k
        end
        evt_w[dict_keys_w[k]] = child_range_w; # create dictionary for evt_w
    end
end

################################################
# Dictionary for extinction (evt_e)
################################################
# key is range before extinction event, value is all possible ranges following an extinction. 
# [(0,0,0)] indicates species extinction
#
# evt_e = {
#     (1,0,0):[(0,0,0)],
#     (0,1,0):[(0,0,0)],
#     (0,0,1):[(0,0,0)],
#     (1,1,0):[(1,0,0),(0,1,0)],
#     (1,0,1):[(1,0,0),(0,0,1)],
#     (0,1,1):[(0,1,0),(0,0,1)],
#     (1,1,1):[(1,1,0),(1,0,1),(0,1,1)]
# }
evt_e       = OrderedDict(); # Create an empty dictionary for extinction
dict_keys_e = copy(dict_keys); # keys field of the dictionary evt_e

for k in 1:length(dict_keys_e)
    if count(==(1),dict_keys_e[k])==1 # if the range is of size 1 (endemic)
        evt_e[dict_keys_e[k]] = zeros(Int64,num_regions); # species extinction 
    elseif count(==(1),dict_keys_e[k])>1 # if the ancestral range is of size > 1 (widespread)
        idx = findall(dict_keys_e[k].== 1 ); # find locations of 1's in species range k before extinction
        new_range_e = Vector{Int}[]; # initialize the vector of vectors to store all possible new ranges after extinction
        for l in 1:length(idx) # loop each location index with a species in range k
            new_range_idx = copy(dict_keys_e[k])
            new_range_idx[idx[l]] = 0 # indicates that extinction occurs in location idx[l] in species with range k
            push!(new_range_e,new_range_idx) # append to possible range new range after extinction for species with range k
        end
        evt_e[dict_keys_e[k]] = new_range_e; # create dictionary for evt_e
    end
end

################################################
# Dictionary for dispersal (evt_d)
################################################
# # key is range before dispersal, value indicates all possible ranges following a dispersal event 
# evt_d = {
#     (1,0,0):[(1,1,0),(1,0,1)],
#     (0,1,0):[(1,1,0),(0,1,1)],
#     (0,0,1):[(0,1,1),(0,1,1)],
#     (1,1,0):[(1,1,1)],
#     (1,0,1):[(1,1,1)],
#     (0,1,1):[(1,1,1)],
#     (1,1,1):[(,)]
# }

evt_d       = OrderedDict(); # Create an empty dictionary for dispersal
dict_keys_d = copy(dict_keys); # keys field of the dictionary evt_d

for k in 1:length(dict_keys_d)
    if count(==(1),dict_keys_d[k])==num_regions # if the species is already widespread in all regions 
        evt_d[dict_keys_d[k]] = [] # can't disperse
    elseif count(==(1),dict_keys_d[k])!=num_regions # if the range is still < num_regions, can disperse
        idx = findall(dict_keys_d[k].== 0 ); # find locations of all 0's in species with range k before dispersal
        new_range_d = Vector{Int}[]; # intialize vector of vectors to store all possible new ranges following a dispersal of species with range k
        for l in 1:length(idx) # loop each location index with a species in range k
            new_range_idx = copy(dict_keys_d[k])
            new_range_idx[idx[l]] = 1 # indicates that species dispersal to that location idx[l]
            push!(new_range_d,new_range_idx) # append to possible range new range after extinction for species with range k
        end
        evt_d[dict_keys_d[k]] = new_range_d; # create dictionary for evt_d
    end
end

################################################
# Dictionary for between-region speciation (evt_b)
################################################
# key indicate ancestral species range, value indicates all possible child ranges for one of the two children 
# following the speciation event. List of all possible ranges of the other child can be computed using:
# ancestral range - child range. E.g. ancestral = (1,1,0), child = (0,1,0). Then the other child is: 
# (1,1,0) - (0,1,0) = (1,0,0)

# evt_b = {
#     (1,0,0):[(,)],
#     (0,1,0):[(,)],
#     (0,0,1):[(,)],
#     (1,1,0):[(1,0,0)],
#     (1,0,1):[(1,0,0)],
#     (0,1,1):[(0,1,0)],
#     (1,1,1):[(1,0,0),(0,1,0),(0,0,1)]
# }

evt_b       = OrderedDict(); # Create an empty dictionary for between-region speciation
dict_keys_b = copy(dict_keys); # keys field of the dictionary evt_b

for k in 1:length(dict_keys_b) # Looping over each ancestral species range in keys field 
    if count(==(1),dict_keys_b[k])==1 # if the ancestral species is endemic, can't speciate 
        evt_b[dict_keys_b[k]] = [] # then assign an empty vector for child range
    else # if the ancestral species is widespread, can speciate
        anc_range_size = count(==(1),dict_keys_b[k]) # ancestral range size for index k in key field
        # looping for each different size of child range after speciation.
        # Divide by 2 because the range size of the other child can be derived from the range of the other one. 
        idx_anc = findall(dict_keys_b[k].== 1 ); # find locations of 1's in ancestral range k
        #
        child_range_b = Vector{Int}[]; # initialize vector of vectors to store all child ranges from ancestral range k 
        for l in 1:floor.(Int64,anc_range_size/2) # loop over each possible range split size for the child range e.g. for k = 4, we can have split of size 1&2 (not including its compelement of size 3)
                child_init = zeros(Int64,num_regions); # Initialize vector of zeros 
                child_init[1:l] .= 1; # Assign 1's equal to the range split size 
                child_comb = unique(permutations(child_init)); # permute all possibilities of combinations of 1's for that split size 
                if l != anc_range_size - l # case: not equal split between left and right child range size e.g. l|m where l!=m
                    for m in 1:length(child_comb) # loop over this permutation 
                        idx_child = findall(child_comb[m].== 1 );
                        test_idx = [x in idx_anc for x = idx_child] # Check if that split of size l is valid (i.e. if its contained in ancestral range)
                        if all(test_idx) == true 
                            push!(child_range_b,child_comb[m]) # If yes, then we accept that split of size l as a valid split
                        end
                    end
                else # case: equal split between left and right child range size e.g. l|m where l = m 
                    track_range = Vector{Int}[] # initialize empty vector of vectors 
                        for n in 1:length(child_comb) # loop over this permutation 
                            idx_child = findall(child_comb[n].== 1 );
                            test_idx = [x in idx_anc for x = idx_child] # Check if that split of size l is valid (i.e. if its contained in ancestral range)
                            if all(test_idx) == true 
                                push!(child_range_b,child_comb[n]) # If yes, then we accept that split of size l as a valid split
                                push!(track_range,child_comb[n]) # If yes, add that to track_range
                                break
                            end
                        end
                    child_comb_reduce = setdiff(child_comb,track_range) # remove the newly added child range from the pool in child_comb
                    for p in 1:length(child_comb_reduce) # loop over the remaining
                        idx_child = findall(child_comb_reduce[p].== 1 );
                        test_idx = [x in idx_anc for x = idx_child] # Check if that split of size l is valid (i.e. if its contained in ancestral range)
                        if (all(test_idx) && child_comb_reduce[p] ∉ [dict_keys_b[k]].-track_range) == true
                            push!(child_range_b,child_comb_reduce[p])
                            push!(track_range,child_comb_reduce[p])
                        end
                    end
                end
        end
        evt_b[dict_keys_b[k]] = child_range_b; # append keys and value fields to dictionary 
    end
end 


###########################################
# # Define rates for event space in each model variant 
###########################################

# log-normal distribution (mean = m, std = std) in Julia 
# p.s. use lognormal to ensure rate is always positive
function myLogNormal(m,std)
    γ = 1+std^2/m^2
    μ = log(m/sqrt(γ))
    σ = sqrt(log(γ))
    return LogNormal(μ,σ)
end

function redraw_dist()
    function flip_coin(min_val,max_val)
        #if rand(Bool)
            val = rand(Uniform(min_val, max_val))
        #else
        #    val = rand(LogUniform(min_val, max_val))
        #end
        return val
    end
    #
    w_d      = 0 #-rand(Uniform(0.0001,1.0));
    #e_d      = rand(Uniform(0.0001,2.0));
    # Remove density dependent extinction
    e_d      = 0;
    d_d_home = 0;
    #d_d_away = -rand(Uniform(0.0001,2.0));
    d_d_away = 0; #-rand(Uniform(0.0001,2.0));
    b_d      = 0; 
    #
    rho_w = flip_coin(0.01,1);
    rho_d = flip_coin(0.0001,1);
    rho_b = flip_coin(0.0001,1);
    rho_e = flip_coin(0.0001,rho_w);
    ans = [rho_w,rho_e,rho_d,rho_b,w_d,e_d,d_d_home,d_d_away,b_d]
    return(ans)
end

#########################
# Rate function for within-region speciation 
#########################

# note: use lognormal distribution with mean 1 (for no effect from non-density)
# use flag (true or false) for using density-independent or not (rate == 1)

function rate_w(base,region_ind,species_count_ind,rate_scalar,only_density,indep_model,model)
    ## Error checking
    @assert base >= 0 "base rate must be positive" # base rate must be non negative
    @assert region_ind <= num_regions && region_ind > 0 "region index must be in {1,2,...,num_regions}" # region ind must be <= num_regions
    n_i   = species_count_ind # Number of species currently in region_ind
    if n_i == 0
        ans = 0 # rate equals 0 if no species 
    else
        ## Define rate parameters 
        rho_w = base # base rate 
        w_d   = rate_scalar # rate scalar for within-region speciation
        if only_density == 1 
            m_i_w = 1 # means no density-independent factors 
        elseif only_density == 0
            if indep_model == 1
                    m_i_w = rand(myLogNormal(1,1)) # draw density-independent rate from lognormal with mean 1, std 1 (same distribution across all regions) 
            elseif indep_model == 0 && num_regions == 2 # standard geosse with 2 regions
                if region_ind == 1
                    m_i_w = 1.0 
                    # m_i_w = 0.8 
                elseif region_ind == 2
                    m_i_w = 0.5
                    # m_i_w = 0.3
                end
            elseif indep_model == 0 && num_regions == 3 # standard geosse with 3 regions
                if region_ind == 1
                    m_i_w = 1.0 
                elseif region_ind == 2
                    m_i_w = 0.5
                elseif region_ind == 3
                    m_i_w = 0.5
                end
            end
        end 
        ## Compute the overall rates
        if model == 0
            m_d_w = (log(n_i)+1)^(w_d)
        elseif model == 1
            m_d_w = n_i^w_d
        end
        #
        if only_density == 1 || (only_density == 0 && indep_model == 1)
            ans   = rho_w*m_i_w*m_d_w
        elseif only_density == 0 && indep_model == 0
            ans   = rho_w*m_i_w
        end
    end
    return(ans)
end

#########################
# Rate function for extinction
#########################

function rate_e(base,region_ind,species_count_ind,rate_scalar,only_density,indep_model,model)
    ## Error checking
    @assert base >= 0 "base rate must be positive" # base rate must be non negative
    @assert region_ind <= num_regions && region_ind > 0 "region index must be in {1,2,...,num_regions}" # region ind must be <= num_regions
    n_i   = species_count_ind # Number of species currently in region_ind
    if n_i == 0
        ans = 0 # rate equals 0 if no species 
    else # if there is species 
        ## Define rate parameters 
        rho_e = base # base rate 
        e_d   = rate_scalar # rate scalar for extinction
        if only_density == 1 
            m_i_e = 1 # means no density-independent factors 
        elseif only_density == 0
            if indep_model == 1
                    m_i_e = rand(myLogNormal(1,1)) # draw density-independent rate from lognormal with mean 1, std 1 (same distribution across all regions)
            elseif indep_model == 0 && num_regions == 2 # standard geosse with 2 regions
                if region_ind == 1
                    m_i_e = 0.2 
                    # m_i_e = 0.5
                elseif region_ind == 2
                    m_i_e = 0.4
                    # m_i_e = 0.2
                end
            elseif indep_model == 0 && num_regions == 3 # standard geosse with 3 regions 
                if region_ind == 1
                    m_i_e = 0.2 
                elseif region_ind == 2
                    m_i_e = 0.4
                elseif region_ind == 3
                    m_i_e = 0.5
                end
            end
        end 
        ## Compute the overall rates
        if model == 0
            m_d_e = (log(n_i)+1)^(e_d)
        elseif model == 1
            m_d_e = n_i^e_d
        end
        #
        if only_density == 1 || (only_density == 0 && indep_model == 1)
            ans   = rho_e*m_i_e*m_d_e
        elseif only_density == 0 && indep_model == 0
            ans   = rho_e*m_i_e
        end
    end
    return(ans)
end

#########################
# Rate function for dispersal 
#########################

function rate_d(base,num_species_per_region,home_ind,away_ind,home_scalar,away_scalar,only_density,indep_model,model)
    ## Error checking
    @assert base >= 0 "base rate must be positive" # base rate must be non negative
    @assert home_ind <= num_regions && home_ind > 0 "home region index must be in {1,2,...,num_regions}" # home_ind must be <= num_regions
    @assert away_ind <= num_regions && away_ind > 0 "away region index must be in {1,2,...,num_regions}" # away_ind must be <= num_regions
    @assert home_ind != away_ind "home region index cannot be equal to away region index"
    n_i   = num_species_per_region[home_ind]
    n_j   = num_species_per_region[away_ind]; @assert n_j >= 0 "number of species in away_ind must be >= 0"
    #
    if n_i == 0
        ans = 0 # rate equals 0 if no species in home region
    else
        ## Define rate parameters 
        rho_d      = base # base rate 
        d_d_home   = home_scalar # rate scalar in home region
        d_d_away   = away_scalar # rate scalar in away region
        if only_density == 1 
            m_i_d = 1 # means no density-independent factors 
        elseif only_density == 0
            if indep_model == 1
                    m_i_d = rand(myLogNormal(1,1)) # draw density-independent rate from lognormal with mean 1, std 1 
            elseif indep_model == 0 && num_regions == 2 # standard geosse with 2 regions
                if home_ind == 1 && away_ind == 2
                    m_i_d = 1.0
                elseif home_ind == 2 && away_ind == 1
                    m_i_d = 0.5
                end 
            elseif indep_model == 0 && num_regions == 3 # standard geosse with 3 regions
                if home_ind == 1 && away_ind == 2
                    m_i_d = 2.5
                elseif home_ind == 2 && away_ind == 1
                    m_i_d = 0.5
                elseif home_ind == 1 && away_ind == 3
                    m_i_d = 0.5
                elseif home_ind == 3 && away_ind == 1
                    m_i_d = 0.5
                elseif home_ind == 2 && away_ind == 3
                    m_i_d = 0.5
                elseif home_ind == 3 && away_ind == 2
                    m_i_d = 0.5
                end 
            end
        end 
        ## Compute the overall rates
        if model == 0
            m_d_d  = (log(n_i)+1)^d_d_home * (log(n_j+1)+1)^d_d_away
        elseif model == 1
            m_d_d  = (n_i^d_d_home)*(n_j+1)^d_d_away
        end
        #
        if only_density == 1 || (only_density == 0 && indep_model == 1)
            ans   = rho_d*m_i_d*m_d_d
        elseif only_density == 0 && indep_model == 0
            ans   = rho_d*m_i_d
        end
    end
    return(ans)
end

#########################
# Rate function for between-region speciation
#########################

# P.S. the range score function already has a factor of 1/2 to take into account
# left or right children orientation. 

# format for left_range/right_range follows the keys field in dict_keys_b (e.g. [1,0,1] to represent range AC)
function rate_b(base,num_species_per_region,left_range,right_range,rate_scalar,only_density,indep_model,model)
    ## Error checking
    @assert base >= 0 "base rate must be positive" # base rate must be positive 
    @assert left_range + right_range in dict_keys_b "Not a valid split" # left and right ranges must be derived from an ancestral range in dict_keys_b
    @assert left_range != right_range "left and right ranges cannot be equal" # left and right child ranges must be different
    #
    idx_leftrange  = findall(left_range .==1)  # find index of every region that are occupied by the species that forms left range 
    idx_rightrange = findall(right_range .==1) # find index of every region that are occupied by the species that forms right range 
    #
    if any(num_species_per_region[idx_leftrange].==0) == true || any(num_species_per_region[idx_rightrange].==0) == true 
        ans = 0 # if no species in the regions that form either left range or right range, rate = 0
    else
        ## Define rate parameters 
        rho_b      = base # base rate 
        b_d        = rate_scalar # rate scalar 
        if only_density == 1 
            m_i_b = 1 # means no density-independent factors 
        elseif only_density == 0
            if indep_model == 1
                m_i_b = rand(myLogNormal(1,1)) # draw density-independent rate from lognormal with mean 1, std 1 
            elseif indep_model == 0 && num_regions == 2 # standard geosse with 2 regions
                m_i_b = 0.5
            elseif indep_model == 0 && num_regions == 3 # standard geosse with 3 regions
                if (left_range == [1,0,0] && right_range == [0,1,0]) || (left_range == [0,1,0] && right_range == [1,0,0]) # A|B or B|A split
                    m_i_b = 0.5
                elseif (left_range == [1,0,0] && right_range == [0,0,1]) || (left_range == [0,0,1] && right_range == [1,0,0]) # A|C or C|A split
                    m_i_b = 0.5
                elseif (left_range == [0,1,0] && right_range == [0,0,1]) || (left_range == [0,0,1] && right_range == [0,1,0]) # B|C or C|B split
                    m_i_b = 0.5
                elseif (left_range == [1,0,0] && right_range == [0,1,1]) || (left_range == [0,1,1] && right_range == [1,0,0]) # A|BC or BC|A split
                    m_i_b = 0.5
                elseif (left_range == [0,1,0] && right_range == [1,0,1]) || (left_range == [1,0,1] && right_range == [0,1,0]) # B|AC or AC|B split
                    m_i_b = 0.5
                elseif (left_range == [0,0,1] && right_range == [1,1,0]) || (left_range == [1,1,0] && right_range == [0,0,1]) # C|AB or AB|C split
                    m_i_b = 0.5
                end
            end
        end 
        ## Compute the range split score f
        f = 0 # initialize the range score for left_range|right_range split 
        for i in 1:length(idx_leftrange)
            for j in 1:length(idx_rightrange)
                if model == 0
                    m_b_d_ij = (log(num_species_per_region[idx_leftrange[i]]*num_species_per_region[idx_rightrange[j]])+1)^b_d
                    m_b_d_ji = (log(num_species_per_region[idx_rightrange[j]]*num_species_per_region[idx_leftrange[i]])+1)^b_d
                elseif model == 1
                    m_b_d_ij = (num_species_per_region[idx_leftrange[i]]+num_species_per_region[idx_rightrange[j]]-1)^b_d
                    m_b_d_ji = (num_species_per_region[idx_rightrange[j]]+num_species_per_region[idx_leftrange[i]]-1)^b_d
                end
                f += 2/(m_b_d_ij + m_b_d_ji) # update f
            end
        end
        f = 1/f
        ## Compute the overall rates
        if only_density == 1 || (only_density == 0 && indep_model == 1)
            ans   = rho_b*m_i_b*f
        elseif only_density == 0 && indep_model == 0
            ans = rho_b*m_i_b
        end
    end
    return(ans)
end

#########################
# Function to get character extant tip data from history_df
#########################

function get_extant_tip_data(history_df,num_regions)
    tip_range_df      = filter([:extinct_status, :split_status] => (l, v) -> ismissing(l) && ismissing(v), history_df) # get tip state information
    tip_character_dat = Matrix{Any}(undef, nrow(tip_range_df), num_regions+1) # empty matrix to store tip state and tip label 
    for i in 1:nrow(tip_range_df)
        region_arr = deepcopy([k for (k,v) in ranges_inv if v==tip_range_df.final_state[i]][1]) # region presence/absence that correspond to range in range_inv
        tip_character_dat[i,1] = tip_range_df.idx[i] # assign tip index 
        for j in 1:num_regions
            tip_character_dat[i,(j+1)] = region_arr[j] # assign region presence/absence
        end
    end
    ans = DataFrame(tip_character_dat,:auto) #convert to dataframe
    return(ans)
end

#########################
# Function to get character of all tip data from history_df
#########################

function get_all_tip_data(history_df,num_regions)
    tip_range_df = filter(row -> 
    (ismissing(row.extinct_status) && ismissing(row.split_status)) ||
    (row.extinct_status == 1 && row.split_status == 0), 
    history_df)
    #
    tip_character_dat = Matrix{Any}(undef, nrow(tip_range_df), num_regions+1) # empty matrix to store tip state and tip label 
    for i in 1:nrow(tip_range_df)
        region_arr = deepcopy([k for (k,v) in ranges_inv if v==tip_range_df.final_state[i]][1]) # region presence/absence that correspond to range in range_inv
        tip_character_dat[i,1] = tip_range_df.idx[i] # assign tip index 
        for j in 1:num_regions
            tip_character_dat[i,(j+1)] = region_arr[j] # assign region presence/absence
        end
    end
    ans = DataFrame(tip_character_dat,:auto) #convert to dataframe
    return(ans)
end


#########################
# Write Newick string from history_df
#########################

function write_newick(history_df)
    newick_str   = "" # initialize an empty string
    #
    # starting from the root branch 
    if ismissing(history_df.split_status[1]) # if the root branch survive till present without any speciation
	    newick_str = string(history_df.idx[1]) * "node" * string(1) * ":" * string(history_df.branch_length[1])*")"
    else
        if history_df.split_status[1] == 1 # if the root branch ends with a speciation event, do (branch2., branch3.):length_1
		newick_str = "(branch2.,branch3.)" * "node" *  string(1) *  ":" * string(history_df.branch_length[1])
        elseif history_df.split_status[1] == 0 # if the root branch ends with a extinction event
            newick_str = string(history_df.idx[1]) * ":" * string(history_df.branch_length[1])*")"
        end
    end

    # loop over branches descending from the root branch
    for i in 2:nrow(history_df)
        if ismissing(history_df.split_status[i]) # if branch i is a tip branch
            prev_label  = "branch"*string(i)*"."
            next_label  = string(i) * ":" * string(history_df.branch_length[i])
            newick_str = replace(newick_str, prev_label => next_label)
        else # if branch i is not a tip branch
            if history_df.split_status[i] == 1 # if branch i ends with a speciation event, replace string i with its child indexes j,k e.g. "branchi." => "(branchj.,branchk.):length_i"
                prev_label = "branch"*string(i)*"."
		next_label = "(" * "branch"*string(history_df.child_1_idx[i]) * "." * "," * "branch"*string(history_df.child_2_idx[i]) * "." * ")" * "node" * string(i) *  ":" * string(history_df.branch_length[i])
                #next_label = "(" * "branch"*string(history_df.child_1_idx[i]) * "." * "," * "branch"*string(history_df.child_2_idx[i]) * "." * ") :" * string(history_df.branch_length[i])
                newick_str = replace(newick_str,prev_label => next_label)
            elseif history_df.split_status[i] == 0 # if branch i ends with an extinction event
                prev_label = "branch"*string(i)*"."
                next_label = "node" * string(i) * ":" * string(history_df.branch_length[i])
                newick_str = replace(newick_str,prev_label => next_label)
            end
        end 
    end
    ans = newick_str * ";"
    return(ans)
end

#########################
# Simulate a single tree
#########################

function simulate(sim_num)
    # print info
    println("calling simulate(" * string(sim_num) * ")")
    # look-up table of ranges -> integers
    global ranges_inv
    # dataset dimensions
    global num_regions, num_ranges
    # simulation identifiers
    global out_dir, out_path, out_prefix
    # # system variables
    global timebin
    # # Whether to resimulate tree if it does not meet min_taxa condition (for "numtaxa" and both conditions)
    global resimulate_tree = true
    #
    global num_trial = 1
    while resimulate_tree == true
        # set up containers
        num_species_per_range  = fill(0, num_ranges); # Empty vector to store number of species in each range
        num_species_per_region = fill(0, num_regions); # Empty vector to store number of species in each region
        # P.S. the entries for num_species_per_region and num_species_per_range will be updated at each time step
        # since rates only depend on n_i at the current time, not past times. 

        ###########################################
        # # Simulate events on tree (set up)
        ###########################################
        
        do_simulation                      = true; # Initial flag to continue simulating events 
        start_time                         = 0.0; # Start time of the simulation (0 = past time)
        start_range                        = floor(Int64,rand(Uniform(1,num_ranges))); # Starting range state at root node e.g. 1,2,...,7 for 3 region
        # start_range                        = 1;
        num_alive_lineage                  = 1; # Start with one species (stem branch)
        
        param_val = redraw_dist() # draw the parameter values 
        rho_w =  param_val[1];
        rho_e =  param_val[2];
        rho_d =  param_val[3];
        rho_b =  param_val[4];
        w_d   =  param_val[5];
        e_d   =  param_val[6];
        d_d_home =  param_val[7];
        d_d_away =  param_val[8];
        b_d     =   param_val[9];

        ### sim_geosse (standard geosse)
        # rho_w = 1;
        # rho_e = 1;
        # rho_d = 1;
        # rho_b = 1;

        ### sim7 (lemma 3)
        # rho_w = 1.75746;
        # rho_e = 0.01;
        # rho_d = 0.01;
        # rho_b = 1;

        ### sim_geosse (standard geosse)
        # w_d      = 0;
        # e_d      = 0;
        # d_d_home = 0;
        # d_d_away = 0;
        # b_d      = 0;

        # w_d      = rand(Normal(-1,0.4));
        # e_d      = rand(Normal(1,0.4));
        # d_d_home = rand(Normal(-1,0.4));
        # d_d_away = rand(Normal(-1,0.4));
        # b_d      = rand(Normal(-1,0.4));

        ### sim7 (lemma 3)
        # w_d      = -1;
        # e_d      = 2; # d_d_home = -1;
        # d_d_away = -1;
        # b_d      = -1;

        model        = 0; # using log DD model
        # model        = 1; # using standard dd model
        # only_density = 0; # including other density-independent factors when computing rates 
        only_density = 1; # exluding other density-independent factors when computing rates 

        # indep_model  = 0; # time-constant density-independent rates (standard geosse)
        indep_model  = 1 # time-varying density-independent rates
        
        
        ##########################
        # Update event at root node
        ##########################
        num_species_per_range[start_range]      = 1;
        range_key_root                          = [k for (k,v) in ranges_inv if v==start_range]; #range key that correspond to range label at the root e.g. (1,1,0) for label 4
        reg_idx_root                            = findall(range_key_root[1] .==1 ); # find the region indexes that correspond to that range at the root e.g. region_ind = 1,2 for range (1,1,0)
        num_species_per_region[reg_idx_root]    .= 1; # update the number of species currently alive in each region e.g. 1 species with range {A,B} means 1 species in A and 1 species in B
        #
        range_count_history  = [] # record num_species_per_range across time  
        region_count_history = [] # record num_species_per_region across time 
        #

        push!(range_count_history,deepcopy(num_species_per_range))
        range_count_history[1] = insert!(Float64[x for x in range_count_history[1]],1,start_time)

        push!(region_count_history,deepcopy(num_species_per_region))
        region_count_history[1] = insert!(Float64[x for x in region_count_history[1]],1,start_time)

        # create array of vectors containing information about num_species across regions in cts time
        timebin_int             = deepcopy(parse(Int,timebin)) # convert string (from args) to int
        cts_time_pts            = collect(range(0,max_time,length=timebin_int)) # continuos time points according to timebin

        species_cts_time_region = [] # initialize array of vectors containing information about num_species across regions in cts time
        for i in 1:timebin_int
            added_vec     = fill(missing,num_regions)
            added_vec_mod = Union{Any,Missing}[added_vec...] # allow another value to be added to a vector that contains missing type
            insert!(added_vec_mod,1,cts_time_pts[i]) # insert each continuous time point to each vector
            #
            push!(species_cts_time_region,added_vec_mod)
        end 
        species_cts_time_region[1][2:(num_regions+1)]=deepcopy(num_species_per_region)
        #
        species_cts_time_range = [] # initialize array of vectors containing information about num_species across ranges in cts time
        for i in 1:timebin_int
            added_vec     = fill(missing,num_ranges)
            added_vec_mod = Union{Any,Missing}[added_vec...] # allow another value to be added to a vector that contains missing type
            insert!(added_vec_mod,1,cts_time_pts[i]) # insert each continuous time point to each vector
            #
            push!(species_cts_time_range,added_vec_mod)
        end 
        species_cts_time_range[1][2:(num_ranges+1)]=deepcopy(num_species_per_range)
    
        
        curr_time = start_time; # Initialize current time of simulation p.s. time starts from 0 (past)
        ###############################################
        # Dataframe to record parent-child relationships 
        ###############################################
        # idx = branch label, time_start = starting time of the branch, time_end = end time of the branch (missing if still continues)
        # final_state = range state on the branch just before an event on that branch, extinct_status = whether its already extinct (1) or not (0) and "missing" if no event on the branch yet
        # split_state = whether it corresponds to a splitting event (1) or not (0) and "missing" if no event occuring on the branch yet,
        # branch_length = length of branch label idx (time_end - time_start), child_1_idx = idx of 1st child (missing if no speciation), child_2_idx = idx of 2nd child (missing if no speciation)
        history_df = DataFrame(idx = [],time_start = [],time_end = [],final_state = [],extinct_status = [],split_status = [],branch_length = [],child_1_idx = [],child_2_idx = [],event_type = [],parent_idx = []) 
        push!(history_df,(1,start_time,missing,start_range,missing,missing,missing,missing,missing,missing,-1)) # add first event at the root node (root node)

        
        ###############################################
        # Dataframe to record anagenetic transitions (local extinction on widespread range or dispersal) 
        ###############################################
        # idx = branch label associated with the event, time = time when the event occurs 
        # from_state = range state before the event, to_state = range state after the event 
        anagen_df = DataFrame(idx = [],time = [], from_state = [], to_state = [])
        #
        while do_simulation == true 
            if num_trial > 1
                param_val_rep = redraw_dist() #redraw the parameter values
                rho_w = param_val_rep[1];
                rho_e = param_val_rep[2];
                rho_d = param_val_rep[3];
                rho_b = param_val_rep[4];
                w_d   = param_val_rep[5];
                e_d   = param_val_rep[6];
                d_d_home = param_val_rep[7];
                d_d_away = param_val_rep[8];
                b_d     =  param_val_rep[9];
                num_trial = 1
            end
            #############################################
            # # get total sum of rates given the number of currently alive lineages
            #############################################
            w_total = 0 # initialize total sum of rates from within-region speciation
            e_total = 0 # initialize total sum of rates from extinction
            d_total = 0 # initialize total sum of rates from dispersal
            b_total = 0 # initialize total sum of rates from between-region speciation
            ############################################
            # # Create vector for state-event relationship 
            ############################################
            # P.S. 1st idx = range idx (e.g. A,B,AB,etc), 2nd idx = event type (w,e,d,b), 3rd idx = left range (for b) or home location for (for d) or region index (which region it occurs for w and e)
            # 4th idx = right range (for b) or away location (for d) or 0 (for w and e)
            state_event_col_name      = Tuple{Any,Any,Any,Any}[] 
            state_event_tot_rate      = zeros(0) # initialize vector to store values from each index in state_event_col_name
            ##
            for i in 1:num_ranges # loop over each range (total sum of rates = N_A*r_A + ...+ N_ABC*r_ABC)
                range_key = [k for (k,v) in ranges_inv if v==i] # find range that correspond to range label i e.g. [[1,0,0]] for i = 1
                reg_idx   = findall(range_key[1] .==1 ) # find region indexes correspond to that range e.g. reg_idx = 1 for range [1,0,0]
                #
                reg_away  = setdiff([1:1:num_regions;],reg_idx) # find all possible regions for dispersing for species with range i
                for j in 1:length(reg_idx) # for region-specific rates
                    n_j = num_species_per_region[reg_idx[j]] # find the number of currently living species in region reg_idx[j]
                    # compute N_i * w_j and N_i * e_j where i = range state and j in i  
                    prod_ni_wj = num_species_per_range[i]*rate_w(rho_w,reg_idx[j],n_j,w_d,only_density,indep_model,model) # e.g. for range {A,B}, w_{A,B} = n_AB * (w_A + w_B) with n_j species in A or B 
                    prod_ni_ej = num_species_per_range[i]*rate_e(rho_e,reg_idx[j],n_j,e_d,only_density,indep_model,model)
                    w_total += prod_ni_wj # e.g. w_total = N_A*(w_A) + N_B*(w_B) + N_AB*(w_A+w_B)
                    e_total += prod_ni_ej # e.g. e_total = N_A*(e_A) + N_B*(e_B) + N_AB*(e_A+e_B)
                    # append each state-event type to state_event, along with its column name
                    push!(state_event_col_name,(i,"w",reg_idx[j],0),(i,"e",reg_idx[j],0)) # e.g. (1,"w",1,0) denotes N_A * w_A or (2,"e",2,0) denotes N_B * w_B
                    append!(state_event_tot_rate,prod_ni_wj,prod_ni_ej)
                    #
                    if reg_away != Int64[] # if the species in range i is not already widespread in all regions 
                        for k in 1:length(reg_away)
                            # if  all(!,idx_match_d) # if that n_j and n_k does not exist yet, compute the rate and store in lookup_d
                                disp_rate = rate_d(rho_d,num_species_per_region,reg_idx[j],reg_away[k],d_d_home,d_d_away,only_density,indep_model,model) # compute dispersal rate from home to away
                                # compute N_i * d_jk where i = range state, and j in i, k not in i
                                prod_ni_dik = num_species_per_range[i]*disp_rate 
                                d_total += prod_ni_dik  # e.g. for range {A,B}, d_AB_j = n_AB*(d_A_C + d_B_C) or d_A_j = n_A*(d_A_B + d_A_C)
                                # append each state-event type to state_event, along with its column name
                                push!(state_event_col_name,(i,"d",reg_idx[j],reg_away[k])) # e.g. (1,"d",1,2) denotes N_A * d_AB or (4,"d",1,3) denotes N_AB * d_AC
                                append!(state_event_tot_rate,prod_ni_dik)
                        end
                    else # if already widepread in all regions
                        prod_ni_dik = num_species_per_range[i]*0  
                        d_total += prod_ni_dik
                        # Then we do not append to state_event for that range state and that particular event type 
                    end
                end
                #
                if count(==(1),range_key[1])==1 # e.g. if range_key = [[1,0,0]], then can't have between-region speciation
                    # compute N_i * b_leftrange_right_range where leftrange U rightrange = i
                    prod_ni_b =  num_species_per_range[i]*0
                    b_total += prod_ni_b
                    # Then we do not append to state_event 
                else
                    num_left_splits = length(evt_b[range_key[1]])
                    for k in 1:num_left_splits # e.g. for range {A,B,C}, can have left range = {A,B} or left_range = {A} (but not considering their change in orientation)
                        left_range  = evt_b[range_key[1]][k] 
                        right_range = range_key[1] - left_range
                        #
                        # compute N_i * b_leftrange_right_range where leftrange U rightrange = i
                        prod_ni_b = num_species_per_range[i]*rate_b(rho_b,num_species_per_region,left_range,right_range,b_d,only_density,indep_model,model)
                        b_total +=  prod_ni_b # e.g. i = {A,B}, then its N_AB*(b_A_B).agreed: no left-right orientation considered.
                        # append each state-event type to state_event, along with its column name
                        left_range_id = [v for (k,v) in ranges_inv if k==left_range] # e.g. 1 for range [1,0,0]
                        right_range_id = [v for (k,v) in ranges_inv if k==right_range]
                        push!(state_event_col_name,(i,"b",left_range_id[1],right_range_id[1])) # e.g. (4,"b",1,2) denotes N_AB * b_A_B
                        append!(state_event_tot_rate,prod_ni_b)
                    end
                end
            end
            # total sum of rates across ranges 
            # e.g. w_total = N_A*w_A + N_B*w_B + N_AB*(w_A+w_B) (loop over each range, and all event associated with each range)
            # e.g. e_total = N_A*e_A + N_B*e_B + N_AB*(e_A+e_B)
            # e.g. d_total = N_A*(d_AB) + N_B*(d_BA)
            # e.g. 
            total_rate = w_total+e_total+d_total+b_total # Note: this should match with sum(state_event_tot_rate) - checked 
            # Check these result total_rate and individual component by hand (Aug 28. Checked for correctness)
            #############################################
            # # draw time to next event
            #############################################
            wait_time = rand(Exponential(1/total_rate),1)[1] # draw waiting time till a next event 
            #
            curr_time += wait_time # update current time 
            #
            #############################################
            # # sampling next event type and associated state at the updated curr_time 
            #############################################
            # PS. each cell in state_event_tot_prob corresponds with a product of N_i and an event rate associated with the range i.
            if (curr_time < max_time && stop_condition == "time") || (stop_condition == "numtaxa") || (curr_time < max_time && stop_condition  == "both")# If conditioned on time, and curr_time > max_time, do not draw next event. 
                state_event_tot_prob  =  state_event_tot_rate/sum(state_event_tot_rate) # convert rates to probabilities (e.g. N_A*w_A | N_AB*w_A | ....| ...)
                rand_u = rand(Uniform(0,1)) # draw random number between [0,1]
                next_event = Tuple{Any,Any,Any}[] 
                next_state = zeros(Int64,0)
                for move in 1:length(state_event_tot_prob) # sampling from probability distribution
                    rand_u = rand_u - state_event_tot_prob[move]
                    if rand_u <0 # we have found over next event type and state associated with the event type 
                        push!(next_event,state_event_col_name[move][2:4]) # find the next event type
                        append!(next_state,state_event_col_name[move][1]) # find the associated state
                        break
                    end
                end 
                #############################################
                # # assign the next event type to a particular alive lineage that are available 
                #############################################
                #############################################
                # # Enable these prints if want to check to see lineage history at each iteration 
                ############################################
                # find row indexes from dataframe that satisfy the conditions (no event occuring yet and has the same range state as next_state)
                # println("parent-child history = ",history_df)
                # println("anagenetic history = ", anagen_df)
                # println("num species per range = ",num_species_per_range)
                # println("num species per region",num_species_per_region)
                # println("chosen state = ",next_state[1])
                # println("next event = ",next_event[1])
                # xxx
        
                row_idx    = intersect(findall(ismissing,history_df.extinct_status),findall(ismissing,history_df.split_status),
                                    findall(history_df.final_state .== next_state[1])) 
                row_chosen = rand(row_idx) # randomly sampling which branch in history_df[row_idx,:].idx that undergoes the next event 
        
                # # Example
                # push!(history_df,(2,start_time,missing,8,missing,missing))
                # push!(history_df,(3,start_time,missing,8,missing,missing))
                # push!(history_df,(4,start_time,missing,8,0,missing))
                # push!(history_df,(5,start_time,missing,8,missing,0))
                # push!(history_df,(6,start_time,missing,next_state[1],missing,missing))
                # #
                # row_idx = intersect(findall(ismissing,history_df.extinct_status),findall(ismissing,history_df.split_status),findall(history_df.final_state .== next_state[1]))
                # println(history_df)
                # #
                # history_df[row_idx,:extinct_status] .= 1 # replace a column entry from those row indexes with a new value
                # history_df[row_idx,:split_status] .= 0
                # println(history_df)
                # println(next_state[1])
                # println("row index = ",row_idx)
                # println(rand(row_idx))
                # xxx
                # 
        
                if next_event[1][1] == "w" # if the next event is within-region speciation (Sep 5: already checked)
                    max_idx = maximum(history_df.idx) # find the maximum idx from current recorded branch labels in history_df
                    push!(history_df,(max_idx+1,curr_time,missing,next_state[1],missing,missing,missing,missing,missing,missing,history_df[[row_chosen],:idx][1])) # child 1 inherits the parent's range
                    push!(history_df,(max_idx+2,curr_time,missing,next_event[1][2],missing,missing,missing,missing,missing,missing,history_df[[row_chosen],:idx][1])) # child 2 inherits the other range in case parent is widespread
                    # update the parent branch's status
                    # P.S. row_chosen must be transformed into a vector instead of scalar
                    history_df[[row_chosen],:time_end]       .= curr_time
                    history_df[[row_chosen],:extinct_status] .= 0
                    history_df[[row_chosen],:split_status]   .= 1
                    history_df[[row_chosen],:branch_length]  .= curr_time - deepcopy(history_df[[row_chosen],:time_start][1])
                    history_df[[row_chosen],:child_1_idx]    .= deepcopy(max_idx+1)
                    history_df[[row_chosen],:child_2_idx]    .= deepcopy(max_idx+2)
                    history_df[[row_chosen],:event_type]     .= "w"
                elseif next_event[1][1] == "e" && next_state[1] > num_regions # if the next event is extinction on a branch with widespread range (Sep 5: already checked)
                    # update the chosen branch's status 
                    original_range  = [k for (k,v) in ranges_inv if v==next_state[1]] # the original range of the branch before the extinction e.g. [[1,1,0]]
                    possible_ranges = [v for (k,v) in evt_e if k == original_range[1]][1] # all the possible ranges after the extinction from evt_e dict
                    lost_range = [k for (k,v) in ranges_inv if v==next_event[1][2]][1] # extinction occurs in this region next_event[1][2], causing the original range to lose this region
                    new_range = [v for v in possible_ranges if v != lost_range] # the new range for this branch after the event 
                    new_label = [v for (k,v) in ranges_inv if k==new_range[1]][1] # the new label e.g. 1 for range [1,0,0] or 7 for [1,1,1]
                    #
                    # assign the new range label to that branch idx (the branch is still in pool for either speciation or extinction)
                    history_df[[row_chosen],:final_state]          .= new_label 
                    # record the anagenetic change on that branch
                    push!(anagen_df,(history_df[[row_chosen],:idx][1],curr_time,next_state[1],new_label))
                    #
                elseif next_event[1][1] == "e" && next_state[1] <= num_regions # if the next event is extinction on a branch with endemic range (Sep 5: already checked)
                    # update the chosen branch's status 
                    history_df[[row_chosen],:time_end]       .= curr_time
                    history_df[[row_chosen],:extinct_status] .= 1
                    history_df[[row_chosen],:split_status]   .= 0
                    history_df[[row_chosen],:branch_length]  .= curr_time - deepcopy(history_df[[row_chosen],:time_start][1])
                    history_df[[row_chosen],:event_type]     .= "e"
                elseif next_event[1][1] == "d" # if the next event is dispersal (Sep 6: already checked)
                    # update the chosen branch's status 
                    original_range  = [k for (k,v) in ranges_inv if v==next_state[1]] # original range before dispersal event 
                    possible_ranges =  [v for (k,v) in evt_d if k==original_range[1]][1] # all possible ranges following a dispersal event in that original_range
                    new_range = [] # initialize empty range
                    for i in 1:length(possible_ranges)
                        if possible_ranges[i][next_event[1][3]] == 1 # find the new range following a dispersal to a new region next_event[1][3] 
                            push!(new_range,possible_ranges[i])
                        end
                    end
                    new_label = [v for (k,v) in ranges_inv if k==new_range[1]][1] # label for the new range e.g. 5 for [1,0,1]
                    # assign the new range label to that branch idx (the branch is still in pool for either speciation or extinction)
                    history_df[[row_chosen],:final_state]          .= new_label 
                    # record the anagenetic change on that branch
                    push!(anagen_df,(history_df[[row_chosen],:idx][1],curr_time,next_state[1],new_label))
                    #
                elseif next_event[1][1] == "b" # if the next event is between-region speciation (Sep 6: already checked)
                    # println("history before = ",history_df)
                    max_idx = maximum(history_df.idx) # find the maximum idx from current recoded branch labels in history_df
                    push!(history_df,(max_idx+1,curr_time,missing,next_event[1][2],missing,missing,missing,missing,missing,missing,history_df[[row_chosen],:idx][1])) # left range for child 1
                    push!(history_df,(max_idx+2,curr_time,missing,next_event[1][3],missing,missing,missing,missing,missing,missing,history_df[[row_chosen],:idx][1])) # right range for child 2
                    # update the parent branch's status
                    # P.S. row_chosen must be transformed into a vector instead of scalar (by doing [row_chosen])
                    history_df[[row_chosen],:time_end]       .= curr_time
                    history_df[[row_chosen],:extinct_status] .= 0
                    history_df[[row_chosen],:split_status]   .= 1
                    history_df[[row_chosen],:branch_length]  .= curr_time - deepcopy(history_df[[row_chosen],:time_start][1])
                    history_df[[row_chosen],:child_1_idx]    .= deepcopy(max_idx+1)
                    history_df[[row_chosen],:child_2_idx]    .= deepcopy(max_idx+2)
                    history_df[[row_chosen],:event_type]     .= "b"
                    #
                end
                
                ############################
                # update species count across regions in continuous time 
                ############################
                chosen_time = findall(x -> x <= curr_time, [v[1] for v in species_cts_time_region])
                    for i in 1:num_regions
                        chosen_time = intersect(chosen_time,findall(x -> ismissing(x), [v[i+1] for v in species_cts_time_region]))
                    end 
                    
                    for j in 1:length(chosen_time)
                        species_cts_time_region[chosen_time[j]][2:(num_regions+1)] = deepcopy(num_species_per_region)
                    end
                
                ############################
                # update species count across ranges in continuous time 
                ############################
                chosen_time = findall(x -> x <= curr_time, [v[1] for v in species_cts_time_range])
                    for i in 1:num_ranges
                        chosen_time = intersect(chosen_time,findall(x -> ismissing(x), [v[i+1] for v in species_cts_time_range]))
                    end 
                                    
                    for j in 1:length(chosen_time)
                        species_cts_time_range[chosen_time[j]][2:(num_ranges+1)] = deepcopy(num_species_per_range)
                    end

                ############################
                # # update num_species_per_region and num_species_per_range following an event
                ############################
                # update num_species_per_range following the previous event 
                for i in 1:num_ranges
                    count_row = count(row -> row.final_state .==i && ismissing(row.extinct_status) && ismissing(row.split_status), eachrow(history_df))
                    num_species_per_range[i] = copy(count_row)
                end
                # update num_species_per_region following the previous event 
                new_species_count = fill(0, num_regions)
                for i in 1:length(num_species_per_range) # loop over each range in current num_species_per_range vector
                    if num_species_per_range[i] != 0 # if there is a species (or more) with that range
                        range_key = [k for (k,v) in ranges_inv if v==i][1] # find the range key correspon to that range label
                        for j in 1:length(range_key) # loop over each region for that given range 
                            if range_key[j] !=0 # if there is a species (or more) in that region
                                global new_species_count[j] += num_species_per_range[i]*1 # record those species in that region
                            end
                        end
                    end
                end

                num_species_per_region = copy(new_species_count) # update species count in each region for the next event iteration
                # record num_species_per_region and num_species_per_range for this time step
                push!(range_count_history,deepcopy(num_species_per_range))
                range_count_history[end] = insert!(Float64[x for x in range_count_history[end]] ,1,curr_time)
                
                #
                push!(region_count_history,deepcopy(num_species_per_region))
                region_count_history[end] = insert!(Float64[x for x in region_count_history[end]] ,1,curr_time)
                #
                # update the number of currently alive lineages 
                num_alive_lineage    = sum(num_species_per_range)
                #print("num_alive lineage = ",num_alive_lineage)
            end
        #################################
        # # check if we stop simulation
        #################################
        
            ### Non tree extinction case
            if stop_condition == "time" && curr_time >= max_time && num_alive_lineage > 0  # if condition on max_time and we don't have tree extinction
                do_simulation = false
            elseif stop_condition == "numtaxa" && num_alive_lineage >= max_taxa && num_alive_lineage > 0  # if condition on max_taxa and we don't have tree extinction
                do_simulation = false
            elseif stop_condition == "both" && (num_alive_lineage >= max_taxa || curr_time >= max_time) && num_alive_lineage > 0  # if condition on time and taxa count and we don't have tree extinction
                do_simulation = false
            ### Tree extinction case
            elseif num_alive_lineage == 0 && resimulate == false  # we have tree extinction and resimulate == false
                do_simulation = false 
            elseif num_alive_lineage == 0 && resimulate == true # we have tree extinction and resimulate == true 
                do_simulation = true
                num_trial     = num_trial + 1
            end

        #################################
        # # post hoc tree treatment 
        #################################
        
            # stop simulation if condition met 
            if do_simulation == false && num_alive_lineage > 0 && stop_condition == "time" # if we don't have tree extinction and condition = "time"
                # # extend the branch length of all currently living lineages to the max_time 
                which_lineages = filter(row -> ismissing(row.extinct_status) && ismissing(row.split_status),history_df).idx
                # println(history_df)
                # println(which_lineages)
                for i in 1:length(which_lineages)
                    history_df.time_end[which_lineages[i]]      = copy(max_time) # extend their external branches to max time
                    history_df.branch_length[which_lineages[i]] = max_time - deepcopy(history_df.time_start[which_lineages[i]][1]) # compute the length of external branch
                end
                # Fill the remaining rows with missing value of species_cts_time_region with num_spec_per_region from previous entry (since we just extend branch lengths here)
                missing_rows = findall(x -> ismissing(x), [v[2] for v in species_cts_time_region])
                for i in 1:length(missing_rows)
                    species_cts_time_region[missing_rows[i]][2:(num_regions+1)] = deepcopy(num_species_per_region)
                end
                # Fill the remaining rows with missing value of  species_cts_time_range with num_spec_per_range from previous entry (since we just extend branch lengths here)
                missing_rows = findall(x -> ismissing(x), [v[2] for v in species_cts_time_range])
                for i in 1:length(missing_rows)
                    species_cts_time_range[missing_rows[i]][2:(num_ranges+1)] = deepcopy(num_species_per_range)
                end
            elseif do_simulation == false && num_alive_lineage > 0 && stop_condition == "numtaxa" # if we don't have tree extinction and condition = "numtaxa"
                # # extend the branch length of all currently living lineages to the curr_time 
                which_lineages = filter(row -> ismissing(row.extinct_status) && ismissing(row.split_status),history_df).idx
                if num_alive_lineage >= max_taxa && curr_time < max_time 
                    added_time = curr_time + 0.001
                    for i in 1:length(which_lineages)
                        history_df.time_end[which_lineages[i]]      = copy(added_time) # extend their external branches to added_time
                        history_df.branch_length[which_lineages[i]] = added_time - deepcopy(history_df.time_start[which_lineages[i]][1]) # compute the length of external branch
                    end
                elseif num_alive_lineage >= max_taxa && curr_time >= max_time 
                    added_time = (curr_time-wait_time) + 0.001
                    for i in 1:length(which_lineages)
                        history_df.time_end[which_lineages[i]]      = copy(added_time) # extend their external branches to added_time
                        history_df.branch_length[which_lineages[i]] = added_time - deepcopy(history_df.time_start[which_lineages[i]][1]) # compute the length of external branch
                    end
                end
                # Fill the remaining rows with missing value of species_cts_time_region with num_spec_per_region from previous entry (since we just extend branch lengths here)
                missing_rows = findall(x -> ismissing(x), [v[2] for v in species_cts_time_region])
                for i in 1:length(missing_rows)
                    species_cts_time_region[missing_rows[i]][2:(num_regions+1)] = deepcopy(num_species_per_region)
                end
                # Fill the remaining rows with missing value of  species_cts_time_range with num_spec_per_range from previous entry (since we just extend branch lengths here)
                missing_rows = findall(x -> ismissing(x), [v[2] for v in species_cts_time_range])
                for i in 1:length(missing_rows)
                    species_cts_time_range[missing_rows[i]][2:(num_ranges+1)] = deepcopy(num_species_per_range)
                end
            elseif do_simulation == false && num_alive_lineage > 0 && stop_condition == "both" # if we don't have tree extinction and condition = "numtaxa"
                which_lineages = filter(row -> ismissing(row.extinct_status) && ismissing(row.split_status),history_df).idx
                if num_alive_lineage >= max_taxa && curr_time < max_time #meeting max taxa requirement, and not meet time requirement
                    added_time = curr_time + 0.001
                    for i in 1:length(which_lineages)
                        history_df.time_end[which_lineages[i]]      = copy(added_time) # extend their external branches to added_time
                        history_df.branch_length[which_lineages[i]] = added_time - deepcopy(history_df.time_start[which_lineages[i]][1]) # compute the length of external branch
                    end
                elseif  num_alive_lineage >= max_taxa && curr_time >= max_time #meeting max taxa requirement, and time requirement
                    added_time = (curr_time-wait_time) + 0.001
                    for i in 1:length(which_lineages)
                        history_df.time_end[which_lineages[i]]      = copy(added_time) # extend their external branches to added_time
                        history_df.branch_length[which_lineages[i]] = added_time - deepcopy(history_df.time_start[which_lineages[i]][1]) # compute the length of external branch
                    end
                elseif num_alive_lineage < max_taxa && curr_time >= max_time # meeting time requirement, and not meet max taxa requirement
                    for i in 1:length(which_lineages)
                        history_df.time_end[which_lineages[i]]      = copy(max_time) # extend their external branches to max time
                        history_df.branch_length[which_lineages[i]] = max_time - deepcopy(history_df.time_start[which_lineages[i]][1]) # compute the length of external branch
                    end
                end
                # Fill the remaining rows with missing value of species_cts_time_region with num_spec_per_region from previous entry (since we just extend branch lengths here)
                missing_rows = findall(x -> ismissing(x), [v[2] for v in species_cts_time_region])
                for i in 1:length(missing_rows)
                    species_cts_time_region[missing_rows[i]][2:(num_regions+1)] = deepcopy(num_species_per_region)
                end
                # Fill the remaining rows with missing value of  species_cts_time_range with num_spec_per_range from previous entry (since we just extend branch lengths here)
                missing_rows = findall(x -> ismissing(x), [v[2] for v in species_cts_time_range])
                for i in 1:length(missing_rows)
                    species_cts_time_range[missing_rows[i]][2:(num_ranges+1)] = deepcopy(num_species_per_range)
                end
            elseif do_simulation == false && num_alive_lineage == 0 # if we have tree extinction and resimulate == false
                println("tree extinction")
            elseif do_simulation == true && num_alive_lineage == 0  # if we have tree extinction and resimulate == true
                println("tree extinction, " * "resimulate tree " * string(sim_num))
                # restart everything
                num_species_per_range  = fill(0, num_ranges); # Empty vector to store number of species in each range
                num_species_per_region = fill(0, num_regions); # Empty vector to store number of species in each region
                start_time                         = 0.0; # Start time of the simulation (0 = past time)
                # global start_range                        = floor(Int64,rand(Uniform(1,num_ranges))); # Starting range state at root node e.g. 1,2,...,7 for 3 region
                start_range                        = 1;
                num_alive_lineage                  = 1; # Start with one species (stem branch)
                num_species_per_range[start_range] = 1;
                #
                range_key_root                          = [k for (k,v) in ranges_inv if v==start_range]; #range key that correspond to range label at the root e.g. (1,1,0) for label 4
                reg_idx_root                            = findall(range_key_root[1] .==1 ); # find the region indexes that correspond to that range at the root e.g. region_ind = 1,2 for range (1,1,0)
                num_species_per_region[reg_idx_root]    .= 1; # update the number of species currently alive in each region e.g. 1 species with range {A,B} means 1 species in A and 1 species in B
                #
                range_count_history  = [] # record num_species_per_range across time  
                region_count_history = [] # record num_species_per_region across time 
                #
                push!(range_count_history,deepcopy(num_species_per_range))
                range_count_history[1] = insert!(Float64[x for x in range_count_history[1]] ,1,start_time)
                #
                push!(region_count_history,deepcopy(num_species_per_region))
                region_count_history[1] = insert!(Float64[x for x in region_count_history[1]] ,1,start_time)
                #
                species_cts_time_region = [] # initialize array of vectors containing information about num_species across regions in cts time
                for i in 1:timebin_int
                    added_vec     = fill(missing,num_regions)
                    added_vec_mod = Union{Any,Missing}[added_vec...] # allow another value to be added to a vector that contains missing type
                    insert!(added_vec_mod,1,cts_time_pts[i]) # insert each continuous time point to each vector
                    #
                    push!(species_cts_time_region,added_vec_mod)
                end 
                species_cts_time_region[1][2:(num_regions+1)]= deepcopy(num_species_per_region)
                #
                species_cts_time_range = [] # initialize array of vectors containing information about num_species across ranges in cts time
                for i in 1:timebin_int
                    added_vec     = fill(missing,num_ranges)
                    added_vec_mod = Union{Any,Missing}[added_vec...] # allow another value to be added to a vector that contains missing type
                    insert!(added_vec_mod,1,cts_time_pts[i]) # insert each continuous time point to each vector
                    #
                    push!(species_cts_time_range,added_vec_mod)
                end 
                species_cts_time_range[1][2:(num_ranges+1)]=deepcopy(num_species_per_range)        
                #
                curr_time = start_time; # Initialize current time of simulation p.s. time starts from 0 (past)
                history_df = DataFrame(idx = [],time_start = [],time_end = [],final_state = [],extinct_status = [],split_status = [],branch_length = [],child_1_idx = [],child_2_idx = [],event_type = [],parent_idx = []) 
                push!(history_df,(1,start_time,missing,start_range,missing,missing,missing,missing,missing,missing,-1)) # add first event at the root node (root node)
                anagen_df = DataFrame(idx = [],time = [], from_state = [], to_state = [])
            end
        end

        ## Enable this if want to validate the simulation
        # println("final history = ",history_df)
        # println("final anagenetic history = ", anagen_df)
        # println("region count over time = ",region_count_history)
        # println("range count over time",range_count_history)
        # println("final num alive lineage = ",num_alive_lineage)
        # # xxx
        # println("length history_df = ",nrow(history_df))

        #########################
        # Character data for extant tips 
        #########################
        extant_tip_character_df = get_extant_tip_data(history_df,num_regions)


        #########################
        # Character data for all tips 
        #########################
        all_tip_character_df = get_all_tip_data(history_df,num_regions)

        #########################
        # Add header names for character matrix
        #########################
        if num_regions == 2
            extant_tip_character_df = rename!(extant_tip_character_df, [:taxa, :region_1, :region_2])
            all_tip_character_df    = rename!(all_tip_character_df, [:taxa, :region_1, :region_2])
        elseif num_regions == 3
            extant_tip_character_df = rename!(extant_tip_character_df, [:taxa, :region_1, :region_2, :region_3])
            all_tip_character_df    = rename!(all_tip_character_df, [:taxa, :region_1, :region_2, :region_3])
        elseif num_regions == 4
            extant_tip_character_df = rename!(extant_tip_character_df, [:taxa, :region_1, :region_2, :region_3, :region_4])
            all_tip_character_df    = rename!(all_tip_character_df, [:taxa, :region_1, :region_2, :region_3, :region_4])
        elseif num_regions == 5
            extant_tip_character_df = rename!(extant_tip_character_df, [:taxa, :region_1, :region_2, :region_3, :region_4, :region_5])
            all_tip_character_df    = rename!(all_tip_character_df, [:taxa, :region_1, :region_2, :region_3, :region_4, :region_5])
        end 

        # #########################
        # # Convert character data at extant tips (multiple binary characters -> single multi-state characters)
        # #########################
        # extant_trans_tip_df     = DataFrame(taxa = [],state = [])
        # for i in 1:nrow(extant_tip_character_df)
        #     key_values = vec(Int.(reduce(vcat, eachcol(extant_tip_character_df[i:i, 2:end]))))
        #     state = [v for (k,v) in ranges_inv if k==key_values][1]
        #     push!(extant_trans_tip_df,(extant_tip_character_df[i,1],state))
        # end

        # #########################
        # # Convert character data at all tips (multiple binary characters -> single multi-state characters)
        # #########################
        # all_trans_tip_df     = DataFrame(taxa = [],state = [])
        # for i in 1:nrow(all_tip_character_df)
        #     key_values = vec(Int.(reduce(vcat, eachcol(all_tip_character_df[i:i, 2:end]))))
        #     state = [v for (k,v) in ranges_inv if k==key_values][1]
        #     push!(all_trans_tip_df,(all_tip_character_df[i,1],state))
        # end
        
        #########################
        # Rate labels and other summary statistics
        #########################
        if indep_model == 0 # if simulated under standard constant geosse model
            if stop_condition == "time"
                param_label_df =  DataFrame(log_base_w = log(rho_w),log_base_e = log(rho_e),log_base_d = log(rho_d),log_base_b = log(rho_b),exp_scalar_w = exp(0),
                                       exp_scalar_e = exp(0), exp_scalar_d_i = exp(0), exp_scalar_d_j = exp(0), exp_scalar_b = exp(0), stop_condition = "time");
            elseif stop_condition == "numtaxa"
                param_label_df =  DataFrame(log_base_w = log(rho_w),log_base_e = log(rho_e),log_base_d = log(rho_d),log_base_b = log(rho_b),exp_scalar_w = exp(0),
                                       exp_scalar_e = exp(0), exp_scalar_d_i = exp(0), exp_scalar_d_j = exp(0), exp_scalar_b = exp(0), stop_condition = "numtaxa");
            elseif stop_condition == "both"
                if num_alive_lineage >= max_taxa
                    param_label_df =  DataFrame(log_base_w = log(rho_w),log_base_e = log(rho_e),log_base_d = log(rho_d),log_base_b = log(rho_b),exp_scalar_w = exp(0),
                                       exp_scalar_e = exp(0), exp_scalar_d_i = exp(0), exp_scalar_d_j = exp(0), exp_scalar_b = exp(0), stop_condition = "numtaxa");
                else
                    param_label_df =  DataFrame(log_base_w = log(rho_w),log_base_e = log(rho_e),log_base_d = log(rho_d),log_base_b = log(rho_b),exp_scalar_w = exp(0),
                                       exp_scalar_e = exp(0), exp_scalar_d_i = exp(0), exp_scalar_d_j = exp(0), exp_scalar_b = exp(0), stop_condition = "time");
                end 
            end
        elseif indep_model == 1 || only_density == 1 # if simulated under dd-geosse with or without non-density dependent factors
            if stop_condition == "time"
                param_label_df =  DataFrame(log_base_w = log(rho_w),log_base_e = log(rho_e),log_base_d = log(rho_d),log_base_b = log(rho_b),exp_scalar_w = exp(w_d),
                exp_scalar_e = exp(e_d), exp_scalar_d_i = exp(d_d_home), exp_scalar_d_j = exp(d_d_away), exp_scalar_b = exp(b_d), stop_condition = "time");
            elseif stop_condition == "numtaxa"
                param_label_df =  DataFrame(log_base_w = log(rho_w),log_base_e = log(rho_e),log_base_d = log(rho_d),log_base_b = log(rho_b),exp_scalar_w = exp(w_d),
                exp_scalar_e = exp(e_d), exp_scalar_d_i = exp(d_d_home), exp_scalar_d_j = exp(d_d_away), exp_scalar_b = exp(b_d), stop_condition = "numtaxa");
            elseif stop_condition == "both"
                if num_alive_lineage >= max_taxa
                    param_label_df =  DataFrame(log_base_w = log(rho_w),log_base_e = log(rho_e),log_base_d = log(rho_d),log_base_b = log(rho_b),exp_scalar_w = exp(w_d),
                    exp_scalar_e = exp(e_d), exp_scalar_d_i = exp(d_d_home), exp_scalar_d_j = exp(d_d_away), exp_scalar_b = exp(b_d), stop_condition = "numtaxa");
                else
                    param_label_df =  DataFrame(log_base_w = log(rho_w),log_base_e = log(rho_e),log_base_d = log(rho_d),log_base_b = log(rho_b),exp_scalar_w = exp(w_d),
                    exp_scalar_e = exp(e_d), exp_scalar_d_i = exp(d_d_home), exp_scalar_d_j = exp(d_d_away), exp_scalar_b = exp(b_d), stop_condition = "time");
                end 
            end
        end
        
        #########################
        # Write Newick string 
        #########################
        newick_str = write_newick(history_df)

        ################## 
        # # Write outputs
        #################

        tree_file = joinpath(out_dir,out_prefix *"."*string(sim_num)* ".tre")

        # Ensure the output directory exists
        isdir(out_dir) || mkpath(out_dir)

        println("done simulating tree "*  string(sim_num))
        open(tree_file,"w") do file
            write(file,newick_str)
        end
        
        # # write tip data to file (as csv)
        # CSV.write(joinpath(out_dir,out_prefix*"."*string(sim_num)*".dat.csv"),extant_tip_character_df;header=true)
        CSV.write(joinpath(out_dir,out_prefix*"."*string(sim_num)*".extant"*".dat.csv"),extant_tip_character_df;header=true)
        CSV.write(joinpath(out_dir,out_prefix*"."*string(sim_num)*".dat.csv"),all_tip_character_df;header=true)

        if num_regions == 2
        # # convert each summary statistics to data frame 
            dat_reg_cts   = DataFrame(time = [x[1] for x in species_cts_time_region],region_1 = [x[2] for x in species_cts_time_region],
                                    region_2 = [x[3] for x in species_cts_time_region])

            dat_range_cts = DataFrame(time = [x[1] for x in species_cts_time_range],range_1 = [x[2] for x in species_cts_time_range],
                                    range_2 = [x[3] for x in species_cts_time_range],range_3 = [x[4] for x in species_cts_time_range])

            dat_reg   = DataFrame(time = [x[1] for x in region_count_history],region_1 = [x[2] for x in region_count_history],
                                region_2 = [x[3] for x in region_count_history])

            dat_range = DataFrame(time = [x[1] for x in range_count_history],range_1 = [x[2] for x in range_count_history],
                                range_2 = [x[3] for x in range_count_history],range_3 = [x[4] for x in range_count_history])   
        elseif num_regions == 3
        # convert each summary statistics to data frame 
            dat_reg_cts   = DataFrame(time = [x[1] for x in species_cts_time_region],region_1 = [x[2] for x in species_cts_time_region],
                                    region_2 = [x[3] for x in species_cts_time_region],region_3 = [x[4] for x in species_cts_time_region])

            dat_range_cts = DataFrame(time = [x[1] for x in species_cts_time_range],range_1 = [x[2] for x in species_cts_time_range],
                                    range_2 = [x[3] for x in species_cts_time_range],range_3 = [x[4] for x in species_cts_time_range],
                                    range_4 = [x[5] for x in species_cts_time_range],range_5 = [x[6] for x in species_cts_time_range],
                                    range_6 = [x[7] for x in species_cts_time_range],range_7 = [x[8] for x in species_cts_time_range])

            dat_reg   = DataFrame(time = [x[1] for x in region_count_history],region_1 = [x[2] for x in region_count_history],
                                    region_2 = [x[3] for x in region_count_history],region_3 = [x[4] for x in region_count_history])

            dat_range = DataFrame(time = [x[1] for x in range_count_history],range_1 = [x[2] for x in range_count_history],
                                    range_2 = [x[3] for x in range_count_history],range_3 = [x[4] for x in range_count_history],
                                    range_4 = [x[5] for x in range_count_history],range_5 = [x[6] for x in range_count_history],
                                    range_6 = [x[7] for x in range_count_history],range_7 = [x[8] for x in range_count_history])
        elseif num_regions == 4
            # convert each summary statistics to data frame 
            dat_reg_cts   = DataFrame(time = [x[1] for x in species_cts_time_region],region_1 = [x[2] for x in species_cts_time_region],
                                    region_2 = [x[3] for x in species_cts_time_region],region_3 = [x[4] for x in species_cts_time_region],
                                    region_4 = [x[5] for x in species_cts_time_region])
                                
            dat_range_cts = DataFrame(time = [x[1] for x in species_cts_time_range],range_1 = [x[2] for x in species_cts_time_range],
                                    range_2 = [x[3] for x in species_cts_time_range],range_3 = [x[4] for x in species_cts_time_range],
                                    range_4 = [x[5] for x in species_cts_time_range],range_5 = [x[6] for x in species_cts_time_range],
                                    range_6 = [x[7] for x in species_cts_time_range],range_7 = [x[8] for x in species_cts_time_range],
                                    range_8 = [x[9] for x in species_cts_time_range],range_9 = [x[10] for x in species_cts_time_range],
                                    range_10 = [x[11] for x in species_cts_time_range],range_11 = [x[12] for x in species_cts_time_range],
                                    range_12 = [x[13] for x in species_cts_time_range],range_13 = [x[14] for x in species_cts_time_range],
                                    range_14 = [x[15] for x in species_cts_time_range],range_15 = [x[16] for x in species_cts_time_range])
                                
            dat_reg   = DataFrame(time = [x[1] for x in region_count_history],region_1 = [x[2] for x in region_count_history],
                                    region_2 = [x[3] for x in region_count_history],region_3 = [x[4] for x in region_count_history],
                                    region_4 = [x[5] for x in region_count_history])
                                
            dat_range = DataFrame(time = [x[1] for x in range_count_history],range_1 = [x[2] for x in range_count_history],
                                    range_2 = [x[3] for x in range_count_history],range_3 = [x[4] for x in range_count_history],
                                    range_4 = [x[5] for x in range_count_history],range_5 = [x[6] for x in range_count_history],
                                    range_6 = [x[7] for x in range_count_history],range_7 = [x[8] for x in range_count_history],
                                    range_8 = [x[9] for x in range_count_history],range_9 = [x[10] for x in range_count_history],
                                    range_10 = [x[11] for x in range_count_history],range_11 = [x[12] for x in range_count_history],
                                    range_12 = [x[13] for x in range_count_history],range_13 = [x[14] for x in range_count_history],
                                    range_14 = [x[15] for x in range_count_history],range_15 = [x[16] for x in range_count_history])
        elseif num_regions == 5
            # convert each summary statistics to data frame 
            dat_reg_cts   = DataFrame(time = [x[1] for x in species_cts_time_region],region_1 = [x[2] for x in species_cts_time_region],
                                    region_2 = [x[3] for x in species_cts_time_region],region_3 = [x[4] for x in species_cts_time_region],
                                    region_4 = [x[5] for x in species_cts_time_region],region_5 = [x[6] for x in species_cts_time_region])
                                
            dat_range_cts = DataFrame(time = [x[1] for x in species_cts_time_range],range_1 = [x[2] for x in species_cts_time_range],
                                    range_2 = [x[3] for x in species_cts_time_range],range_3 = [x[4] for x in species_cts_time_range],
                                    range_4 = [x[5] for x in species_cts_time_range],range_5 = [x[6] for x in species_cts_time_range],
                                    range_6 = [x[7] for x in species_cts_time_range],range_7 = [x[8] for x in species_cts_time_range],
                                    range_8 = [x[9] for x in species_cts_time_range],range_9 = [x[10] for x in species_cts_time_range],
                                    range_10 = [x[11] for x in species_cts_time_range],range_11 = [x[12] for x in species_cts_time_range],
                                    range_12 = [x[13] for x in species_cts_time_range],range_13 = [x[14] for x in species_cts_time_range],
                                    range_14 = [x[15] for x in species_cts_time_range],range_15 = [x[16] for x in species_cts_time_range],
                                    range_16 = [x[17] for x in species_cts_time_range],range_17 = [x[18] for x in species_cts_time_range],
                                    range_18 = [x[19] for x in species_cts_time_range],range_19 = [x[20] for x in species_cts_time_range],
                                    range_20 = [x[21] for x in species_cts_time_range],range_21 = [x[22] for x in species_cts_time_range],
                                    range_22 = [x[23] for x in species_cts_time_range],range_23 = [x[24] for x in species_cts_time_range],
                                    range_24 = [x[25] for x in species_cts_time_range],range_25 = [x[26] for x in species_cts_time_range],
                                    range_26 = [x[27] for x in species_cts_time_range],range_27 = [x[28] for x in species_cts_time_range],
                                    range_28 = [x[29] for x in species_cts_time_range],range_29 = [x[30] for x in species_cts_time_range],
                                    range_30 = [x[31] for x in species_cts_time_range],range_31 = [x[32] for x in species_cts_time_range])
                                
            dat_reg   = DataFrame(time = [x[1] for x in region_count_history],region_1 = [x[2] for x in region_count_history],
                                    region_2 = [x[3] for x in region_count_history],region_3 = [x[4] for x in region_count_history],
                                    region_4 = [x[5] for x in region_count_history],region_5 = [x[6] for x in region_count_history])
                                
            dat_range = DataFrame(time = [x[1] for x in range_count_history],range_1 = [x[2] for x in range_count_history],
                                    range_2 = [x[3] for x in range_count_history],range_3 = [x[4] for x in range_count_history],
                                    range_4 = [x[5] for x in range_count_history],range_5 = [x[6] for x in range_count_history],
                                    range_6 = [x[7] for x in range_count_history],range_7 = [x[8] for x in range_count_history],
                                    range_8 = [x[9] for x in range_count_history],range_9 = [x[10] for x in range_count_history],
                                    range_10 = [x[11] for x in range_count_history],range_11 = [x[12] for x in range_count_history],
                                    range_12 = [x[13] for x in range_count_history],range_13 = [x[14] for x in range_count_history],
                                    range_14 = [x[15] for x in range_count_history],range_15 = [x[16] for x in range_count_history],
                                    range_16 = [x[17] for x in range_count_history],range_17 = [x[18] for x in range_count_history],
                                    range_18 = [x[19] for x in range_count_history],range_19 = [x[20] for x in range_count_history],
                                    range_20 = [x[21] for x in range_count_history],range_21 = [x[22] for x in range_count_history],
                                    range_22 = [x[23] for x in range_count_history],range_23 = [x[24] for x in range_count_history],
                                    range_24 = [x[25] for x in range_count_history],range_25 = [x[26] for x in range_count_history],
                                    range_26 = [x[27] for x in range_count_history],range_27 = [x[28] for x in range_count_history],
                                    range_28 = [x[29] for x in range_count_history],range_29 = [x[30] for x in range_count_history],
                                    range_30 = [x[31] for x in range_count_history],range_31 = [x[32] for x in range_count_history])
        end
        # # write summary statistics and rates to file
        CSV.write(joinpath(out_dir,out_prefix*"."*string(sim_num)*".labels.csv"),param_label_df) # write rate parameters and number of alive tips
        CSV.write(joinpath(out_dir,"anagenetic_changes_tree_"*string(sim_num)*".csv"),anagen_df) # write anagenetic history along branches over event times 
        CSV.write(joinpath(out_dir,"species_count_region_tree_"*string(sim_num)*".csv"),dat_reg) # write species count across regions over event times
        CSV.write(joinpath(out_dir,"species_count_range_tree_"*string(sim_num)*".csv"),dat_range) # write species count across ranges over event times
        CSV.write(joinpath(out_dir,"cts_species_count_region_tree_"*string(sim_num)*".csv"),dat_reg_cts) # write species count across regions over cts times
        CSV.write(joinpath(out_dir,"cts_species_count_range_tree_"*string(sim_num)*".csv"),dat_range_cts) # write species count across ranges over cts times
        CSV.write(joinpath(out_dir,"cladogenetic event history_"*string(sim_num)*".csv"),history_df) # write cladogenetic event history over event times 
        #
        if num_alive_lineage < min_taxa && (stop_condition == "both" || stop_condition == "numtaxa")
            resimulate_tree = true 
            println("Final number of alive lineage = ",num_alive_lineage)
            println("resimulate tree = true")
        else
            resimulate_tree = false
        end
    end
end

#########################
# Simulate a set of trees
#########################

for i in start_idx:(start_idx + batch_size - 1)
    println("Simulation ", i)
    Random.seed!(i)
    #global max_taxa       = rand(Uniform(100,150)); # Maximum number of extant taxa to simulate
    #global max_time       = rand(Uniform(70,100)); # Maximum simulation time 
    simulate(i)
end

exit()


