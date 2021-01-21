using Distributed
using ArgParse
Distributed.@everywhere using VIDA
using CSV
using DataFrames
using Random


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "arg1"
            help = "file list of fits images to read"
            arg_type = String
            required = true
        "--stride"
             help = "Checkpointing stride, i.e. number of steps."
             arg_type = Int
             default = 200
        "--out"
            help = "name of output files with extractions"
            default = "fit_summaries.txt"
        "--restart"
            help = "Tells the sampler to read in the old files and restart the run."
            action = :store_true
        "-c"
            help = "lower percentage of flux to clip from image"
            arg_type = Float64
            default = 0.0
        "--div"
            help = "Divergence type to be used in image reconstruction"
            arg_type = String
            default = "Bh"
        "--filter"
            help = "Cosine filter to use for feature extraction. Expects two numbers"
            action = :store_arg
            nargs = 2
            required=true
        "--seed"
            help = "Random seed for initial positions in extract"
            arg_type = Int
            default = 42
        "--istart"
            help = "start index of fits file"
            arg_type = Int
            default = 1
        "--iend"
            help = "end index of fits file"
            arg_type = Int
            default = -1
    end
    return parse_args(s)
end

function main()
    #Parse the command line arguments
    parsed_args = parse_commandline()
    #Assign the arguments to specific variables
    fitsfiles = parsed_args["arg1"]
    out_name = parsed_args["out"]
    clip_percent = parsed_args["c"]
    startindx = parsed_args["istart"]
    endindx = parsed_args["iend"]
    seed = parsed_args["seed"]
    div_type = parsed_args["div"]
    stride = parsed_args["stride"]
    if (parsed_args["div"]=="KL")
      div_type = "KL"
    elseif (parsed_args["div"]=="Bh")
      div_type = "Bh"
    else
      error("$div_type not found! Must be Bh or KL")
    end
    filter_type = parse.(Int, parsed_args["filter"])
    println("Filter type $filter_type")
    restart = parsed_args["restart"]
    println("Using options: ")
    println("list of files: $fitsfiles, ")
    println("output name: $out_name, ")
    println("Clipping $(clip_percent*100) percent")
    println("random seed $seed")
    println("divergence type $div_type")
    println("Checkpoint stride $stride")

    #Read in a file and create list of images to filter
    #the last line is the termination of the file
    files = split(read(fitsfiles,String),"\n")

    #check if the last entry of files is an empty string
    if files[end] == ""
        files = files[1:end-1]
    end

    if startindx == 0
      startindx = 1
    end

    if endindx == -1 || endindx > length(files)
      endindx = length(files)
    end

    println("starting index $startindx")
    println("ending index $endindx")
    println("Starting fit")
    files = files[startindx:endindx]


    #Now run on the files for real
    main_sub(files, out_name,
             div_type,CosineRing{filter_type[1],filter_type[2]},
             clip_percent, seed,
             restart, stride)
    println("Done! Check $out_name for summary")
    return 0
end

function make_initial_filter(::Type{CosineRing{N,M}}) where {N,M}
    @show N,M
    lower_σ = [0.01, [ -2.0 for i in 1:N]... ]
    upper_σ = [15.0, [ 2.0 for i in 1:N]... ]
    lower_ξσ = Float64[]
    upper_ξσ = Float64[]
    if N > 0
        lower_ξσ = Float64[-π for i in 1:N ]
        upper_ξσ = Float64[ π for i in 1:N ]
    end
    lower_s = [0.001, [-0.99 for i in 2:M]...]
    upper_s = [0.999, [0.99 for i in 2:M]...]
    lower_ξs = [-π for i in 1:M ]
    upper_ξs = [π for i in 1:M ]

    lower = [5.0 ,
             lower_σ..., lower_ξσ...,
             0.001, 0.0,
             lower_s..., lower_ξs...,
             -60.0, -60.0
            ]
    upper = [ 30.0 ,
              upper_σ..., upper_ξσ...,
              0.999, π,
              upper_s..., upper_ξs...,
              60.0, 60.0
            ]
    filter = CosineRing{N,M}(lower)
    filt_lower = CosineRing{N,M}(lower)
    filt_upper = CosineRing{N,M}(upper)
    return (filter, filt_lower, filt_upper)
end

function create_initial_df!(start_indx::Int, fitsfiles, filter, restart, out_name)
    df = DataFrame()
    nfiles = length(fitsfiles)
    if !restart
      #we want the keynames to match the model parameters
      key_names = fieldnames(typeof(filter))
      n = length(filter.θ1.σ)
      m = length(filter.θ1.s)

      key_names = Symbol[:r0,
                    Symbol[:σ for i in 1:n]...,
                    Symbol[:ξσ for i in 2:n]...,
                    :τ,
                    :ξτ,
                    Symbol[:s for i in 1:m]...,
                    Symbol[:ξs for i in 1:m]...,
                    :x0,
                    :y0,
                    :Irel
                ]
      for i in 1:length(key_names)
        insertcols!(df, ncol(df)+1, Symbol(key_names[i]) => zeros(nfiles); makeunique=true)
      end
      @show key_names
      #fill the data frame with some likely pertinent information
      insertcols!(df, ncol(df)+1, :divmin=>zeros(nfiles))
      insertcols!(df, ncol(df)+1, :fitsfiles=>fitsfiles)
    else
      df = DataFrame(CSV.File(out_name, delim=";"))
      start_indx = convert(Int, findfirst(isequal(0.0), df[:,1]) )::Int
      println("Restarting run for $out_name at index $start_indx")
    end
    return df
end

@everywhere function make_div(image::EHTImage{S,T},clip) where {S,T}
    cimage = VIDA.clipimage(clip,image)
    div = VIDA.Bhattacharyya(cimage)
    return div
end

@everywhere function fit_func(filter,lower, upper, clip_percent)
    function (file)
        println("Extracting $file")
        image = load_fits(String(file))
        div = make_div(image, clip_percent)
        prob = ExtractProblem(div, filter, lower, upper)
        θ,divmin = extractor(prob, BBO(maxevals=5_000))
        prob_new = ExtractProblem(div, θ, lower, upper)
        θ2,divmin2 = extractor(prob_new, CMAES(ftol=1e-6, cov_scale=0.01, verbosity=0))
        return θ2,divmin2
    end
end


function main_sub(fitsfiles, out_name,
                  div_type, ::Type{CosineRing{N,M}},
                  clip_percent, seed,
                  restart, stride) where {N,M}

    #"Define the filter I want to use and the var bounds"
    matom = make_initial_filter(CosineRing{N,M})
    model = matom[1] + 1.0*Constant()
    lower = matom[2] + 1e-8*Constant()
    upper = matom[3] + 1.0*Constant()


    #Need to make sure all the procs know this information
    @everywhere model = $(model)
    @everywhere lower = $(lower)
    @everywhere upper = $(upper)
    @everywhere div_type = $(div_type)
    @everywhere clip_percent = $(clip_percent)

    #Set up the data frame to hold the optimizer output that
    #will be saved
    start_indx = 1
    df = create_initial_df!(start_indx,fitsfiles, model, restart, out_name)
    rng = MersenneTwister(seed) #set the rng

    #Now fit the files!
    @everywhere fit = fit_func(model,lower,upper, clip_percent)
    indexpart = Iterators.partition(start_indx:length(fitsfiles), stride)
    for ii in indexpart
      results = pmap(fit, fitsfiles[ii])
      #return df, results
      #println(hcat(unpack.(first.(results))...)')
      df[ii,1:VIDA.size(typeof(lower))] = hcat(unpack.(first.(results))...)'
      df[ii,VIDA.size(typeof(lower))+1] = last.(results)
      df[ii,end] = fitsfiles[ii]
      #save the file
      println("Checkpointing $(ii)")
      CSV.write(out_name, df, delim=',')
    end

    return df
end




main()
