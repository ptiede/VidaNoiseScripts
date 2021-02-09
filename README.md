# Using VIDA
To use VIDA on your image reconstructions the first step is to install the VIDA package using Julia.
To install Julia I recommend either of the following:

1. Download the Julia binary from [julialang](https://julialang.org/downloads/)
2. If on Linux or MacOS you can use the jill installer [here](https://github.com/abelsiqueira/jill)

This will install the most recent version of Julia, 1.5.3 and should work.

To install VIDA please:

1. Run julia from the terminal by typing `julia` to launch the REPL.
2. Type `]` in the REPL to drop into the Pkg environment. 
3. Type `instantiate` to install VIDA and some additional dependencies needed. The `instantiate` will set up your local Julia environment to match the versions so that everything will work. Note that this will only install the VIDA packages in the local VIDANoiseScripts environment. To install is to the global Julia installation open a new Julia REPL and type 
```julia-repl
]
add VIDA, ArgParse, CSV, DataFrames
```
4. Finally, before running VIDA I would recommend to precompile the installed packages so lower latency. To do this enter the Pkg environment in the REPL
```julia-repl
]
precompile
```

Now to run VIDA on your images first create a txt document with all the **Fits** image file paths.
To do this I typically use 
```bash
readlink -f PATH/TO/IMAGE/DIRECTORY/*.fits > imagelist
```

To run VIDA on the image list then use
With an asymmetric Gaussian:
```bash
julia -p ncores bbimage_extractor.jl outdir/imagelist -c 0.0 --filter AsymG --stride 250
```
or with a general cosine ring with a first and fourth order expansion in thickness and brightness:
```bash
julia -p ncores cosine_extractor.jl outdir/imagelist -c 0.0 --filter 1 4 --stride 250
```
which will make a `vida_summary.csv` csv file with the optimal filter parameters.

# Using VIDA and Themis together
To run VIDA on a Themis chain I have also created a number of helper scripts. These section describes
how to use them.  

The TLDR is that to create a vida and noise chain file run the following

To create set of fits images and imagelist from a Themis chain

```bash
./gen_images.sh chainfile nimages outdir
```
 which also creates a `subchain.dat` file containing the parameters of the image.
 
 **NOTE: You need to ensure you have the model_image.tag file in the directory you run this in.**

To run VIDA
With an asymmetric Gaussian:
```bash
julia -p ncores bbimage_extractor.jl outdir/imagelist -c 0.0 --filter AsymG --stride 250
```
or with a general cosine ring with a first and fourth order expansion in thickness and brightness:
```bash
julia -p ncores cosine_extractor.jl outdir/imagelist -c 0.0 --filter 1 4 --stride 250
```
which will make a `vida_summary.csv` file with the optimal filter parameters

If your chain included noise model parameters we can then merge the VIDA results with the 
chain to create a merged feature, noise model chain.
To create a merged chain list:
```bash
python merge_chains.py -c Test/subchain.dat -v vida_summary.csv -o output_chain
```

**Again make sure you have a model_image.tag file in this directory.**


For a more detailed breakdown please read below


This is done in three steps

1. Run `./gen_images.sh "CHAINFILE" "NIMAGES" "OUTDIRECTORY"` where
    - `CHAINFILE` is a Themis chainfile
    - `NIMAGES` is the number of images you want to generate from the chain
    - `OUTDIRECTORY` is where you want to save the fits images and imagelist to be used with VIDA

    This script will result in a set of images and a  file that contains the paths to each image called `imagelist`.

2. Run either 
    - `julia -p #NCORES bbimage_extractor.jl imagelist -c 0.0 --filter TYPEOFFILTER --stride STRIDE`
    - `julia -p #NCORES cosine_extractor.jl imagelist -c 0.0 --filter N M --stride STRIDE`
    
    The options mean:
    - `imagelist` is the file with the list of images set in step 1.
    - `-c` Clips the pixel intensities to 0 if they have a pixel intensity less than the specified model. Basically always set this to 0.0 to ensure there are no negative intensity pixels.
    - `--stride` which sets the checkpoint frequency. This should be largish. Checkpointing is expensive so set this number high. For instance if your imagelist has 1000 images I would checkpoint every 250 images.
    - `--filter` This is the filter function you want to use depending on your image structure. The available filters for the bbimage_extractor are:
        1. Gen: This is a GeneralGaussianRing filter that has a radius, width, asymmetry, slash, position angle, and center.
        2. Slash: This is a SlashedGaussianRing filter that only has a linear slash and no asymmetry
        3. Ellip: This is a EllipticalGaussianRing filter that only has asymmetry no slash.
        4. Circ: This is a GaussianRing filter that is is a Gaussian ring with some radius and width.
        5. AsymG: This is an AsymGaussian filter that is an asymmetric Gaussian
  
      Note that to add multiple filter together you would just do e.g. `--filter Gen AsymG AsymG`.

3. To create a merged chain file run:
    - `python merge_chains.py -c OUTDIRECTORY/subchain.dat -v vida_summary.csv` 

which will merge the subchain.dat file from the Themis run and the fit_summaries file from VIDA. Note that it removes the Themis model components not related to the uncertainty, and replaces them with the VIDA optimal filter parameters. 
