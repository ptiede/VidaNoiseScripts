import themispy as tpy
import pandas as pd
import argparse as ap
import numpy as np


parser = ap.ArgumentParser()
parser.add_argument("-c", type=str, action="store", required=True, help="subchain file used to make imagelist")
parser.add_argument("-v", type=str, action="store", required=True, help="VIDA fit_summary from imagelist")
parser.add_argument("-o", type=str, action="store", help="Output file name", default="vida_themis_noise_chain.dat")
parser.add_argument("-t", type=str, action="store", help="model_image tag file", default="model_image.tag")
args = parser.parse_args()

chainfile = args.c
vidafile = args.v
tagfile = args.t
outfile = args.o


#construct the model and extract the uncertainty params
model_image = tpy.vis.construct_model_image_from_glob(tagfile, 'tagvers-1.0')
chain = np.loadtxt(chainfile, skiprows=2)
chain_noise = chain[:,model_image.size:]
with open(chainfile) as f:
    cheader = f.readline()
cheader = cheader.split("#")[-1]
#load the VIDA fits summary file
df = pd.read_csv(vidafile, sep=";")
vida = np.array(df)[:,:-2]
image = np.array(df)[:,-1]

#Now stack the two chains together
chain_tot = np.zeros((len(vida), vida.shape[1]+chain_noise.shape[1]))
chain_tot[:,:vida.shape[1]] = vida
chain_tot[:,vida.shape[1]:] = chain_noise
headerv = ["         vida_"+str(i) for i in range(vida.shape[1])]
headern = ["         noise_"+str(i) for i in range(chain_noise.shape[1])]
header = "          ".join(headerv+headern)#+["imagename"]
np.savetxt(outfile, chain_tot, header=cheader+header, fmt="%24s")

