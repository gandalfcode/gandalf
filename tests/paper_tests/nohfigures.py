from gandalf.analysis.facade import *
import matplotlib.pyplot as plt
import argparse

def plot_rrho():
    window()
    plotanalytical("r","rho")
    addplot("r","rho")
    plt.xlim(0,1)
    plt.yscale('log')
    plt.ylim(1,1e2)
    plt.savefig("noh_rho.png")

def plot_render():
    window()
    renderslice("x","y","rho",0,res=1024)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.savefig("noh_render.png")



switch_nongui()

parser=argparse.ArgumentParser()
parser.add_argument("snap",type=int)
parser.add_argument("--sim",default='NOH1')
args=parser.parse_args()


loadsim(args.sim)
snap(args.snap)

plot_rrho()

plot_render()