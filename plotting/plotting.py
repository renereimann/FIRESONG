import matplotlib
matplotlib.use("Agg")
import os, glob, re, argparse
import numpy as np
import matplotlib.pyplot as plt

def plot_firesong(infile, savepath=None, 
                  densitiy=1e-9, fluxnorm=1e-8, index=2.0, 
                  evolution = "HB2006SFR", zmax=10., 
                  LF="LG", sigma=1.0, L=4e52):
                  # HB2006SFR, NoEvolution
                  # SC, LG, PL
    indata = []
    with open(infile, "r") as open_file:
        for line in open_file:
            if line.startswith("#"): continue
            if len(line.split()) != 3: continue
            indata.append( line.split() )
    
    data = np.zeros(len(indata), dtype=[("dec", ">f4"), ("z", ">f4"), ("flux"">f4")])
    data["dec"] = [d[0] for d in indata]
    data["z"] = [d[1] for d in indata]
    data["flux"] = [d[2] for d in indata]
            
    nbins_sindec = 40
    nbins_z = 40
    nbins_flux = 40
    
    label1 = r"$\rho=%.1e/\mathrm{Mpc}^3$"%densitiy+"\n"+r"$N_\mathrm{sources}=%d$"%len(data)+"\n"+"$\gamma=-%.2f$"%index+"\n"+"$F_\mathrm{diff}$=%.2e"%fluxnorm
    label2 = "Evolution: {evolution}\nz-range: 0.0-{zmax}".format(**locals()) 
    label3 = "Lognormal dist.\n" if LF=="LG" else "Power-law dist.\n" if LF=="PL" else "Standart candle\n"
    label3 += r"$L_\nu=${L} erg/yr".format(**locals())
    if LF != "SC":
        label3 += "\n$\sigma={sigma}$".format(**locals())
    
    fig = plt.figure(figsize=(10,8))
    
    plt.subplot(221)
    plt.hist(np.sin(np.radians(data["dec"])), histtype="step", bins=np.linspace(-1., 1., nbins_sindec), label=label1)
    plt.xlabel(r"$\sin(\delta)$")
    plt.ylabel(r"counts")
    plt.legend(loc="lower center")
    
    plt.subplot(222)
    plt.hist(data["z"], histtype="step", bins=np.linspace(0., 10., nbins_z), label=label2)
    plt.xlabel(r"z")
    plt.ylabel(r"counts")
    plt.legend(loc="best")
    
    plt.subplot(223)
    plt.hist(np.log10(data["flux"]), histtype="step", bins=np.linspace(-17, -7, nbins_flux), label=label3)
    plt.yscale("log", nonposy="clip")
    plt.xlabel(r"$\log_{10}(\phi_0 / (\mathrm{GeV}\, \mathrm{cm}^{-2}\, \mathrm{s}^{-1}))$")
    plt.ylabel(r"counts")
    plt.legend(loc="lower center")
    
    ax1 = plt.subplot(224)
    hi, xedges, yedges, p0 = ax1.hist2d(data["z"], np.log10(data["flux"]), bins=(np.linspace(0., 10., nbins_z),np.linspace(-17, -7, nbins_flux)), cmin=1e-5) 
    
    plt.xlabel(r"z")
    plt.ylabel(r"$\log_{10}(\phi / (\mathrm{GeV}\, \mathrm{cm}^{-2}\, \mathrm{s}^{-1}))$")
    
    cbaxes = fig.add_axes([0.65, 0.45, 0.3, 0.02])
    plt.colorbar(mappable=p0, cax=cbaxes, orientation="horizontal").set_label("counts")
    
    plt.tight_layout()
    
    if savepath is not None:
        plt.savefig(savepath)
        plt.close()
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infiles", type=str, nargs="+")
    args = parser.parse_args()
    
    for f in args.infiles:
        try:
            density = float(re.search('_density_[0-9\.e\-]*_', os.path.basename(f)).group(0).split("_")[-2])
            evolution = re.search('_evolution_[0-9a-zA-Z]*_', os.path.basename(f)).group(0).split("_")[-2]
            zmax = float(re.search('_zmax_[0-9\.]*_', os.path.basename(f)).group(0).split("_")[-2])
            fluxnorm = float(re.search('_fluxnorm_[0-9\.e\-]*_', os.path.basename(f)).group(0).split("_")[-2])
            index = float(re.search('_index_[0-9\.]*_', os.path.basename(f)).group(0).split("_")[-2])
            LF = re.search('_LF_[A-Z]*_', os.path.basename(f)).group(0).split("_")[-2]
            sigma = float(re.search('_sigma_[0-9\.]*_', os.path.basename(f)).group(0).split("_")[-2])
            L = float(re.search('_L_[0-9\.e\-]*\.', os.path.basename(f)).group(0).split("_")[-1][:-1])
            savepath = f.replace(".out", ".png")
            if density == 1e-5: continue
            plot_firesong(f, savepath=savepath, densitiy=density, fluxnorm=fluxnorm, index=index, evolution=evolution, zmax=zmax, LF=LF, sigma=sigma, L=L)
        except:
            continue
