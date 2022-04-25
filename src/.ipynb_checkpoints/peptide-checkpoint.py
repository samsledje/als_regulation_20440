import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class Peptide:
    def __init__(self, mass, charge, scans, rt, seq, spectrum):
        self.mass = float(mass)
        self.charge = charge
        self.scans = scans
        self.rt = float(rt)
        self.seq = seq
        self.spectrum = spectrum
        
    def plot_spectrum(self, **kwargs):
        title = f"Peptide Sequence: {self.seq} Mass: {self.mass} ({self.charge})"
        self.spectrum.plot(title=title, **kwargs)
        
class Spectrum:
    def __init__(self, peaks):
        self.peaks = pd.DataFrame(np.array([[float(mz), float(y)] for (mz,y) in peaks]),columns=['MZ','Y'])
        
    def plot(self, title='', annot_thresh=-1, figsize=(15,10)):
        plt.figure(figsize=figsize)
        plt.vlines(self.peaks['MZ'], 0, self.peaks['Y'])
        plt.title(title)
        if annot_thresh > 0:
            for _, (mz, y) in self.peaks.iterrows():
                if y > annot_thresh:
                    plt.text(mz, y, f"{mz:.3f}")
        sns.despine()
        plt.show()
        
def parse_peptide_lines(line_list):
    attributes = {}
    peaks = []
    for line in line_list:
        if '=' in line:
            k,v = line.split('=')
            attributes[k] = v
        else:
            peaks.append(line.split())
            
    return Peptide(
        attributes['PEPMASS'],
        attributes['CHARGE'],
        attributes['SCANS'],
        attributes['RTINSECONDS'],
        attributes['SEQ'],
        Spectrum(peaks)
    )

def parse_peptide_file(file_path):
    peptides = []

    with open(file_path,"r") as f:
        inpep = False
        peplines = []
        for line in f:
            line = line.strip()
            if line == 'BEGIN IONS':
                inpep = True
            elif line == 'END IONS':
                inpep = False
                peptides.append(parse_peptide_lines(peplines))
                peplines = []
            elif inpep:
                peplines.append(line)
    
    return peptides