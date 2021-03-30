#
import pickle

def read_sed_pickle(filename):
    infile = open(filename, 'rb')
    data = pickle.load(infile)
    infile.close()
    return data['wavelengths'], data['fluxes']
