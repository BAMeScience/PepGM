import pandas as pd 
import time
import argparse

parser = argparse.ArgumentParser(description = 'Filter host spectra from mgf or mzML file')
#filepath = '/home/tholstei/repos/PepGM_all/PepGM/results/PXD002936_avian_bronchitis/chicken_refseq_Default_PSM_Report.txt'
#SpectrumFile ='/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD002936_avian_bronchitis/PXD002936_avian_bronchitis.mgf'
#SpectrumFile = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD003013_Cowpox_BR/PXD003013_Cowpox_BR.mzML'
#output = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD002936_avian_bronchitis/PXD002936_avian_bronchitis.mgf'

parser.add_argument('--SpectrumFile',type = str, required = True, help = 'path to your spectrum file')
parser.add_argument('--PSMReport', type = str, required = True, help = 'path to PeptideShaker simple psm report file')
parser.add_argument('--out',type = str, required = True, help= 'output filepath')


args = parser.parse_args()

def getSpectraNames(filepath):

    '''
    get spectrum titles from PeptideShaker output
    input: path to spectrum file
    returns: list of spectra tiltes
    '''

    raw = pd.read_csv(filepath, sep = '\t', error_bad_lines=False, usecols=['Spectrum Title'])
    spectrumNames = raw['Spectrum Title'].to_list()
    return spectrumNames


def FilterMGF(SpectraToFilter,SpectrumFile,output):
    
    '''
    filter spectra with title provided in a list from .mgf file
    input: list of spectrum titles, .mgf file of spectra
    output: filtered .mgf file
    '''

    with open(SpectrumFile) as Spectra:

        with open(output, mode ='w') as cleanedSpectra:

            lastline = None
            writeLine = False
            
            for line in Spectra:
                
                if any(title in  line for title in SpectraToFilter):
                    writeLine = False
                elif line != "END IONS":
                    writeLine = False
                else:
                    writeLine = True

                if writeLine:
                    cleanedSpectra.write(lastline)

                lastline= line

def filterMzML(SpectraToFilter,SpectrumFile,output):

    '''
    filter spectra with title provided in a list from .mgf file
    input: list of spectrum titles, .mzML file of spectra
    output: filtered .mzML file
    '''

    with open(SpectrumFile) as Spectra:

        with open(output, mode ='w') as cleanedSpectra:

            lastline = None
            writeLine = False
            
            for line in Spectra:
                
                if any(title in  line for title in SpectraToFilter):
                    writeLine = False
                elif line != "</spectrum>":
                    writeLine = False
                else:
                    writeLine = True

                if writeLine:
                    cleanedSpectra.write(lastline)

                lastline= line


FilterList = getSpectraNames(args.PSMReport)

if args.filepath[-1] == 'L':
    filterMzML(FilterList,args.SpectrumFile,args.out)
else:
    FilterMGF(FilterList,args.SpectrumFile,args.out)




