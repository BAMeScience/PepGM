import pandas as pd 
import time
import argparse
import difflib

parser = argparse.ArgumentParser(description = 'Filter host spectra from mgf or mzML file')
#PSMReport = '/home/tholstei/repos/PepGM_all/PepGM/results/PXD002936_avian_bronchitis/chicken_refseq_Default_PSM_Report.txt'
#PSMReport ='/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD003013_Cowpox_BR/human_refseq_Default_PSM_Report.txt' 
#SpectrumFile ='/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD002936_avian_bronchitis/BeauR2.mgf'
#SpectrumFile = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD003013_Cowpox_BR/PXD003013_Cowpox_BR.mzML'
#out = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD002936_avian_bronchitis/PXD002936_avian_bronchitis_filtered_test.mgf'
#out = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD003013_Cowpox_BR/PXD003013_Cowpox_BR_filtertest.mzML'

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



def filterMzML(SpectraToFilter,SpectrumFile,output):

    '''
    filter spectra with title provided in a list from .mgf file
    input: list of spectrum titles, .mzML file of spectra
    output: filtered .mzML file
    '''
    SpectraToFilter = set(SpectraToFilter)
    LinesToWrite = []
    Spectra = open(SpectrumFile,'r')
    cleanedSpectra = open(output,'w')
    lastline = ''

    writeLine = True
    
    for line in Spectra:
            
        if line.find('<spectrum') != -1:

            if any(title in line for title in SpectraToFilter):
                writeLine = False

            else:
                writeLine = True        
                              
        if writeLine:
            LinesToWrite.append(lastline)
        
        lastline = line
        
    print('filtered '+str(len(SpectraToFilter))+ 'host or crap spectra')
    cleanedSpectra.writelines(LinesToWrite)

    Spectra.close()
    cleanedSpectra.close()


def FilterMGF(SpectraToFilter,SpectrumFile,output):
    
    '''
    filter spectra with title provided in a list from .mgf file
    input: list of spectrum titles, .mgf file of spectra
    output: filtered .mgf file
    '''
    
    SpectraToFilter = set(SpectraToFilter)
    LinesToWrite = []
    Spectra = open(SpectrumFile,'r', newline='\n')
    cleanedSpectra = open(output,'w')
    lastline = ''

    writeLine = True

    for line in Spectra:
         
        lineNoN = line.rstrip()
            
        if line.find('TITLE') != -1:

            if lineNoN[6:] in SpectraToFilter:
                writeLine = False
            else:
                writeLine = True        
                              
        if writeLine:
            LinesToWrite.append(lastline)
        
        lastline = line
        
    print('filtered '+str(len(SpectraToFilter))+ 'host or crap spectra')
    cleanedSpectra.writelines(LinesToWrite)

    Spectra.close()
    cleanedSpectra.close()



FilterList = getSpectraNames(args.PSMReport)

if args.SpectrumFile[-1] == 'L':
    filterMzML(FilterList,args.SpectrumFile,args.out)
else:
    FilterMGF(FilterList,args.SpectrumFile,args.out)




