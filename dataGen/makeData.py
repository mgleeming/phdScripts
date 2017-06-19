import os
import re
import sys
import pymzml
import random
import argparse
import numpy as np
from pyteomics import mass, fasta, parser, mgf
from subprocess import Popen, PIPE, STDOUT

try:
    import better_exceptions
    print 'better_exceptions imported'
except:
    print 'better_exceptions import failed'

argParser = argparse.ArgumentParser(description = 'Generate LC-MS/MS data for arbitrary CRM adducted peptides')

argParser.add_argument('--fasta',
					help = 'Input FATA file',
					nargs = '+')
argParser.add_argument('--CRM',
					help = 'Input Formula for desired CRM',
					type = str)
argParser.add_argument('--targetResidue',
					help = 'Amino acid site modified by input CRM',
					type = str)
argParser.add_argument('--targetProteins',
                                        help = 'FASTA sequence file names to which CRM is applied',
                                        nargs = '+')
argParser.add_argument('--logFile',
					help = 'Output log file name',
					default = 'logFile.dat',
					type = str)
argParser.add_argument('--mzmlOut',
					help = 'Output mzML file name',
					default = 'output.mzML',
					type = str)
argParser.add_argument('--mgfOut',
					help = 'Output mgf file name',
					default = 'output.mgf',
					type = str)
argParser.add_argument('--mzDelta',
					help = 'Mass difference (Da) between heavy and light CRMs',
					default = 6.0201,
					type = float)


'''
Note:

For some reason, running MSSimulator using either subprocess.Popen or os.system crashes if pyopenms
has been previously imported.

Error msg is:
Warning: Unable to fetch subsection parameters! Addressing subsection parameters will not work for this tool.
Error: File not found (the file 'CHEMISTRY/Enzymes.xml' could not be found)

Not sure what's going on here.

Workaround - process all files and save temporary mzml files.
'''

class Spectrum(object):
	def __init__(self, rt, mzs, ints, original):
		self.rt = rt
		self.mzs = mzs
		self.ints = ints
		self.original = original

class Peptide(object):
    def __init__(self, sequence, protein):
        self.sequence = sequence
        self.protein = protein

def getPeptides(inFastas, minLength=5):
	''' Perform digestion on input FASTA sequence and return peptides
	'''

	peptides = set()

        def digest (fastaFile):
            entries = []

            for description, sequence in fasta.read(fastaFile):
		new_peptides = parser.cleave(sequence, parser.expasy_rules['trypsin'])
	        new_peptides = [Peptide(p, fastaFile.split('/')[-1]) for p in new_peptides if len(p) > minLength and p.count('K') + p.count('R') == 1]
                entries.extend(new_peptides)

            return entries

        if isinstance(inFastas, str):
            # single fast file
            new_peptides = digest(inFastas)
            peptides.update(new_peptides)

        elif isinstance(inFastas, list):
            for inFasta in inFastas:
                new_peptides = digest(inFasta)
                peptides.update(new_peptides)

        else:
            print 'Unknown fasta file specification'
            sys.exit()
        for p in peptides:
            print p.sequence, p.protein
        return peptides

def produceMGFEntry(peptide, pepRT, pepMass, CRMmass = None, target = None):
	frags = []

	if CRMmass is not None:
		pepMass += float(CRMmass) / 2

	for i in range(len(peptide.sequence)):
		yi = peptide.sequence[i:]
		mz = mass.fast_mass(yi, ion_type='y', charge=1)
		if all([CRMmass, target]) and target in yi:
			mz += CRMmass

		frags.append(mz)

	ints = [random.randint(1000, 5000) for _ in range(len(frags))]

	charges = [1 for _ in ints]

	pepHeader = '%s_%s_MOD' %(peptide.sequence, peptide.protein) if CRMmass else '%s_%s' %(peptide.sequence, peptide.protein)

        mgf_entry = {
		'm/z array': frags,
		'intensity array': ints,
		'charge array': charges,
		'params': {
			'TITLE': '%s_%s_%s' % (pepMass, pepRT, pepHeader),
			'PEPMASS': pepMass,
			'RTINSECONDS': pepRT,
			'CHARGE': '2+'
		}
	}
	return mgf_entry

def writeMGF(spectra, output):
	headers = {
		'COM': 'OpenMS_search',
		'USERNAME': 'OpenMS',
		'FORMAT': 'Mascot generic',
		'TOLU': 'Da',
		'ITOLU': 'Da',
		'FORMVER': '1.01',
		'DB': 'MSDB',
		'SEARCH': 'MIS',
		'REPORT': 'AUTO',
		'CLE': 'Trypsin',
		'MASS': 'monoisotopic',
		'INSTRUMENT': 'Default',
		'PFA': '1',
		'TOL': '3',
		'ITOL': '0.3',
		'TAXONOMY': 'All entries',
		'CHARGE': '1,2,3'
	}

	mgf.write(spectra=spectra, output=output)

	return

def combineSpectra(mzMLFiles, targetFiles, mzDelta, CRMmass):

	data = []

	for i, filei in enumerate(mzMLFiles):
		spectra = MZMLtoSpectrum(filei)

		if i == 0:
			for spectrum in spectra:
				data = [s for s in spectra]
				spectra = MZMLtoSpectrum(filei)

		for spectrum in spectra:
			for d in data:
				if d.rt == spectrum.rt:
					if i != 0:
						d.mzs = np.concatenate((d.mzs, spectrum.mzs))
						d.ints = np.concatenate((d.ints, spectrum.ints))

					if filei in targetFiles:
						# add in CRM_modified light and heavy peptides

						# light
						d.mzs = np.concatenate((d.mzs, spectrum.mzs + CRMmass / 2))
						d.ints = np.concatenate((d.ints, spectrum.ints))
						# heavy
						d.mzs = np.concatenate((d.mzs, spectrum.mzs + CRMmass / 2 + mzDelta))
						d.ints = np.concatenate((d.ints, spectrum.ints))

	return data

def writeCombinedMzML(data, outputFile):
	import pyopenms, copy

	output_file = pyopenms.MzMLFile()
	output_experiment = pyopenms.MSExperiment()

	for d in data:

		new_spectrum = copy.deepcopy(d.original)
		sortmap = np.argsort(d.mzs)
		new_spectrum.set_peaks((d.mzs[sortmap], d.ints[sortmap]))
		output_experiment.addSpectrum(new_spectrum)

	output_file.store(outputFile, output_experiment)
	return

def MZMLtoSpectrum(filename):
	spectra = []
	import pyopenms

	mzml_file = pyopenms.MzMLFile()
	experiment = pyopenms.MSExperiment()
	mzml_file.load(filename, experiment)
	for n,spectrum in enumerate(experiment):
		(mzData, intData) = spectrum.get_peaks()
		if mzData.shape[0] == 0:
			mzData = np.empty(0, dtype = 'float32')
			intData = np.empty(0, dtype = 'float32')
		try:
			time = spectrum.getRT()
		except KeyError, e:
			time_prev = time
			if delta_time > 0:
				time += delta_time
			else:
				time += 1.0
		#spectra.append(Spectrum(time, mzData, intData, spectrum))
		yield Spectrum(time, mzData, intData, spectrum)
	#return spectra

def main(options):

	'''
	NOTE:
	sys.path[0] = absolute path to top level script being executed
	os.getcwd() = absolute path to directory the script is called from
	'''

	wd = os.getcwd()
	paramDirectory = os.path.join(sys.path[0], 'params')

	logF = open(os.path.join(wd, options.logFile),'wt')
	logF.write('# Fragment, entry, RT, m/z, modified?, sequence\n')

	print 'Parsing peptides'
	outputMGF = options.mgfOut
	outputMzML = options.mzmlOut
	iniFile = os.path.join(paramDirectory, 'dataGen.ini')
	targetResidue = options.targetResidue
	crmMod = options.CRM
	mzDelta = options.mzDelta
	CRMmass = mass.calculate_mass(formula=crmMod)
        targetProteins = options.targetProteins

	peptides = getPeptides(options.fasta)
	print 'Number of peptides is: %s' %len(peptides)

	# data for later
	mgf_specs = []
	data = []
	fastaFiles = []
	mzMLFiles = []
	targetFiles = []
	entry = 0
	native_peptides = 0
	modified_peptides = 0
	total_peptides = 0

	# create temp file if needed
	tmpFiles = os.path.join(sys.path[0], 'tempFiles')
	if not os.path.exists(tmpFiles):
		os.makedirs(tmpFiles)

	# remove mgf file if present - for testing
	try:
		os.remove(outputMGF)
		print 'Removed existing MGF file'
	except OSError:
		pass

	#rtStep = 1500 / len(peptides) 

	print 'Geneating MS features'
	for pn, peptide in enumerate(peptides):

		pepMass = mass.fast_mass(peptide.sequence, charge=2)
		pepRT = random.randint(500, 2000)
		#pepRT = 500 + rtStep * pn

		pepInt =  100000 # random.randint(10000, 100000)

		fasta = os.path.join(tmpFiles, 'temp_%s.fasta' %pn)
		mzml = os.path.join(tmpFiles, 'temp_%s.mzML' %pn)

		# write fasta files
		of1 = open(fasta, 'wt')
		of1.write('>seq %s [# intensity= %s, RT=%s #]\n' % (pn, pepInt, pepRT))
		of1.write('%s' % (peptide.sequence))
		of1.close()

		fastaFiles.append(fasta)
		mzMLFiles.append(mzml)

		# produce mzML file for fasta fragment
		p = Popen( ['MSSimulator' + ' -ini %s' %iniFile + ' -in %s' %fasta + ' -out %s' % mzml], stdout=PIPE, stderr=STDOUT, shell = True, cwd = paramDirectory)
		stdout, nothing = p.communicate()

		# generate fragments
		mgf_specs.append(produceMGFEntry(peptide, pepRT, pepMass))

		logF.write('%s, %s, %s, %s, %s, %s, %s\n' %(pn, entry, pepRT, pepMass, 0, peptide.sequence, peptide.protein))
		entry += 1
		total_peptides += 1

		# add modified peptide MGF if needed
		if targetResidue in peptide.sequence and peptide.protein in targetProteins:
			mgf_specs.append(produceMGFEntry(peptide, pepRT, pepMass, CRMmass = CRMmass, target = targetResidue))
			targetFiles.append(mzml)
			logF.write('%s, %s, %s, %s, %s, %s\n' %(pn, entry, pepRT, pepMass + CRMmass/2, 1, '%s_%s_MOD'%(peptide.sequence, peptide.protein)))
			entry += 1
			modified_peptides += 1
			total_peptides += 1
		else: native_peptides += 1

	# write MGF file
	print 'Writing MGF file'
	writeMGF(mgf_specs, os.path.join(wd, outputMGF))

	# combine mzML file data and write output
	print 'Combining mzML files'
	data = combineSpectra(mzMLFiles, targetFiles, mzDelta, CRMmass)

	print 'Writing final mzML file'
	outputmzML = os.path.join(wd, outputMzML)
	writeCombinedMzML(data, os.path.join(wd, outputmzML))

	# clean up
	rubbish = fastaFiles + mzMLFiles
	for f in rubbish:
		os.remove(f)
	os.rmdir(tmpFiles)

	logF.write('\n# Summary:\n')
	logF.write('#-----------------------\n')
	logF.write('# native peptides: %s\n' %native_peptides)
	logF.write('# modified peptides: %s\n' %modified_peptides)
	logF.write('# total peptides: %s\n' %total_peptides)
	logF.close()

	return

if __name__ == '__main__':
	options = argParser.parse_args()
	main(options)
