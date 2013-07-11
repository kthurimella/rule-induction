#!/usr/bin/env python

__author__ = "Kumar Thurimella"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Kumar Thurimella"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Kumar Thurimella"
__email__ = "kthurimella@gmail.com"
__status__ = "Release"

import subprocess
from os import remove, path, devnull
from os.path import join, split, splitext
from sys import stdout
from time import sleep
from tempfile import mkdtemp
from numpy import array, set_printoptions, nan, asarray
from cogent.app.util import CommandLineApplication, CommandLineAppResult, FilePath, ResultPath, ApplicationError
from qiime.util import get_tmp_filename
from qiime.util import get_qiime_project_dir
from qiime.format import format_otu_table
from cogent.app.parameters import Parameters
from qiime.parse import parse_otu_table, parse_mapping_file
from biom.parse import convert_biom_to_table, parse_biom_table
from qiime.util import parse_command_line_parameters, make_option
from biom.table import *

def delimitedSelf(self, delim='\t', header_key=None, header_value=None, 
    metadata_formatter=str):
    """Return self as a string in a delimited form
    
    Default str output for the Table is just row/col ids and table data
    without any metadata

    Including observation metadata in output: If ``header_key`` is not 
    ``None``, the observation metadata with that name will be included 
    in the delimited output. If ``header_value`` is also not ``None``, the 
    observation metadata will use the provided ``header_value`` as the 
    observation metadata name (i.e., the column header) in the delimited 
    output. 
    
    ``metadata_formatter``: a function which takes a metadata entry and 
    returns a formatted version that should be written to file
    """
    if self.isEmpty():
        raise TableException, "Cannot delimit self if I don't have data..."

    samp_ids = delim.join(map(str, self.SampleIds))
    # 17hrs of programing straight later...
    if header_key is not None:
        if header_value is None:
            raise TableException, "You need to specify both header_key and header_value"
    if header_value is not None:
        if header_key is None:
            raise TableException, "You need to specify both header_key and header_value"

    if header_value:
        output = ['# Constructed from biom file','#OTU ID%s%s' % (delim, 
            samp_ids)]
    else:
        output = ['# Constructed from biom file','#OTU ID%s%s' % (delim, 
            samp_ids)]

    for obs_id, obs_values in zip(self.ObservationIds, self._iter_obs()):
        str_obs_vals = delim.join(map(str, self._conv_to_np(obs_values)))

        if header_key and self.ObservationMetadata is not None:
            md = self.ObservationMetadata[self._obs_index[obs_id]]
            md_out = metadata_formatter(md.get(header_key,None))
            output.append('%s%s%s%s' % (obs_id, md_out, delim, str_obs_vals))
        else:
            output.append('%s%s%s' % (obs_id, delim, str_obs_vals))
            
    return '\n'.join(output)



def taxonomy_biom_to_table(biom_f, header_key=None, header_value=None, \
        md_format=None):
    """Convert a biom table to a contigency table"""
    table = parse_biom_table(biom_f)

    if md_format is None:
        md_format = lambda x: '; '.join(x)
    
    if table.ObservationMetadata is None:
        return delimitedSelf(table)
    
    if header_key in table.ObservationMetadata[0]:
        return delimitedSelf(table, header_key=header_key, 
                                       header_value=header_value,
                                       metadata_formatter=md_format)
        
    else:
        return delimitedSelf(table)

class RruleInduction(CommandLineApplication):
    """R Supervised Learner application controller
       Runs R with a source script (from qiime/support_files/R), and 
       passes in an OTU table and mapping file. Causes R to run a supervised
       classifier to predict labels from a given category from the mapping file
       using the provided OTUs.
    """
    _input_handler = '_input_as_path'
    _command = "R"
    _options ={}

    _R_parameters = {
        'flags': '--slave --vanilla'
        }

    # The name of the R script (located under qiime/support_files/R/)
    _R_script = 'print_rules.r'

    _parameters = {}
    _parameters.update(_options)
    _parameters.update(_R_parameters)

    def __call__(self, predictor_fp, confidence, support, sort = "lift", taxonomy = "no", output_dir='/scratch/thurimellak/', verbose=False):
        """Run the application with the specified kwargs on data

            remove_tmp: if True, removes tmp files
            
            returns a CommandLineAppResult object
        """
        suppress_stdout = self.SuppressStdout
        suppress_stderr = self.SuppressStderr
        if suppress_stdout:
            outfile = devnull
        else:
            outfilepath = FilePath(self.getTmpFilename(self.TmpDir))
            outfile = open(outfilepath,'w')
        if suppress_stderr:
            errfile = devnull
        else:
            errfilepath = FilePath(self.getTmpFilename(self.TmpDir))
            errfile = open(errfilepath, 'w')
        if output_dir is None:
            output_dir = mkdtemp(prefix='R_output_')
        
        ## temporary hack: this converts a biom file to classic otu table
        ##  format for use within R
        if verbose:
            print 'Converting BIOM format to tab-delimited...'

        temp_predictor_fp = join(output_dir,
                                 splitext(split(predictor_fp)[1])[0]+'.txt')
        temp_predictor_f = open(temp_predictor_fp,'w')
        if taxonomy == "no":
      temp_predictor_f.write(convert_biom_to_table(open(predictor_fp,'U')))
        else:
            temp_predictor_f.write(taxonomy_biom_to_table(open(predictor_fp, 'U'), 'taxonomy', 'taxonomy'))
        temp_predictor_f.close()
        predictor_fp = temp_predictor_fp

        rflags = self.RParameters['flags']
        rscript = self._get_R_script_path()
        base_command = self._get_base_command()
        cd_command, base_command = base_command.split(';')
        cd_command += ';'
        source_dir = self._get_R_script_dir()
        
        # Build up the command, consisting of a BaseCommand followed byp

        # input and output (file) specifications
        args = ['-i', predictor_fp, 
		'-c', confidence,
		'-s', support,
		'-o', sort,
		'--source_dir', source_dir]
        if verbose:
                args += ['-v']

        command = self._commandline_join(
            [   cd_command, base_command,
                '--args'
            ] + args + [' < %s ' %(rscript)]
            ) 
	print "Command: ", command
        if self.HaltExec: 
            raise AssertionError, "Halted exec with command:\n" + command

        # run command, wait for output, get exit status
        proc = subprocess.Popen(command, shell=True, stdout=outfile, stderr=errfile)

        if verbose:
            print '\nR output\n'
            tmpoutfile = open(outfilepath,'U')
            while proc.poll() is None:
                stdout.write(tmpoutfile.readline())
                sleep(0.00001)
            tmpoutfile.close()
        proc.wait()
        exit_status = proc.returncode

        # Determine if error should be raised due to exit status of 
        # appliciation
        if not self._accept_exit_status(exit_status):
            if exit_status == 2:
                raise ApplicationError, \
                    'R library not installed: \n' + \
                    ''.join(open(errfilepath,'r').readlines()) + '\n'
            else:
                raise ApplicationError, \
                    'Unacceptable application exit status: %s, command: %s'\
                    % (str(exit_status),command) +\
                    ' Program output: \n\n%s\n'\
                     %(''.join(open(errfilepath,'r').readlines()))

    def _get_result_paths(self, output_dir):
        """Returns the filepaths for all result files"""
        files = {
            'Transpose_File': 'min_sparse_otu_table.txt.txt',
            'Eclat_Rules': 'eclat_rules.txt',
            'Apriori_Rules': 'apriori_rules.txt',
            'logfile': 'logfile.txt',
        }
        result_paths = {}
        for name, file in files.iteritems():
            result_paths[name] = ResultPath(
                Path=path.join(output_dir, file), IsWritten=True)
        return result_paths

    def _get_R_script_dir(self):
        """Returns the path to the qiime R source directory
        """
        qiime_dir = get_qiime_project_dir()
        script_dir = path.join(qiime_dir,'qiime','support_files','R')
        return script_dir

    def _get_R_script_path(self):
        """Returns the path to the R script to be executed
        """
        return path.join(self._get_R_script_dir(), self._R_script)

    def _commandline_join(self, tokens):
        """Formats a list of tokens as a shell command
        """
        commands = filter(None, map(str, tokens))
        return self._command_delimiter.join(commands).strip()

    def _accept_exit_status(self,exit_status):
        """ Return False to raise an error due to exit_status !=0 of application
        """
        if exit_status != 0:
            return False
        return True

    @property
    def RParameters(self):
        return self.__extract_parameters('R')

    def __extract_parameters(self, name):
        """Extracts parameters in self._<name>_parameters from self.Parameters

        Allows the program to conveniently access a subset of user-
        adjusted parameters, which are stored in the Parameters
        attribute.
        
        Relies on the convention of providing dicts named according to
        "_<name>_parameters" and "_<name>_synonyms".  The main
        parameters object is expected to be initialized with the
        contents of these dicts.  This method will throw an exception
        if either convention is not adhered to.
        """
        parameters = getattr(self, '_' + name + '_parameters')
        result = Parameters(parameters)
        for key in result.keys():
            result[key] = self.Parameters[key]
        return result
    
script_info = {}
script_info['brief_description'] = "Rule Induction"
script_info['script_description'] = "This script uses the arules package in R to induce rules on datasets"
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [\
make_option('-i', '--input_data', type='existing_filepath',
            help='Input data file containing predictors (e.g. otu table)',),
make_option('-s', '--support', type='float', help='Support for the given OTU Table'),
make_option('-c', '--confidence', type='float', help='Confidence for the given OTU Table')]

script_info['optional_options'] = [\
make_option('-o', '--sort', type= 'string', default ="lift", help='Sort rules based on measure of choice'),
make_option('-t', '--taxonomy', type='string', default ="no", help='Return Taxonomy String with Rules.')]
script_info['version'] = __version__



def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # run the rule induction algorithm
    learner = RruleInduction()
    learner(opts.input_data, opts.confidence, opts.support, sort = opts.sort, taxonomy = opts.taxonomy)
if __name__ == "__main__":
    main()
