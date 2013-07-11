#!/usr/bin/env python
# File created on 14 Mar 2013
from __future__ import division

__author__ = "Kumar Thurimella"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Kumar Thurimella"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Kumar Thurimella"
__email__ = "kthurimella@gmail.com"
__status__ = "Development"


from qiime.util import parse_command_line_parameters, make_option
from biom.parse import *
from biom.table import *
import numpy as np

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [make_option('-i','--input',type="string", help='the input file in .biom format'),
                                   make_option('-r','--iterations', type="int", help='the number of iterations'),
                                   make_option('-o','--output', type="string", help='the output of this file')]
script_info['version'] = __version__

def randomlyShuffle(bt):

    def itself(value, ID, MetaData):
        
        np.random.shuffle(value)
        
        return value

    return bt.transformSamples(itself)

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    input_f = open(opts.input, 'U')
    x = parse_biom_table(input_f)  
    iterations = opts.iterations


    i = 0
    while i < iterations:
        output_f = open(opts.output+"_"+str(i)+".biom", 'w')
        new_bt = randomlyShuffle(x)
        output_f.write(new_bt.getBiomFormatJsonString(generatedby()))
        i = i+1

    input_f.close()
    output_f.close()


if __name__ == "__main__":
    main()
