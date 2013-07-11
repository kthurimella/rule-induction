#!/usr/bin/env python
# File created on 19 Sep 2012
from __future__ import division

__author__ = "Kumar Thurimella"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Kumar Thurimella"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Kumar Thurimella"
__email__ = "kthurimella@gmail.com"
__status__ = "Release"
 


from qiime.util import parse_command_line_parameters, make_option
from biom.table import *
from biom.parse import *
import numpy as np

script_info = {}
script_info['brief_description'] = "This program calculates the relative abundances of different OTU's."
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= "A .biom file that will contain relative abundances."
script_info['required_options'] = [make_option('-i','--input',type="string", help='the input file in .biom format'),
                                   make_option('-o','--output',type="string",help='the output filepath with a designated name'),
                                   make_option('-d','--density', type="float", help='the density of the table you want')]
script_info['version'] = __version__


def randomlyPopulate(bt, density):

    def itself(value, ID, MetaData):
        
        for i in np.nditer(value, op_flags=['readwrite']):
            randNum = np.random.random_sample()
            if randNum < density:
                i[...] = 1
            else:
                i[...] = 0
 
        return value

    return bt.transformSamples(itself)

def main():
    
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    input_f = open(opts.input, 'U')
    output_f = open(opts.output, 'w')
    density = opts.density
    x = parse_biom_table(input_f)  
    
    new_bt = randomlyPopulate(x, density)

    output_f.write(new_bt.getBiomFormatJsonString(generatedby()))

    input_f.close()
    output_f.close()

if __name__ == "__main__":
    main()
