#!/usr/bin/env python
# File created on 09 Mar 2013
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
from biom.table import *
from biom.parse import *
import numpy as np

script_info = {}
script_info['brief_description'] = "This program filters out OTU's based on its standard deviation and average."
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [make_option('-i', '--input', type = "string", help = 'the input file in .biom format'), 
make_option('-o', '--output', type = "string", help = 'the output filepith with a designated name'),
make_option('-t', '--thres', type = "float", help = "the threshold at which to filter by")
]
script_info['optional_options'] = [\
 # Example optional option
 #make_option('-o','--output_dir',type="new_dirpath",help='the output directory [default: %default]'),\
]
script_info['version'] = __version__

def filterBySample(bt, filterThreshold):

    def itself(o_val, o_id, o_md):
        row_sum = o_val.sum()
        filterMeasure = filterThreshold * row_sum
        for i in np.nditer(o_val, op_flags=['readwrite']):
            if i > filterMeasure:
                i[...] = 1
            else:
                i[...] = 0
        return o_val

    return bt.transformSamples(itself)


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    input_f = open(opts.input, 'U')
    output_f = open(opts.output, 'w')
    threshold = opts.thres

    x = parse_biom_table(input_f)

    new_bt = filterBySample(x, threshold)

    output_f.write(new_bt.getBiomFormatJsonString(generatedby()))

    input_f.close()
    output_f.close()

if __name__ == "__main__":
    main()
