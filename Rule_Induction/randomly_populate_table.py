#!/usr/bin/env python
# File created on 19 Sep 2013
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
import random as rd

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

    return bt.transformObservations(itself)

def random(percent, width, height):
    n = width * height
    t = int(round(percent*n))
    y = [1 for x in xrange(t)]
    y = y + [0 for x in xrange(n-t)]
    rd.shuffle(y)
    matrix = []
    for i in xrange(height):
        temp = []
        for j in xrange(width):
            temp.append(y[width*i+j])
        matrix.append(temp)
    return matrix

def main():
    
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    input_f = open(opts.input, 'U')
    output_f = open(opts.output, 'w')
    density = opts.density


    obs_ids = []
    samp_ids = []

    bt = parse_biom_table(input_f)

    for obs_v, obs_id, obs_md in bt.iterObservations():
  obs_ids.append(obs_id)

    for samp_v, samp_id, samp_md in bt.iterSamples():
	samp_ids.append(samp_id)

    matrix_rand = random(density, 20, 101)
    numpy_matrix = np.array(matrix_rand)

    new_bt = table_factory(numpy_matrix, samp_ids, obs_ids)

    print(new_bt.getTableDensity())
    
#    new_bt = randomlyPopulate(x, density)


    output_f.write(new_bt.getBiomFormatJsonString(generatedby()))

    input_f.close()
    output_f.close()

if __name__ == "__main__":
    main()
