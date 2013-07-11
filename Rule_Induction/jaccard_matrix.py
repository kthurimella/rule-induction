#!/usr/bin/env python
# File created on 11 Jan 2013
from __future__ import division

__author__ = "Kumar Thurimella"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Kumar Thurimella"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Kumar Thurimella"
__email__ = "kthurimella@gmail.com"
__status__ = "Release"
 


from qiime.util import parse_command_line_parameters, make_option
import re
import itertools
from biom.parse import parse_biom_table
from biom.table import *
import numpy as np

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
    make_option('-r','--rules_input',type="string",help='the rules input'),
    make_option('-t','--rules_input_2',type="string",help='the rules input 2'),
    make_option('-o','--output',type="string",help='the sample output')
]
script_info['optional_options'] = [\
 # Example optional option
 #make_option('-o','--output_dir',type="new_dirpath",help='the output directory [default: %default]'),\
]
script_info['version'] = __version__

def parseRules(rf):
    """Returns OTU ID's for sample comparison (LHS, RHS)"""
    line = rf.read().splitlines()
    k = []
    
    # Get rid of all white space
    for i in line:
        j = i.replace(' ', '')
        k.append(j)
    newline = ''.join(k)
    
    # Split up each rule based on decimals present at the end of each rule
    # Ex: 8  {59Root; k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__}                => {56Root; k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__}            0.7777778  1.0000000 1.000000
    # This would split once it reached the 0.7777778 and keep string before then.
    rules = re.findall(r'\{(.+?)\}', newline)
    #rules = re.split("\d\.[\d]+", newline)

    new_rules = []
    for rule in rules:
  new_rules.append(re.findall(r'\d+', rule))

    final_rules = map (lambda x: new_rules[2*x:(x+1)*2], range (int(len(new_rules)/2)))

    explode_list = [item for t in final_rules for item in t]

    newlist = []

    for i in xrange(len(explode_list)):
        if i % 2 == 0:
            list1 = explode_list[i]
            list2 = explode_list[i + 1]
	    newlist.append(list1 + list2)

    intlist = []

    for i in newlist:
	intlist.append(map(int, i))

    return intlist

def compute_jaccard_index(list_1, list_2):
    set_1 = set(list_1)
    set_2 = set(list_2)
    n = len(set_1.intersection(set_2))
    return n / float(len(set_1) + len(set_2) - n) 

def compute_jaccard_matrix(setlist1, setlist2):
    m = len(setlist1)
    n = len(setlist2)

    answer = []

    for i in xrange(m):
	answer.append([])
	for j in xrange(n):
	    answer[i].append(compute_jaccard_index(setlist1[i], setlist2[j]))

#    for i, j in itertools.product(xrange(m), xrange(n)):
#    for x in np.nditer(jaccard_matrix, op_flags=['readwrite']):
#	x[...] = answer

    return answer

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    input_rf = open(opts.rules_input, 'r')
    input_rf1 = open(opts.rules_input_2, 'r')
    output_sf = open(opts.output, 'w')

    OtuIds1 = parseRules(input_rf)
    OtuIds2 = parseRules(input_rf1)

    jaccard_matrix = compute_jaccard_matrix(OtuIds1, OtuIds2)

    with output_sf as file:
	file.writelines('\t'.join(str(j) for j in i) + '\n' for i in jaccard_matrix)

    input_rf.close()
    output_sf.close()

if __name__ == "__main__":
    main()
