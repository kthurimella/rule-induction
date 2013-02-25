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

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
    make_option('-r','--temp_input',type="existing_filepath",help='the rules input'),
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
    for i in line:
        j = i.replace(' ', '')
        k.append(j)
    newline = ''.join(k)
    rules = re.split("\d\.[\d]+", newline)
    listOfOtuIds = []
    for rule in rules:
        temp = re.findall(r"(\d+)R", rule)
        if temp != []:
            listOfOtuIds.append(temp)
    return listOfOtuIds

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    input_rf = open(opts.temp_input, 'r+')

    listOfIds = parseRules(input_rf)
    print listOfIds

if __name__ == "__main__":
    main()
