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
import re, operator
from biom.parse import parse_biom_table
from biom.table import *
from qiime.parse import parse_mapping_file_to_dict

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
    make_option('-r','--rules_input',type="string",help='the rules input'),
    make_option('-b','--biom_input',type="string",help='the biom input'),
    make_option('-o','--output',type="string",help='the sample output'),
    make_option('-m','--mapping_input',type="string",help='the mapping file'),
    make_option('-c','--category',type="string",help='the category of interest for the mapping file')
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
    rules = re.split("\d\.[\d]+", newline)



    listOfOtuIds = []
    for rule in rules:
        # This goes inside each split rule and finds all the numbers within it per rule using findall
	# NOTE: The anchor for this regex seems to be fishy within data sets. Need to find
	# a universal one. Either it's Root or k__
        temp = re.findall(r"(\d+)k", rule)
        if temp != []:
            listOfOtuIds.append(temp)
    return listOfOtuIds

def interpretBiom(bf, mf, c, OTUIds):
    """Returns a list comprising of rules and the samples with each rule"""
    biom_file = parse_biom_table(bf)
    mapping_file = parse_mapping_file_to_dict(mf)
    mapping_file = mapping_file[0]

    category_dict = dict( [ ( key, val [ c ] ) for ( key, val ) in mapping_file.iteritems() ] )
    sorted_category_dict = sorted(category_dict.iteritems(), key = operator.itemgetter(1))

    print sorted_category_dict

    samp_ids = []
    for vals, ids, md in biom_file.iterSamples():
        samp_ids.append(ids)
    
    samples_present = []
    final_list = []
    count = 0
    counter = 0

    # This takes in the list of OTU ID's and matches them with ID from
    # the biom file using the getValueByIds. If it isn't 0 then keep track
    # of it (i.e. the OTU is present in that sample) and do an intersection
    # between all of the said OTU's within each sample. Unfortunately, matching
    # is O(n^2) no matter what.
    for j in OTUIds:
	for id in samp_ids:
	    for k in j:
	    	if int(biom_file.getValueByIds(k, id)) != 0:
                    count = count + 1
            if count == len(j):
                samples_present.append(id)
            count = 0
	    if id == 
        counter = counter + 1
        final_list.append(counter)
	# temporary hack: used the set function to make a unique list 
	# I should clear the list after each iteration through the OTUIds
	# but it somehow clears the entire list even if I append it before.
        final_list.append(set(samples_present))

    # this overcomes the temporary hack and converts from set to list
    for i in xrange(len(final_list)):
	if i % 2 != 0:
	    final_list[i] = list(final_list[i])


    return final_list


		    
#for key in mf:
#    new_mf[key] = mf[key]['OBESITYCAT']
# colors = dict( [ ( key, val [ 'OBESITYCAT' ] ) for ( key, val ) in original.iteritems() ] )


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    input_rf = open(opts.rules_input, 'r')
    input_bf = open(opts.biom_input, 'U')
    output_sf = open(opts.output, 'w')
    input_mf = open(opts.mapping_input, 'U')
    input_category = opts.category

    OtuIds = parseRules(input_rf)
    sample_ids = interpretBiom(input_bf, input_mf, input_category, OtuIds)
	
    output_sf.write('\t Rule Number \n \t Samples Present \n')

    for item in sample_ids:
	output_sf.write("%s\n" % item)

    input_rf.close()
    input_bf.close()
    output_sf.close()

if __name__ == "__main__":
    main()
