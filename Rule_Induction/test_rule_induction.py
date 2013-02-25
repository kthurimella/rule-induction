#!/usr/bin/env python
# File created on 02 Jan 2013
from __future__ import division

__author__ = "Kumar Thurimella"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Kumar Thurimella"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Kumar Thurimella"
__email__ = "kthurimella@gmail.com"
__status__ = "Release"
 

from cogent.util.unit_test import TestCase, main

from rule_induction import transpose_biom_to_table
from biom.parse import parse_biom_table

"""Biom table 1 is a dense table """
biom1 = '{"rows": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_4", "metadata": null}, {"id": "GG_OTU_5", "metadata": null}], "format": "Biological Observation Matrix 1.0.0", "data": [[0,0,1,0,0,0], [5,1,0,2,3,1], [0,0,1,4,2,0], [2,1,1,0,0,1], [0,1,1,0,0,0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}, {"id": "Sample5", "metadata": null}, {"id": "Sample6", "metadata": null}], "generated_by": "BIOM-Format 1.0.0c", "matrix_type": "dense", "shape": [5, 6], "format_url": "http://biom-format.org", "date": "2012-09-27T12:33:33.854501", "type": "OTU table", "id": null, "matrix_element_type": "float"}'
"""Biom table 2 is a sparse table """
biom2 = '{"rows": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_4", "metadata": null}, {"id": "GG_OTU_5", "metadata": null}], "format": "Biological Observation Matrix 1.0.0", "data": [[0,2,1], [1,0,5], [1,1,1], [1,3,2], [1,4,3], [1,5,1], [2,2,1], [2,3,4], [2,5,2], [3,0,2], [3,1,1], [3,2,1], [3,5,1], [4,1,1], [4,2,1]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}, {"id": "Sample5", "metadata": null}, {"id": "Sample6", "metadata": null}], "generated_by": "BIOM-Format 1.0.0c", "matrix_type": "sparse", "shape": [5, 6], "format_url": "http://biom-format.org", "date": "2012-09-25T16:04:34.723468", "type": "OTU table", "id": null, "matrix_element_type": "float"}'
"""Biom table 4 is a table full of digits (sparse)"""
biom4 = '{"rows": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_4", "metadata": null}, {"id": "GG_OTU_5", "metadata": null}], "format": "Biological Observation Matrix 0.9.3", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [0, 3, 4.0], [0, 4, 5.0], [0, 5, 6.0], [1, 0, 5.0], [1, 1, 1.0], [1, 2, 4.0], [1, 3, 2.0], [1, 4, 3.0], [1, 5, 1.0], [2, 0, 3.0], [2, 1, 2.0], [2, 2, 1.0], [2, 3, 4.0], [2, 4, 2.0], [2, 5, 5.0], [3, 0, 2.0], [3, 1, 1.0], [3, 2, 1.0], [3, 3, 3.0], [3, 4, 5.0], [3, 5, 1.0], [4, 0, 2.0], [4, 1, 1.0], [4, 2, 1.0], [4, 3, 5.0], [4, 4, 7.0], [4, 5, 8.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}, {"id": "Sample5", "metadata": null}, {"id": "Sample6", "metadata": null}], "generated_by": "BIOM-Format 0.9.3", "matrix_type": "sparse", "shape": [5, 6], "format_url": "http://biom-format.org", "date": "2012-09-27T13:10:46.201528", "type": "OTU table", "id": null, "matrix_element_type": "float"}'



class ruleInductionTests(TestCase):
    
    def setUp(self):
        pass
          
    def test_ruleInduction(self):
        """tests that rule induction runs correctly"""
        biom_output1 = transpose_biom_to_table(biom1)
        self.assertEqual(biom_output1, '# Constructed from biom file\n#OTU ID GG_OTU_1 GG_OTU_2 GG_OTU_3 GG_OTU_4 GG_OTU_5\nSample1\t\t0.0\t5.0\t0.0\t2.0\t0.0\nSample2\t\t0.0\t1.0\t0.0\t1.0\t1.0\nSample3\t\t1.0\t0.0\t1.0\t1.0\t1.0\nSample4\t\t0.0\t2.0\t4.0\t0.0\t0.0\nSample5\t\t0.0\t3.0\t2.0\t0.0\t0.0\nSample6\t\t0.0\t1.0\t0.0\t1.0\t0.0', but expected '# Constructed from biom file #OTU ID GG_OTU_1 GG_OTU_2 GG_OTU_3 GG_OTU_4 GG_OTU_5 Sample1\t\t0.0\t5.0\t0.0\t2.0\t0.0 Sample2\t\t0.0\t1.0\t0.0\t1.0\t1.0 Sample3\t\t1.0\t0.0\t1.0\t1.0\t1.0 Sample4\t\t0.0\t2.0\t4.0\t0.0\t0.0 Sample5\t\t0.0\t3.0\t0.0\t0.0\t0.0 Sample6\t\t0.0\t1.0\t2.0\t1.0\t0.0')
        
#        biom_output2 = relativeAbundances(biom2).getBiomFormatJsonString("a")
#        self.assertEqual(biom_output2, '{"rows": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_4", "metadata": null}, {"id": "GG_OTU_5", "metadata": null}], "format": "Biological Observation Matrix 1.0.0", "data": [[0.0, 0.0, 0.25, 0.0, 0.0, 0.0], [0.714, 0.333, 0.0, 0.333, 0.6, 0.5], [0.0, 0.0, 0.25, 0.667, 0.4, 0.0], [0.286, 0.333, 0.25, 0.0, 0.0, 0.5], [0.0, 0.333, 0.25, 0.0, 0.0, 0.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}, {"id": "Sample5", "metadata": null}, {"id": "Sample6", "metadata": null}], "generated_by": "BIOM-Format 1.0.0c", "matrix_type": "dense", "shape": [5, 6], "format_url": "http://biom-format.org", "date": "2012-09-27T12:33:33.854501", "type": "OTU table", "id": null, "matrix_element_type": "float"}')
#        
#        #biom_output3 = relativeAbundances(biom3)
#        #self.assertRaises(TypeError, biom_output3)
#            
#        biom_output4 = relativeAbundances(biom4).getBiomFormatJsonString("a")
#        self.assertEqual(biom_output4, '{"rows": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_4", "metadata": null}, {"id": "GG_OTU_5", "metadata": null}], "format": "Biological Observation Matrix 1.0.0", "data": [[0, 0, 0.077], [0, 1, 0.286], [0, 2, 0.3], [0, 3, 0.222], [0, 4, 0.227], [0, 5, 0.286], [1, 0, 0.385], [1, 1, 0.143], [1, 2, 0.4], [1, 3, 0.111], [1, 4, 0.136], [1, 5, 0.048], [2, 0, 0.231], [2, 1, 0.286], [2, 2, 0.1], [2, 3, 0.222], [2, 4, 0.091], [2, 5, 0.238], [3, 0, 0.154], [3, 1, 0.143], [3, 2, 0.1], [3, 3, 0.167], [3, 4, 0.227], [3, 5, 0.048], [4, 0, 0.154], [4, 1, 0.143], [4, 2, 0.1], [4, 3, 0.278], [4, 4, 0.318], [4, 5, 0.381]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}, {"id": "Sample5", "metadata": null}, {"id": "Sample6", "metadata": null}], "generated_by": "BIOM-Format 1.0.0c", "matrix_type": "sparse", "shape": [5, 6], "format_url": "http://biom-format.org", "date": "2012-09-27T13:31:34.758279", "type": "OTU table", "id": null, "matrix_element_type": "float"}')
    


if __name__ == "__main__":
    main()
