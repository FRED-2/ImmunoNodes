# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
from tempfile import NamedTemporaryFile

__author__ = 'mohr,schubert'

import os
import subprocess
import logging
import itertools
import pandas

from Fred2.Core import Allele, AExternal
import DistanceMatrices
from DistanceMatrix import DistanceMatrix


class Distance2Self(object):
    """
        Implements calulcation routine of distance to (self) peptides
        Calculate k closest distances of peptide to peptide set represented as trie

        All our matrices have the same ordering of letters.
        If you use a new matrix, pleas make sure to use the same ordering! Otherwise the tries have to be recomputed!
    """

    def __init__(self, _matrix, trie=None, saveTrieFile=False):
        self.__saveTrieFile = saveTrieFile
        self.__matrix = _matrix
        self.__trie = trie

        this_dir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        self.__externalPathDistanceCalculator = os.path.join(this_dir, 'compute_distances_ivac')
        self.__externalPathTrieGenerator =  os.path.join(this_dir, 'get_TrieArray')


    def __del__(self):
        if not self.__saveTrieFile:
            pass

    def generate_trie(self, fastaFile, outfile='peptideTrie', peptideLength=9):

        cmd = self.__externalPathTrieGenerator + " %s %s %s %s"
        specifiedTrie = outfile
        self.__trie = specifiedTrie

        subprocess.check_output(cmd%(fastaFile, self.__matrix.path_to_matrix_file, peptideLength, specifiedTrie), shell=True)

    def calculate_distances(self, peptides, pep_header="neopeptide", specifiedTrie="uniprot_proteome_l9", n=10):
      def __load_trie(trieSource):
         current = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
         return os.path.join(current,"data","tries","{}.trie".format(trieSource))

      # create temporary file with peptides for distance computation
      tmpFile = NamedTemporaryFile(delete=False)

      with open(tmpFile.name, "w") as peptidesFile:
         for pep in peptides:
            peptidesFile.write('%s\n' % pep)

      cmd = self.__externalPathDistanceCalculator + " %s %s %s %s"
      results = {}

      trie = specifiedTrie if os.path.isfile(specifiedTrie) else __load_trie(specifiedTrie)
      method = os.path.basename(specifiedTrie).split('.')[0] if os.path.isfile(specifiedTrie) else specifiedTrie
      try:
         re = self.parse_external_result(
            subprocess.check_output(cmd % (self.__matrix.path_to_matrix_file, trie, tmpFile.name, n),shell=True))

         for k, vs in re.iteritems():
            results.setdefault(pep_header, []).append(k)
            results.setdefault("trie", []).append(method)
            for i,v in enumerate(vs):
               if i > 0:
                  results.setdefault("distance_{i}".format(i=i),[]).append(float(v))
               else:
                  results.setdefault("distance", []).append(float(v))
      except:
         logging.warning("Could not make distance calculation for trie {}".format(trie))

      os.remove(tmpFile.name)
      return pandas.DataFrame.from_dict(results)

    def parse_external_result(self, result):

        """

        :rtype : DataFrame
        """
        parsedResult = {}

        for line in result.strip().split('\n'):
            splitted = line.strip().split(" ")[-1].split(";")
            distanceValues = []
            peptide = splitted[0].split(":")[0]

            for s in splitted[:-1]:
                distanceValues.append(float(s.split(",")[-1])/float(len(peptide)))

            parsedResult[peptide] = distanceValues
        return parsedResult
