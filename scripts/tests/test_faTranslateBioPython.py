import unittest
from Bio import SeqIO
from faTranslateBioPython import pad_seq, translate_records
try:
    from StringIO import StringIO # Python 2
except ImportError:
    from io import StringIO # Python 3

# Unit testing class
class TestFunctions(unittest.TestCase):
    """Unit test functions"""

    def test_pad_seq(self):
        """Test pad string function"""
        pad_seq_test_cases = [
        ('ATG', 'ATG'),
        ('ATGA', 'ATGANN'),
        ('ATGAA', 'ATGAAN')]
        for in_test, out_test in pad_seq_test_cases:
            self.assertEqual(pad_seq(in_test), out_test)

    def test_translate_records(self):
        """Test translate records"""
        # Basic test cases
        self.assertEqual(str(next(translate_records(StringIO('>1\nAATGGGC'), ['M'], 11)).seq), 'MG')
        self.assertEqual(str(next(translate_records(StringIO('>1\nAATGGGC'), ['M', 'V', 'L'], 11)).seq), 'MG')

        # More complex test case
        test_fasta = """
        >mph(D)_1_AB048591
        CTCCTGTAACCAAGCCAATTGCTACATGCGCTCTTCATCAACAATCTCCTGAAACTCCTC
        GTCCCAGACATCGCAGATTCGCTTGAGAACAGCAATGTGGAAGATGTAAGTCTTTGAAAT
        CAGCGCATCCCGGGCCGCGCAGACATGTTGGCGGATCTCCTCCCAGAGTGGCTATTGGTG
        GCGTGTCTCAAGGGTTGGCTGAGGGTGTGGCTAGATTTGTTGATCGNCTCAAACGCGGTG
        CTCCTTCTGTATTCGATGCTCGTCGCCTCGGCGGCGGGTGTGTTCTTGCTTGGCTCGTAA
        ACCCCCGCCTCGGTCGTGCCGCCGGCTTCGGCAGATTTCGGCGAGAGCTTTTCTACCTTC
        CTCGCAGCCCTTCATTCGGCAACAGCCTGTGCAACCATTGACAACTGCATGTTGAGACGT
        CCAGAAACGGTAGCCACACAACATTCGAGCAATCTTCATCATCGCGTGCGTCAAGAGCTG
        CTGAATCGCTACACGCTACTGCTGCGAAGAGCACGTCGGCCCGGGATCGATTCCGGCTGC
        GGTGACTCCTTGGTCTCGCCCCGTGGCGCCTCCAATGCTTCCGAAGTGGCTATTGGCATC
        ACGGATCGTCGCTGCGGAAGCTGCGCATAGAAACAGTCGTCGTCTCACGCCGAAGACATC
        GCGATCGTCGCCAGCCCGGTGAGCTGTGCCTCAAACGCGAAAGTGCCCCTTCTCGTATCC
        GATGTCGCGCCTCGGCGGCGGGTGTGTTCTTGCTTGGCTCGTAAACCCCCGCTGGGTCGT
        GCGCAGGACTCGGAGGTCTTCGCGGAGAGCTCCCTCGCCTAATCTGGTCGGGGTTGATAA
        """.replace(' ', '')
        out1 = """
        >mph(D)_1_AB048591
        MMKIARMLCGYRFWTSQHAVVNGCTGCCRMKGCEEGRKALAEICRSRRHDRGGGLRAKQEHTRRRGDEHRIQKEHRV
        """.replace(' ', '')
        out2 = """
        >mph(D)_1_AB048591
        VVPPASADFGESFSTFLAALHSATACATIDNCMLRRPETVATQHSSNLHHRVRQELLNRYTLLLRRARRPGIDSGCGDSLVSPRGASNASEVAIGITDRRCGSCA
        """.replace(' ', '')
        # Note, this tests it behaves as it is programmed to do, but it's not actually clear what the correct protein is for this case!
        # The top blastx hit starts: RIPGRADMLAD
        self.assertEqual(str(next(translate_records(StringIO(test_fasta), ['M'], 11)).seq), str(SeqIO.read(StringIO(out1), 'fasta').seq))
        self.assertEqual(str(next(translate_records(StringIO(test_fasta), ['M', 'V', 'L'], 11)).seq), str(SeqIO.read(StringIO(out2), 'fasta').seq))
        # TODO: expand test cases

def main():
    unittest.main()

if __name__ == '__main__':
    main()
