import unittest
import tempfile
import pandas as pd
import os
import sys
from unittest.mock import patch, MagicMock

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.utils.blast_utils import run_blast, stitch_blast
from src.utils.seq_utils import write_temp_fasta, check_input
from src.database.blastdb import blast_db_exists, db_list

class TestBlastPipelineUnit(unittest.TestCase):
# In TestBlastPipelineUnit.setUpClass, fix the temporary file creation:
    @classmethod
    def setUpClass(cls):
        # Create test data
        cls.test_header = "test_sequence"
        cls.test_sequence = "ATGCATGCATGCATGC"
        cls.test_fasta = f">{cls.test_header}\n{cls.test_sequence}\n"
        
        # Use write_temp_fasta instead of manual file creation
        cls.tmp_query_fasta = write_temp_fasta(cls.test_header, cls.test_sequence)
        with open(cls.tmp_query_fasta.name, "w") as f:
            f.write(cls.test_fasta)
            
        # Mock database paths
        cls.db_list = {
            "ship": {
                "nucl": "test_db/ships.fa",
                "prot": "test_db/ships.faa"
            }
        }

    def setUp(self):
        # Create temporary output file for each test
        self.tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast", delete=False)

    def tearDown(self):
        # Clean up temporary files after each test
        if os.path.exists(self.tmp_blast.name):
            os.unlink(self.tmp_blast.name)

    @classmethod
    def tearDownClass(cls):
        # Clean up all temporary files
        if os.path.exists(cls.tmp_query_fasta.name):
            os.unlink(cls.tmp_query_fasta.name)

    def test_input_validation(self):
        # Test empty input
        input_type, header, query = check_input("", None)
        self.assertIsNone(input_type)
        self.assertIsNone(header)
        self.assertIsNone(query)

        # Test valid FASTA input
        input_type, header, query = check_input(self.test_fasta, None)
        self.assertEqual(input_type, "text")
        self.assertEqual(header, self.test_header)
        self.assertEqual(query, self.test_sequence)

    @patch('src.utils.blast_utils.subprocess.run')
    @patch('src.database.blastdb.blast_db_exists')
    def test_blast_execution(self, mock_db_exists, mock_run):
        # Mock database existence check
        mock_db_exists.return_value = True
        
        # Mock successful BLAST execution
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_run.return_value = mock_process
        
        try:
            result_df = run_blast(
                db_list=self.db_list,
                query_type="nucl",
                query_fasta=self.tmp_query_fasta,
                tmp_blast=self.tmp_blast.name,
                input_eval=0.01,
                threads=2
            )
            
            # Verify BLAST was called
            mock_run.assert_called_once()
            
        except Exception as e:
            self.fail(f"run_blast raised an exception: {str(e)}")
    def test_blast_results_processing(self):
        # Create mock BLAST results with proper column headers
        mock_results = pd.DataFrame({
            'qseqid': ['test1', 'test1'],
            'sseqid': ['hit1', 'hit2'],
            'pident': [95.0, 85.0],
            'length': [100, 90],
            'mismatch': [5, 15],
            'gapopen': [0, 2],
            'qstart': [1, 1],
            'qend': [100, 90],
            'sstart': [1, 1],
            'send': [100, 90],  # Make sure these are integers
            'evalue': [1e-50, 1e-40],
            'bitscore': [200, 180],
            'qseq': ['ATGC'*25, 'ATGC'*22],
            'sseq': ['ATGC'*25, 'ATGC'*22]
        })
        
        # Write mock results with header
        mock_results.to_csv(self.tmp_blast.name, sep='\t', index=False)
        
        # Test results processing
        processed_df = stitch_blast(self.tmp_blast.name, self.tmp_blast.name)
        
        self.assertIsInstance(processed_df, pd.DataFrame)
        self.assertTrue(len(processed_df) > 0)
        self.assertIn('pident', processed_df.columns)
        self.assertIn('evalue', processed_df.columns)

class TestBlastPipelineIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Use a real Starship sequence
        cls.test_header = "aspfum3_s01168"
        cls.test_sequence = """
        ctagaatatggggacctgtgaaatcgacgccgaaaatgacagctactgaggccataaacg
        gccgcatttccacctgtattttccctctgcgagccgccgctttagcactcaagcacgtgc
        tacacaccactattctgcacccgccaagcattgttttagtcaccgtcctgaggagacaat
        cctcagcgagaactcgctcggtttcaattcacgaccattagtgacgccgacaatggcctc
        ctcagcagccagatccagctccctcgaattcaagcagattgcccatgtgctctataaaat
        gccggccatcctcactcctttctttcctccctcgcactctcgcgccttcctcccaccttg
        acgccttccaaatcaatccgccaggcgatcgctgctcgaatcactgcagaaaacaatcag
        ctttctcttatatttggatcagcttatcttgatgtcaccgcctcctctccctcctgcaga
        tcccaatgagcatcaatactggaccgacccgatcctctgtgaagagacccgagccagact
        cgagcattttcgcagcattggatggctacccccgaacttcaaacccaagacattggaagg
        tcttgcagtggtcgaaagatactggcgaaggtacgtacgcgctttctcagattgctgcag
        ggatcgtacatccttaacgtaccgtcgacagattttgcatccattcaaaggaggactacg
        tggagtatttactcctggaggaccaaacaatttacatgaatttcttcgattggatgtaca
        ggacctcacgacagaagcttcttcagtcgtacgatgagtactggcgacggctttgccaat
        actttggtctgtttgcacggcgccgtctgaatggcaacgtctatgaacagatgcgacgag
        tacgttgtttcgttcttggtcccagcagctccaccattaacatcggcgttagttcctgaa
        taaagtattccccgcagagcgaaagatccccagacgtatgaagaagaaggacacccttga
        cgttaatgtcttctgcatcttatatcgccagcactggatatattcgagatttttccgtca
        cgggagtatgattatccaatttgccgtagtgcaactctggtctgccatcactggaacgcg
        accgggtgtactgcttcctcagaaagcttcccagactggccataagtcatccctgggcaa
        acgcaaacgcgagcagacttttcagagtgaccttcccaaacacgtttcggccgacgatct
        tccagattcagtctgctaccgcgatattgagctcttcattcttaaagaccccgatagcaa
        gcgtgatgtcctttgtgccatcattgagtttcgcaacctcaaaggccggccagaaggcgc
        tgacgggtaagcggatccatagcagatccatagacagctagttgctgacctgctgaggct
        aactttattttccctcagaactaaattcttcatgcatggtgattaccaattggcttactg
        cccaatcaaccaaattatctcactagctttccgcgatggggcttttgtcaatcccttgac
        tccagagctgatatggcggcttcgagttccaaagcgccgacaggatcttcccctccggtg
        gaaggaggagatccttgatacccctctccttcgacgtattgagcgcactccgcacggcta
        cgagctccataagtcattgccgatgacgtacgactcaagtcgacaagcgctcaaagaact
        aggtcgagatgctaagtttgaggacgatattggccattataactaccggcgctgggttgc
        caatgaagttaatcgtatgtcctcccgcgtggaagatggctctgcaaaacgaaagaacct
        gaaactgtatctaacttgaccttacagggaacttcacaagccaggagcggcagagagtgc
        taggccagtccggcgacgcggttttcgagaggcattaccagtcacagtttatcggccgcg
        atctccagcatgtcgtgcttctccgaccctctcaggaaggtctgcttcgggttgctggaa
        gcatgctcaggaagagggatccccttgccccgtcaggtcttaccgaccagcagaaacgtg
        ctatttgtcgagaccccagaatccttgagcttagacgagagcagagggagctaaaggagg
        agatgcgatcgctagctggtacagtggcgaaagctcgggagtcattccctcatctgcatg
        aaagacatgaggcagttaaaaaagagctctctcgggtaaggaagaacctgacaaaggata
        ctcgcgagacagctagaaaggaatattttcacaatgagcctgtgctcgaagttgacagac
        agatcaagcaacttcttggccagtctgatgctgaaagttgtgaagatctcgactccgacg
        atgaagactgggatattcctattccagaatacgtctttcccgagcgagcgcgactggtgg
        agaatttctatggccctgaagcggaatgttttgatgcggataagctacttgcacgacgta
        tccaggttaccaaggacatggttgcgctctcgaagctctgcgaaccaagccggcgaggca
        accaaatcaactgggatgtgagtgatggcgctaaggctgatgaaactgaacgtgaagagt
        cggaagaatcaccaccttctgaagaagaagccctggactgcccaacagatgtttgcatta
        tctgttgcggtttgtcacgtcgttcagcttctaatccgcctccacataaatttccttcca
        agcgaaaagattcactccgtcgccatctcattgattttcatctcgcccgtgcgcacgagg
        ggattagctgtacctggtcggcctgctgcaatgtccctaagttcgcaaaggctactgaat
        tcttggcgcatgccgttgatgtgcactcatatgatatcggaatcaaattaaaacatctcc
        caagcagacctctgttaacttgcagcgatacttcctctatcgacagcacagacgcttctc
        """
        
        # Create temporary query file - store the path
        cls.tmp_query_fasta = write_temp_fasta(cls.test_header, cls.test_sequence.strip())
        
        # Verify database exists
        if not blast_db_exists(db_list["ship"]["nucl"]):
            raise unittest.SkipTest("BLAST database not found - skipping integration tests")

    def setUp(self):
        self.tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast", delete=False)

    def tearDown(self):
        if os.path.exists(self.tmp_blast.name):
            os.unlink(self.tmp_blast.name)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.tmp_query_fasta):
            os.unlink(cls.tmp_query_fasta)

    def test_real_blast_search(self):
        """Test BLAST search against real database"""
        try:
            result_df = run_blast(
                db_list=db_list,
                query_type="nucl",
                query_fasta=self.tmp_query_fasta,
                tmp_blast=self.tmp_blast.name,
                input_eval=0.01,
                threads=2
            )
            
            if result_df is None:
                self.fail("BLAST search returned None")            
            # Verify results
            self.assertIsInstance(result_df, pd.DataFrame)
            self.assertTrue(len(result_df) > 0)
            
            # Check expected columns
            expected_columns = ['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore']
            for col in expected_columns:
                self.assertIn(col, result_df.columns)
            
            # Verify reasonable e-values
            self.assertTrue(all(result_df['evalue'] <= 1.0))
            
            # Check alignment quality
            self.assertTrue(any(result_df['pident'] > 80.0))
            
        except Exception as e:
            self.fail(f"Real BLAST search failed: {str(e)}")

    def test_end_to_end_pipeline(self):
        """Test entire BLAST pipeline with real data"""
        try:
            # Run BLAST
            blast_df = run_blast(
                db_list=db_list,
                query_type="nucl",
                query_fasta=self.tmp_query_fasta,
                tmp_blast=self.tmp_blast.name,
                input_eval=0.01,
                threads=2
            )
            
            # Process results
            processed_df = stitch_blast(self.tmp_blast.name, self.tmp_blast.name)
            
            # Verify pipeline output
            self.assertTrue(len(processed_df) > 0)
            self.assertIn('familyName', processed_df.columns)
            
            # Check if we found expected Starship families
            self.assertTrue(any(processed_df['familyName'].notna()))
            
        except Exception as e:
            self.fail(f"End-to-end pipeline test failed: {str(e)}")

if __name__ == '__main__':
    unittest.main() 