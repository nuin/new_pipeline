import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile
import shutil

from src.pipeline import find_fastq


class TestPipeline(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for testing
        self.test_dir = Path(tempfile.mkdtemp())

        # Create some mock FASTQ files
        (self.test_dir / "sample1_R1.fastq.gz").touch()
        (self.test_dir / "sample1_R2.fastq.gz").touch()
        (self.test_dir / "sample2_R1.fastq.gz").touch()
        (self.test_dir / "sample2_R2.fastq.gz").touch()

    def tearDown(self):
        # Remove the temporary directory after the test
        shutil.rmtree(self.test_dir)

    def test_find_fastq(self):
        fastqs = find_fastq(self.test_dir)
        self.assertEqual(len(fastqs), 4)
        self.assertTrue(all(fastq.suffix == '.gz' for fastq in fastqs))


if __name__ == '__main__':
    unittest.main()