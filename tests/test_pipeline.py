import sys
from pathlib import Path
import unittest
from pathlib import Path
# Add the project root directory to Python's module search path
project_root = Path(__file__).parents[1]
sys.path.insert(0, str(project_root))

from src.pipeline import find_fastq, get_sample_ids, create_directories, align_files, process_sample

class TestMain(unittest.TestCase):

    def setUp(self):
        self.test_dir = Path('/tmp/test_pipeline')
        self.test_dir.mkdir(exist_ok=True)
        (self.test_dir / 'sample1_R1.fastq.gz').touch()
        (self.test_dir / 'sample1_R2.fastq.gz').touch()
        (self.test_dir / 'sample2_R1.fastq.gz').touch()
        (self.test_dir / 'sample2_R2.fastq.gz').touch()

    def tearDown(self):
        for file in self.test_dir.glob('*'):
            file.unlink()
        self.test_dir.rmdir()

    def test_find_fastq(self):
        fastqs = find_fastq(self.test_dir)
        self.assertEqual(len(fastqs), 4)
        self.assertTrue(all(fastq.suffix == '.gz' for fastq in fastqs))

    def test_get_sample_ids(self):
        fastqs = list(self.test_dir.glob('*.fastq.gz'))
        sample_ids = get_sample_ids(fastqs)
        self.assertEqual(sample_ids, ['sample1', 'sample2'])

    # Add more tests here...


if __name__ == '__main__':
    unittest.main()