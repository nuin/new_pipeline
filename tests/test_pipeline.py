import unittest
from pathlib import Path
import tempfile
import shutil
import gzip
import random
from unittest import mock

from src.pipeline import find_fastq, create_directories, align_files, process_dir, generate_analysis


def generate_fastq_entry(read_id, length=100):
    bases = 'ACGT'
    sequence = ''.join(random.choice(bases) for _ in range(length))
    quality = ''.join(chr(random.randint(33, 74)) for _ in range(length))  # Phred+33 encoding
    return f"@{read_id}\n{sequence}\n+\n{quality}\n"


class TestPipeline(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for testing
        self.test_dir = Path(tempfile.mkdtemp())

        # Create some mock FASTQ files
        self.create_test_fastq(self.test_dir / "sample1_R1.fastq.gz")
        self.create_test_fastq(self.test_dir / "sample1_R2.fastq.gz")
        self.create_test_fastq(self.test_dir / "sample2_R1.fastq.gz")
        self.create_test_fastq(self.test_dir / "sample2_R2.fastq.gz")

        # Define panel samples and panel for testing
        self.panel_samples = ["sample1", "sample2"]
        self.panel = "test_panel"

    def tearDown(self):
        # Remove the temporary directory after the test
        shutil.rmtree(self.test_dir)

    def create_test_fastq(self, filename, num_reads=1000, read_length=100):
        with gzip.open(filename, 'wt') as f:
            for i in range(num_reads):
                f.write(generate_fastq_entry(f"READ_{i + 1}", read_length))

    def test_find_fastq_basic(self):
        fastqs = find_fastq(self.test_dir, self.panel_samples, self.panel)
        self.assertEqual(len(fastqs), 4)
        self.assertTrue(all(fastq.suffix == '.gz' for fastq in fastqs))
        self.assertTrue(all(any(sample in fastq.name for sample in self.panel_samples) for fastq in fastqs))

    def test_find_fastq_with_subset_of_samples(self):
        subset_samples = ["sample1"]
        fastqs = find_fastq(self.test_dir, subset_samples, self.panel)
        self.assertEqual(len(fastqs), 2)  # Assuming it only returns files for sample1
        self.assertTrue(all('sample1' in fastq.name for fastq in fastqs))

    def test_find_fastq_with_nonexistent_sample(self):
        nonexistent_samples = ["sample3"]
        fastqs = find_fastq(self.test_dir, nonexistent_samples, self.panel)
        self.assertEqual(len(fastqs), 0)  # Assuming it returns an empty list for non-existent samples

    def test_find_fastq_empty_directory(self):
        empty_dir = self.test_dir / "empty"
        empty_dir.mkdir()
        fastqs = find_fastq(empty_dir, self.panel_samples, self.panel)
        self.assertEqual(len(fastqs), 0)

    def test_find_fastq_non_fastq_files(self):
        (self.test_dir / "not_a_fastq.txt").touch()
        fastqs = find_fastq(self.test_dir, self.panel_samples, self.panel)
        self.assertEqual(len(fastqs), 4)  # Should still only find the 4 FASTQ files

    def test_create_directories(self):
        sample_ids = ['sample1', 'sample2']
        create_directories(self.test_dir, sample_ids, self.panel)
        for sample in sample_ids:
            for sub_dir in ["", "BAM", "VCF", "QC", "Metrics"]:
                self.assertTrue((self.test_dir / "BAM" / sample / sub_dir).exists())

    def test_create_directories_existing_dirs(self):
        sample_ids = ['sample1']
        (self.test_dir / "BAM" / "sample1").mkdir(parents=True)
        create_directories(self.test_dir, sample_ids, self.panel)
        self.assertTrue((self.test_dir / "BAM" / "sample1" / "VCF").exists())

    @unittest.mock.patch('src.pipeline.run_bwa')
    def test_align_files(self, mock_run_bwa):
        samples = ['sample1', 'sample2']
        fastqs = list(self.test_dir.glob("*.fastq.gz"))
        align_files(self.test_dir, samples, fastqs)
        self.assertEqual(mock_run_bwa.call_count, 2)

    @unittest.mock.patch('src.pipeline.find_fastq')
    @unittest.mock.patch('src.pipeline.create_directories')
    @unittest.mock.patch('src.pipeline.align_files')
    @unittest.mock.patch('src.pipeline.process_sample')
    def test_process_dir(self, mock_process_sample, mock_align_files, mock_create_directories, mock_find_fastq):
        mock_find_fastq.return_value = ['sample1_R1.fastq.gz', 'sample1_R2.fastq.gz']
        mock_process_sample.return_value = {'dedup': 'success', 'recalibration1': 'success', 'recalibrate': 'success'}

        results = process_dir(self.test_dir, ['sample1'], self.panel, False)

        self.assertEqual(len(results), 1)
        self.assertIn('sample1', results)
        self.assertEqual(results['sample1'],
                         {'dedup': 'success', 'recalibration1': 'success', 'recalibrate': 'success'})

    @unittest.mock.patch('src.pipeline.process_dir')
    def test_generate_analysis(self, mock_process_dir):
        mock_process_dir.return_value = {
            'sample1': {'dedup': 'success', 'recalibration1': 'success', 'recalibrate': 'success'}}

        results = generate_analysis('config.yaml', self.test_dir, ['sample1'], self.panel, False)

        self.assertEqual(results,
                         {'sample1': {'dedup': 'success', 'recalibration1': 'success', 'recalibrate': 'success'}})

    def test_generate_analysis_invalid_config(self):
        with self.assertRaises(FileNotFoundError):
            generate_analysis('nonexistent_config.yaml', self.test_dir, ['sample1'], self.panel, False)

    def test_generate_analysis_no_samples(self):
        results = generate_analysis('config.yaml', self.test_dir, [], self.panel, False)
        self.assertEqual(results, {})  # Assuming it returns an empty dict when no samples are provided


if __name__ == '__main__':
    unittest.main()