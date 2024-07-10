import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile
import shutil

from src.pipeline import (
    find_fastq,
    # get_sample_ids,  # Commented out as it doesn't exist
    create_directories,
    align_files,
    process_sample,
    process_dir,
    generate_analysis
)


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

    # def test_get_sample_ids(self):
    #     fastqs = list(self.test_dir.glob("*.fastq.gz"))
    #     sample_ids = get_sample_ids(fastqs)
    #     self.assertEqual(sample_ids, ['sample1', 'sample2'])

    def test_create_directories(self):
        sample_ids = ['sample1', 'sample2']
        create_directories(self.test_dir, sample_ids)
        for sample in sample_ids:
            for sub_dir in ["", "BAM", "VCF", "QC", "Metrics"]:
                self.assertTrue((self.test_dir / "BAM" / sample / sub_dir).exists())

    @patch('src.pipeline.bwa_align.run_bwa')
    def test_align_files(self, mock_run_bwa):
        samples = ['sample1', 'sample2']
        fastqs = list(self.test_dir.glob("*.fastq.gz"))
        align_files(self.test_dir, samples, fastqs)
        self.assertEqual(mock_run_bwa.call_count, 2)

    @patch('src.pipeline.dup_indels.remove_duplicates')
    @patch('src.pipeline.utils.move_bam')
    @patch('src.pipeline.recalibration.base_recal1')
    @patch('src.pipeline.recalibration.recalibrate')
    def test_process_sample(self, mock_recalibrate, mock_base_recal1, mock_move_bam, mock_remove_duplicates):
        mock_remove_duplicates.return_value = 'success'
        mock_base_recal1.return_value = 'success'
        mock_recalibrate.return_value = 'success'

        results = process_sample(self.test_dir, 'sample1', 'test_panel', False)

        self.assertEqual(results['dedup'], 'success')
        self.assertEqual(results['recalibration1'], 'success')
        self.assertEqual(results['recalibrate'], 'success')
        self.assertEqual(mock_move_bam.call_count, 2)

    @patch('src.pipeline.find_fastq')
    @patch('src.pipeline.create_directories')
    @patch('src.pipeline.align_files')
    @patch('src.pipeline.process_sample')
    def test_process_dir(self, mock_process_sample, mock_align_files, mock_create_directories, mock_find_fastq):
        mock_find_fastq.return_value = ['sample1_R1.fastq.gz', 'sample1_R2.fastq.gz']
        mock_process_sample.return_value = {'dedup': 'success', 'recalibration1': 'success', 'recalibrate': 'success'}

        results = process_dir(self.test_dir, ['sample1'], 'test_panel', False)

        self.assertEqual(len(results), 1)
        self.assertIn('sample1', results)
        self.assertEqual(results['sample1'],
                         {'dedup': 'success', 'recalibration1': 'success', 'recalibrate': 'success'})

    @patch('src.pipeline.process_dir')
    def test_generate_analysis(self, mock_process_dir):
        mock_process_dir.return_value = {
            'sample1': {'dedup': 'success', 'recalibration1': 'success', 'recalibrate': 'success'}}

        results = generate_analysis('config.yaml', self.test_dir, ['sample1'], 'test_panel', False)

        self.assertEqual(results,
                         {'sample1': {'dedup': 'success', 'recalibration1': 'success', 'recalibrate': 'success'}})


if __name__ == '__main__':
    unittest.main()