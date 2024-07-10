"""
Unit tests for the pipeline module
"""

import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
from pipeline import find_fastq, get_sample_ids, create_directories, align_files, process_sample


class TestPipeline(unittest.TestCase):

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

    @patch('pipeline.log_to_api')
    def test_create_directories(self, mock_log):
        sample_ids = ['sample1', 'sample2']
        create_directories(self.test_dir, sample_ids)
        for sample in sample_ids:
            for sub_dir in ["", "BAM", "VCF", "QC", "Metrics"]:
                self.assertTrue((self.test_dir / 'BAM' / sample / sub_dir).exists())
        self.assertEqual(mock_log.call_count, 10)  # 5 directories * 2 samples

    @patch('pipeline.bwa_align.run_bwa')
    def test_align_files(self, mock_run_bwa):
        samples = ['sample1', 'sample2']
        fastqs = list(self.test_dir.glob('*.fastq.gz'))
        align_files(self.test_dir, samples, fastqs)
        self.assertEqual(mock_run_bwa.call_count, 2)

    @patch('pipeline.dup_indels.remove_duplicates')
    @patch('pipeline.utils.move_bam')
    @patch('pipeline.recalibration.base_recal1')
    @patch('pipeline.recalibration.recalibrate')
    def test_process_sample(self, mock_recalibrate, mock_base_recal1, mock_move_bam, mock_remove_duplicates):
        mock_remove_duplicates.return_value = 'success'
        mock_base_recal1.return_value = 'success'
        mock_recalibrate.return_value = 'success'

        results = process_sample(self.test_dir, 'sample1', 'test_panel', False)

        self.assertEqual(results['dedup'], 'success')
        self.assertEqual(results['recalibration1'], 'success')
        self.assertEqual(results['recalibrate'], 'success')
        self.assertEqual(mock_move_bam.call_count, 2)


if __name__ == '__main__':
    unittest.main()