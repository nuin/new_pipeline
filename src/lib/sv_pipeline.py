#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path
import logging
import concurrent.futures
from typing import List
import click
import shutil

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class SVDetectionPipeline:
    def __init__(self, bam_file: Path, reference: Path, output_dir: Path, threads: int = 4, force: bool = False):
        self.bam_file = str(bam_file)
        self.reference = str(reference)
        self.output_dir = str(output_dir)
        self.threads = threads
        self.force = force
        os.makedirs(self.output_dir, exist_ok=True)

    def log_vcf_content(self, vcf_file: str):
        cmd = f"bcftools view -H {vcf_file} | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        count = result.stdout.strip()
        logging.info(f"Number of variants in {vcf_file}: {count}")

    def run_gridss(self):
        logging.info("Running GRIDSS")

        self.bam_file = os.path.abspath(self.bam_file)
        self.reference = os.path.abspath(self.reference)
        self.output_dir = os.path.abspath(self.output_dir)

        gridss_out = os.path.join(self.output_dir, "gridss_output.vcf.gz")
        gridss_work_dir = os.path.join(self.output_dir, "gridss_work")
        gridss_assembly = os.path.join(self.output_dir, "gridss_assembly.bam")

        if os.path.exists(gridss_out) and not self.force:
            logging.info(f"GRIDSS output file already exists: {gridss_out}")
            self.log_vcf_content(gridss_out)
            return gridss_out

        if self.force:
            shutil.rmtree(gridss_work_dir, ignore_errors=True)

        os.makedirs(gridss_work_dir, exist_ok=True)

        cmd = [
            "gridss",
            f"OUTPUT={gridss_out}",
            f"REFERENCE_SEQUENCE={self.reference}",
            f"INPUT={self.bam_file}",
            f"ASSEMBLY={gridss_assembly}",
            f"THREADS={self.threads}",
            f"WORKING_DIR={gridss_work_dir}",
            "TMP_DIR=/tmp"
        ]

        cmd_str = " ".join(cmd)

        try:
            logging.info(f"GRIDSS command: {cmd_str}")
            result = subprocess.run(cmd_str, shell=True, check=True, capture_output=True, text=True)
            logging.info(f"GRIDSS stdout: {result.stdout}")
            logging.info(f"GRIDSS stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            logging.error(f"GRIDSS failed with exit code {e.returncode}")
            logging.error(f"GRIDSS stdout: {e.stdout}")
            logging.error(f"GRIDSS stderr: {e.stderr}")
            raise

        self.log_vcf_content(gridss_out)
        return gridss_out

    def run_smoove(self) -> str:
        logging.info("Running Smoove")
        smoove_out_dir = os.path.join(self.output_dir, "smoove_output")
        os.makedirs(smoove_out_dir, exist_ok=True)

        bam_file = os.path.abspath(self.bam_file)
        reference = os.path.abspath(self.reference)
        smoove_out_dir = os.path.abspath(smoove_out_dir)

        bam_dir = os.path.dirname(bam_file)

        output_vcf = os.path.join(smoove_out_dir,
                                  f"{os.path.splitext(os.path.basename(bam_file))[0]}-smoove.genotyped.vcf.gz")

        if os.path.exists(output_vcf) and not self.force:
            logging.info(f"Smoove output file already exists: {output_vcf}")
            self.log_vcf_content(output_vcf)
            return output_vcf

        if self.force:
            shutil.rmtree(smoove_out_dir, ignore_errors=True)
            os.makedirs(smoove_out_dir, exist_ok=True)

        cmd = [
            "sudo docker run --rm",
            f"-v {bam_dir}:/data",
            f"-v {os.path.dirname(reference)}:/ref",
            f"-v {smoove_out_dir}:/out",
            "brentp/smoove",
            "smoove call",
            f"--outdir /out",
            f"--name {os.path.splitext(os.path.basename(bam_file))[0]}",
            "--fasta /ref/" + os.path.basename(reference),
            f"-p {self.threads}",
            "--genotype",
            "--duphold",
            "--removepr",
            "/data/" + os.path.basename(bam_file)
        ]

        cmd_str = " ".join(cmd)
        logging.info(f"Smoove command: {cmd_str}")

        try:
            result = subprocess.run(cmd_str, shell=True, check=True, capture_output=True, text=True)
            logging.info(f"Smoove stdout: {result.stdout}")
            logging.info(f"Smoove stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Smoove failed with exit code {e.returncode}")
            logging.error(f"Smoove stdout: {e.stdout}")
            logging.error(f"Smoove stderr: {e.stderr}")
            raise

        self.log_vcf_content(output_vcf)
        return output_vcf

    def run_svaba(self) -> str:
        logging.info("Running SVABA")
        svaba_out = os.path.join(self.output_dir, "svaba_output")
        output_vcf = f"{svaba_out}.svaba.sv.vcf"

        if os.path.exists(output_vcf) and not self.force:
            logging.info(f"SVABA output file already exists: {output_vcf}")
            self.log_vcf_content(output_vcf)
            return output_vcf

        if self.force:
            for file in Path(self.output_dir).glob("svaba_output*"):
                file.unlink()

        cmd = [
            "svaba run",
            f"-t {self.bam_file}",
            f"-p {self.threads}",
            f"-a {svaba_out}",
            f"-G {self.reference}"
        ]
        subprocess.run(" ".join(cmd), shell=True, check=True)

        self.log_vcf_content(output_vcf)
        return output_vcf

    def merge_sv_calls(self, vcf_files: List[str]) -> str:
        logging.info("Merging SV calls")
        merged_vcf = os.path.join(self.output_dir, "merged_svs.vcf")

        if os.path.exists(merged_vcf) and not self.force:
            logging.info(f"Merged VCF file already exists: {merged_vcf}")
            self.log_vcf_content(merged_vcf)
            return merged_vcf

        with open(os.path.join(self.output_dir, "vcf_list.txt"), "w") as f:
            for vcf in vcf_files:
                f.write(f"{vcf}\n")

        logging.info(f"Content of VCF list file:")
        with open(os.path.join(self.output_dir, 'vcf_list.txt'), 'r') as f:
            logging.info(f.read())

        cmd = [
            "SURVIVOR merge",
            os.path.join(self.output_dir, 'vcf_list.txt'),
            "1000",  # Max distance between breakpoints
            "1",  # Minimum number of supporting callers (reduced from 2 to 1)
            "1",  # Take the type into account
            "1",  # Take the strand into account
            "0",  # Disabled
            "30",  # Minimum size of SVs to be taken into account
            merged_vcf
        ]
        subprocess.run(" ".join(cmd), shell=True, check=True)

        self.log_vcf_content(merged_vcf)
        return merged_vcf

    def run_pipeline(self) -> Path:
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
            future_gridss = executor.submit(self.run_gridss)
            future_smoove = executor.submit(self.run_smoove)
            future_svaba = executor.submit(self.run_svaba)

            gridss_out = future_gridss.result()
            smoove_out = future_smoove.result()
            svaba_out = future_svaba.result()

        merged_vcf = self.merge_sv_calls([gridss_out, smoove_out, svaba_out])

        logging.info(f"Pipeline completed. Final merged VCF: {merged_vcf}")
        return merged_vcf


@click.command()
@click.option('-b', '--bam', required=True, type=click.Path(exists=True), help='Input BAM file')
@click.option('-r', '--reference', required=True, type=click.Path(exists=True), help='Reference genome FASTA file')
@click.option('-o', '--output', required=True, type=click.Path(), help='Output directory')
@click.option('-t', '--threads', default=4, help='Number of threads to use')
@click.option('-f', '--force', is_flag=True, help='Force overwrite of existing results')
def main(bam: str, reference: str, output: str, threads: int, force: bool):
    """Run the Structural Variant Detection Pipeline."""
    pipeline = SVDetectionPipeline(Path(bam), Path(reference), Path(output), threads, force)
    final_vcf = pipeline.run_pipeline()

    click.echo(f"SV detection pipeline completed. Final VCF: {final_vcf}")


if __name__ == "__main__":
    main()