import os
import subprocess
from pathlib import Path
import logging
import concurrent.futures

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class SVDetectionPipeline:
    def __init__(self, bam_file, reference, output_dir, threads=4):
        self.bam_file = Path(bam_file)
        self.reference = Path(reference)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run_gridss(self):
        logging.info("Running GRIDSS")
        gridss_out = self.output_dir / "gridss_output.vcf.gz"
        cmd = [
            "gridss",
            f"OUTPUT={gridss_out}",
            f"INPUT={self.bam_file}",
            f"THREADS={self.threads}",
            f"WORKING_DIR={self.output_dir / 'gridss_work'}",
        ]

        try:
            logging.info(f"GRIDSS command: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logging.info(f"GRIDSS stdout: {result.stdout}")
            logging.info(f"GRIDSS stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            logging.error(f"GRIDSS failed with exit code {e.returncode}")
            logging.error(f"GRIDSS stdout: {e.stdout}")
            logging.error(f"GRIDSS stderr: {e.stderr}")
            raise
        return gridss_out

    def run_smoove(self):
        logging.info("Running Smoove")
        smoove_out_dir = self.output_dir / "smoove_output"
        smoove_out_dir.mkdir(exist_ok=True)
        cmd = [
            "smoove call",
            f"--outdir {smoove_out_dir}",
            f"--name {self.bam_file.stem}",
            f"--fasta {self.reference}",
            f"-p {self.threads}",
            self.bam_file
        ]
        subprocess.run(" ".join(cmd), shell=True, check=True)
        return smoove_out_dir / f"{self.bam_file.stem}-smoove.genotyped.vcf.gz"

    def run_svaba(self):
        logging.info("Running SVABA")
        svaba_out = self.output_dir / "svaba_output"
        cmd = [
            "svaba run",
            f"-t {self.bam_file}",
            f"-p {self.threads}",
            f"-a {svaba_out}",
            f"-G {self.reference}"
        ]
        subprocess.run(" ".join(cmd), shell=True, check=True)
        return svaba_out.with_suffix(".svaba.sv.vcf")

    def run_cnvnator(self):
        logging.info("Running CNVnator")
        cnvnator_out = self.output_dir / "cnvnator_output"
        bin_size = 100  # Adjust as needed
        cmd = [
            f"cnvnator -root {cnvnator_out}.root -tree {self.bam_file} -chrom $(seq 1 22) X Y",
            f"cnvnator -root {cnvnator_out}.root -his {bin_size} -d {self.reference.parent}",
            f"cnvnator -root {cnvnator_out}.root -stat {bin_size}",
            f"cnvnator -root {cnvnator_out}.root -partition {bin_size}",
            f"cnvnator -root {cnvnator_out}.root -call {bin_size} > {cnvnator_out}.txt"
        ]
        for c in cmd:
            subprocess.run(c, shell=True, check=True)
        return cnvnator_out.with_suffix(".txt")

    def merge_sv_calls(self, vcf_files):
        logging.info("Merging SV calls")
        merged_vcf = self.output_dir / "merged_svs.vcf"
        with open(self.output_dir / "vcf_list.txt", "w") as f:
            for vcf in vcf_files:
                f.write(f"{vcf}\n")
        cmd = [
            "SURVIVOR merge",
            f"{self.output_dir / 'vcf_list.txt'}",
            "1000",  # Max distance between breakpoints
            "2",  # Minimum number of supporting callers
            "1",  # Take the type into account
            "1",  # Take the strand into account
            "0",  # Disabled
            "30",  # Minimum size of SVs to be taken into account
            merged_vcf
        ]
        subprocess.run(" ".join(cmd), shell=True, check=True)
        return merged_vcf

    def annotate_svs(self, merged_vcf):
        logging.info("Annotating SVs")
        annotated_vcf = self.output_dir / "annotated_svs.vcf"
        cmd = [
            "AnnotSV",
            "-SVinputFile", merged_vcf,
            "-outputFile", annotated_vcf,
            "-genomeBuild", "GRCh38"  # Adjust as needed
        ]
        subprocess.run(" ".join(cmd), shell=True, check=True)
        return annotated_vcf

    def run_pipeline(self):
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            future_gridss = executor.submit(self.run_gridss)
            future_smoove = executor.submit(self.run_smoove)
            future_svaba = executor.submit(self.run_svaba)
            future_cnvnator = executor.submit(self.run_cnvnator)

            gridss_out = future_gridss.result()
            smoove_out = future_smoove.result()
            svaba_out = future_svaba.result()
            cnvnator_out = future_cnvnator.result()

        merged_vcf = self.merge_sv_calls([gridss_out, smoove_out, svaba_out])
        annotated_vcf = self.annotate_svs(merged_vcf)

        logging.info(f"Pipeline completed. Final annotated VCF: {annotated_vcf}")
        return annotated_vcf


def main():
    bam_file = "/path/to/your/sample.bam"
    reference = "/path/to/reference/genome.fa"
    output_dir = "/path/to/output/directory"

    pipeline = SVDetectionPipeline(bam_file, reference, output_dir)
    final_vcf = pipeline.run_pipeline()

    logging.info(f"SV detection pipeline completed. Final VCF: {final_vcf}")


if __name__ == "__main__":
    main()