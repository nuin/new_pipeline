"""
.. module:: snpEff_ann
    :platform: any
    :synopsis: This module calls snpEff to annotate variants from the merged VCF file, including HGVS notation
.. moduleauthor:: Paulo Nuin, June 2016, Updated September 2024

"""

import subprocess
from pathlib import Path
from typing import Dict, Union, List
import logging
import shlex
from rich.console import Console
from rich.progress import Progress
from rich.panel import Panel
from rich.logging import RichHandler
import os
import gzip

# Set up logging
logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)]
)
logger = logging.getLogger("snpEff_ann")

console = Console()

# Configuration
SNPEFF_CONFIG = {
    "database": "GRCh37.75",  
    "transcript_list": "/apps/data/src/bundle/transcripts_only.txt",
    "annotation_params": [
        "-noLog",       
        "-quiet",       
        "-canon",
    ],
    "fallback_data_dir": "/apps/data/src/bin/snpEff/data"  # Fallback database directory
}

def get_snpeff_info(snpEff: str, config: Dict) -> Dict[str, str]:
    """Get snpEff version and database directory information."""
    version_cmd = shlex.split(snpEff) + ["-version"]
    version_result = subprocess.run(version_cmd, capture_output=True, text=True)
    
    config_cmd = shlex.split(snpEff) + ["config", "-v"]
    config_result = subprocess.run(config_cmd, capture_output=True, text=True)
    
    info = {
        "version": version_result.stdout.strip() if version_result.returncode == 0 else "Unknown",
        "database_dir": ""
    }
    
    for line in config_result.stdout.split('\n'):
        if line.startswith("data_dir"):
            info["database_dir"] = line.split('=')[1].strip()
            break
    
    if not info["database_dir"]:
        logger.warning("Database directory not found in snpEff config. Using fallback directory.")
        info["database_dir"] = config["fallback_data_dir"]
    
    return info

def check_database_files(database_dir: Path, database: str) -> Dict[str, bool]:
    """Check if the necessary database files exist."""
    db_path = database_dir / database
    return {
        "snpEffectPredictor.bin": (db_path / "snpEffectPredictor.bin").exists(),
        "sequences.fa": (db_path / "sequences.fa").exists(),
        "genes.gtf": (db_path / "genes.gtf").exists(),
    }

def download_database(snpEff: str, database: str) -> bool:
    """Download the specified database."""
    cmd = shlex.split(snpEff) + ["download", "-v", database]
    logger.info(f"Attempting to download database: {database}")
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Database download output:\n{result.stdout}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to download database {database}: {e.stderr}")
        return False

def read_transcript_list(file_path: str) -> List[str]:
    """Read the transcript list file and return a list of transcript IDs."""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def run_snpeff(input_vcf: Path, output_vcf: Path, snpEff: str, config: Dict, database_dir: Path, transcripts: List[str]) -> subprocess.CompletedProcess:
    """Run snpEff with the specified parameters."""
    cmd = shlex.split(snpEff) + [
        "ann",
        *config["annotation_params"],
        f"-dataDir", str(database_dir),
    ]

    if transcripts:
        transcript_arg = ','.join(transcripts)
        cmd.append(f"-onlyTr")
        cmd.append(transcript_arg)
    
    cmd.extend([
        config["database"],
        str(input_vcf)
    ])
    
    logger.info(f"Running snpEff command: {' '.join(cmd)}")
    
    with open(output_vcf, 'w') as out_file:
        return subprocess.run(cmd, stdout=out_file, stderr=subprocess.PIPE, text=True, check=True)

def post_process_vcf(input_vcf: Path, output_vcf: Path):
    """Post-process the VCF file to extract HGVS notations if needed."""
    with gzip.open(input_vcf, 'rt') if input_vcf.suffix == '.gz' else open(input_vcf, 'r') as infile, \
         gzip.open(output_vcf, 'wt') if output_vcf.suffix == '.gz' else open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                fields = line.strip().split('\t')
                info = dict(item.split('=') for item in fields[7].split(';') if '=' in item)
                
                # Extract HGVS notations
                hgvs_c = info.get('HGVS.c', '')
                hgvs_p = info.get('HGVS.p', '')
                
                # Add HGVS notations to the INFO field if not already present
                if 'HGVS.c' not in fields[7]:
                    fields[7] += f";HGVS.c={hgvs_c}"
                if 'HGVS.p' not in fields[7]:
                    fields[7] += f";HGVS.p={hgvs_p}"
                
                outfile.write('\t'.join(fields) + '\n')

def annotate_merged(sample_id: str, datadir: Union[str, Path], snpEff: str, config: Dict = SNPEFF_CONFIG) -> str:
    """
    Function that calls snpEff to annotate the merged VCF file.

    :param sample_id: ID of the patient/sample being analysed
    :param datadir: Location of the data files
    :param snpEff: snpEff command string from .env file
    :param config: Configuration dictionary for snpEff
    :return: Status string for compatibility with existing pipeline

    :rtype: str
    """
    datadir = Path(datadir)
    vcf_dir = datadir / "BAM" / sample_id / "VCF"
    input_vcf = vcf_dir / f"{sample_id}_merged.vcf"
    output_vcf = vcf_dir / f"{sample_id}_merged.ann.vcf"
    final_vcf = vcf_dir / f"{sample_id}_merged.ann.hgvs.vcf"
    
    if final_vcf.exists():
        logger.info(f"Annotated VCF file with HGVS already exists for {sample_id}")
        return "exists"
    
    if not input_vcf.exists():
        logger.error(f"Input VCF file not found for {sample_id}")
        return "error"
    
    try:
        # Get snpEff information
        snpeff_info = get_snpeff_info(snpEff, config)
        logger.info(f"snpEff version: {snpeff_info['version']}")
        logger.info(f"snpEff database directory: {snpeff_info['database_dir']}")
        
        database_dir = Path(snpeff_info['database_dir'])
        
        # Check database files
        db_files = check_database_files(database_dir, config["database"])
        logger.info(f"Database files status: {db_files}")
        
        if not all(db_files.values()):
            logger.warning(f"Some database files are missing for {config['database']}. Attempting to download...")
            if not download_database(snpEff, config["database"]):
                logger.error("Failed to download the database. Please download it manually or check your internet connection.")
                return "error"
            db_files = check_database_files(database_dir, config["database"])
            if not all(db_files.values()):
                logger.error("Database files are still missing after download attempt. Please check your snpEff installation.")
                return "error"
        
        # Read transcript list
        transcripts = read_transcript_list(config["transcript_list"])
        logger.info(f"Number of transcripts in list: {len(transcripts)}")
        
        if not transcripts:
            logger.warning("No transcripts found in the transcript list. Will use all available transcripts.")
        
        with Progress() as progress:
            task = progress.add_task("[cyan]Running snpEff annotation...", total=100)
            
            # Run snpEff
            result = run_snpeff(input_vcf, output_vcf, snpEff, config, database_dir, transcripts)
            progress.update(task, advance=75)
            
            # Post-process VCF to ensure HGVS notations are present
            post_process_vcf(output_vcf, final_vcf)
            progress.update(task, advance=25)
        
        # Log snpEff output for debugging
        logger.debug(f"snpEff stderr output:\n{result.stderr}")
        
        if final_vcf.exists() and final_vcf.stat().st_size > 0:
            console.print(Panel(f"[bold green]Variant annotation with HGVS completed successfully for {sample_id}[/bold green]"))
            return "success"
        else:
            logger.error(f"snpEff ran without errors, but output file is missing or empty for {sample_id}")
            return "error"
        
    except subprocess.CalledProcessError as e:
        logger.error(f"snpEff annotation failed for {sample_id}: {e}")
        logger.error(f"snpEff stderr output:\n{e.stderr}")
        return "error"
    
    except Exception as e:
        logger.exception(f"Unexpected error during snpEff annotation for {sample_id}")
        return "error"

if __name__ == "__main__":
    # Example usage and testing
    from dotenv import load_dotenv
    
    load_dotenv()
    
    sample_id = "SAMPLE001"
    datadir = Path("/path/to/data/directory")
    snpEff_cmd = os.getenv('SNPEFF', 'java -Xmx8g -jar /apps/data/src/bin/snpEff/snpEff.jar')

    result = annotate_merged(sample_id, datadir, snpEff_cmd)
    print(f"Annotation result: {result}")