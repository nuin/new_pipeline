import subprocess
from pathlib import Path
from typing import Dict, Optional


def get_vep_version(vep: str) -> str:
    cmd = [vep, "--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    for line in result.stdout.split("\n"):
        if "ensembl-vep" in line:
            return line.strip()
    return "Unknown"


def vep_annotate(
        sample_id: str,
        datadir: Path,
        vep: str,
        reference: Path,
        db: Dict,
        transcript_list: Optional[str] = None,
        max_retries: int = 3,
) -> str:
    vcf_dir = datadir / "BAM" / sample_id / "VCF"
    input_vcf = vcf_dir / f"{sample_id}_merged.vcf"
    output_vcf = vcf_dir / f"{sample_id}_merged.vep.vcf"
    vep_dir = Path("/apps/data/src/bin/vep")

    if not input_vcf.exists():
        return "error"

    if output_vcf.exists():
        return "exists"

    vep_input = vep_dir / input_vcf.name
    vep_output = vep_dir / output_vcf.name
    subprocess.run(["sudo", "cp", str(input_vcf), str(vep_input)], check=True)
    subprocess.run(["sudo", "chown", "systemd-coredump:input", str(vep_input)], check=True)
    subprocess.run(["sudo", "chmod", "644", str(vep_input)], check=True)

    vep_cmd = f"""
    sudo docker run --rm \
    -v {vep_dir}:/opt/vep/.vep:Z \
    ensemblorg/ensembl-vep \
    vep \
    -i /opt/vep/.vep/{input_vcf.name} \
    -o /opt/vep/.vep/{output_vcf.name} \
    --offline --cache --dir_cache /opt/vep/.vep \
    --assembly GRCh37 --species homo_sapiens \
    --fasta /opt/vep/.vep/homo_sapiens/112_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
    --format vcf --vcf --force_overwrite --hgvs \
    --dir /opt/vep/.vep --config /opt/vep/.vep/vep_config.ini
    """

    if transcript_list:
        vep_cmd = vep_cmd.replace(
            "vep \\",
            f"vep --transcript_filter file=/opt/vep/.vep/{transcript_list} \\",
        )

    for attempt in range(max_retries):
        process = subprocess.run(vep_cmd, shell=True, capture_output=True, text=True)

        if process.returncode == 0:
            break
        elif attempt == max_retries - 1:
            return "error"

    if vep_output.exists():
        subprocess.run(["sudo", "chown", f"{Path.owner()}:{Path.group()}", str(vep_output)], check=True)
        subprocess.run(["sudo", "cp", str(vep_output), str(output_vcf)], check=True)
        subprocess.run(["sudo", "rm", str(vep_input)], check=True)
        subprocess.run(["sudo", "rm", str(vep_output)], check=True)
        return "success"
    else:
        return "error"