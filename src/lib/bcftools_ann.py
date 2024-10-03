import subprocess
from pathlib import Path
from typing import Dict, Optional

def get_bcftools_version(bcftools: str) -> str:
    result = subprocess.run([bcftools, "--version"], capture_output=True, text=True)
    return result.stdout.split("\n")[0]

def annotate_merged(
    sample_id: str,
    datadir: Path,
    bcftools: str,
    reference: Path,
    gff: Path,
    transcript_list: Optional[Path] = None,
    max_retries: int = 3,
) -> str:
    vcf_dir = Path(datadir) / "BAM" / sample_id / "VCF"
    input_vcf = vcf_dir / f"{sample_id}_merged.vcf"
    output_vcf = vcf_dir / f"{sample_id}_merged.csq.vcf"
    final_vcf = vcf_dir / f"{sample_id}_merged.csq.filtered.vcf"

    if final_vcf.exists():
        return "exists"

    bcftools_cmd = [
        bcftools,
        "csq",
        "-f",
        str(reference),
        "-g",
        str(gff),
        "--local-csq",
        "--ncsq",
        "20",
        "-Ov",
        "-o",
        str(output_vcf),
        str(input_vcf),
    ]

    for attempt in range(max_retries):
        process = subprocess.run(bcftools_cmd, capture_output=True, text=True)
        
        if process.returncode == 0:
            break
        elif attempt == max_retries - 1:
            return "error"

    if output_vcf.exists():
        if transcript_list:
            try:
                with open(transcript_list, "r") as f:
                    desired_transcripts = set(line.strip() for line in f)

                with open(output_vcf, "r") as infile, open(final_vcf, "w") as outfile:
                    for line in infile:
                        if line.startswith("#"):
                            outfile.write(line)
                            continue

                        fields = line.strip().split("\t")
                        info = fields[7]
                        csq_field = next((f for f in info.split(";") if f.startswith("BCSQ=")), None)

                        if csq_field:
                            csq_parts = csq_field.split("|")
                            transcript = csq_parts[4] if len(csq_parts) > 4 else None

                            if transcript in desired_transcripts:
                                outfile.write(line)
                        else:
                            outfile.write(line)

                final_vcf.replace(output_vcf)
            except Exception:
                return "error"
        else:
            final_vcf = output_vcf

        return "success"
    else:
        return "error"