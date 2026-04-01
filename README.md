# prostate-str-pipeline

Automated AR CAG repeat genotyping pipeline for prostate cancer WES samples using STRipy and Seven Bridges CGC

# Prostate Cancer AR STR Analysis Pipeline
 
Automated pipeline for genotyping androgen receptor (AR) CAG repeat expansions across whole-exome sequencing (WES) samples from prostate cancer patients.
 
Built as part of a research project at Georgetown University investigating prostate cancer disparities.
 
---
 
## Overview
 
This pipeline integrates STRipy, ExpansionHunter, REViewer, samtools, and the Seven Bridges Cancer Genomics Cloud (CGC) API to automate batch STR genotyping from cloud-stored BAM files — without requiring large local storage.
 
```
CGC BAM Files → Download → Index → STRipy/ExpansionHunter → REViewer → Pileup → HTML Report → Delete BAM
```
 
---
 
## Features
 
- **Cloud-integrated** — reads BAM files directly from Seven Bridges CGC via API
- **Storage-efficient** — downloads and deletes one BAM at a time; only one file on disk at any point
- **Parallel downloading** — next BAM downloads in background while current sample is being analyzed
- **AR CAG repeat genotyping** — uses STRipy + ExpansionHunter with male sex parameter for prostate samples
- **Smart REViewer** — only runs read alignment visualization for flagged (non-normal) samples
- **Interactive pileup** — Plotly-based zoomable read depth graph per sample
- **Merged HTML reports** — genotyping summary, STRipy report, pileup, and REViewer all in one file
- **Flask web dashboard** — browser-based UI with folder browser, loci selection, progress bar, and download buttons
- **One-click launcher** — Windows Desktop `.bat` shortcut to start the app
 
---
 
## Requirements
 
### System
- Windows with WSL2 (Ubuntu)
- Docker Desktop (WSL2 integration enabled)
- Python 3.12+
 
### Python Packages
```bash
pip install sevenbridges-python flask plotly matplotlib --break-system-packages
```
 
### Tools (in WSL)
```bash
sudo apt-get install samtools -y
```
 
### Other
- STRipy Docker image (`stripy:v2.5`) — built automatically on first run via `runDocker.sh`
- Reference genome: `hg38.fa` with `.fai` index stored at `/home/<user>/ref/hg38.fa`
- Seven Bridges CGC account with API token
 
---
 
## File Structure
 
```
~/
├── stripy-pipeline/              # This repo (cloned from GitLab)
│   ├── runDocker.sh              # Docker launcher script
│   ├── config.json               # STRipy configuration
│   └── stripy_results/           # Intermediate output (not committed)
├── ref/
│   └── hg38.fa                   # Reference genome (not committed)
├── bam/                          # Temporary BAM storage (auto-deleted)
└── stripy_webapp/
    ├── app.py                    # Flask backend
    └── templates/
        └── index.html            # Web dashboard UI
```
 
---
 
## Setup
 
**1. Clone STRipy pipeline:**
```bash
git clone https://gitlab.com/andreassh/stripy-pipeline.git ~/stripy-pipeline
chmod +x ~/stripy-pipeline/runDocker.sh
```
 
**2. Install Python dependencies:**
```bash
pip install sevenbridges-python flask plotly matplotlib --break-system-packages
```
 
**3. Set up reference genome:**
```bash
mkdir -p ~/ref
# Copy your hg38.fa and hg38.fa.fai to ~/ref/
```
 
**4. Create required folders:**
```bash
mkdir -p ~/bam
mkdir -p ~/stripy-pipeline/stripy_results
```
 
**5. Set up the web app:**
```bash
mkdir -p ~/stripy_webapp/templates
# Copy app.py to ~/stripy_webapp/app.py
# Copy index.html to ~/stripy_webapp/templates/index.html
```
 
**6. Create Desktop launcher (optional):**
```bash
echo 'wsl -e bash -c "cd ~/stripy_webapp && python3 app.py; exec bash"' > "/mnt/c/Users/<username>/OneDrive/Desktop/Start_STRipy.bat"
```
 
---
 
## Usage
 
### Web App (Recommended)
 
1. Open **Docker Desktop** and wait for it to load
2. Double-click **Start_STRipy.bat** on your Desktop (or run `cd ~/stripy_webapp && python3 app.py` in WSL)
3. Open **http://localhost:5000** in your browser
4. Enter your CGC API token
5. Select a folder from the dropdown (or click Refresh to fetch latest)
6. Select loci (default: AR)
7. Choose reference genome and sample sex
8. Click **▶ Run Pipeline**
9. Monitor progress in the live log
10. Download reports from the Results section when complete
 
### Command Line (Single Sample)
 
```bash
cd ~/stripy-pipeline
 
# Download and index BAM
samtools index /path/to/sample.bam
 
# Run STRipy
./runDocker.sh \
  -g hg38 \
  -r /home/<user>/ref/hg38.fa \
  -o /home/<user>/stripy-pipeline/stripy_results/ \
  -l AR \
  -s male \
  -i /path/to/sample.bam
```
 
### Batch Processing via Python
 
```python
import sevenbridges as sbg
 
api = sbg.Api(url='https://cgc-api.sbgenomics.com/v2', token='YOUR_TOKEN')
folder = api.files.get('YOUR_FOLDER_ID')
files = list(api.files.query(parent=folder.id))
bam_files = [f for f in files if f.name.endswith('.sorted.bam')]
 
for bam in bam_files:
    bam.download('/home/<user>/bam/sample.bam')
    # ... run STRipy, generate report, delete BAM
```
 
---
 
## Output
 
Each sample produces a single `_full_report.html` file containing:
 
| Section | Description |
|---------|-------------|
| Genotyping Summary | AR repeat count, confidence interval, disease range, outlier status |
| STRipy Full Report | Complete STRipy HTML output with read visualizations |
| Read Depth Pileup | Interactive Plotly chart of coverage across AR locus |
| REViewer Alignment | Read alignment SVG (only for flagged/expanded samples) |
 
### AR CAG Repeat Ranges
 
| Classification | Repeat Count | Notes |
|----------------|-------------|-------|
| Normal | 9–36 | No disease association |
| Intermediate | 37–38 | Uncertain significance |
| Pathogenic | ≥38 | Spinal and bulbar muscular atrophy (SBMA) |
 
---
 
## Configuration
 
Key settings in `config.json`:
 
```json
{
  "num_threads": "all",
  "output_html": true,
  "output_tsv": true,
  "output_svg_flag_threshold": 0
}
```
 
| Parameter | Default | Description |
|-----------|---------|-------------|
| `num_threads` | `"all"` | CPU threads for analysis |
| `output_html` | `true` | Generate HTML report |
| `output_tsv` | `true` | Generate TSV summary |
| `output_svg_flag_threshold` | `0` | Generate read images for all loci (0) or flagged only (1-3) |
 
---
 
## Tools & References
 
- [STRipy-pipeline](https://gitlab.com/andreassh/stripy-pipeline) — Halman et al., 2022
- [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) — Dolzhenko et al., 2017
- [REViewer](https://github.com/Illumina/REViewer) — Dolzhenko et al., 2022
- [Seven Bridges CGC](https://cgc.sbgenomics.com)
- [Protocol reference](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.70010) — Halman, Lonsdale & Oshlack, 2024
 
---
 
## Notes
 
- **Token security**: Never commit your CGC API token. Enter it at runtime in the web app.
- **Patient data**: BAM files are never permanently stored locally — they are downloaded, processed, and immediately deleted.
- **Sex parameter**: Always use `-s male` for prostate cancer samples. AR is X-linked; male samples are hemizygous (one allele).
- **Sleep prevention**: Disable laptop sleep before running large batches to prevent interruptions.
