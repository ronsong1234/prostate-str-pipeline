from flask import Flask, render_template, request, jsonify, send_file
import subprocess, os, threading
import sevenbridges as sbg
import plotly.graph_objects as go
from concurrent.futures import ThreadPoolExecutor

app = Flask(__name__)

RESULTS_DIR = '/home/USER/stripy-pipeline/stripy_results'  # Update with your WSL username
OUTPUT_DIR = '/mnt/c/Users/USER/OneDrive/Desktop/STRipy Results'  # Update with your Windows username

import shutil
BAM_DIR = '/home/USER/bam'  # Update with your WSL username
REF = '/home/USER/ref/hg38.fa'  # Update with your WSL username
PIPELINE = '/home/USER/stripy-pipeline/runDocker.sh'  # Update with your WSL username
CATALOG = '/usr/local/bin/stripy-pipeline/catalog.json'

LOCI = [
    'AR', 'HTT', 'DMPK', 'ATXN1', 'ATXN2', 'ATXN3', 'ATXN7', 'ATXN10',
    'FMR1', 'FXN', 'C9ORF72', 'CNBP', 'NOP56', 'TBP', 'CACNA1A',
    'JPH3', 'PPP2R2B', 'RFC1', 'NOTCH2NLC', 'PABPN1'
]

FOLDERS = [
    # Add your CGC folder IDs here after running list_folders()
    # Example: {'name': 'MyProject / OutputFolder', 'id': 'your_folder_id_here'},
]

log_messages = []
is_running = False
progress = {'current': 0, 'total': 0, 'sample': ''}


def log(msg):
    log_messages.append(msg)
    print(msg)


def download_and_index(bam, dest, all_files=None):
    log(f'Downloading {bam.name}...')
    bam.download(dest)

    # Try to find and download matching .bai file
    bai_downloaded = False
    if all_files:
        bai_name = bam.name + '.bai'
        bai = next((f for f in all_files if f.name == bai_name), None)
        if bai:
            log(f'Downloading index for {bam.name}...')
            bai.download(dest + '.bai')
            bai_downloaded = True

    # Fall back to samtools index if no .bai found
    if not bai_downloaded:
        log(f'Indexing {bam.name}...')
        result = subprocess.run(['samtools', 'index', dest])
        if result.returncode != 0:
            log(f'WARNING: Indexing failed for {bam.name} - file may be corrupted')
            return False
    return True


def is_flagged(bam_name):
    """Check TSV to see if any locus is outside normal range."""
    tsv_path = f'{RESULTS_DIR}/sample.bam.tsv'
    if not os.path.exists(tsv_path):
        return False
    with open(tsv_path) as f:
        lines = f.readlines()
    if len(lines) < 2:
        return False
    headers = lines[0].strip().split('\t')
    for line in lines[1:]:
        values = line.strip().split('\t')
        row = dict(zip(headers, values))
        flag = row.get('Flag', '0')
        if flag and flag != '0' and flag != 'normal':
            return True
    return False


def process_sample(bam_name, loci_str, genome, sex='male'):
    log(f'Running STRipy on {bam_name} (sex: {sex})...')
    subprocess.run([
        'bash', PIPELINE,
        '-g', genome,
        '-r', REF,
        '-o', f'{RESULTS_DIR}/',
        '-l', loci_str,
        '-s', sex,
        '-i', f'{BAM_DIR}/sample.bam'
    ])

    # Only run REViewer if sample is flagged (outside normal range)
    if is_flagged(bam_name):
        log(f'Sample flagged - running REViewer on {bam_name}...')
        subprocess.run([
            'docker', 'run', '--rm',
            '-v', f'{BAM_DIR}:/mnt/data',
            '-v', '/home/ronak/ref:/mnt/ref',
            '-v', f'{RESULTS_DIR}:/mnt/results',
            'stripy:v2.5', 'bash', '-c',
            f'REViewer '
            f'--reads /mnt/data/sample.bam '
            f'--vcf /mnt/results/sample.bam.vcf '
            f'--reference /mnt/ref/hg38.fa '
            f'--catalog {CATALOG} '
            f'--locus {loci_str} '
            f'--output-prefix /mnt/results/sample_reviewer'
        ])
    else:
        log(f'Normal range - skipping REViewer for {bam_name}')

    log(f'Generating pileup for {bam_name}...')
    pileup_file = f'{RESULTS_DIR}/AR_pileup.txt'
    subprocess.run([
        'samtools', 'mpileup',
        '-r', 'chrX:67545316-67545395',
        '-f', REF,
        f'{BAM_DIR}/sample.bam'
    ], stdout=open(pileup_file, 'w'))

    positions, depths = [], []
    with open(pileup_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                positions.append(int(parts[1]))
                depths.append(int(parts[3]))

    pileup_html = ''
    if positions:
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=positions, y=depths,
            fill='tozeroy',
            fillcolor='rgba(70,130,180,0.4)',
            line=dict(color='steelblue'),
            hovertemplate='Position: %{x}<br>Depth: %{y}<extra></extra>'
        ))
        fig.update_layout(
            title='AR CAG Repeat Locus - Read Depth',
            xaxis_title='Genomic Position (chrX)',
            yaxis_title='Read Depth',
            hovermode='x unified',
            template='plotly_white'
        )
        pileup_html = fig.to_html(full_html=False, include_plotlyjs='cdn')

    stripy_html = ''
    stripy_path = f'{RESULTS_DIR}/sample.bam.html'
    if os.path.exists(stripy_path):
        with open(stripy_path) as f:
            stripy_html = f.read()

    reviewer_svg = ''
    reviewer_path = f'{RESULTS_DIR}/sample_reviewer.AR.svg'
    if os.path.exists(reviewer_path):
        with open(reviewer_path) as f:
            reviewer_svg = f.read()

    tsv_table = ''
    tsv_path = f'{RESULTS_DIR}/sample.bam.tsv'
    if os.path.exists(tsv_path):
        with open(tsv_path) as f:
            lines = f.readlines()
        if len(lines) >= 2:
            headers = lines[0].strip().split('\t')
            values = lines[1].strip().split('\t')
            tsv_table = '<table style="border-collapse:collapse;width:100%">'
            for h, v in zip(headers, values):
                tsv_table += (
                    f'<tr>'
                    f'<td style="padding:6px;border:1px solid #ddd;font-weight:bold;background:#f5f5f5">{h}</td>'
                    f'<td style="padding:6px;border:1px solid #ddd">{v}</td>'
                    f'</tr>'
                )
            tsv_table += '</table>'

    merged_html = f'''<!DOCTYPE html>
<html>
<head>
    <title>STRipy Report - {bam_name}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 30px; background: #fafafa; }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #2980b9; padding-bottom: 10px; }}
        h2 {{ color: #2980b9; margin-top: 40px; }}
        .section {{ background: white; padding: 20px; margin-bottom: 30px; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,0.1); }}
        .sample-name {{ color: #7f8c8d; font-size: 14px; }}
    </style>
</head>
<body>
    <h1>STRipy Analysis Report</h1>
    <p class="sample-name">Sample: {bam_name}</p>
    <div class="section"><h2>Genotyping Summary</h2>{tsv_table or "<p>No data</p>"}</div>
    <div class="section"><h2>STRipy Full Report</h2>{stripy_html}</div>
    <div class="section"><h2>Read Depth Pileup (Interactive)</h2>{pileup_html or "<p>No pileup data</p>"}</div>
    <div class="section"><h2>REViewer Read Alignment</h2>{reviewer_svg or "<p>No REViewer image (normal range)</p>"}</div>
</body>
</html>'''

    # Save to WSL results dir first (where Docker can access)
    merged_path = f'{RESULTS_DIR}/{bam_name}_full_report.html'
    with open(merged_path, 'w') as f:
        f.write(merged_html)
    # Copy final report to Windows Desktop folder
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        shutil.copy2(merged_path, f'{OUTPUT_DIR}/{bam_name}_full_report.html')
    except Exception as e:
        log(f'Warning: Could not copy to Desktop: {str(e)}')

    for temp in [stripy_path, reviewer_path, pileup_file, tsv_path]:
        if os.path.exists(temp):
            os.remove(temp)


def run_pipeline(token, folder_id, loci, genome, sex='male'):
    global is_running
    is_running = True
    log_messages.clear()

    try:
        log('Connecting to CGC...')
        api = sbg.Api(url='https://cgc-api.sbgenomics.com/v2', token=token)
        folder = api.files.get(folder_id)
        # Fetch all files with pagination
        files = []
        offset = 0
        while True:
            batch = list(api.files.query(parent=folder.id, offset=offset, limit=100))
            files.extend(batch)
            if len(batch) < 100:
                break
            offset += 100
        # Pick best BAM per sample
        # Priority: sorted.indexed.bam > sorted.bam > .bam
        all_bams = [f for f in files if f.name.endswith('.bam') and f.type == 'file']
        indexed_bams = [f for f in all_bams if f.name.endswith('.sorted.indexed.bam')]
        sorted_bams = [f for f in all_bams if f.name.endswith('.sorted.bam') and not f.name.endswith('.sorted.indexed.bam')]

        if indexed_bams:
            bam_files = indexed_bams
            log(f'Using sorted.indexed.bam files')
        elif sorted_bams:
            bam_files = sorted_bams
            log(f'Using sorted.bam files')
        else:
            bam_files = all_bams
            log(f'Using all BAM files')

        bam_files.sort(key=lambda f: f.name)
        log(f'Found {len(bam_files)} BAM files to process.')
        for b in bam_files:
            log(f'  - {b.name}')

        loci_str = ','.join(loci)

        with ThreadPoolExecutor(max_workers=1) as executor:
            next_future = None

            for i, bam in enumerate(bam_files):
                # Use pre-downloaded file or download first one
                if next_future:
                    success = next_future.result()
                    if not success:
                        log(f'Skipping corrupted file: {bam.name}')
                        for f in [f'{BAM_DIR}/next.bam', f'{BAM_DIR}/next.bam.bai']:
                            if os.path.exists(f):
                                os.remove(f)
                        continue
                    os.rename(f'{BAM_DIR}/next.bam', f'{BAM_DIR}/sample.bam')
                    if os.path.exists(f'{BAM_DIR}/next.bam.bai'):
                        os.rename(f'{BAM_DIR}/next.bam.bai', f'{BAM_DIR}/sample.bam.bai')
                else:
                    success = download_and_index(bam, f'{BAM_DIR}/sample.bam', files)
                    if not success:
                        log(f'Skipping {bam.name} - corrupted')
                        continue

                # Start pre-downloading next sample in background
                if i + 1 < len(bam_files):
                    next_bam = bam_files[i + 1]
                    for f in [f'{BAM_DIR}/next.bam', f'{BAM_DIR}/next.bam.bai']:
                        if os.path.exists(f):
                            os.remove(f)
                    next_future = executor.submit(download_and_index, next_bam, f'{BAM_DIR}/next.bam', files)
                else:
                    next_future = None

                # Process current sample
                progress['current'] = i + 1
                progress['sample'] = bam.name
                process_sample(bam.name, loci_str, genome, sex)

                # Cleanup current BAM
                for f in [f'{BAM_DIR}/sample.bam', f'{BAM_DIR}/sample.bam.bai']:
                    if os.path.exists(f):
                        os.remove(f)

                log(f'Done with {bam.name}! ({i+1}/{len(bam_files)})')

        log('All samples complete!')

    except Exception as e:
        log(f'ERROR: {str(e)}')
    finally:
        is_running = False


@app.route('/')
def index():
    return render_template('index.html', loci=LOCI, folders=FOLDERS)

@app.route('/get_folders')
def get_folders():
    return jsonify({'folders': FOLDERS})


@app.route('/refresh_folders', methods=['POST'])
def refresh_folders():
    try:
        data = request.json
        token = data.get('token')
        if not token:
            return jsonify({'error': 'Token required'})

        api = sbg.Api(url='https://cgc-api.sbgenomics.com/v2', token=token)

        def get_all_files(project=None, parent=None):
            all_files = []
            offset = 0
            limit = 100
            while True:
                if parent:
                    batch = list(api.files.query(parent=parent, offset=offset, limit=limit))
                else:
                    batch = list(api.files.query(project=project, offset=offset, limit=limit))
                all_files.extend(batch)
                if len(batch) < limit:
                    break
                offset += limit
            return all_files

        def get_all_folders(project=None, parent=None, prefix=''):
            files = get_all_files(project=project, parent=parent)
            result = []
            for f in files:
                if f.type == 'folder':
                    display_name = f'{prefix}{f.name}' if prefix else f.name
                    result.append({'id': f.id, 'name': display_name})
                    subfolders = get_all_folders(parent=f.id, prefix=f'{display_name} / ')
                    result.extend(subfolders)
            return result

        # Query all accessible projects
        projects = list(api.projects.query())
        all_folders = []
        for p in projects:
            try:
                folders = get_all_folders(project=p.id)
                all_folders.extend(folders)
            except Exception:
                pass

        return jsonify({'folders': all_folders, 'count': len(all_folders)})
    except Exception as e:
        return jsonify({'error': str(e)})


@app.route('/projects', methods=['POST'])
def projects():
    try:
        data = request.json
        token = data.get('token')
        api = sbg.Api(url='https://cgc-api.sbgenomics.com/v2', token=token)
        project_list = list(api.projects.query())
        return jsonify({'projects': [{'id': p.id, 'name': p.name} for p in project_list]})
    except Exception as e:
        return jsonify({'error': str(e)})


@app.route('/folders', methods=['POST'])
def folders():
    try:
        data = request.json
        token = data.get('token')
        project_id = data.get('project_id')
        parent_id = data.get('parent_id')
        api = sbg.Api(url='https://cgc-api.sbgenomics.com/v2', token=token)

        def get_all_folders(parent, prefix=''):
            try:
                files = list(api.files.query(parent=parent))
            except Exception:
                return []
            result = []
            for f in files:
                if f.type == 'folder':
                    display_name = f'{prefix}{f.name}' if prefix else f.name
                    result.append({'id': f.id, 'name': display_name})
                    subfolders = get_all_folders(f.id, prefix=f'{display_name}/')
                    result.extend(subfolders)
            return result

        if parent_id:
            folder_list = get_all_folders(parent_id)
        else:
            top_files = list(api.files.query(project=project_id))
            folder_list = []
            for f in top_files:
                if f.type == 'folder':
                    folder_list.append({'id': f.id, 'name': f.name})
                    subfolders = get_all_folders(f.id, prefix=f'{f.name}/')
                    folder_list.extend(subfolders)

        return jsonify({'folders': folder_list})
    except Exception as e:
        return jsonify({'error': str(e)})


@app.route('/run', methods=['POST'])
def run():
    global is_running
    if is_running:
        return jsonify({'status': 'error', 'message': 'Pipeline already running!'})

    data = request.json
    token = data.get('token')
    folder_id = data.get('folder_id')
    loci = data.get('loci', ['AR'])
    genome = data.get('genome', 'hg38')
    sex = data.get('sex', 'male')

    if not token or not folder_id:
        return jsonify({'status': 'error', 'message': 'Token and Folder ID are required!'})

    thread = threading.Thread(target=run_pipeline, args=(token, folder_id, loci, genome, sex))
    thread.daemon = True
    thread.start()

    return jsonify({'status': 'ok', 'message': 'Pipeline started!'})


@app.route('/logs')
def logs():
    return jsonify({'logs': log_messages, 'running': is_running})


@app.route('/progress')
def get_progress():
    pct = 0
    if progress['total'] > 0:
        pct = int((progress['current'] / progress['total']) * 100)
    return jsonify({
        'current': progress['current'],
        'total': progress['total'],
        'sample': progress['sample'],
        'percent': pct,
        'running': is_running
    })


@app.route('/results')
def results():
    files = []
    if os.path.exists(OUTPUT_DIR):
        files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('_full_report.html')]
    return jsonify({'files': sorted(files)})


@app.route('/download/<filename>')
def download(filename):
    path = os.path.join(OUTPUT_DIR, filename)
    if os.path.exists(path):
        return send_file(path, as_attachment=True)
    return 'File not found', 404


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
