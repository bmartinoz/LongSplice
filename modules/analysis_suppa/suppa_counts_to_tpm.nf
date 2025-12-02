// modules/analysis_suppa/suppa_counts_to_tpm.nf
nextflow.enable.dsl=2

process SUPPA_COUNTS_TO_TPM {
  label 'process_low'

  publishDir "${params.outdir ?: 'results'}/suppa",
             mode: 'copy', overwrite: true

  input:
  path counts_tx
  path gtf

  output:
  path "bambu.tpm.tsv", emit: tpm_tsv

  shell:
  '''
  set -euo pipefail

  GTF_IN="!{gtf}"
  CTS_IN="!{counts_tx}"

  if [[ -z "$GTF_IN" || ! -s "$GTF_IN" ]]; then
      GTF_IN="$(ls -1 *.gtf 2>/dev/null | head -n1 || true)"
  fi
  if [[ -z "$CTS_IN" || ! -s "$CTS_IN" ]]; then
      CTS_IN="$(ls -1 *counts_transcript*.txt 2>/dev/null | head -n1 || true)"
  fi

  [[ -s "$GTF_IN" ]] || { echo '[SUPPA] ERROR: no se encontró GTF para SUPPA_COUNTS_TO_TPM.' >&2; ls -la || true; exit 1; }
  [[ -s "$CTS_IN" ]] || { echo '[SUPPA] ERROR: no se encontró counts_transcript.txt para SUPPA_COUNTS_TO_TPM.' >&2; ls -la || true; exit 1; }

  python3 - "$GTF_IN" "$CTS_IN" << 'PY'
import sys, os

TAB = chr(9)
NL  = chr(10)

gtf_path = sys.argv[1]
cts_path = sys.argv[2]

# ---------- Longitudes por transcript ----------
lens = {}
with open(gtf_path, 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.rstrip().split(TAB)
        if len(parts) < 9 or parts[2] != 'exon':
            continue
        start, end = int(parts[3]), int(parts[4])
        info = parts[8]
        key = 'transcript_id "'
        i = info.find(key)
        if i == -1:
            continue
        i += len(key)
        j = info.find('"', i)
        if j == -1:
            continue
        tid = info[i:j]
        L = end - start + 1
        if L > 0:
            lens[tid] = lens.get(tid, 0) + L

# ---------- Lee counts y normaliza a TPM ----------
rows = []
with open(cts_path, 'r', encoding='utf-8', errors='ignore') as f:
    header = f.readline().rstrip().split(TAB)
    if len(header) < 3 or header[0].upper() not in ('TXNAME','TRANSCRIPT','TRANSCRIPT_ID'):
        sys.exit("counts_transcript.txt: cabecera inesperada: " + str(header[:3]))
    samples = header[2:]  # tal como viene de bambu
    for line in f:
        if not line.strip():
            continue
        parts = line.rstrip().split(TAB)
        tid = parts[0]
        vals = [float(x) if x not in ('','NA') else 0.0 for x in parts[2:]]
        rows.append((tid, vals))

def order_key(sample_name):
    if sample_name.endswith(".sorted"):
        head = sample_name[:-7]
        if "_" in head:
            grp, num = head.split("_", 1)
            if grp in ("KO","WT") and num.isdigit():
                return (0, 0 if grp=="KO" else 1, int(num), sample_name)
    return (1, 9, 999999, sample_name)

order_idx = sorted(range(len(samples)), key=lambda i: order_key(samples[i]))
samples_sorted = [samples[i] for i in order_idx]

# ---------- TPM ----------
denom = [0.0]*len(samples)
num   = []
for tid, vec in rows:
    L = float(lens.get(tid, 1.0))
    lk = max(1.0, L) / 1000.0   # kb
    w  = [x / lk for x in vec]  # RPK
    num.append((tid, w))
    for i, x in enumerate(w):
        denom[i] += x

def to_tpm(vec):
    return [0.0 if denom[i]==0 else (vec[i]/denom[i]*1e6) for i in range(len(vec))]

# ---- Salida SUPPA-friendly (sin ID en header) ----
with open('bambu.tpm.tsv','w', encoding='utf-8') as o:
    o.write(TAB.join(samples_sorted) + NL)
    for tid, w in num:
        tpm = to_tpm(w)
        tpm_sorted = [tpm[i] for i in order_idx]
        o.write(tid + TAB + TAB.join(str(x) for x in tpm_sorted) + NL)

# ---- Copia con ID ----
with open('bambu.tpm.with_id.tsv','w', encoding='utf-8') as o:
    o.write('transcript_id' + TAB + TAB.join(samples_sorted) + NL)
    for tid, w in num:
        tpm = to_tpm(w)
        tpm_sorted = [tpm[i] for i in order_idx]
        o.write(tid + TAB + TAB.join(str(x) for x in tpm_sorted) + NL)
PY
  '''
}
