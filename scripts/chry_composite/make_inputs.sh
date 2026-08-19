#!/bin/bash
# make_inputs.sh — regenerates data/ for chry_composite.R from public data.
#
# Needs samtools, mosdepth, bwa, python3 and about 4 GB of scratch. Everything
# it fetches is public: BAM slices come out of the Platinum Pedigree Consortium
# bucket over HTTP range requests, so only chrY and a 5 Mb chr20 window are
# transferred, not whole genomes.
#
#   usage:  ./make_inputs.sh [workdir]        (default: ./work)
#
# Writes work/data/{windows.tsv,counts.tsv,cytoband_chrY.tsv,*.parwin.regions.bed.gz}
# and the pileups analysis.py reads.

set -euo pipefail
WORK=${1:-work}
DATA=$WORK/data
BUCKET=https://platinum-pedigree-data.s3.amazonaws.com/data/illumina/mapped/GRCh38
mkdir -p "$WORK" "$DATA"
cd "$WORK"

# CEPH 1463: secondary id -> who they are. NA12889 is NA12877's father, so his Y
# is the one NA12878's sons carry; NA12891 is NA12878's father, an unrelated Y.
SAMPLES="2188:NA12878_mother 2212:NA12882_son 2189:NA12886_son2 2209:NA12877_father
         2216:NA12879_daughter 2281:NA12889_unrel_male 2214:NA12891_her_father"

# ---------------------------------------------------------------- reference
# The analysis-set chrY: UCSC hg38 chrY with both PARs hard-masked, which
# reproduces the 1000 Genomes chrY (md5 ce3e31103314a704255f3cd90369ecce).
if [ ! -f chrY.analysis.fa ]; then
  curl -sO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrY.fa.gz
  python3 - <<'PY'
import gzip
seq = []
with gzip.open('chrY.fa.gz', 'rt') as f:
    f.readline()
    for line in f:
        seq.append(line.strip())
s = list(''.join(seq).upper())
for a, b in [(10001, 2781479), (56887903, 57217415)]:
    for i in range(a - 1, b):
        s[i] = 'N'
s = ''.join(s)
with open('chrY.analysis.fa', 'w') as o:
    o.write('>chrY\n')
    for i in range(0, len(s), 60):
        o.write(s[i:i + 60] + '\n')
PY
  samtools faidx chrY.analysis.fa
fi

# chrX + chrY with PAR left in place, for the realignment in panel B
if [ ! -f XY.fa.bwt ]; then
  curl -sO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrX.fa.gz
  { zcat chrX.fa.gz; zcat chrY.fa.gz; } | awk '{print toupper($0)}' \
    | sed 's/^>CHR/>chr/' > XY.fa
  samtools faidx XY.fa
  bwa index XY.fa 2> bwa_index.log
fi

curl -s "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cytoBand;chrom=chrY" \
  | python3 -c '
import json, sys
print("start\tend\tname\tstain")
for b in json.load(sys.stdin)["cytoBand"]:
    print(b["chromStart"], b["chromEnd"], b["name"], b["gieStain"], sep="\t")' \
  > "$DATA/cytoband_chrY.tsv"

# ---------------------------------------------------------------- slices
: > counts.tsv
for pair in $SAMPLES; do
  id=${pair%%:*}; s=${pair##*:}
  [ -f "$s.chrY.bam" ] || {
    samtools view -b -@2 -o "$s.chrY.bam" "$BUCKET/$id.v4.2.4.grc38.bam" chrY
    samtools index "$s.chrY.bam"
  }
  [ -f "$s.chr20.bam" ] || {
    samtools view -b -@2 -o "$s.chr20.bam" "$BUCKET/$id.v4.2.4.grc38.bam" chr20:30000000-35000000
    samtools index "$s.chr20.bam"
  }
  mosdepth -n -x --by 10000 -c chrY "$s.win" "$s.chrY.bam"
  mosdepth -n -x --by 10000 -c chrY -Q 30 "$s.winq30" "$s.chrY.bam"
  tot=$(samtools view -c -F 0x900 "$s.chrY.bam")
  q30=$(samtools view -c -F 0x900 -q 30 "$s.chrY.bam")
  auto=$(samtools depth -a -r chr20:30000000-35000000 "$s.chr20.bam" \
         | awk '{t += $3} END {printf "%.4f", t / NR}')
  printf '%s\t%s\t%s\t%s\n' "$s" "$tot" "$q30" "$auto" >> counts.tsv
done
cp counts.tsv "$DATA/counts.tsv"

# ---------------------------------------------------------------- 10 kb windows
python3 - <<'PY'
import gzip
W = 10000
seq = []
with open('chrY.analysis.fa') as f:
    f.readline()
    for line in f:
        seq.append(line.strip())
seq = ''.join(seq)
nfrac = {s: seq[s:s + W].count('N') / max(1, len(seq[s:s + W]))
         for s in range(0, len(seq), W)}

XTR = (3_050_000, 6_234_000); CEN = (10_000_000, 12_500_000)
SAT2 = (56_670_000, 56_890_000); PAR1 = (10_000, 2_781_479)
PAR2 = (56_887_902, 57_217_415); MSY_END = 26_673_214

def klass(start):
    mid = start + W // 2
    if PAR1[0] <= mid < PAR1[1]: return 'PAR1(masked)'
    if PAR2[0] <= mid < PAR2[1]: return 'PAR2(masked)'
    if XTR[0] <= mid < XTR[1]: return 'XTR'
    if CEN[0] <= mid < CEN[1]: return 'centromeric satellite'
    if SAT2[0] <= mid < SAT2[1]: return 'proximal-PAR2 satellite'
    return 'MSY unique' if mid < MSY_END else 'Yq12 gap'

samples = [l.split('\t')[0] for l in open('counts.tsv')]
auto = {l.split('\t')[0]: float(l.split('\t')[3]) for l in open('counts.tsv')}
depth = {}
for s in samples:
    d = {}
    with gzip.open(f'{s}.win.regions.bed.gz', 'rt') as f:
        for line in f:
            c, a, b, v = line.split()
            if c == 'chrY':
                d[int(a)] = float(v)
    depth[s] = d

with open('data/windows.tsv', 'w') as o:
    o.write('start\tclass\tnfrac\t' + '\t'.join(samples) + '\n')
    for start in sorted(nfrac):
        row = [str(start), klass(start), f'{nfrac[start]:.3f}']
        row += [f'{depth[s].get(start, 0) / auto[s]:.5g}' for s in samples]
        o.write('\t'.join(row) + '\n')
PY

# ---------------------------------------------------------------- panel B overlay
# PAR is hard-masked on chrY in these BAMs, so the pale tracks come from pulling
# each sample's chrX PAR1 reads and remapping them to a reference that keeps both
# copies. Roughly half land on chrY, nearly all at MAPQ 0.
for s in NA12878_mother NA12882_son; do
  id=$([ "$s" = NA12878_mother ] && echo 2188 || echo 2212)
  [ -f "$s.par1.bam" ] || samtools view -b -@2 -o "$s.par1.bam" \
    "$BUCKET/$id.v4.2.4.grc38.bam" chrX:10001-2781479
  samtools collate -u -O -@2 "$s.par1.bam" \
    | samtools fastq -n -0 /dev/null -s /dev/null - 2>/dev/null \
    | bwa mem -p -t 8 XY.fa - 2>/dev/null \
    | samtools sort -@2 -o "$s.par1.xy.bam" -
  samtools index "$s.par1.xy.bam"
  mosdepth -n -x --by 10000 -c chrY "$s.parwin" "$s.par1.xy.bam"
  cp "$s.parwin.regions.bed.gz" "$DATA/"
done

# ---------------------------------------------------------------- pileups
# analysis.py turns these into the counts quoted in panels C and D.
for s in NA12878_mother NA12879_daughter NA12882_son NA12886_son2 \
         NA12889_unrel_male NA12891_her_father; do
  samtools mpileup -f chrY.analysis.fa -r chrY:13000000-20000000 -q 0 -Q 20 -B \
    "$s.chrY.bam" 2>/dev/null | awk '{print $2"\t"$3"\t"$4"\t"$5}' > "$s.msy.pile"
  samtools mpileup -f chrY.analysis.fa -r chrY:3050000-6234000 -q 0 -Q 20 -B \
    "$s.chrY.bam" 2>/dev/null \
    | awk '{n = split($5, a, ""); alt = 0;
            for (i = 1; i <= n; i++) { c = toupper(a[i]);
              if (c == "A" || c == "C" || c == "G" || c == "T") alt++ }
            print $2"\t"$4"\t"alt}' > "$s.xtr.pile"
done

echo "data/ written — now: python3 analysis.py $WORK && Rscript chry_composite.R $DATA"
