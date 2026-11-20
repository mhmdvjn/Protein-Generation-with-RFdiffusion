#!/usr/bin/env bash
set -euo pipefail

CONDITIONS=("cond_a_canonical" "cond_b_alternative" "cond_c_control")
MPNN_CKPT="${HOME}/.foundry/checkpoints/ligandmpnn_v_32_010_25.pt"
RF3_CKPT="${HOME}/.foundry/checkpoints/rf3_foundry_01_24_latest_remapped.ckpt"

# Stage 1: RFD3 All-Atom Generative Diffusion
# (Matched budget: adjust n_batches & diffusion_batch_size so all 3 folders have equal structure counts)
for COND in "${CONDITIONS[@]}"; do
  echo "=== [Stage 1] Running RFD3 for ${COND} ==="
  rfd3 design \
    out_dir="test/results/structures/${COND}" \
    inputs="test/configs/${COND}.yaml" \
    n_batches=2 \
    diffusion_batch_size=2 \
    prevalidate_inputs=True \
    skip_existing=True
done

Stage 2: Inverse Folding with LigandMPNN (4 sequences per backbone, target Chain A fixed)
for COND in "${CONDITIONS[@]}"; do
  echo "=== [Stage 2] Running MPNN for ${COND} ==="
  mkdir -p "test/results/mpnn/${COND}"
  find "test/results/structures/${COND}" -type f \( -name "*.cif.gz" -o -name "*.cif" -o -name "*.pdb" \) | sort | while read -r STRUCT; do
    mpnn \
      --model_type ligand_mpnn \
      --checkpoint_path "${MPNN_CKPT}" \
      --is_legacy_weights True \
      --structure_path "$STRUCT" \
      --out_directory "test/results/mpnn/${COND}" \
      --fixed_chains "A" \
      --batch_size 4 \
      --number_of_batches 1 \
      --temperature 0.1 \
      --write_fasta True \
      --write_structures True \
      --seed 42
  done
done

export CCD_MIRROR_PATH="/tmp"
export PDB_MIRROR_PATH="/tmp"

# Stage 3: Independent Fold-Back Validation with RoseTTAFold3 (RF3) on CPU
for COND in "${CONDITIONS[@]}"; do
  echo "=== [Stage 3] Running RF3 Fold-Back Validation for ${COND} ==="
  mkdir -p "test/results/refold/${COND}" "test/results/mpnn_fasta/${COND}"

  # Copy all FASTA files (.fa or .fasta) recursively from Stage 2 output
  find "test/results/mpnn/${COND}" -type f \( -name "*.fa" -o -name "*.fasta" \) -exec cp {} "test/results/mpnn_fasta/${COND}/" \;

  FASTA_COUNT=$(find "test/results/mpnn_fasta/${COND}" -type f | wc -l | tr -d ' ')
  echo "[INFO] Found ${FASTA_COUNT} FASTA sequence files for ${COND}"

  python -c "
import concurrent.futures, functools, sys, torch

# 1. Disable MPS cleanly while preserving __wrapped__ for torch._dynamo
orig_mps = torch.backends.mps.is_available
@functools.wraps(orig_mps)
def no_mps():
    return False
torch.backends.mps.is_available = no_mps

# 2. Prevent ProcessPoolExecutor(max_workers=0) crash on macOS
OrigPool = concurrent.futures.ProcessPoolExecutor
class SafeProcessPoolExecutor(OrigPool):
    def __init__(self, max_workers=None, *args, **kwargs):
        if max_workers is not None and max_workers < 1:
            max_workers = 1
        super().__init__(max_workers=max_workers, *args, **kwargs)
concurrent.futures.ProcessPoolExecutor = SafeProcessPoolExecutor

from rf3.cli import app
sys.exit(app())
" fold \
    inputs="test/results/mpnn_fasta/${COND}" \
    out_dir="test/results/refold/${COND}" \
    ckpt_path="${RF3_CKPT}"
done

# Stage 4: Quantitative Metric Aggregation & Plotting
echo "=== [Stage 4] Evaluating Metrics ==="
python test/analyze_metrics.py

# and others python codes