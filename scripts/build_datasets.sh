#!/usr/bin/env bash
# Rebuild the OrthoPan test datasets under scratch/datasets/.
#
# Usage:
#   scripts/build_datasets.sh                 # all four datasets
#   scripts/build_datasets.sh demo hog        # a subset
#
# Datasets:
#   demo    OrthoFinder demo results (24 genomes) + synthetic genome metadata
#   hog     tests/data Format A + real descriptions + demo tree + --fasta
#   large   resources/large-dataset (2,193 genomes) + UPGMA tree + synthetic files
#   sim10k  scripts/simulate_pangenome.py at 10,000 x 50,000, then import
#
# Idempotent: each target (and its generated input directory) is removed first.
# Inputs live under resources/ (read-only, gitignored) or tests/data, which is
# downloaded here when missing. Run from the repo root or any subdirectory.
set -euo pipefail

cd "$(dirname "$0")/.."
OUT=scratch/datasets

need() {
    if [[ ! -e "$1" ]]; then
        echo "error: missing input '$1'" >&2
        exit 1
    fi
}

time_build() {
    local name=$1 start=$SECONDS
    shift
    "$@"
    echo "built $name in $((SECONDS - start))s"
}

# tests/data is gitignored and no test requests it, so download it here when
# missing. Mirrors tests/conftest.py REMOTE_FILES.
TESTS_DATA_FILES=(
    "orthologs.tsv|https://demo.abrinca.com/files_html/orthologs/orthologs.tsv?download"
    "orthologs-descriptions.tsv|https://demo.abrinca.com/files_html/orthologs/HOG_best_name.tsv?download"
)

ensure_tests_data() {
    mkdir -p tests/data
    local entry name url dest tmp
    for entry in "${TESTS_DATA_FILES[@]}"; do
        name=${entry%%|*}
        url=${entry#*|}
        dest=tests/data/$name
        [[ -s "$dest" ]] && continue
        echo "downloading tests/data/$name ..." >&2
        tmp=$dest.part
        if ! curl -fsSL -A OrthoPan-Build "$url" -o "$tmp"; then
            rm -f "$tmp"
            echo "error: failed to download $url" >&2
            exit 1
        fi
        mv "$tmp" "$dest"
    done
}

build_demo() {
    local results=resources/orthofinder-demo/results/Results_Oct04_1
    local input=$OUT/demo-input target=$OUT/demo
    need "$results/Orthogroups/Orthogroups.tsv"
    need "$results/Species_Tree/SpeciesTree_rooted.txt"
    rm -rf "$target" "$input"
    uv run python scripts/make_test_metadata.py \
        --orthologs "$results/Orthogroups/Orthogroups.tsv" \
        --tree "$results/Species_Tree/SpeciesTree_rooted.txt" \
        --out "$input"
    uv run orthopan import --orthofinder "$results" \
        --genome-metadata "$input/genome-metadata.tsv" -o "$target"
}

build_hog() {
    local results=resources/orthofinder-demo/results/Results_Oct04_1
    local target=$OUT/hog
    ensure_tests_data
    need "$results/Species_Tree/SpeciesTree_rooted.txt"
    need resources/orthofinder-demo/faa-core
    need resources/orthofinder-demo/faa-add
    rm -rf "$target"
    uv run orthopan import tests/data/orthologs.tsv \
        --tree "$results/Species_Tree/SpeciesTree_rooted.txt" \
        --annotation tests/data/orthologs-descriptions.tsv \
        --fasta resources/orthofinder-demo/faa-core resources/orthofinder-demo/faa-add \
        -o "$target"
}

build_large() {
    local orthologs=resources/large-dataset/orthologs.tsv.gz
    local input=$OUT/large-input target=$OUT/large
    need "$orthologs"
    rm -rf "$target" "$input"
    mkdir -p "$input"
    # Build the tree first (it needs the whole matrix), then derive the
    # synthetic metadata/annotations from it and import with --tree.
    uv run orthopan-tree "$orthologs" -o "$input/tree.nwk"
    uv run python scripts/make_test_metadata.py \
        --orthologs "$orthologs" --tree "$input/tree.nwk" \
        --out "$input" --annotations
    uv run orthopan import "$orthologs" --tree "$input/tree.nwk" \
        --annotation "$input/descriptions.tsv" \
        --annotation "$input/cog.tsv" \
        --genome-metadata "$input/genome-metadata.tsv" \
        -o "$target"
}

build_sim10k() {
    local input=$OUT/sim10k-input target=$OUT/sim10k
    rm -rf "$target" "$input"
    mkdir -p "$input"
    uv run python scripts/simulate_pangenome.py \
        --genomes 10000 --orthogroups 50000 --out "$input"
    uv run orthopan import "$input/orthologs.tsv.gz" --tree "$input/tree.nwk" \
        --annotation "$input/annotations.descriptions.tsv" \
        --annotation "$input/annotations.cog.tsv" \
        --genome-metadata "$input/genomes.meta.tsv" -o "$target"
}

targets=("$@")
if [[ ${#targets[@]} -eq 0 ]]; then
    targets=(demo hog large sim10k)
fi
# Validate every requested name before starting any build.
for name in "${targets[@]}"; do
    case "$name" in
        demo|hog|large|sim10k) ;;
        *) echo "error: unknown dataset '$name' (choose: demo hog large sim10k)" >&2
           exit 1 ;;
    esac
done
for name in "${targets[@]}"; do
    case "$name" in
        demo)   time_build demo build_demo ;;
        hog)    time_build hog build_hog ;;
        large)  time_build large build_large ;;
        sim10k) time_build sim10k build_sim10k ;;
    esac
done
