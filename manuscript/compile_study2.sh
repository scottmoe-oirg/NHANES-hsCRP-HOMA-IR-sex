#!/usr/bin/env bash
set -euo pipefail

# Build from the directory containing this script so that the command works
# whether invoked from the repository root or from manuscript/.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MAIN="NHANES_sex_interaction_manuscript"
SUPP="NHANES_sex_interaction_supplement"

echo "Building main manuscript..."
pdflatex -interaction=nonstopmode -halt-on-error "$MAIN.tex"
bibtex "$MAIN"
pdflatex -interaction=nonstopmode -halt-on-error "$MAIN.tex"
pdflatex -interaction=nonstopmode -halt-on-error "$MAIN.tex"

echo "Building supplementary material..."
pdflatex -interaction=nonstopmode -halt-on-error "$SUPP.tex"
pdflatex -interaction=nonstopmode -halt-on-error "$SUPP.tex"

echo
echo "Build complete:"
echo "  manuscript/$MAIN.pdf"
echo "  manuscript/$SUPP.pdf"
