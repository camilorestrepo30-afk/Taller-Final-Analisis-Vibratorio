#!/usr/bin/env bash
# Genera el PDF del informe a partir de las fuentes en Markdown.
# Requiere: pandoc y una distribución de LaTeX con xelatex (TeX Live o MiKTeX).
set -e
cd "$(dirname "$0")"
cat 00_inicio.md 01_nf1.md 02_nf3.md 03_nf4.md 04_nf5.md 05_nf9.md \
    06_f1.md 07_f2.md 08_f3.md 09_f7.md 10_f10.md 11_final.md > informe.md
pandoc informe.md -o informe.tex --standalone
xelatex -interaction=nonstopmode informe.tex > /dev/null
xelatex -interaction=nonstopmode informe.tex > /dev/null
cp informe.pdf ../Informe_Taller_Final_Vibraciones_corregido.pdf
echo "Listo: informe/Informe_Taller_Final_Vibraciones_corregido.pdf"
