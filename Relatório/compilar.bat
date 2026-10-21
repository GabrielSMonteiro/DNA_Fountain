@echo off
title Compilando TCC - Synthetic DNA Data Storage
echo ======================================================
echo   Compilando TCC: Synthetic DNA Data Storage (UnB)
echo ======================================================
cd /d "%~dp0"

echo [1/4] Executando pdflatex (1a passada)...
pdflatex -interaction=nonstopmode -synctex=1 main.tex > nul

echo [2/4] Executando bibtex...
bibtex main > nul

echo [3/4] Executando pdflatex (2a passada)...
pdflatex -interaction=nonstopmode -synctex=1 main.tex > nul

echo [4/4] Executando pdflatex (3a passada final)...
pdflatex -interaction=nonstopmode -synctex=1 main.tex > nul

if exist main.pdf (
    echo.
    echo ======================================================
    echo   SUCESSO! PDF compilado com exito: main.pdf
    echo ======================================================
) else (
    echo.
    echo ERRO na compilacao. Verifique o arquivo main.log.
)

pause
