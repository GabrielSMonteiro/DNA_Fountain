# Script de compilacao automatica do TCC (UnB / ABNT)
# Executa pdflatex + bibtex + pdflatex + pdflatex para resolver referencias cruzadas, sumario e bibliografia.

$ErrorActionPreference = "Stop"

Write-Host "======================================================" -ForegroundColor Cyan
Write-Host "  Compilando TCC: Synthetic DNA Data Storage (UnB)   " -ForegroundColor Cyan
Write-Host "======================================================" -ForegroundColor Cyan

$CurrentDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $CurrentDir

# 1. Primeira passada do PDFLaTeX
Write-Host "[1/4] Executando pdflatex (1a passada)..." -ForegroundColor Yellow
pdflatex -interaction=nonstopmode -synctex=1 main.tex | Out-Null
if ($LASTEXITCODE -ne 0) {
    Write-Warning "Aviso na 1a passada do pdflatex. Verifique main.log se necessario."
}

# 2. BibTeX para processar a bibliografia
Write-Host "[2/4] Executando bibtex..." -ForegroundColor Yellow
bibtex main | Out-Null

# 3. MakeIndex se arquivo .idx existir com entradas
if (Test-Path "main.idx") {
    makeindex -q main.idx 2>$null
}

# 4. Segunda passada para atualizar referencias e citacoes
Write-Host "[3/4] Executando pdflatex (2a passada)..." -ForegroundColor Yellow
pdflatex -interaction=nonstopmode -synctex=1 main.tex | Out-Null

# 5. Terceira passada para estabilizar sumario e listas
Write-Host "[4/4] Executando pdflatex (3a passada final)..." -ForegroundColor Yellow
pdflatex -interaction=nonstopmode -synctex=1 main.tex | Out-Null

if (Test-Path "main.pdf") {
    $pdf = Get-Item "main.pdf"
    Write-Host "`nSUCESSO: PDF compilado com exito!" -ForegroundColor Green
    Write-Host "Arquivo: $($pdf.FullName) ($([math]::Round($pdf.Length / 1KB, 1)) KB)" -ForegroundColor Green
} else {
    Write-Error "Falha na geracao de main.pdf. Verifique main.log."
}
