@echo off
setlocal DisableDelayedExpansion
rem Export only this revised normalized-betting publication package.
rem Double-click: create a NEW copy in this folder's results\exports directory.
rem Optional first argument: new folder name, or absolute path under results\exports.
rem Explicit allowlist only. Never deletes, overwrites, changes research files, or uploads.
set "COPA_EXPORT_BAT=%~f0"
set "COPA_EXPORT_SOURCE=%~dp0"
set "COPA_EXPORT_DEST=%~1"
powershell.exe -NoLogo -NoProfile -Command "$ErrorActionPreference='Stop'; try { $text=[System.IO.File]::ReadAllText($env:COPA_EXPORT_BAT); $parts=$text -split '(?m)^# POWERSHELL_PAYLOAD\r?$',2; if ($parts.Count -ne 2) { throw 'Embedded export script was not found.' }; & ([scriptblock]::Create($parts[1])) } catch { Write-Host ('EXPORT FAILED: ' + $_.Exception.Message) -ForegroundColor Red; exit 1 }"
set "COPA_EXPORT_EXIT=%ERRORLEVEL%"
if not defined COPA_EXPORT_NO_PAUSE pause
exit /b %COPA_EXPORT_EXIT%
# POWERSHELL_PAYLOAD
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version 2.0

function Get-PackageHash([string] $path) {
    $hasher = [System.Security.Cryptography.SHA256]::Create()
    try {
        $stream = [System.IO.File]::OpenRead($path)
        try {
            return [System.BitConverter]::ToString($hasher.ComputeHash($stream)).Replace('-','').ToLowerInvariant()
        } finally { $stream.Dispose() }
    } finally { $hasher.Dispose() }
}

try {
    $source = (Resolve-Path -LiteralPath $env:COPA_EXPORT_SOURCE).ProviderPath.TrimEnd('\','/')
    if (-not (Test-Path -LiteralPath $source -PathType Container)) { throw 'Source must be a folder.' }
    $exportRoot = [System.IO.Path]::GetFullPath((Join-Path $source 'results\exports')).TrimEnd('\','/')
    if ([string]::IsNullOrWhiteSpace($env:COPA_EXPORT_DEST)) {
        $suffix = (Get-Date -Format 'yyyyMMdd_HHmmss') + '_' + [guid]::NewGuid().ToString('N').Substring(0,6)
        $destination = Join-Path $exportRoot ('github_copa2026_upload_' + $suffix)
    } elseif ([System.IO.Path]::IsPathRooted($env:COPA_EXPORT_DEST)) {
        $destination = [System.IO.Path]::GetFullPath($env:COPA_EXPORT_DEST).TrimEnd('\','/')
    } else {
        $destination = [System.IO.Path]::GetFullPath((Join-Path $exportRoot $env:COPA_EXPORT_DEST)).TrimEnd('\','/')
    }
    # Restrict writes to a new direct child of this package's ignored export root.
    if (-not $destination.StartsWith($exportRoot + '\',[System.StringComparison]::OrdinalIgnoreCase) -or
        -not (Split-Path -Parent $destination).Equals($exportRoot,[System.StringComparison]::OrdinalIgnoreCase)) {
        throw 'Destination must be a NEW direct child of this publication folder''s results\exports directory.'
    }
    if (Test-Path -LiteralPath $destination) { throw "Destination exists; nothing changed. Choose a NEW name: $destination" }
    # Check existing ancestors before creating any directory; never follow a junction for writes.
    $existingParent = $exportRoot
    while (-not (Test-Path -LiteralPath $existingParent)) { $existingParent = Split-Path -Parent $existingParent }
    $ancestor = Get-Item -LiteralPath $existingParent -Force
    if (-not $ancestor.PSIsContainer) { throw 'The export parent must be a directory.' }
    while ($null -ne $ancestor) {
        if (($ancestor.Attributes -band [System.IO.FileAttributes]::ReparsePoint) -ne 0) {
            throw "Export parent uses a junction/symlink: $($ancestor.FullName)"
        }
        $ancestor = $ancestor.Parent
    }

    $matlabFiles = @('combination.m','knn_betting_prepare.m','knn_ellipsoid_pdf.m','pvalues.m',
        'run_knn_cm_experiment.m','summarize_results.m','run_benchmarks.m',
        'run_full_figure_sweeps.m','reconstruct_copa_figures.m',
        'prepare_airquality_input.m','run_real.m')
    $dataFiles = @('CASP.mat','MiniBooNEPID.mat','Sensorlessdrivediagnosis.mat','creditcard_X.mat',
        'README.md','SOURCES_AND_LICENSES.md','PREPROCESSING.md','CITATIONS.bib','SHA256SUMS.txt','manifest.json')
    $requiredFiles = $matlabFiles + @('README.md','EXPERIMENTS.md','LICENSE','LICENSE-NOTICE.md','DATA_ATTRIBUTION.md',
        'REPRODUCTION_SCOPE.md','PUBLICATION_VERIFICATION.md','.gitignore','PACKAGE_MANIFEST.json',
        'Prepare_GitHub_Package.bat','results\README.md')
    foreach ($name in $dataFiles) { $requiredFiles += 'datasets\' + $name }
    Write-Host "Checking this publication package: $source"
    $hashes = @{}
    foreach ($name in $requiredFiles) {
        $inputPath = Join-Path $source $name
        if (-not (Test-Path -LiteralPath $inputPath -PathType Leaf)) {
            throw "Required file missing: $name. No output folder was created."
        }
        $hashes[$name] = Get-PackageHash $inputPath
    }
    $manifest = Get-Content -LiteralPath (Join-Path $source 'PACKAGE_MANIFEST.json') -Raw | ConvertFrom-Json
    if ($manifest.DefaultBetting -ne 'normalized' -or $manifest.Epsilon -ne 100 -or
        $manifest.LookBackWindow -ne 5000 -or $manifest.AlarmThreshold -ne 100 -or
        $manifest.MatlabFiles -ne 11 -or $manifest.DatasetCount -ne 4) {
        throw 'Unexpected publication manifest settings.'
    }
    $manifestNames = @($manifest.Files | ForEach-Object { $_.File.Replace('/','\') } | Sort-Object)
    $expectedNames = @($requiredFiles | Where-Object { $_ -ne 'PACKAGE_MANIFEST.json' } | Sort-Object)
    if ($manifestNames.Count -ne $expectedNames.Count -or
        [string]::Join('|',$manifestNames) -ne [string]::Join('|',$expectedNames)) {
        throw 'Publication manifest and export allowlist differ.'
    }
    foreach ($entry in $manifest.Files) {
        if ($hashes[$entry.File.Replace('/','\')] -ne $entry.SHA256) {
            throw "Publication file changed since its manifest was recorded: $($entry.File). Review and update the manifest first."
        }
    }
    $dataManifest = Get-Content -LiteralPath (Join-Path $source 'datasets\manifest.json') -Raw | ConvertFrom-Json
    if (@($dataManifest.Files).Count -ne 4) { throw 'Dataset manifest must contain exactly four inputs.' }
    foreach ($spec in $dataManifest.Files) {
        if (-not $hashes.ContainsKey('datasets\' + $spec.File) -or
            $hashes['datasets\' + $spec.File] -ne $spec.SHA256) {
            throw "Dataset checksum/allowlist mismatch: $($spec.File)"
        }
    }

    # Preflight is complete. All mutations remain under this publication folder.
    if (-not (Test-Path -LiteralPath $exportRoot)) {
        [System.IO.Directory]::CreateDirectory($exportRoot) | Out-Null
    }
    New-Item -ItemType Directory -Path $destination -ErrorAction Stop | Out-Null
    New-Item -ItemType Directory -Path (Join-Path $destination 'datasets') -ErrorAction Stop | Out-Null
    New-Item -ItemType Directory -Path (Join-Path $destination 'results') -ErrorAction Stop | Out-Null
    foreach ($name in $requiredFiles) {
        $outputPath = Join-Path $destination $name
        [System.IO.File]::Copy((Join-Path $source $name),$outputPath,$false)
        if ((Get-PackageHash $outputPath) -ne $hashes[$name]) { throw "Copy verification failed: $name" }
    }
    Write-Host ''
    Write-Host 'PUBLICATION EXPORT READY' -ForegroundColor Green
    Write-Host "Folder to publish: $destination"
    Write-Host "Files: $($requiredFiles.Count), including 11 MATLAB files and four benchmark datasets."
    Write-Host 'Includes proposed-method benchmarks, Figures 1/2 and the local-input Figure 3 workflow.'
    Write-Host 'Read REPRODUCTION_SCOPE.md for missing baseline pipelines and historical differences.'
    Write-Host 'Existing results, logs, caches and Air Quality input were not copied. Nothing was uploaded.'
    exit 0
} catch {
    Write-Host ('EXPORT FAILED: ' + $_.Exception.Message) -ForegroundColor Red
    Write-Host 'Source publication files were not changed. A partial NEW export, if any, is not ready to publish.'
    exit 1
}
