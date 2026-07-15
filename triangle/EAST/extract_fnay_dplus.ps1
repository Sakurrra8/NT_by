param(
    [Parameter(Mandatory = $true)]
    [string[]] $CasePath,
    [switch] $Replace
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$utf8NoBom = [System.Text.UTF8Encoding]::new($false)

foreach ($case in $CasePath) {
    $exports = @(Get-ChildItem -LiteralPath $case -Directory -Filter 'nt_export_*')
    if ($exports.Count -ne 1) {
        throw "Expected one nt_export directory in $case, found $($exports.Count)"
    }

    $source = Join-Path $exports[0].FullName 'b2plot/fnay_all_species.dat'
    $lines = [System.IO.File]::ReadAllLines($source)
    $starts = @(for ($index = 0; $index -lt $lines.Length; ++$index) {
        if ($lines[$index].StartsWith('# Radial particle flux') -and
            $lines[$index].Contains('D^{+1}')) {
            $index
        }
    })
    if ($starts.Count -ne 1) {
        throw "Expected one D+ FNIY block in $source, found $($starts.Count)"
    }

    $start = $starts[0]
    $end = $lines.Length
    for ($index = $start + 1; $index -lt $lines.Length; ++$index) {
        if ($lines[$index].StartsWith('#')) {
            $end = $index
            break
        }
    }
    $sourceBlock = @($lines[$start..($end - 1)])
    if ($sourceBlock.Count -ne 1 + 98 * 38) {
        throw "D+ FNIY block in $source has $($sourceBlock.Count - 1) rows; expected 3724"
    }

    $block = [System.Collections.Generic.List[string]]::new()
    $block.Add($sourceBlock[0])
    for ($index = 1; $index -lt $sourceBlock.Count; ++$index) {
        $columns = @($sourceBlock[$index].Split(
            [char[]] " `t", [System.StringSplitOptions]::RemoveEmptyEntries))
        if ($columns.Count -ne 4 -or [int] $columns[2] -ne 1) {
            throw "Unexpected all-species FNIY row: $($sourceBlock[$index])"
        }
        $block.Add("$($columns[0]) $($columns[1]) 0 $($columns[3])")
    }

    $destination = Join-Path $case '2D_data/fnay_Dplus.dat'
    if (Test-Path -LiteralPath $destination) {
        $existing = [System.IO.File]::ReadAllLines($destination)
        $matches = $existing.Count -eq $block.Count
        for ($index = 1; $matches -and $index -lt $block.Count; ++$index) {
            $old = @($existing[$index].Split(
                [char[]] " `t", [System.StringSplitOptions]::RemoveEmptyEntries))
            $new = @($block[$index].Split(
                [char[]] " `t", [System.StringSplitOptions]::RemoveEmptyEntries))
            $matches = $old.Count -eq 4 -and $new.Count -eq 4
            for ($column = 0; $matches -and $column -lt 3; ++$column) {
                $matches = [int] $old[$column] -eq [int] $new[$column]
            }
            if ($matches) {
                $oldValue = [double]::Parse(
                    $old[3], [System.Globalization.CultureInfo]::InvariantCulture)
                $newValue = [double]::Parse(
                    $new[3], [System.Globalization.CultureInfo]::InvariantCulture)
                $scale = [Math]::Max(1.0, [Math]::Max([Math]::Abs($oldValue), [Math]::Abs($newValue)))
                $matches = [Math]::Abs($oldValue - $newValue) -le 1.0e-12 * $scale
            }
        }
        if (-not $matches -and -not $Replace) {
            throw "Existing $destination does not match the exported D+ block"
        }
        if (-not $matches) {
            [System.IO.File]::WriteAllLines($destination, $block, $utf8NoBom)
            Write-Host "replaced $destination"
            continue
        }
        Write-Host "verified $destination"
        continue
    }

    [System.IO.Directory]::CreateDirectory((Split-Path $destination)) | Out-Null
    [System.IO.File]::WriteAllLines($destination, $block, $utf8NoBom)
    Write-Host "wrote $destination"
}
