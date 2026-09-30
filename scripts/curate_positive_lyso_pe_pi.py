"""Curate bundled positive LPE-O/LPE-P and PI ammonium MSP entries.

Back up the previous archive in the caller's test directory. The operation is
idempotent: an archive already containing LPE-P entries is left unchanged.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import re
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

from lipidgate.ms2.library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C,
    HYDROGEN_MONOISOTOPIC_MASS as H,
    NITROGEN_MONOISOTOPIC_MASS as N,
    OXYGEN_MONOISOTOPIC_MASS as O,
    PHOSPHORUS_MONOISOTOPIC_MASS as P,
)


PE_HEADGROUP_LOSS = 2 * C + 8 * H + N + 4 * O + P
PI_CHAIN_LOSS = re.compile(r'^\[M-\((?:ROOH|R=O)\)\+H\]\+')
ETHER_NAME = re.compile(r'^PE\(O-(\d+):(\d+)\)$')


def blocks(handle):
    current = []
    for line in handle:
        if line.strip():
            current.append(line)
        elif current:
            yield current
            current = []
    if current:
        yield current


def header(lines, field):
    prefix = field + ':'
    return next((line.split(':', 1)[1].strip() for line in lines if line.startswith(prefix)), '')


def peak_lines(lines):
    return [line for line in lines if line[:1].isdigit()]


def peak_name(line):
    return line.split('"')[1] if '"' in line else ''


def rendered(lines, peaks, replacements=None):
    replacements = replacements or {}
    first_peak = next(i for i, line in enumerate(lines) if line[:1].isdigit())
    heading = lines[:first_peak]
    rewritten = []
    for line in heading:
        field = line.split(':', 1)[0]
        if field == 'Num Peaks':
            rewritten.append(f'Num Peaks: {len(peaks)}\n')
        elif field in replacements:
            rewritten.append(f'{field}: {replacements[field]}\n')
        else:
            rewritten.append(line)
    return rewritten + peaks


def peak(mz, name, role):
    return f'{float(mz):.4f} 100.00 "{name}" "{role}"\n'


def curate(source: Path, backup_dir: Path):
    source = Path(source)
    backup_dir = Path(backup_dir)
    backup_dir.mkdir(parents=True, exist_ok=True)
    backup = backup_dir / 'positive_before_lyso_pe_pi.msp.gz'
    with gzip.open(source, 'rt', encoding='utf-8') as check:
        if any(header(block, 'CompoundClass') == 'LPE-P' for block in blocks(check)):
            print('Already curated: LPE-P records are present')
            return
    if not backup.exists():
        shutil.copy2(source, backup)
    temporary = source.with_suffix('.tmp')
    changed_o = created_p = changed_pi = removed_pi = 0
    with gzip.open(source, 'rt', encoding='utf-8') as incoming, gzip.open(
        temporary, 'wt', encoding='utf-8', compresslevel=6
    ) as outgoing:
        for block in blocks(incoming):
            cls = header(block, 'CompoundClass')
            adduct = header(block, 'PrecursorType')
            raw = peak_lines(block)
            if cls == 'LPE-O' and adduct == '[M+H]+':
                by_name = {peak_name(line): float(line.split()[0]) for line in raw}
                mz = float(header(block, 'PrecursorMZ'))
                o_peaks = [
                    peak(mz - PE_HEADGROUP_LOSS, '[M+H-141]+', 'Diagnostic_HG'),
                    peak(by_name['[M-H2O+H]+'], '[M+H-H2O]+', 'Common'),
                    peak(mz, '[M+H]+', 'Common'),
                ]
                outgoing.writelines(rendered(block, o_peaks))
                outgoing.write('\n')
                changed_o += 1
                matched = ETHER_NAME.fullmatch(header(block, 'Name'))
                if matched and int(matched.group(2)) >= 1:
                    name = f'LPE(P-{matched.group(1)}:{int(matched.group(2)) - 1})'
                    p_peaks = [
                        peak(by_name['[M-C3H11NO5P+H]+'], 'LPE-P diagnostic 294', 'Diagnostic_HG'),
                        peak(by_name['[M-C3H9NO4P+H]+'], 'LPE-P diagnostic 312', 'Diagnostic_HG'),
                        peak(by_name['[M-C3H6O2+H]+'], '[P-chain+PEtn]+', 'Common'),
                        peak(by_name['[M-H2O+H]+'], '[M+H-H2O]+', 'Common'),
                        peak(mz, '[M+H]+', 'Common'),
                    ]
                    outgoing.writelines(rendered(block, p_peaks, {
                        'Name': name,
                        'CompoundClass': 'LPE-P',
                        'Comment': f'MS1_name={name};polarity=+',
                    }))
                    outgoing.write('\n')
                    created_p += 1
                continue
            if cls == 'PI' and adduct == '[M+NH4]+':
                kept = [line for line in raw if not PI_CHAIN_LOSS.match(peak_name(line))]
                if len(kept) != len(raw):
                    outgoing.writelines(rendered(block, kept))
                    outgoing.write('\n')
                    changed_pi += 1
                    removed_pi += len(raw) - len(kept)
                    continue
            outgoing.writelines(block)
            outgoing.write('\n')
    if changed_o != 138 or created_p == 0 or changed_pi != 3907:
        temporary.unlink(missing_ok=True)
        raise ValueError(f'Unexpected source counts: O={changed_o}, P={created_p}, PI={changed_pi}')
    temporary.replace(source)
    print(f'LPE-O revised: {changed_o}; LPE-P added: {created_p}; '
          f'PI ammonium revised: {changed_pi}; PI fragments removed: {removed_pi}; backup: {backup}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1] /
                        'libraries/ms2/current_positive.msp.gz')
    parser.add_argument('--backup-dir', type=Path, required=True)
    args = parser.parse_args()
    curate(args.source, args.backup_dir)
