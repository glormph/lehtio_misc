#!/usr/bin/env python

import re
import sys


def main():
    psms = sys.argv[1]
    outfile = sys.argv[2]
    pepcol = int(sys.argv[3]) - 1
    amount_spectra_files = sys.argv[4]
    platepatterns = sys.argv[5:]
    with open(outfile, 'w') as outfp:
        newpsms = parse_psms(psms, pepcol, platepatterns)
        outfp.write('\t'.join(next(newpsms)))
        for outline in newpsms:
            outfp.write('\n{}'.format('\t'.join(outline)))
    with open('amount_spectra_plates', 'w') as wfp, open(amount_spectra_files) as fp:
        wfp.write('Plate_ID\tamount_ms2\n')
        for line in fp:
            line = line.strip('\n').split('|')
            wfp.write('{}\t{}\n'.format(get_plate_id(line[0], line[1],
                                                     platepatterns), line[2]))


def parse_header(oldheader):
    header = ['Biological set', 'Retention time(min)', 'PrecursorError(ppm)',
              'Peptide', 'percolator svm-score', 'MS1 area', 'Fractions',
              'Delta pI',
              'tmt10plex_126',
              'tmt10plex_127N',
              'tmt10plex_127C',
              'tmt10plex_128N',
              'tmt10plex_128C',
              'tmt10plex_129N',
              'tmt10plex_129C',
              'tmt10plex_130N',
              'tmt10plex_130C',
              'tmt10plex_131',
              'tmt10plex_131N',
              'tmt10plex_131C',
              'itraq8plex_113',
              'itraq8plex_114',
              'itraq8plex_115',
              'itraq8plex_116',
              'itraq8plex_117',
              'itraq8plex_118',
              'itraq8plex_119',
              'itraq8plex_120',
              'itraq8plex_121',
              ]
    newheader, colnrs = [], []
    for field in header:
        try:
            colnrs.append(oldheader.index(field))
        except ValueError:
            # happens when not running hirief, TMT, e.g.
            pass
        else:
            newheader.append(field)
    return newheader + ['Plate_ID', 'missed_cleavage'], colnrs


def parse_psms(infile, pepseqcol, platepatterns):
    with open(infile) as fp:
        oldheader = next(fp).strip('\n').split('\t')
        newheader, colnrs = parse_header(oldheader)
        yield newheader
        biosetcol = oldheader.index('Biological set')
        fncol = oldheader.index('#SpecFile')
        for line in fp:
            line = line.strip('\n').split('\t')
            plate_id = get_plate_id(line[biosetcol], line[fncol], platepatterns)
            yield [line[x] for x in colnrs] + [
                plate_id, str(count_missed_cleavage(line[pepseqcol]))]


def get_plate_id(bioset, fn, patterns):
    for pattern in patterns:
        if pattern in fn:
            return '{}_{}'.format(bioset, pattern)
    raise RuntimeError('Could not match patterns {} to filename {} to detect '
                       'name of plate'.format(patterns, fn))


def count_missed_cleavage(full_pepseq, count=0):
    '''Regex .*[KR][^P] matches until the end and checks if there is a final
    charachter so this will not match the tryptic residue'''
    pepseq = re.sub('[\+\-]\d*.\d*', '', full_pepseq)
    match = re.match('.*[KR][^P]', pepseq)
    if match:
        count += 1
        return count_missed_cleavage(match.group()[:-1], count)
    else:
        return count


if __name__ == '__main__':
    main()
