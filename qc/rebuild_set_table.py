import sys
import re
import argparse
from collections import OrderedDict


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-t', dest='table', help='Table to transform',
                        required=True)
    parser.add_argument('-o', dest='outtable', help='Output file name',
                        required=True)
    parser.add_argument('-i', dest='identifier', help='Identifier colnr in '
                        'table, e.g. of peptide, protein ID. First column is '
                        '1', required=True, nargs='+', type=int)
    parser.add_argument('-s', dest='setnames', help='Setnames as in table.'
                        'Header fields should be $setname$_$fieldname$.',
                        nargs='+', required=True)
    parser.add_argument('-f', dest='fieldpats', nargs='+', help='Field name '
                        'patterns to extract from table. Header fields should '
                        'be $setname$_$fieldname$. Patterns are python regex',
                        required=True)
    args = parser.parse_args(sys.argv[1:])
    with open(args.outtable, 'w') as outfp:
        newtablelines = parse_input(args.table, args.setnames, args.fieldpats,
                                    [x - 1 for x in args.identifier])
        outfp.write('\t'.join(next(newtablelines)))
        for line in newtablelines:
            outfp.write('\n{}'.format('\t'.join(line)))


def parse_input(table, setnames, outfieldpats, idcols):
    with open(table) as fp:
        header = next(fp).strip('\n').split('\t')
        setindices = get_setindices(header, setnames)
        outfields = get_outfields(header, outfieldpats, setnames)
        newheader = [header[x] for x in idcols] + ['Set'] + outfields
        newheader = [x.replace('MS1 area (highest of all PSMs)',
                               'MS1 precursor area') for x in newheader]
        newheader = [x.replace('Protein accession', 'Accession') for x in newheader]
        newheader = [x.replace('Peptide sequence', 'Accession') for x in newheader]
        yield [x.replace('#', 'nr_') for x in newheader]
        for line in fp:
            line = line.strip('\n').split('\t')
            baseline = [line[x] for x in idcols]
            for setname in setnames:
                # FIXME this does not address the other fields than set fields
                # e.g coverage
                yield baseline + [setname] + [line[setindices[f][setname]]
                                              for f in outfields]


def get_outfields(header, patterns, setnames):
    outfields = OrderedDict()
    for pat in patterns:
        for field in header:
            for setname in setnames:
                if field.startswith('{}_'.format(setname)):
                    fieldname = field[len(setname) + 1:]
                    if re.match(pat, fieldname):
                        outfields[fieldname] = 1
                    break
    return [k for k in outfields]


def get_setindices(header, setnames):
    """From header like ---ID, coverage, set1_q-value set2_q-value---
    this returns indices for different sets {'q-value': {'set1': 2, 'set2: 3}}
    """
    setindices = OrderedDict()
    for index, field in enumerate(header):
        for setname in setnames:
            if field.startswith('{}_'.format(setname)):
                fieldname = field[len(setname) + 1:]
                try:
                    setindices[fieldname][setname] = index
                except KeyError:
                    setindices[fieldname] = {setname: index}
    return setindices


if __name__ == '__main__':
    main()
