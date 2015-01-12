from collections import defaultdict, namedtuple

SAMRecordLite = namedtuple('SAMRecord', 'queryid flags refid startpos mapqual cigar pairid pairstart pairfraglen querylen AS XS YS XN XM XO XG NM YF YT MD')

def gather_by_query(sam_fp, aligned=True):
    """
    Given a path to a Bowtie2 SAM file, collect the records into a
    dictionary keyed on the Query ID.
    
    :type aligned: bool
    :param aligned: If True, gather only aligned records from the SAM file
    """
    recs = defaultdict(list)
    opt_fields = ['AS', 'XS', 'YS', 'XN', 'XM', 'XO', 
                  'XG', 'NM', 'YF', 'YT', 'MD']

    with open(sam_fp, 'rU') as sf:
        split_opt = lambda split_entry: (split_entry[0], split_entry[-1])
        for line in sf:
            line = line.split('\t')
            file_opt = dict([split_opt(o.split(':')) for o in line[11:]])
            full_opt = [file_opt.get(opt, '') for opt in opt_fields]
            line = list(line[:9]) + [len(line[9])]
            line.extend(full_opt)
            rec = SAMRecordLite(*line)
            recs[rec.queryid].append(rec)

    return recs