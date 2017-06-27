import pybedtools
from tqdm import trange
import gffutils
from collections import defaultdict
from collections import OrderedDict
import copy

GENE_PRIORITY = [
    ['protein_coding','CDS'],
    ['protein_coding','start_codon'],
    ['protein_coding','stop_codon'],
    ['protein_coding','THREE_AND_FIVE_PRIME_UTR'],
    ['protein_coding','5UTR'],
    ['protein_coding','3UTR'],
    ['protein_coding','intron'],
    ['protein_coding','Selenocysteine'],
    ['non_coding','exon'],
    ['non_coding','intron'],
    ['non_coding','transcript'],
    ['non_coding','gene'],
    ['non_coding','THREE_AND_FIVE_PRIME_UTR'],
    ['non_coding','3UTR'],
    ['non_coding','5UTR'],
    ['non_coding','Selenocysteine'],
    ['non_coding', 'CDS'],  # shouldn't occur?
    ['non_coding','start_codon'],  # shouldn't occur?
    ['non_coding','stop_codon'],  # shouldn't occur?
    ['protein_coding', 'exon'],  # shouldn't occur?
    ['protein_coding', 'transcript'],  # shouldn't occur?
    ['protein_coding', 'gene'],  # shouldn't occur?

]

TRANSCRIPT_PRIORITY = [
    ['protein_coding','CDS'],
    ['protein_coding','start_codon'],
    ['protein_coding','stop_codon'],
    ['protein_coding','THREE_AND_FIVE_PRIME_UTR'],
    ['protein_coding','5UTR'],
    ['protein_coding','3UTR'],
    ['protein_coding','intron'],
    ['protein_coding','Selenocysteine'],
    ['non_coding','exon'],
    ['non_coding','intron'],
    ['non_coding','transcript'],
    ['non_coding','gene'],
    ['non_coding','THREE_AND_FIVE_PRIME_UTR'],
    ['non_coding','3UTR'],
    ['non_coding','5UTR'],
    ['non_coding','Selenocysteine'],
    ['non_coding','CDS'], # shouldn't occur?
    ['non_coding','start_codon'],  # shouldn't occur?
    ['non_coding','stop_codon'],  # shouldn't occur?
    ['protein_coding','exon'], # shouldn't occur?
    ['protein_coding','transcript'], # shouldn't occur?
    ['protein_coding','gene'], # shouldn't occur
]

HASH_VAL = 1000000
MAXVAL = 1000000000
MINVAL = 0

class Annotator():
    """

    class to get genomic features from gffutils _db

    """

    def __init__(self, db_file, chroms=[]):
        """
        
        :param db_file: 
        :param chroms: 
        """
        self.num_features = 0

        progress = trange(8, desc='Initializing/creating defs', leave=False)

        self._db = gffutils.FeatureDB(db_file)
        progress.update(1)
        self._featuretypes = self._db.featuretypes()
        progress.update(1)
        self._geneid_to_name_dict = self._gene_id_to_name()
        progress.update(1)
        if len(chroms) != 0:  # if use specific chromosomes, otherwise hash all
            self._chromosomes = set(chroms)
        else:
            self._chromosomes = self._chromosome_set()

        self._hash_features()

        progress.update(1)
        self.exons_dict = self._get_all_exons_dict()
        progress.update(1)
        self.transcripts_dict = self._get_all_transcripts_dict()
        progress.update(1)
        self._update_introns()
        progress.update(1)
        self.cds_dict = self._get_all_cds_dict()
        progress.update(1)

    def _chromosome_set(self):
        """
        Returns the set of chromosomes that exist in a database.
        
        :return: 
        """
        ret = self._db.execute("SELECT seqid FROM features").fetchall()
        all_chromosomes = [r['seqid'] for r in ret]
        return set(all_chromosomes)

    def _hash_features(self):
        """
        hashes features by position.
        :return features_dict : collections.defaultdict()
            dictionary of features{[chrom, pos/HASH_VAL, strand] : feature_list}
        """
        num_features = 0

        features_dict = defaultdict(list)
        progress = trange(
            len(self._chromosomes),
            leave=False,
            desc='Build Location Index'
        )
        for chrom in self._chromosomes:
            for element in self._db.region(seqid=chrom):
                start = int(element.start / HASH_VAL)
                end = int(element.end / HASH_VAL)
                for i in range(start, end+1):
                    features_dict[chrom, i, element.strand].append(element)
                num_features+=1
            progress.update(1)

        self.features_dict = features_dict
        self.num_features = num_features

    def _update_introns(self):
        # progress = trange(self.num_features, leave=False, desc='inferring introns from exon/genes')
        for hash_val, features in self.features_dict.iteritems():
            for feature in features:
                if feature.featuretype == 'transcript':
                    for transcript_id in feature.attributes['transcript_id']:
                        exons = self.exons_dict[transcript_id]
                        transcript = self.transcripts_dict[transcript_id]
                        introns = self.find_introns(transcript, exons)
                        for intron in introns:
                            intron_feature = copy.deepcopy(feature)
                            intron_feature.start = intron['start']
                            intron_feature.end = intron['end']
                            intron_feature.featuretype = 'intron'
                            self.features_dict[hash_val].append(
                                intron_feature
                            )
                # progress.update(1)

    def _get_all_cds_dict(self):
        """
        For every cds-annotated transcript id (ENST), return a 
        dictionary containing the lowest and highest
        cds start and end vals for that transcript.

        :return cds_dict : defaultdict{transcript:{'start':START, 'end':END}} 
        """
        cds_dict = defaultdict(lambda: {'low': MAXVAL, 'hi': MINVAL})
        for cds_feature in self._db.features_of_type('CDS'):
            for transcript_id in cds_feature.attributes['transcript_id']:

                if cds_feature.start <= cds_dict[transcript_id]['low']:
                    cds_dict[transcript_id]['low'] = cds_feature.start
                if cds_feature.end >= cds_dict[transcript_id]['hi']:
                    cds_dict[transcript_id]['hi'] = cds_feature.end
        return cds_dict

    def _get_all_exons_dict(self):
        """
        :return:
        """
        exons_dict = defaultdict(list)
        for exon_feature in self._db.features_of_type('exon'):
            for transcript_id in exon_feature.attributes['transcript_id']:
                exons_dict[transcript_id].append(
                    {
                        'start': exon_feature.start,
                        'end': exon_feature.end
                    }
                )
        return exons_dict

    def _get_all_transcripts_dict(self):
        """

        :return:
        """
        transcripts_dict = defaultdict(list)
        for transcript_feature in self._db.features_of_type('transcript'):
            for transcript_id in transcript_feature.attributes['transcript_id']:
                transcripts_dict[transcript_id] = {
                    'start': transcript_feature.start,
                    'end': transcript_feature.end
                }
        return transcripts_dict

    def _classify_utr(self, utr_feature):
        """
        Given a list of features, return a dictionary of 
        gene_id : {start: cds_start, end: cds_end} sites

        :param features : list[gffutils.Feature]
        :return cds_ends : dict
        """
        three_prime_utr = False
        five_prime_utr = False

        for transcript_id in utr_feature.attributes['transcript_id']:
            if utr_feature.strand == '+':
                if self.cds_dict[transcript_id]['low'] > utr_feature.end:
                    five_prime_utr = True
                if self.cds_dict[transcript_id]['hi'] < utr_feature.start:
                    three_prime_utr = True
            elif utr_feature.strand == '-':
                if self.cds_dict[transcript_id]['low'] > utr_feature.end:
                    three_prime_utr = True
                if self.cds_dict[transcript_id]['hi'] < utr_feature.start:
                    five_prime_utr = True

        if five_prime_utr and three_prime_utr:
            return 'THREE_AND_FIVE_PRIME_UTR'
        elif five_prime_utr:
            return '5UTR'
            # return 'five_prime_utr'
        elif three_prime_utr:
            return '3UTR'
            # return 'three_prime_utr'
        else:
            return 'UNCLASSIFIED_UTR'

    def _gene_id_to_name(self):
        """
        Returns a dictionary containing a gene_id:name translation
        Note: may be different if the 'gene_id' or 'gene_name' 
        keys are not in the source GTF file
        (taken from gscripts.region_helpers)
        
        :return gene_name_dict : dict
            dict of {gene_id : gene_name}
        """

        genes = self._db.features_of_type('gene')
        gene_name_dict = {}

        for gene in genes:
            gene_id = gene.attributes['gene_id'][0] if type(
                gene.attributes['gene_id']
            ) == list else gene.attributes['gene_id']
            try:
                gene_name_dict[gene_id] = gene.attributes['gene_name'][0]
            except KeyError:
                print(gene.attributes.keys())
                print("Warning. Key not found for {}".format(gene))
                return 1
        return gene_name_dict

    def find_introns(self, transcript, exons):
        positions = []
        introns = []
        for exon in exons:
            positions.append(exon['start'] - 1)
            positions.append(exon['end'] + 1)
        positions = sorted(positions)

        if positions[0] < transcript[
            'start']:  # there is no intron at the start of the feature
            positions.pop(0)
        else:
            positions.insert(transcript['start'])
        if positions[-1] > transcript['end']:
            positions.pop(-1)
        else:
            positions.append(transcript['end'])
        for i in range(0, len(positions) - 1, 2):
            introns.append({'start': positions[i], 'end': positions[i + 1]})
        return introns

    def get_all_overlapping_features_from_query(self, chrom, qstart, qend,
                                                strand):
        """
        Given a query location (chr, start, end), return all features that
        overlap by at least one base. Functions similarly to gffutils db.region(),
        but uses the pre-hashed self.features_dict to greatly speed things up.

        :param chrom : string
        :param qstart : int
        :param qend : int
        :param strand : string
        :return features : list
            list of gffutils.Feature objects.
        """
        features = []
        start_key = int(qstart / HASH_VAL)
        end_key = int(qend / HASH_VAL)
        qstart = qstart + 1 # change 0-based bed to 1-based gff

        for i in range(start_key, end_key + 1):
            for feature in self.features_dict[chrom, i, strand]:
                if qstart <= feature.start and qend >= feature.end:  # feature completely contains query
                    features.append(feature)
                elif qstart >= feature.start and qend <= feature.end:  # query completely contains feature
                    features.append(feature)
                elif qstart <= feature.start and qend >= feature.start:  # feature partially overlaps (qstart < fstart < qend)
                    features.append(feature)
                elif qstart <= feature.end and qend >= feature.end:  # feature partially overlaps (qstart < fend < qend)
                    features.append(feature)

        return features

    def annotate(
            self, interval,
            transcript_priority=TRANSCRIPT_PRIORITY,
            gene_priority=GENE_PRIORITY
    ):
        """
        Given an interval, annotates using the priority

        :param interval:
        :return:
        """
        overlapping_features = self.get_all_overlapping_features_from_query(
            interval.chrom,
            interval.start,
            interval.end,
            interval.strand
        )
        if len(overlapping_features) == 0:
            return 'INTERGENIC', 'INTERGENIC'
        to_append = ''  # full list of genes overlapping features
        transcript = defaultdict(list)
        for feature in overlapping_features:  # for each overlapping feature
            for transcript_id in feature.attributes[
                'transcript_id'
            ]:  # multiple genes can be associated with one feature
                transcript[transcript_id].append(
                    feature)  # append features to their respective genes
        for transcript, features in transcript.iteritems():
            for feature in features:
                # if 'protein_coding' not in feature.attributes['transcript_type']:
                #     if feature.featuretype == 'exon' or feature.featuretype == 'UTR':
                #         feature.featuretype = 'noncoding_exon'
                if feature.featuretype == 'UTR':
                    feature.featuretype = self._classify_utr(feature)
                to_append += "{}:{}:{}:{}:{}:".format(
                    transcript,
                    feature.start,
                    feature.end,
                    feature.strand,
                    feature.featuretype,
                )
                for t in feature.attributes['gene_id']:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
                for t in feature.attributes['gene_name']:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
                for t in feature.attributes['transcript_type']:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + '|'
        to_append = to_append[:-1]
        priority = self.prioritize_transcript_then_gene(
            self.parse_annotation_string(to_append),
            transcript_priority,
            gene_priority
        )
        # print(feature.start, feature.end, feature.strand, to_append)
        return priority, to_append

    def parse_annotation_string(self, features_string):
        """
        Splits a feature string into a list of feature strings

        :param features_string:
        :return:
        """
        features = features_string.split('|')
        return features

    def is_protein_coding(self, transcript_type):
        """
        if defined protein coding, return True else False

        :param transcript_type:
        :return:
        """
        if transcript_type == 'protein_coding':
            return True
        return False

    def return_highest_priority_feature(self, formatted_features, priority):
        # Build dict
        combined_dict = defaultdict(list)
        for feature_string in formatted_features:
            transcript, start, end, strand, feature_type, gene_id, gene_name, transcript_type_list = feature_string.split(
                ':')
            transcript_type_list = transcript_type_list.split(',')
            for transcript_type in transcript_type_list:
                if self.is_protein_coding(
                        transcript_type):  # simplify all the types at first
                    combined_dict['protein_coding', feature_type].append(
                        feature_string)
                else:
                    combined_dict['non_coding',  feature_type].append(
                        feature_string)
        # return the highest one
        combined_dict = OrderedDict(
            combined_dict)  # turn into ordered dict, is that ok?
        combined_dict = sorted(  # sort based on priority list
            combined_dict.iteritems(),
            key=lambda x: priority.index([x[0][0], x[0][1]])
        )

        return combined_dict[0]

    def prioritize_transcript_then_gene(self, formatted_features,
                                        transcript_priority, gene_priority):


        unique_transcript_features = defaultdict(list)
        unique_transcripts = defaultdict(list)
        unique_genes = defaultdict(list)
        final = []
        for feature_string in formatted_features:
            if feature_string.split(':')[4] != 'gene':
                transcript = feature_string.split(':')[0]
                unique_transcript_features[transcript].append(
                    feature_string)

        if len(unique_transcript_features.keys()) == 0:
            return 'INTERGENIC'

        ### PRIORITIZE TRANSCRIPT ###
        for transcript in unique_transcript_features.keys(): # For each unique transcript
            top_transcript = self.return_highest_priority_feature(
                unique_transcript_features[transcript],
                transcript_priority
            )[1][0]  # [0] contains the dictionary key
            unique_transcripts[transcript].append(
                top_transcript
            )
            # add gene key
            gene_list = top_transcript.split(':')[5].split(',')
            for gene in gene_list:
                unique_genes[gene].append(top_transcript)

        ### PRIORITIZE GENE ###
        for gene, transcripts in unique_genes.iteritems():
            for transcript in transcripts:
                final.append(transcript)
        feature_type, final = self.return_highest_priority_feature(
            final, gene_priority
        )

        if feature_type[0] == 'non_coding':  # TODO: fix. kind of hacky
            return final[0].replace('exon', 'noncoding_exon').replace('intron', 'noncoding_intron')
        return final[0]


def annotate(db_file, bed_file, out_file, chroms):
    """
    Given a bed6 file, return the file with an extra column containing
    '|' delimited gene annotations

    :param db_file:
    :param bed_file:
    :param out_file:
    :param chroms:
    :return:
    """
    annotator = Annotator(db_file, chroms)
    bed_tool = pybedtools.BedTool(bed_file)
    with open(out_file, 'w') as o:
        for interval in bed_tool:  # for each line in bed file
            priority, annotation = annotator.annotate(interval)
            o.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                interval.chrom, interval.start,
                interval.end, interval.name,
                interval.score, interval.strand,
                priority,
                annotation
            ))