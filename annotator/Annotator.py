import pybedtools
from tqdm import trange
import gffutils
from collections import defaultdict
from collections import OrderedDict
import copy
import sys

HASH_VAL = 1000000
MAXVAL = 1000000000
MINVAL = 0

class Annotator():
    """

    class to get genomic features from gffutils _db

    """

    def __init__(self, db_file, chroms=[], species='hg19', append_chr=False, fuzzy=0):
        """
        
        :param db_file: 
        :param chroms: 
        """

        if species == 'ce11':
            self.cds_key = 'CDS'
            self.utr3_key = 'three_prime_UTR'
            self.utr5_key = 'five_prime_UTR'
            self.utr_key = None
            self.gene_name_key = 'gene_id'
            self.trancript_id_key = 'transcript_id'
            self.type_key = 'transcript_biotype'
        else:
            self.cds_key = 'CDS'
            self.utr3_key = None #
            self.utr5_key = None # in human/mice, this key doesn't really exist
            self.utr_key = 'UTR'
            self.gene_name_key = 'gene_name'
            self.trancript_id_key = 'transcript_id'
            self.type_key = 'transcript_type'

        self.num_features = 0
        self._db = gffutils.FeatureDB(db_file)

        progress = trange(7, desc='Initializing/creating defs', leave=False)
        self._featuretypes = [x.lower() for x in list(self._db.featuretypes())]

        progress.update(1)

        print("featuretypes = {}".format(self._featuretypes))

        progress.set_description("Building gene ids -> gene names dictionary")
        self._geneid_to_name_dict = self._gene_id_to_name()
        progress.update(1)

        if len(chroms) != 0:  # if use specific chromosomes, otherwise hash all
            progress.set_description("Adding all features in chroms {} to dictionary".format(chroms))
            self._chromosomes = set(chroms)
        else:
            progress.set_description("Adding all features in DB to dictionary")
            self._chromosomes = self._chromosome_set()

        ### Store all features into dictionary
        self._hash_features(append_chr, fuzzy)
        progress.update(1)

        ### If 'exon' features exist in db, catalog it.
        if 'exon' in self._featuretypes:
            progress.set_description("Adding all exon boundaries")
            self.exons_dict = self._get_all_exons_dict()
        progress.update(1)
        ### If 'exon' features exist in db, catalog it.
        if 'transcript' in self._featuretypes:
            progress.set_description("Adding all transcript boundaries")
            self.transcripts_dict = self._get_all_transcripts_dict()
        progress.update(1)
        ### If 'exon' features exist in db, catalog it.
        if 'intron' not in self._featuretypes:
            progress.set_description("Inferring introns using exon/transcripts")
            self._update_introns()
        progress.update(1)
        ### If 'exon' features exist in db, catalog it.
        if 'cds' in self._featuretypes:
            progress.set_description("Adding all CDS boundaries")
            self.cds_dict = self._get_all_cds_dict()
        progress.update(1)
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
            try:
                gene_id = gene.attributes['gene_id'][0] if type(
                    gene.attributes['gene_id']
                ) == list else gene.attributes['gene_id']
                gene_name_dict[gene_id] = gene.attributes[self.gene_name_key][0]
            except KeyError:
                print(gene.attributes.keys())
                print("Warning. Key not found for {}".format(gene))
                sys.exit(1)
        return gene_name_dict

    def _chromosome_set(self):
        """
        Returns the set of chromosomes that exist in a database.

        :return: 
        """
        ret = self._db.execute("SELECT seqid FROM features").fetchall()
        all_chromosomes = [r['seqid'] for r in ret]
        return set(all_chromosomes)

    def _hash_features(self, append_chr, fuzzy):
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
            desc='Building feature index'
        )
        for chrom in self._chromosomes:
            fixed_chrom = 'chr' + chrom if append_chr else chrom

            progress.set_description("Adding {}...".format(chrom))
            for element in self._db.region(seqid=chrom):

                if element.featuretype == 'gene':
                    if element.strand == '+':
                        element.start = element.start - fuzzy
                    elif element.strand == '-':
                        element.end = element.end + fuzzy
                if element.featuretype == 'gene':
                    if element.strand == '+':
                        element.end = element.end + fuzzy
                    elif element.strand == '-':
                        element.start = element.start - fuzzy

                start = int(element.start / HASH_VAL)
                end = int(element.end / HASH_VAL)
                for i in range(start, end+1):
                    element.chrom = fixed_chrom
                    features_dict[fixed_chrom, i, element.strand].append(element)
                num_features+=1
            progress.update(1)

        self.features_dict = features_dict
        self.num_features = num_features

    def _get_all_cds_dict(self):
        """
        For every cds-annotated transcript id (ENST), return a 
        dictionary containing the lowest and highest
        cds start and end vals for that transcript.

        :return cds_dict : defaultdict{transcript:{'start':START, 'end':END}} 
        """
        cds_dict = defaultdict(lambda: {'low': MAXVAL, 'hi': MINVAL})

        for cds_feature in self._db.features_of_type(self.cds_key):
            for transcript_id in cds_feature.attributes['transcript_id']:
                if cds_feature.start <= cds_dict[transcript_id]['low']:
                    cds_dict[transcript_id]['low'] = cds_feature.start
                if cds_feature.end >= cds_dict[transcript_id]['hi']:
                    cds_dict[transcript_id]['hi'] = cds_feature.end
        return cds_dict

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
                if self.cds_dict[transcript_id]['hi'] < utr_feature.start + 1:
                    three_prime_utr = True
            elif utr_feature.strand == '-':
                if self.cds_dict[transcript_id]['low'] > utr_feature.end:
                    three_prime_utr = True
                if self.cds_dict[transcript_id]['hi'] < utr_feature.start + 1:
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

    def _get_all_exons_dict(self):
        """
        Returns dictionary of exons as transcript_id:{
            [
                {'start':START, 'end':END},
                {'start':START, 'end':END},
                ...
            ]
        }.
        
        :return:
        """
        exons_dict = defaultdict(list)
        for exon_feature in self._db.features_of_type('exon'):
            for transcript_id in exon_feature.attributes[self.trancript_id_key]:
                exons_dict[transcript_id].append(
                    {
                        'start': exon_feature.start,
                        'end': exon_feature.end
                    }
                )
        return exons_dict

    def _get_all_transcripts_dict(self):
        """
        Returns dictionary of transcript_id:{'start':START, 'end':END}.
        
        :return transcripts_dict: defaultdict(dict)
            hash of transcripts and their start/end coordinates
        """
        transcripts_dict = defaultdict(dict)
        for transcript_feature in self._db.features_of_type('transcript'):
            for transcript_id in transcript_feature.attributes[self.trancript_id_key]:
                transcripts_dict[transcript_id] = {
                    'start': transcript_feature.start,
                    'end': transcript_feature.end
                }
        return transcripts_dict

    def _update_introns(self):
        """
        Iterates over each transcript, uses each transcript_dict
        and exon_dict to group transcript + exon coordinates, 
        and appends to features_dict a list of intron features
        with the same attributes as the original transcript, but
        with the featuretype set as 'intron'.
        
        :return: 
        """
        for hash_val, features in self.features_dict.iteritems():
            for feature in features:
                if feature.featuretype == 'transcript':
                    for transcript_id in feature.attributes[self.trancript_id_key]:
                        exons = self.exons_dict[transcript_id]
                        transcript = self.transcripts_dict[transcript_id]
                        introns = self._find_introns(transcript, exons)
                        for intron in introns:
                            intron_feature = copy.deepcopy(feature)
                            intron_feature.start = intron['start']
                            intron_feature.end = intron['end']
                            intron_feature.featuretype = 'intron'
                            self.features_dict[hash_val].append(
                                intron_feature
                            )

    def _find_introns(self, transcript, exons):
        """
        Given transcript coordinates and a list of exons,
        return a list of introns (inverse coordinates)
        
        :param transcript: dictionary
            {'start':int, 'end':int}
        :param exons: list[dict]
            [
                {'start':int,'end':int},
                {'start':int,'end':int},
                ...
            ]
        :return introns: list[dict]
            [
                {'start':int,'end':int},
                {'start':int,'end':int},
                ...
            ]
        """
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

    def get_all_overlapping_features_from_query_stranded(self, chrom, qstart, qend,
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
                if qstart <= feature.start and qend >= feature.end:  # query completely contains feature
                    feature.attributes['overlap'] = 'query_contains_feature'
                    features.append(feature)
                elif qstart >= feature.start and qend <= feature.end:  # feature completely contains query
                    feature.attributes['overlap'] = 'feature_contains_query'
                    features.append(feature)
                elif qstart <= feature.start and qend > feature.start:  # feature partially overlaps (qstart < fstart < qend)
                    feature.attributes['overlap'] = 'partial_by_{}_bases'.format(qend- (feature.start - 1))
                    features.append(feature)
                elif qstart <= feature.end and qend >= feature.end:  # feature partially overlaps (qstart < fend < qend)
                    feature.attributes['overlap'] = 'partial_by_{}_bases'.format(feature.end - qstart)
                    features.append(feature)
        # for f in features:
        #     f.start = f.start - 1 # fix to make gff-coordinates 0-based

        return features

    def get_all_overlapping_features_from_query_unstranded(self, chrom, qstart, qend, strand=None):
        """
        Returns overlapping features from a query.
        
        :param chrom: 
        :param qstart: 
        :param qend: 
        :param strand: 
        :return: 
        """
        if strand is None:  # If no strand given, return all features from both strands
            positive_features = self.get_all_overlapping_features_from_query_stranded(
                chrom, qstart, qend, '+'
            )
            negative_features = self.get_all_overlapping_features_from_query_stranded(
                chrom, qstart, qend, '-'
            )
            return positive_features + negative_features
        else:  # return features on the same strand, only return features on opposite strand if it can't find anything otherwise.
            priority_features = self.get_all_overlapping_features_from_query_stranded(
                chrom, qstart, qend, strand
            )
            if len(priority_features) == 0: # found nothing on the proper strand
                if strand == '+':
                    return self.get_all_overlapping_features_from_query_stranded(
                        chrom, qstart, qend, '-'
                    )
                elif strand == '-':
                    return self.get_all_overlapping_features_from_query_stranded(
                        chrom, qstart, qend, '+'
                    )
                else:
                    print("Strand not correct: {}".format(strand))
                    return []
            else:
                return priority_features

    def get_all_overlapping_features_from_query(self, chrom, qstart, qend, strand, stranded=True):
        """
        Returns overlapping features from a query.
        If stranded is True, then return only the features that overlap on the same strand.
        Otherwise, features returned need not to be with respect to strand (although if 
        strand is given, it will return same-stranded features first before attempting to 
        find those on the opposite strand).
        
        :param chrom: 
        :param qstart: 
        :param qend: 
        :param strand: 
        :param stranded: 
        :return: 
        """
        if stranded:
            return self.get_all_overlapping_features_from_query_stranded(
                chrom, qstart, qend, strand
            )
        else:
            return self.get_all_overlapping_features_from_query_unstranded(
                chrom, qstart, qend, strand=None  # TODO: implement the priority better. For now, this will return all features for both strands
            )

    def parse_annotation_string(self, features_string):
        """
        Splits a feature string into a list of feature strings

        :param features_string:
        :return:
        """
        features = features_string.split('|')
        return features

    def _is_protein_coding(self, transcript_type):
        """
        if defined protein coding, return True else False

        :param transcript_type:
        :return:
        """
        if transcript_type == 'protein_coding':
            return True

        return False

    def return_highest_priority_feature(self, formatted_features, priority):
        """
        Given a list of formatted features, return the feature with the highest
        priority. If there are ties, return just the first one.

        :param formatted_features: list
            list of feature_strings
        :param priority: list
            list containing feature types in order of priority.
        :return:
        """
        # Build dict
        combined_dict = defaultdict(list)
        for feature_string in formatted_features:
            transcript, start, end, strand, feature_type, gene_id, gene_name, transcript_type_list, overlap = feature_string.split(
                ':')
            transcript_type_list = transcript_type_list.split(',')
            for transcript_type in transcript_type_list:
                if self._is_protein_coding(
                        transcript_type):  # simplify all the types at first
                    combined_dict['protein_coding', feature_type].append(
                        feature_string)
                else:
                    combined_dict['non_coding',  feature_type].append(
                        feature_string)
        # return the highest one
        combined_dict = OrderedDict(
            combined_dict
        )  # turn into ordered dict, is that ok?
        # print("combined dict: {}".format(combined_dict))
        try:
            # print("PRIORITY: ", priority)
            for key, _ in combined_dict.iteritems(): # append nonexisting regions to the end of the priority list
                if key not in priority:
                    priority.append(list(key))
            # print("PRIORITY: ",priority)
            combined_dict = sorted(  # sort based on priority list
                combined_dict.iteritems(),
                key=lambda x: priority.index([x[0][0], x[0][1]])
            )
            # print("PRIORITY: {}".format(priority))
            # print("combined dict (sorted): {}".format(combined_dict))
        except ValueError as e:
            print(e)
            print(combined_dict)
        return combined_dict[0]

    def prioritize_transcript_then_gene(self, formatted_features,
                                        transcript_priority, gene_priority):
        """
        Given a list of features, group them first by transcript
        and return the highest priority feature for each transcript.
        Then group each transcript into associated genes and 
        return the highest priority feature.
        
        :param formatted_features: list[string]
        :param transcript_priority: list of tuples (featuretype, transcripttype)
        :param gene_priority: list of tuples (featuretype, transcripttype
        :return priority: string 
        """

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
        # print("PRIORITIZING TRANSCRIPT!!!!!")
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
        # print("PRIORITIZING GENE!!!!!")
        ### PRIORITIZE GENE ###
        for gene, transcripts in unique_genes.iteritems():
            for transcript in transcripts:
                final.append(transcript)
        feature_type, final = self.return_highest_priority_feature(
            final, gene_priority
        )

        if feature_type[0] == 'non_coding':  # TODO: fix. hacky
            return final[0].replace('exon', 'noncoding_exon').replace('intron', 'noncoding_intron')
        return final[0]

    def annotate(self, interval, stranded, transcript_priority, gene_priority):
        """
        Given an interval, annotates using the priority

        :param interval:
        :return:
        """

        overlapping_features = self.get_all_overlapping_features_from_query(
            interval.chrom,
            interval.start,
            interval.end,
            interval.strand,
            stranded
        )

        #   If we find no overlapping features, return 'intergenic' (no feature seen)
        if len(overlapping_features) == 0:
            return 'INTERGENIC', 'INTERGENIC'
        to_append = ''  # full list of genes overlapping features
        transcript = defaultdict(list)
        for feature in overlapping_features:  # for each overlapping feature
            #   for M10 annotations, I see gene featuretypes without transcript_ids.
            #   That's fine. We can re-construct genes from the transcript features, so ignore
            #   anything that doesn't have a transcript_id in it.
            if self.trancript_id_key in feature.attributes.keys():
                for transcript_id in feature.attributes[
                    self.trancript_id_key
                ]:  # multiple genes can be associated with one feature
                    transcript[transcript_id].append(
                        feature)  # append features to their respective genes

        for transcript, features in transcript.iteritems():
            for feature in features:
                # if 'protein_coding' not in feature.attributes['transcript_type']:
                #     if feature.featuretype == 'exon' or feature.featuretype == 'UTR':
                #         feature.featuretype = 'noncoding_exon'
                if feature.featuretype == 'UTR':  # if we only have UTR, specify what kind of UTR (3' or 5').
                    feature.featuretype = self._classify_utr(feature)
                to_append += "{}:{}:{}:{}:{}:".format(
                    transcript,
                    feature.start - 1,  # report 0 based
                    feature.end,
                    feature.strand,
                    feature.featuretype,
                )
                if 'gene_id' in feature.attributes.keys():
                    for t in feature.attributes['gene_id']:
                        to_append += '{},'.format(t)
                    to_append = to_append[:-1] + ':'
                else:
                    to_append = to_append + '-:'
                if 'gene_name' in feature.attributes.keys():
                    for t in feature.attributes['gene_name']:
                        to_append += '{},'.format(t)
                    to_append = to_append[:-1] + ':'
                elif 'transcript_id' in feature.attributes.keys():
                    for t in feature.attributes['transcript_id']:
                        to_append += '{},'.format(t)
                    to_append = to_append[:-1] + ':'
                else:
                    to_append = to_append + '-:'
                if self.type_key in feature.attributes.keys():
                    for t in feature.attributes[self.type_key]:
                        to_append += '{},'.format(t)
                    to_append = to_append[:-1] + ':'
                else:
                    to_append = to_append + '-:'
                if 'overlap' in feature.attributes.keys():
                    for t in feature.attributes['overlap']:
                        to_append += '{},'.format(t)
                    to_append = to_append[:-1] + '|'
                else:
                    to_append = to_append + '-:'

        to_append = to_append[:-1]
        priority = self.prioritize_transcript_then_gene(
            self.parse_annotation_string(to_append),
            transcript_priority,
            gene_priority
        )
        # print(feature.start, feature.end, feature.strand, to_append)
        return priority, to_append


def annotate(db_file, bed_file, out_file, unstranded, chroms,
             transcript_priority, gene_priority, species, append_chr, fuzzy):
    """
    Given a bed6 file, return the file with an extra column containing
    '|' delimited gene annotations

    :param db_file:
    :param bed_file:
    :param out_file:
    :param chroms:
    :param species:
    :param append_chr: boolean
        True if we need to add 'chr' to the database annotations
        For example, wormbase/Ensembl annotations use I/II/1/2/etc.
        instead of chrI/chr1, so we need to let the database know that.
    :return:
    """
    annotator = Annotator(db_file, chroms, species, append_chr, fuzzy)
    bed_tool = pybedtools.BedTool(bed_file)
    with open(out_file, 'w') as o:
        for interval in bed_tool:  # for each line in bed file
            priority, annotation = annotator.annotate(
                interval, unstranded, transcript_priority, gene_priority,
            )
            o.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                interval.chrom, interval.start,
                interval.end, interval.name,
                interval.score, interval.strand,
                priority,
                annotation
            ))