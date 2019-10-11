from Bio import Align

aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -0.5
aligner.extend_gap_score = -0.1
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0


def score(seqA, seqB):
    return aligner.score(seqA, seqB)
