import warnings


class TranslateSingleSeq:
    def __init__(self, seq, verbose):
        self.seq = seq
        self.verbose = verbose

    def three_frame_translation(self):
        if self.verbose:
            warning_filter = 'default'
        else:
            warning_filter = 'ignore'

        with warnings.catch_warnings():
            warnings.simplefilter(warning_filter)
            for frame in range(3):
                translation = self.seq[frame:].translate()
                yield translation

    def three_frame_translation_no_stops(self, sign):
        """Remove translations with stop codons &
        keep track of reading frame"""
        return {
            sign * (i + 1): t
            for i, t in enumerate(self.three_frame_translation(self.seq))
            if '*' not in t
        }

    def six_frame_translation_no_stops(self):
        forward_translations = self.three_frame_translation_no_stops(
            self.seq, sign=1)

        # Sign=-1 sets the reading frames as negative
        # to make it obvious they are
        # from the reverse strand
        reverse_translations = self.three_frame_translation_no_stops(
            self.seq.reverse_complement(), sign=-1)
        forward_translations.update(reverse_translations)
        return forward_translations
