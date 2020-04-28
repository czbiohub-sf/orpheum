import warnings


class TranslateSingleSeq:
    def __init__(self, seq, verbose):
        self.seq = seq
        self.verbose = verbose
        self.sign = 1

    def three_frame_translation(self):
        from Bio import BiopythonWarning
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            for frame in range(3):
                if self.sign == 1:
                    translation = self.seq[frame:].translate()
                elif self.sign == -1:
                    translation = self.seq.reverse_complement()[
                        frame:].translate()
                yield translation

    def three_frame_translation_no_stops(self, sign):
        """Remove translations with stop codons &
        keep track of reading frame"""
        self.sign = sign
        return {
            self.sign * (i + 1): t
            for i, t in enumerate(self.three_frame_translation())
            if '*' not in t
        }

    def three_frame_translation_stops(self, sign):
        """Remove translations with stop codons &
        keep track of reading frame"""
        self.sign = sign
        return {
            self.sign * (i + 1): t
            for i, t in enumerate(self.three_frame_translation())
        }

    def six_frame_translation_no_stops(self):
        forward_translations = self.three_frame_translation_no_stops(1)

        # Sign=-1 sets the reading frames as negative
        # to make it obvious they are
        # from the reverse strand
        reverse_translations = self.three_frame_translation_no_stops(-1)
        forward_translations.update(reverse_translations)
        return forward_translations

    def six_frame_translation(self):
        forward_translations = self.three_frame_translation_stops(1)

        # Sign=-1 sets the reading frames as negative
        # to make it obvious they are
        # from the reverse strand
        reverse_translations = self.three_frame_translation_stops(-1)
        forward_translations.update(reverse_translations)
        return forward_translations
