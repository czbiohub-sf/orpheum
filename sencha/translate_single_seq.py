from .constants_translate import (
    STANDARD_CODON_TABLE_MAPPING,
    REVERSE_COMPLEMENT_TABLE,
)


class TranslateSingleSeq:
    def __init__(self, seq, verbose):
        self.seq = seq
        self.verbose = verbose
        self.sign = 1

    @staticmethod
    def _translate_codon(codon):
        if len(codon) == 3:
            return STANDARD_CODON_TABLE_MAPPING.get(codon, "X")
        else:
            # only length 3 is a valid codon, return an empty string otherwise
            return ""

    def _single_seq_translation(self, seq, frame=0):
        translated = ""
        for i in range(int(len(seq) / 3)):
            start = frame + (i * 3)
            end = start + 3
            codon = seq[start:end]
            translated += self._translate_codon(codon)
        return translated

    @staticmethod
    def _reverse_complement(seq):
        return seq.translate(REVERSE_COMPLEMENT_TABLE)[::-1]

    def three_frame_translation(self):

        for frame in range(3):
            if self.sign == 1:
                translation = self._single_seq_translation(self.seq, frame)
            elif self.sign == -1:
                reverse_complement = self._reverse_complement(self.seq)
                translation = self._single_seq_translation(reverse_complement, frame)
            yield translation

    def three_frame_translation_no_stops(self, sign):
        """Remove translations with stop codons & keep track of reading frame"""
        self.sign = sign
        return {
            self.sign * (i + 1): t
            for i, t in enumerate(self.three_frame_translation())
            if "*" not in t
        }

    def three_frame_translation_stops(self, sign):
        """Remove translations with stop codons & keep track of reading frame"""
        self.sign = sign
        return {
            self.sign * (i + 1): t for i, t in enumerate(self.three_frame_translation())
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
