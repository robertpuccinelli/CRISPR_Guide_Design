"""crRNA generation utilities."""
import logging

logger = logging.getLogger()


class GuideRNAGeneratorBase():
    """Base class for crRNA generation.
    
    Parameters
    ----------
    seq_repeat : str
        Initial repeat sequence to match the desired Cas application

    Attributes
    ----------
    seq_repeat : str
        Direct repeat sequence for Cas.
    seq_target : str
        Target sequence for hybridization.
    seq_crRNA : str
        Final crRNA sequence to be loaded in the ribonucleoprotein complex.
    seq_crRNA_DNA_template : str
        Reverse transcription of the final crRNA sequence.

    """
    dna_bases = {'A','T','G','C'}
    rna_bases = {'A','U','G','C'}
    dna_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rna_complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    dna_to_rna = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'U'}
    rna_to_dna = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'T'}

    def __init__(self, seq_repeat: str):
        self.seq_repeat = seq_repeat
        self._seq_target = ''
        self.seq_crRNA = ''
        self.seq_crRNA_DNA_template = ''

    @property
    def seq_repeat(self) -> str:
        """Return direct repeat sequence in use."""
        return self._seq_repeat
    
    @seq_repeat.setter
    def seq_repeat(self, input_seq_repeat: str):
        """Set direct repeat sequence."""
        self._sequence_in_processing = input_seq_repeat
        self._validateSequence()
        self._seq_repeat = self._sequence_in_processing
        logger.info('Repeat sequence: {}'.format(self._seq_repeat))

    @property
    def seq_target(self) -> str:
        """Return target sequence in use."""
        return self._seq_target
    
    @seq_target.setter
    def seq_target(self, input_seq_target: str):
        """Set target sequence."""
        self._sequence_in_processing = input_seq_target
        self._validateSequence()
        self._seq_target = self._sequence_in_processing
        logger.info('Target sequence: {}'.format(self._seq_target))

    def crRNAGenerate(self):
        """Synthesize crRNA and ssDNA template from provided repeat and target sequences.
        
        ########################################
        # crRNA Generation Process Definitions #
        ########################################
        RS : Repeat Sequence
        TS : Target Sequence
        T  : Transform bases from DNA to RNA or reverse
        RC : Reverse Complement
        J[a,b] : Join element 'a' and 'b'

        ####################################
        # RS (RNA) + Target Sequence (DNA) #
        ####################################
        1. T^TS
        2. RC^T^TS
        3. J[RS, RC^T^TS] (crRNA)
        4. RC^J[RS, RC^T^TS]
        5. T^RC^J[RS, RC^T^TS] (DNA template for crRNA)
        6. Postprocess (Generic overloaded function)

        ####################################
        # RS (RNA) + Target Sequence (RNA) #
        ####################################

        1. RC^TS
        2. J[RS, RC^TS] (crRNA)
        3. RC^J[RS, RC^TS]
        4. T^RC^J[RS, RC^TS] (DNA template for crRNA)
        5. Postprocess (Generic overloaded function)
        """

        if (len(self.seq_target) & len(self.seq_repeat)) < 1:
            logger.error('seq_target must be defined!')
            return
        
        self._sequence_in_processing = self.seq_target

        # If target sequence is in DNA representation, convert to RNA representation
        if set.issuperset(self.dna_bases, set(self.seq_target)):
            self._transform()

        # Generate crRNA sequence
        self._reverseComplement()
        self._sequence_in_processing = self.seq_repeat + self._sequence_in_processing
        self.seq_crRNA = self._sequence_in_processing

        # Generate ssDNA template
        self._reverseComplement()
        self._transform()
        self.seq_crRNA_DNA_template = self._sequence_in_processing

        self._postprocess()

        logger.info('crRNA sequence: {}'.format(self.seq_crRNA))
        logger.info('ssDNA template: {}'.format(self.seq_crRNA_DNA_template))

    def _postprocess(self):
        pass

    def _reverseComplement(self):
        """Perform reverse complement of RNA or DNA sequence."""

        if set.issuperset(self.dna_bases, set(self._sequence_in_processing)):
            logger.debug('Running reverse complement on DNA sequence. Input sequence: {}'.format(self._sequence_in_processing))
            complement = self.dna_complement

        elif set.issuperset(self.rna_bases, set(self._sequence_in_processing)):
            logger.debug('Running reverse complement on RNA sequence. Input sequence: {}'.format(self._sequence_in_processing))
            complement = self.rna_complement

        else:
            logger.error('Sequence being processed for reverse complement is not strictly RNA or DNA. Input sequence: {}'.format(self._sequence_in_processing))
            return 0
        self._sequence_in_processing = ''.join([complement[base] for base in self._sequence_in_processing[::-1]])
        logger.debug('Reverse complement sequence: {}'.format(self._sequence_in_processing))
        return 1
    
    def _transform(self):
        """Swap U's and T's in the sequence to convert between RNA and DNA."""

        if set.issuperset(self.dna_bases, set(self._sequence_in_processing)):
            logger.debug('Converting DNA to RNA. Input sequence: {}'.format(self._sequence_in_processing))
            look_up_table = self.dna_to_rna

        elif set.issuperset(self.rna_bases, set(self._sequence_in_processing)):
            logger.debug('Converting RNA to DNA. Input sequence: {}'.format(self._sequence_in_processing))
            look_up_table = self.rna_to_dna

        else:
            logger.error('Sequence being processed for conversion is not strictly RNA or DNA. Input sequence: {}'.format(self._sequence_in_processing))
            return 0
        self._sequence_in_processing = ''.join([look_up_table[base] for base in self._sequence_in_processing])
        logger.debug('Converted sequence: {}'.format(self._sequence_in_processing))
        return 1

    def _validateSequence(self):
        """Check for unusual characters and change string to upper-case."""

        valid = 0
        self._sequence_in_processing = self._sequence_in_processing.upper()
        seq_set = set(self._sequence_in_processing)

        if len(seq_set) < 1:
            logger.error('Input string cannot be empty.')
            return valid

        if set.issuperset(self.dna_bases, seq_set):
            logger.debug('Valid DNA input sequence: {}'.format(self._sequence_in_processing))
            valid = 1

        elif set.issuperset(self.rna_bases, seq_set):
            logger.debug('Valid RNA input sequence: {}'.format(self._sequence_in_processing))
            valid = 1

        else:
            logger.error('Invalid input sequence! DNA should only contain {}. RNA should only contain {}. Input contains: {}'.format(self.dna_bases, self.rna_bases, set(self._sequence_in_processing)))

        return valid
        

class T7crRNAGenerator(GuideRNAGeneratorBase):
    """A crRNA generator that appends modifications required to synthesize crRNA with T7.
    
    See GuideRNAGeneratorBase class for documentation.
    """

    _t7_reverse_comp_promoter = 'CTATAGTGAGTCGTATTA'
    _t7_initial_transcription = 'G'

    def _postprocess(self):
        """Appends T7 promoter and inital transcribed nucleotides to crRNA and DNA template."""
        super()._postprocess()
        self.seq_crRNA_DNA_template = self.seq_crRNA_DNA_template + self._t7_reverse_comp_promoter
        self.seq_crRNA = self._t7_initial_transcription + self.seq_crRNA
