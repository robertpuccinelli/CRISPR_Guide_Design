import logging

logger = logging.getLogger()


class GuideRNAGeneratorBase():
    """Base class for crRNA generation.
    

    ########################################
    # crRNA Generation Process Definitions #
    ########################################
    RS : Repeat Sequence
    TS : Target Sequence
    T  : Transform bases from DNA to RNA or reverse
    RC : Reverse Complement
    J[a,b] : Join element 'b' to end of element 'a'

    ####################################
    # RS (RNA) + Target Sequence (DNA) #
    ####################################
    T^TS
    RC^T^TS
    J[RS, RC^T^TS] (crRNA)
    RC^J[RS, RC^T^TS]
    T^RC^J[RS, RC^T^TS] (DNA template for crRNA)

    ####################################
    # RS (RNA) + Target Sequence (RNA) #
    ####################################

    RC^TS
    J[RS, RC^TS] (crRNA)
    RC^J[RS, RC^TS]
    T^RC^J[RS, RC^TS] (DNA template for crRNA)
    """
    dna_bases = {'A','T','G','C'}
    rna_bases = {'A','U','G','C'}
    dna_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rna_complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    dna_to_rna = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'U'}
    rna_to_dna = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'T'}

    _verified_type = {'NONE': 0,'DNA': 1, 'RNA': 2}

    def __init__(self, seq_repeat: str):
        self.seq_repeat = seq_repeat
        self._seq_target = ''
        self.seq_crRNA = ''
        self.seq_crRNA_DNA_template = ''

    @property
    def seq_repeat(self) -> str:
        """Direct repeat sequence in use."""
        return self._seq_repeat
    
    @seq_repeat.setter
    def seq_repeat(self, input_seq_repeat: str):
        self._sequence_in_processing = input_seq_repeat
        self._validateSequence()
        self._seq_repeat = self._sequence_in_processing

    @property
    def seq_target(self) -> str:
        """Target sequence to be appended to guide."""
        return self._seq_target
    
    @seq_target.setter
    def seq_target(self, input_seq_target: str):
        self._sequence_in_processing = input_seq_target
        self._validateSequence()
        self._seq_target = self._sequence_in_processing

    def crRNAGenerate(self):
        if (len(self.seq_target) & len(self.seq_repeat)) < 1:
            logger.error('seq_target must be defined!')
            return
        
        self._sequence_in_processing = self.seq_target

        # If target sequence is in DNA representation, convert to RNA representation
        if set.issuperset(self.dna_bases, set(self.seq_target)):
            self._transcribe()

        # Generate crRNA sequence
        self._reverseComplement()
        self._sequence_in_processing = self.seq_repeat + self._sequence_in_processing
        self.seq_crRNA = self._sequence_in_processing

        # Generate ssDNA template
        self._reverseComplement()
        self._transcribe()
        self.seq_crRNA_DNA_template = self._sequence_in_processing

    def _reverseComplement(self):
        if set.issuperset(self.dna_bases, set(self._sequence_in_processing)):
            logger.debug('Running reverse complement on DNA sequence.\
                         \n\r  Input sequence: {}'.format(self._sequence_in_processing))
            complement = self.dna_complement

        elif set.issuperset(self.rna_bases, set(self._sequence_in_processing)):
            logger.debug('Running reverse complement on RNA sequence.\
                         \n\r  Input sequence: {}'.format(self._sequence_in_processing))
            complement = self.rna_complement

        else:
            logger.error('Sequence being processed for reverse complement is not strictly RNA or DNA.\
                         \n\r  Input sequence: {}'.format(self._sequence_in_processing))
            return 0
        self._sequence_in_processing = ''.join([complement[base] for base in self._sequence_in_processing[::-1]])
        logger.debug('Reverse complement sequence: {}'.format(self._sequence_in_processing))
        return 1
    
    def _transcribe(self):
        if set.issuperset(self.dna_bases, set(self._sequence_in_processing)):
            logger.debug('Converting DNA to RNA.\
                         \n\r  Input sequence: {}'.format(self._sequence_in_processing))
            look_up_table = self.dna_to_rna

        elif set.issuperset(self.rna_bases, set(self._sequence_in_processing)):
            logger.debug('Converting RNA to DNA.\
                         \n\r  Input sequence: {}'.format(self._sequence_in_processing))
            look_up_table = self.rna_to_dna

        else:
            logger.error('Sequence being processed for conversion is not strictly RNA or DNA.\
                         \n\r  Input sequence: {}'.format(self._sequence_in_processing))
            return 0
        self._sequence_in_processing = ''.join([look_up_table[base] for base in self._sequence_in_processing])
        logger.debug('Converted sequence: {}'.format(self._sequence_in_processing))
        return 1

    def _validateSequence(self):
        valid = 0
        self._sequence_in_processing = self._sequence_in_processing.upper()
        seq_set = set(self._sequence_in_processing)

        if len(seq_set) < 1:
            logger.error('Input string cannot be empty.')
            return valid

        if set.issuperset(self.dna_bases, seq_set):
            self._last_verified_NA = self._verified_type['DNA']
            logger.debug('Valid DNA input sequence: {}'.format(self._sequence_in_processing))
            valid = 1

        elif set.issuperset(self.rna_bases, seq_set):
            self._last_verified_NA = self._verified_type['RNA']
            logger.debug('Valid RNA input sequence: {}'.format(self._sequence_in_processing))
            valid = 1

        else:
            logger.error('Invalid input sequence!\
                         \n\r  DNA should only contain {}\
                         \n\r  RNA should only contain {}\
                         \n\r  Input sequence: {}'.format(self.dna_bases, self.rna_bases, self._sequence_in_processing))

        return valid
        