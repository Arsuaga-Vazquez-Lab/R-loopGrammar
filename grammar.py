import region_handler


class GrammarSymbol:
    ALPHA = '\N{greek small letter alpha}'
    OMEGA = '\N{greek small letter omega}'
    DELTA = '\N{greek small letter delta}'
    BETA = '\N{greek small letter beta}'
    SIGMA = '\N{greek small letter sigma}'
    SIGMA_HAT = '\N{greek small letter sigma with macron}'
    TAU = '\N{greek small letter tau}'
    TAU_HAT = '\N{greek small letter tau with macron}'
    RHO = '\N{greek small letter rho}'
    GAMMA = '\N{greek small letter gamma}'


class GrammarDictionary:
    pass

def _partition_rloop(sequence, rloop: bed.GenomicRegion):
	before = sequence[:rloop.start]
	inside = sequence[rloop.start:rloop.end]
	after = sequence[rloop.end:]

	return before, inside, after


class ParsingBlock:
	def __init__(self, block, width):
		self.block = block
		self.sub_blocks = [self.block[i:i + width] for i in range(0, len(block), width)]

def create_dictionary(sequence: str, rloops: List[bed.GenomicRegion],
                      region_summaries: List[region_handler.RegionSummary]):
    for rloop in rloops:
    	before, inside, after = _partition_rloop(sequence, rloop)
    	inside_reverse = inside[::-1]
    	after_reverse = after[::-1]
