# BAM File Format
# ===============

module BAM

using BioGenerics
using GenomicFeatures
using XAM.SAM
import ..XAM: flag, XAMRecord, XAMReader, XAMWriter,
	ismapped, isprimary, ispositivestrand, isnextmapped #TODO: Deprecate import of flag queries. These were imported to preseve existing API.

import BGZFStreams
import BioAlignments
import Indexes
import BioSequences
import BioGenerics: isfilled, header

import GenomicFeatures: eachoverlap


include("bai.jl")
include("auxdata.jl")
include("reader.jl")
include("record.jl")
include("writer.jl")
include("overlap.jl")

end
