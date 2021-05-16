# BAM File Format
# ===============

module BAM

using BioGenerics
using GenomicFeatures
using XAM.SAM

import BGZFStreams
import BioAlignments
import Indexes
import BioSequences
import BioGenerics: isfilled, header

import GenomicFeatures: eachoverlap

# Overwrite virtualoffset method as per https://github.com/BioJulia/BGZFStreams.jl/issues/22#issuecomment-839746122.
function BGZFStreams.virtualoffset(stream::BGZFStreams.BGZFStream)
    if stream.mode == BGZFStreams.READ_MODE
        i = BGZFStreams.ensure_buffered_data(stream)
        if i == 0
            block = stream.blocks[end]
        else
            block = stream.blocks[i]
        end
    else
        block = stream.blocks[1]
    end
    return BGZFStreams.VirtualOffset(block.block_offset, block.position - 1)
end

include("bai.jl")
include("auxdata.jl")
include("reader.jl")
include("record.jl")
include("writer.jl")
include("overlap.jl")

end
