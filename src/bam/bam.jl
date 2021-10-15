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

# Overwrite `ensure_buffered_data` to address https://github.com/BioJulia/BGZFStreams.jl/issues/22.
@inline function BGZFStreams.ensure_buffered_data(stream)
    #@assert stream.mode == READ_MODE
    @label doit
    while stream.block_index ≤ lastindex(stream.blocks)
        @inbounds block = stream.blocks[stream.block_index]
        if BGZFStreams.is_eof_block(block.compressed_block) # Note: `read_blocks!` does not necessarily fill/overwrite blocks till `lastindex(stream.blocks)`, we need to stop incrementing `stream.block_index` when an eof block is encountered.
            break
        end
        if block.position ≤ block.size
            return stream.block_index
        end
        stream.block_index += 1
    end
    if !BGZFStreams.eof(stream.io)
        BGZFStreams.read_blocks!(stream)
        @goto doit
    end
    return 0
end

include("bai.jl")
include("auxdata.jl")
include("reader.jl")
include("record.jl")
include("writer.jl")
include("overlap.jl")

end
