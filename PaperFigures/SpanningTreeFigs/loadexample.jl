# load parameters specific to each example

function loadexample(whichcase)
    if whichcase==1
        n = 9
        numEdges = 12
        nblocks = 218
        adjvec = [0	 1	0	1	0	0	0	0	0	1	0	1	0	1	0	0	0	0	0	1	0	0	0	1	0	0	0	1	0	0	0	1	0	1	0	0	0	1	0	1	0	1	0	1	0	0	0	1	0	1	0	0	0	1	0	0	0	1	0	0	0	1	0	0	0	0	0	1	0	1	0	1	0	0	0	0	0	1	0	1	0]
    elseif whichcase==2
        n = 7
        numEdges = 12
        nblocks = 95
        adjvec = [0 1	1	1	1	1	1	1	0	1	0	0	0	1	1	1	0	1	0	0	0	1	0	1	0	1	0	0	1	0	0	1	0	1	0	1	0	0	0	1	0	1	1	1	0	0	0	1	0]
    elseif whichcase==3
        n = 8
        numEdges = 12
        nblocks = 167
        adjvec = [0	1	1	0	1	0	0	0	1	0	0	1	0	1	0	0	1	0	0	1	0	0	1	0	0	1	1	0	0	0	0	1	1	0	0	0	0	1	1	0	0	1	0	0	1	0	0	1	0	0	1	0	1	0	0	1	0	0	0	1	0	1	1	0]
    elseif whichcase ==4
        n = 6
        numEdges = 15
        nblocks = 63
        adjvec = [0	 1	1	1	1	1	1	0	1	1	1	1	1	1	0	1	1	1	1	1	1	0	1	1	1	1	1	1	0	1	1	1	1	1	1	0]
    end

    return n,numEdges,nblocks,adjvec

end