/*
 * Copyright (c) 2017 NCIC, Institute of Computing Technology, Chinese Academy of Sciences
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller;

import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.SeqGraph;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readthreading.ReadThreadingGraph;

/**
 * Author: wbc
 */
public class AssemblyResult {
    private final Status status;
    private ReadThreadingGraph threadingGraph;
    private final SeqGraph graph;

    /**
     * Create a new assembly result
     * @param status the status, cannot be null
     * @param graph the resulting graph of the assembly, can only be null if result is failed
     */
    public AssemblyResult(final Status status, final SeqGraph graph) {
        if ( status == null ) throw new IllegalArgumentException("status cannot be null");
        if ( status != Status.FAILED && graph == null ) throw new IllegalArgumentException("graph is null but status is " + status);

        this.status = status;
        this.graph = graph;
    }

    /**
     * Returns the threading-graph associated with this assembly-result.
     */
    public void setThreadingGraph(final ReadThreadingGraph threadingGraph) {
        this.threadingGraph = threadingGraph;
    }

    public ReadThreadingGraph getThreadingGraph() {
        return threadingGraph;
    }

    public Status getStatus() { return status; }
    public SeqGraph getGraph() { return graph; }

    public int getKmerSize() {
        return graph.getKmerSize();
    }


    /**
     * Status of the assembly result
     */
    public enum Status {
        /** Something went wrong, and we couldn't produce a meaningful graph */
        FAILED,
        /** Assembly succeeded, but graph degenerated into just the reference sequence */
        JUST_ASSEMBLED_REFERENCE,
        /** Assembly succeeded, and the graph has some meaningful structure */
        ASSEMBLED_SOME_VARIATION
    }
}
