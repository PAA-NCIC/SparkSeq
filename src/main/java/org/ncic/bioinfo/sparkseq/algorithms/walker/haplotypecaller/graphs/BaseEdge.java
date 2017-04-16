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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs;

import java.io.Serializable;
import java.util.Collection;
import java.util.Comparator;

/**
 * simple edge class for connecting nodes in the graph
 *
 * Works equally well for all graph types (kmer or sequence)
 *
 * Author: wbc
 */
public class BaseEdge {
    private int multiplicity;
    private boolean isRef;

    /**
     * Create a new BaseEdge with weight multiplicity and, if isRef == true, indicates a path through the reference
     *
     * @param isRef indicates whether this edge is a path through the reference
     * @param multiplicity the number of observations of this edge
     */
    public BaseEdge(final boolean isRef, final int multiplicity) {
        if ( multiplicity < 0 ) throw new IllegalArgumentException("multiplicity must be >= 0 but got " + multiplicity);

        this.multiplicity = multiplicity;
        this.isRef = isRef;
    }

    /**
     * Create a new copy of this BaseEdge
     */
    public BaseEdge copy() {
        return new BaseEdge(isRef(), getMultiplicity());
    }

    /**
     * Get the number of observations of paths connecting two vertices
     * @return a positive integer >= 0
     */
    public int getMultiplicity() {
        return multiplicity;
    }

    /**
     * Get the DOT format label for this edge, to be displayed when printing this edge to a DOT file
     * @return a non-null string
     */
    public String getDotLabel() {
        return Integer.toString(getMultiplicity());
    }

    /**
     * Increase the multiplicity of this edge by incr
     * @param incr the change in this multiplicity, must be >= 0
     */
    public void incMultiplicity(final int incr) {
        if ( incr < 0 ) throw new IllegalArgumentException("incr must be >= 0 but got " + incr);
        multiplicity += incr;
    }

    /**
     * A special assessor that returns the multiplicity that should be used by pruning algorithm
     *
     * Can be overloaded by subclasses
     *
     * @return the multiplicity value that should be used for pruning
     */
    public int getPruningMultiplicity() {
        return getMultiplicity();
    }

    /**
     * Set the multiplicity of this edge to value
     * @param value an integer >= 0
     */
    public void setMultiplicity( final int value ) {
        if ( multiplicity < 0 ) throw new IllegalArgumentException("multiplicity must be >= 0");
        multiplicity = value;
    }

    /**
     * Does this edge indicate a path through the reference graph?
     * @return true if so
     */
    public boolean isRef() {
        return isRef;
    }

    /**
     * Indicate that this edge follows the reference sequence, or not
     * @param isRef true if this is a reference edge
     */
    public void setIsRef( final boolean isRef ) {
        this.isRef = isRef;
    }

    /**
     * Sorts a collection of BaseEdges in decreasing order of weight, so that the most
     * heavily weighted is at the start of the list
     */
    public static class EdgeWeightComparator implements Comparator<BaseEdge>, Serializable {
        @Override
        public int compare(final BaseEdge edge1, final BaseEdge edge2) {
            return edge2.multiplicity - edge1.multiplicity;
        }
    }

    /**
     * Add edge to this edge, updating isRef and multiplicity as appropriate
     *
     * isRef is simply the or of this and edge
     * multiplicity is the sum
     *
     * @param edge the edge to add
     * @return this
     */
    public BaseEdge add(final BaseEdge edge) {
        if ( edge == null ) throw new IllegalArgumentException("edge cannot be null");
        this.multiplicity += edge.getMultiplicity();
        this.isRef = this.isRef || edge.isRef();
        return this;
    }

    /**
     * Create a new BaseEdge with multiplicity and isRef that's an or of all edges
     *
     * @param edges a collection of edges to or their isRef values
     * @param multiplicity our desired multiplicity
     * @return a newly allocated BaseEdge
     */
    public static BaseEdge orRef(final Collection<BaseEdge> edges, final int multiplicity) {
        for ( final BaseEdge e : edges )
            if ( e.isRef() )
                return new BaseEdge(true, multiplicity);
        return new BaseEdge(false, multiplicity);
    }

    /**
     * Return a new edge whose multiplicity is the max of this and edge, and isRef is or of this and edge
     *
     * isRef is simply the or of this and edge
     * multiplicity is the max
     *
     * @param edge the edge to max
     */
    public BaseEdge max(final BaseEdge edge) {
        if ( edge == null ) throw new IllegalArgumentException("edge cannot be null");
        return new BaseEdge(isRef() || edge.isRef(), Math.max(getMultiplicity(), edge.getMultiplicity()));
    }

    @Override
    public String toString() {
        return "BaseEdge{" +
                "multiplicity=" + multiplicity +
                ", isRef=" + isRef +
                '}';
    }
}
