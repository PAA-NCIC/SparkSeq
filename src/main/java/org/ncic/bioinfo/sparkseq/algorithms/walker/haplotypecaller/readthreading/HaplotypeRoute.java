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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readthreading;

import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.MultiSampleEdge;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.graphs.Route;

import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

/**
 * Author: wbc
 */
public class HaplotypeRoute extends Route<MultiDeBruijnVertex,MultiSampleEdge> {

    protected final Set<MultiDeBruijnVertex> vertexSet;

    protected final Map<MultiDeBruijnVertex,Integer> vertexOrder;

    protected final Set<MultiDeBruijnVertex> forkAndJoins;

    /**
     * Constructs a HaplotypeRoute given its route.
     *
     * @param route the haplotype route.
     */
    public HaplotypeRoute(final Route<MultiDeBruijnVertex, MultiSampleEdge> route) {
        super(route);
        vertexOrder = new LinkedHashMap<>(route.length() + 1);
        int nextOrder = 0;
        vertexOrder.put(getFirstVertex(),nextOrder++);
        for (final MultiSampleEdge edge : edgesInOrder)
            vertexOrder.put(graph.getEdgeTarget(edge), nextOrder++);
        Route<MultiDeBruijnVertex,MultiSampleEdge> currentRoute = this;
        forkAndJoins = new HashSet<>(route.length());
        while (currentRoute != null) {
            if (currentRoute.lastVertexIsForkOrJoin())
                forkAndJoins.add(currentRoute.getLastVertex());
            currentRoute = currentRoute.getPrefixRouteWithLastVertexThatIsForkOrJoin();
        }
        vertexSet = Collections.unmodifiableSet(new HashSet<>(vertexOrder.keySet()));
    }



    @SuppressWarnings("unused")
    public Route<MultiDeBruijnVertex,MultiSampleEdge> subRoute(final MultiDeBruijnVertex start, final MultiDeBruijnVertex end) {
        final Integer startOrder = vertexOrder.get(start);
        final Integer endOrder = vertexOrder.get(end);
        if (startOrder == null || endOrder == null)
            return null;
        else if (startOrder > endOrder)
            return null;
        else {
            Route<MultiDeBruijnVertex,MultiSampleEdge> result = new Route<>(start,graph);
            for (final MultiSampleEdge edge : edgesInOrder.subList(startOrder,endOrder))
                result = new Route(result,edge);
            return result;
        }
    }

    /**
     * Returns the set of vertex on the route.
     * @return read only, never {@code null} vertex set.
     */
    public Set<MultiDeBruijnVertex> vertexSet() {
        return vertexSet;
    }


    /**
     * Returns the position of the vertex in the route.
     *
     * @param vertex the query vertex.
     *
     * @throws NullPointerException if {@code vertex} is {@code null}.
     *
     * @return -1 if there is no such a vertex in the route, otherwise a number between 0 and {@link #length()} - 1.
     */
    public int getVertexPosition(final MultiDeBruijnVertex vertex) {
        final Integer result = vertexOrder.get(vertex);
        return result == null ? -1 : result;
    }
}
