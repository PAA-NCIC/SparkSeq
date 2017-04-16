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

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.data.basic.CountSet;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.EventMap;
import org.ncic.bioinfo.sparkseq.algorithms.utils.haplotype.Haplotype;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.readthreading.ReadThreadingGraph;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Author: wbc
 */
public class AssemblyResultSet {

    private final Map<Integer, AssemblyResult> assemblyResultByKmerSize;
    private final Set<Haplotype> haplotypes;
    private final Map<Haplotype, AssemblyResult> assemblyResultByHaplotype;
    private ActiveRegion regionForGenotyping;
    private byte[] fullReferenceWithPadding;
    private GenomeLoc paddedReferenceLoc;
    private boolean variationPresent;
    private Haplotype refHaplotype;
    private boolean wasTrimmed = false;
    private final CountSet kmerSizes;
    private TreeSet<VariantContext> variationEvents;
    private boolean debug;
    private static Logger logger = Logger.getLogger(AssemblyResultSet.class);

    /**
     * Constructs a new empty assembly result set.
     */
    public AssemblyResultSet() {
        assemblyResultByKmerSize = new LinkedHashMap<>(4);
        haplotypes = new LinkedHashSet<>(10);
        assemblyResultByHaplotype = new LinkedHashMap<>(10);
        kmerSizes = new CountSet(4);
    }


    /**
     * Change the debug status for this assembly-result-set.
     *
     * @param newValue new value for the debug status.
     */
    void setDebug(final boolean newValue) {
        debug = newValue;
    }

    /**
     * Trims an assembly result set down based on a new set of trimmed haplotypes.
     *
     * @param trimmedActiveRegion the trimmed down active region.
     * @return never {@code null}, a new trimmed assembly result set.
     * @throws NullPointerException     if any argument in {@code null} or
     *                                  if there are {@code null} entries in {@code originalByTrimmedHaplotypes} for trimmed haplotype keys.
     * @throws IllegalArgumentException if there is no reference haplotype amongst the trimmed ones.
     */
    public AssemblyResultSet trimTo(final ActiveRegion trimmedActiveRegion) {

        final Map<Haplotype, Haplotype> originalByTrimmedHaplotypes = calculateOriginalByTrimmedHaplotypes(trimmedActiveRegion);
        if (refHaplotype == null) throw new IllegalStateException();
        if (trimmedActiveRegion == null) throw new NullPointerException();
        final AssemblyResultSet result = new AssemblyResultSet();

        for (final Haplotype trimmed : originalByTrimmedHaplotypes.keySet()) {
            final Haplotype original = originalByTrimmedHaplotypes.get(trimmed);
            if (original == null)
                throw new NullPointerException("all trimmed haplotypes must have an original one");
            final AssemblyResult as = assemblyResultByHaplotype.get(original);
            if (as == null) result.add(trimmed);
            else result.add(trimmed, as);
        }

        result.setRegionForGenotyping(trimmedActiveRegion);
        result.setFullReferenceWithPadding(this.fullReferenceWithPadding);
        result.setPaddedReferenceLoc(this.paddedReferenceLoc);
        if (result.refHaplotype == null)
            throw new IllegalStateException("missing reference haplotype in the trimmed set");
        result.wasTrimmed = true;
        return result;
    }

    private Map<Haplotype, Haplotype> calculateOriginalByTrimmedHaplotypes(final ActiveRegion trimmedActiveRegion) {
        if (debug)
            logger.info("Trimming active region " + getRegionForGenotyping() + " with " + getHaplotypeCount() + " haplotypes");

        final List<Haplotype> haplotypeList = getHaplotypeList();

        // trim down the haplotypes
        final Map<Haplotype, Haplotype> originalByTrimmedHaplotypes = new HashMap<>();

        for (final Haplotype h : haplotypeList) {
            final Haplotype trimmed = h.trim(trimmedActiveRegion.getExtendedLoc());

            if (trimmed != null) {
                if (originalByTrimmedHaplotypes.containsKey(trimmed)) {
                    if (trimmed.isReference()) {
                        originalByTrimmedHaplotypes.remove(trimmed);
                        originalByTrimmedHaplotypes.put(trimmed, h);
                    }
                } else
                    originalByTrimmedHaplotypes.put(trimmed, h);
            } else if (h.isReference())
                throw new IllegalStateException("trimming eliminates the reference haplotype");
            else if (debug) {
                logger.info("Throwing out haplotype " + h + " with cigar " + h.getCigar() +
                        " because it starts with or ends with an insertion or deletion when trimmed to " +
                        trimmedActiveRegion.getExtendedLoc());
            }
        }

        // create the final list of trimmed haplotypes
        final List<Haplotype> trimmedHaplotypes = new ArrayList<>(originalByTrimmedHaplotypes.keySet());

        // resort the trimmed haplotypes.
        Collections.sort(trimmedHaplotypes, new HaplotypeSizeAndBaseComparator());
        final Map<Haplotype, Haplotype> sortedOriginalByTrimmedHaplotypes = new LinkedHashMap<>(trimmedHaplotypes.size());
        for (final Haplotype trimmed : trimmedHaplotypes)
            sortedOriginalByTrimmedHaplotypes.put(trimmed, originalByTrimmedHaplotypes.get(trimmed));


        if (debug) logger.info("Trimmed region to " + trimmedActiveRegion.getLocation() + " size " +
                trimmedActiveRegion.getLocation().size() + " reduced number of haplotypes from " +
                haplotypeList.size() + " to only " + trimmedHaplotypes.size());
        if (debug)
            for (final Haplotype remaining : trimmedHaplotypes)
                logger.info("Remains: " + remaining + " cigar " + remaining.getCigar());
        return sortedOriginalByTrimmedHaplotypes;
    }

    /**
     * Query the reference haplotype in the result set.
     *
     * @return {@code null} if none wasn't yet added, otherwise a reference haplotype.
     */
    public Haplotype getReferenceHaplotype() {
        return refHaplotype;
    }

    /**
     * Checks whether there is any variation present in the assembly result set.
     * <p>
     * <p>
     * This is equivalent to whether there is more than one haplotype.
     * </p>
     *
     * @return {@code true} if there is variation present, {@code false} otherwise.
     */
    public boolean isVariationPresent() {
        return variationPresent && haplotypes.size() > 1;
    }

    /**
     * Dumps debugging information into a print-writer.
     *
     * @param pw where to dump the information.
     * @throws NullPointerException if {@code pw} is {@code null}.
     */
    public void debugDump(final PrintWriter pw) {
        if (getHaplotypeList().size() == 0) {
            return;
        }
        pw.println("Active Region " + this.regionForGenotyping.getLocation());
        pw.println("Extended Act Region " + this.getRegionForGenotyping().getExtendedLoc());
        pw.println("Ref haplotype coords " + getHaplotypeList().get(0).getGenomeLocation());
        pw.println("Haplotype count " + haplotypes.size());
        final Map<Integer, Integer> kmerSizeToCount = new HashMap<>();

        for (final Map.Entry<Haplotype, AssemblyResult> e : assemblyResultByHaplotype.entrySet()) {
            final AssemblyResult as = e.getValue();
            final int kmerSize = as.getGraph().getKmerSize();
            if (kmerSizeToCount.containsKey(kmerSize)) {
                kmerSizeToCount.put(kmerSize, kmerSizeToCount.get(kmerSize) + 1);
            } else {
                kmerSizeToCount.put(kmerSize, 1);
            }
        }
        pw.println("Kmer sizes count " + kmerSizeToCount.entrySet().size());
        Integer[] kmerSizes = new Integer[kmerSizeToCount.size()];
        kmerSizes = kmerSizeToCount.keySet().toArray(kmerSizes);
        Arrays.sort(kmerSizes);
        pw.println("Kmer sizes values " + Arrays.toString(kmerSizes));
        for (int size : kmerSizes) {
            pw.println("Kmer size " + size + " count " + kmerSizeToCount.get(size));
        }
    }

    /**
     * Adds a haplotype to the result set without indicating a generating assembly result.
     * <p>
     * <p>
     * It is possible to call this method with the same haplotype several times. In that the second and further
     * calls won't have any effect (thus returning {@code false}).
     * </p>
     *
     * @param h the haplotype to add to the assembly result set.
     * @return {@code true} if the assembly result set has been modified as a result of this call.
     * @throws NullPointerException     if {@code h} is {@code null}
     * @throws IllegalArgumentException if {@code h} does not have a genome location.
     */
    public boolean add(final Haplotype h) {
        if (h == null) throw new NullPointerException("input haplotype cannot be null");
        if (h.getGenomeLocation() == null)
            throw new IllegalArgumentException("the haplotype provided must have a genomic location");
        if (haplotypes.contains(h))
            return false;
        haplotypes.add(h);
        updateReferenceHaplotype(h);
        return true;
    }

    /**
     * Adds simultaneously a haplotype and the generating assembly-result.
     * <p>
     * <p>
     * Haplotypes and their assembly-result can be added multiple times although just the first call will have
     * any effect (return value is {@code true}).
     * </p>
     *
     * @param h  haplotype to add.
     * @param ar assembly-result that is assumed to have given rise to that haplotype.
     * @return {@code true} iff this called changes the assembly result set.
     * @throws NullPointerException     if {@code h} or {@code ar} is {@code null}.
     * @throws IllegalArgumentException if {@code h} has not defined genome location.
     */
    public boolean add(final Haplotype h, final AssemblyResult ar) {
        if (h == null) throw new NullPointerException("input haplotype cannot be null");
        if (ar == null) throw new NullPointerException("input assembly-result cannot be null");
        if (h.getGenomeLocation() == null)
            throw new IllegalArgumentException("the haplotype provided must have a genomic location");

        final boolean assemblyResultAdditionReturn = add(ar);

        if (haplotypes.contains(h)) {
            final AssemblyResult previousAr = assemblyResultByHaplotype.get(h);
            if (previousAr == null) {
                assemblyResultByHaplotype.put(h, ar);
                return true;
            } else if (!previousAr.equals(ar))
                throw new IllegalStateException("there is already a different assembly result for the input haplotype");
            else
                return assemblyResultAdditionReturn;
        } else {
            haplotypes.add(h);
            assemblyResultByHaplotype.put(h, ar);
            updateReferenceHaplotype(h);
            if (h.isNonReference()) variationPresent = true;
            return true;
        }
    }

    /**
     * Add a assembly-result object.
     *
     * @param ar the assembly result to add.
     * @return {@code true} iff this addition changed the assembly result set.
     * @throws NullPointerException  if {@code ar} is {@code null}.
     * @throws IllegalStateException if there is an assembly result with the same kmerSize.
     */
    public boolean add(final AssemblyResult ar) {
        if (ar == null)
            throw new NullPointerException();
        final int kmerSize = ar.getKmerSize();
        if (assemblyResultByKmerSize.containsKey(kmerSize)) {
            if (!assemblyResultByKmerSize.get(kmerSize).equals(ar))
                throw new IllegalStateException("a different assembly result with the same kmerSize was already added");
            return false;
        } else {
            assemblyResultByKmerSize.put(kmerSize, ar);
            kmerSizes.add(kmerSize);
            return true;
        }
    }

    /**
     * Returns the current region for genotyping.
     *
     * @return might be {@code null}.
     */
    public ActiveRegion getRegionForGenotyping() {
        return regionForGenotyping;
    }

    /**
     * Sets the region for genotyping.
     *
     * @param regionForGenotyping the new value.
     */
    public void setRegionForGenotyping(final ActiveRegion regionForGenotyping) {
        this.regionForGenotyping = regionForGenotyping;
    }

    /**
     * Returns the current full reference with padding.
     *
     * @return might be {@code null}.
     */
    public byte[] getFullReferenceWithPadding() {
        return fullReferenceWithPadding;
    }

    /**
     * Sets the full reference with padding base sequence.
     *
     * @param fullReferenceWithPadding the new value.
     */
    public void setFullReferenceWithPadding(final byte[] fullReferenceWithPadding) {
        this.fullReferenceWithPadding = fullReferenceWithPadding;
    }

    /**
     * Returns the padded reference location.
     *
     * @return might be {@code null}
     */
    public GenomeLoc getPaddedReferenceLoc() {
        return paddedReferenceLoc;
    }

    /**
     * Changes the padded reference location.
     *
     * @param paddedReferenceLoc the new value.
     */
    public void setPaddedReferenceLoc(final GenomeLoc paddedReferenceLoc) {
        this.paddedReferenceLoc = paddedReferenceLoc;
    }

    /**
     * Returns the number of haplotypes in the assembly result set.
     *
     * @return {@code 0} or greater.
     */
    public int getHaplotypeCount() {
        return haplotypes.size();
    }

    /**
     * Returns the haplotypes as a list.
     * <p>
     * <p>
     * The result is unmodifiable.
     * </p>
     *
     * @return never {@code null}, but perhaps a empty list if no haplotype was generated during assembly.
     */
    public List<Haplotype> getHaplotypeList() {
        return Arrays.asList(haplotypes.toArray(new Haplotype[haplotypes.size()]));
    }

    /**
     * Returns the maximum kmerSize available.
     *
     * @return greater than 0.
     * @throws IllegalStateException if no assembly-result was added to the set, thus there is no kmerSize.
     */
    public int getMaximumKmerSize() {
        if (kmerSizes.size() == 0)
            throw new IllegalStateException("there is yet no kmerSize in this assembly result set");
        return kmerSizes.max();
    }

    /**
     * Indicates whether there are more than one kmerSize in the set.
     *
     * @return {@code true} iff there is more than one kmerSize assembly in the set.
     */
    public boolean hasMultipleKmerSizes() {
        return kmerSizes.size() > 1;
    }

    /**
     * Returns the minimum kmerSize available.
     *
     * @return greater than 0.
     * @throws IllegalStateException if no assembly-result was added to the set, thus there is no kmerSize.
     */
    public int getMinimumKmerSize() {
        if (kmerSizes.size() == 0)
            throw new IllegalStateException("there is yet no kmerSize in this assembly result set");
        return kmerSizes.min();
    }

    /**
     * Returns a read-threading graph in the assembly set that has a particular kmerSize.
     *
     * @param kmerSize the requested kmerSize.
     * @return {@code null} if there is no read-threading-graph amongst assembly results with that kmerSize.
     */
    public ReadThreadingGraph getUniqueReadThreadingGraph(final int kmerSize) {
        final AssemblyResult assemblyResult = assemblyResultByKmerSize.get(kmerSize);
        if (assemblyResult == null) return null;
        return assemblyResult.getThreadingGraph();
    }

    /**
     * Checks whether this assembly result set was trimmed.
     *
     * @return {@code true} iff this assembly result set was trimmed.
     */
    public boolean wasTrimmed() {
        return wasTrimmed;
    }

    /**
     * Marks the assembly as not having variation even if it has more than one haplotype.
     */
    public void resetVariationPresent() {
        variationPresent = false;
    }

    /**
     * Dumps debugging information into a logger.
     *
     * @param logger where to dump the information.
     * @throws NullPointerException if {@code logger} is {@code null}.
     */
    public void debugDump(final Logger logger) {
        final StringWriter sw = new StringWriter();
        final PrintWriter pw = new PrintWriter(sw);
        debugDump(pw);
        final String str = sw.toString();
        final String[] lines = str.split("\n");
        for (final String line : lines) {
            if (line.isEmpty()) {
                continue;
            }
            logger.debug(line);
        }
    }

    /**
     * Given whether a new haplotype that has been already added to {@link #haplotypes} collection is the
     * reference haplotype and updates {@link #refHaplotype} accordingly.
     * <p>
     * <p>
     * This method assumes that the colling code has verified that the haplotype was not already in {@link #haplotypes}
     * I.e. that it is really a new one. Otherwise it will result in an exception if it happen to be a reference
     * haplotype and this has already be set. This is the case even if the new haplotypes and the current reference
     * are equal.
     * </p>
     *
     * @param newHaplotype the new haplotype.
     * @throws NullPointerException  if {@code newHaplotype} is {@code null}.
     * @throws IllegalStateException if there is already a reference haplotype.
     */
    private void updateReferenceHaplotype(final Haplotype newHaplotype) {
        if (!newHaplotype.isReference()) return;
        if (refHaplotype == null)
            refHaplotype = newHaplotype;
        else // assumes that we have checked wether the haplotype is already in the collection and so is no need to check equality.
            throw new IllegalStateException("the assembly-result-set already have a reference haplotype that is different");
    }

    /**
     * Returns a sorted set of variant events that best explain the haplotypes found by the assembly
     * across kmerSizes.
     * <p>
     * <p/>
     * The result is sorted incrementally by location.
     *
     * @return never {@code null}, but perhaps an empty collection.
     */
    public TreeSet<VariantContext> getVariationEvents() {
        if (variationEvents == null) {
            final List<Haplotype> haplotypeList = getHaplotypeList();
            EventMap.buildEventMapsForHaplotypes(haplotypeList, fullReferenceWithPadding, paddedReferenceLoc, debug);
            variationEvents = EventMap.getAllVariantContexts(haplotypeList);
        }
        return variationEvents;
    }
}
