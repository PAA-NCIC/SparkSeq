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
package org.ncic.bioinfo.sparkseq.algorithms.walker.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.DiploidGenotype;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Author: wbc
 */
public class CandidateMutation {
    private GenomeLoc location;
    private String sequenceContext;
    private char refAllele;
    private boolean dbsnpSite = false;
    private boolean cosmicSite = false;
    private VariantContext panelOfNormalsVC;
    private VariantContext dbsnpVC;
    private boolean covered = false;

    private double power;
    private double tumorPower;
    private double normalPower;
    private double normalPowerWithSNPPrior;
    private double normalPowerNoSNPPrior;

    private char altAllele = 'N';
    private String tumorSampleName = "TUMOR";
    private String normalSampleName = "NORMAL";

    private double contaminationFraction;

    private double contaminantLod;
    private int score;

    private int tumorQ20Count;
    private int normalQ20Count;

    private int totalReads;
    private int mapQ0Reads;
    private int initialTumorRefCounts;
    private int initialTumorAltCounts;
    private int initialTumorRefQualitySum;
    private int initialTumorAltQualitySum;
    private int initialTumorNonRefQualitySum;
    private int initialTumorReadDepth;
    private int initialNormalRefCounts;
    private int initialNormalAltCounts;
    private int initialNormalRefQualitySum;
    private int initialNormalAltQualitySum;
    private int tumorRefMaxMapQ;
    private int tumorAltMaxMapQ;
    private int initialNormalReadDepth;
    private DiploidGenotype initialNormalBestGenotype;

    private double initialTumorLod;
    private double initialNormalLod;

    private double tumorF;
    private double tumorLodFStar;
    private double tumorLodFStarForward;
    private double tumorLodFStarReverse;

    private double normalF;

    private double powerToDetectPositiveStrandArtifact;
    private double powerToDetectNegativeStrandArtifact;

    private int[] strandContingencyTable;

    private List<Integer> tumorAltForwardOffsetsInRead;
    private List<Integer> tumorAltReverseOffsetsInRead;

    private Double tumorForwardOffsetsInReadMedian;
    private Double tumorForwardOffsetsInReadMad;
    private Double tumorReverseOffsetsInReadMedian;
    private Double tumorReverseOffsetsInReadMad;

    private int tumorInsertionCount;
    private int tumorDeletionCount;

    // quality score info
    private List<Integer> tumorRefQualityScores;
    private List<Integer> tumorAltQualityScores;
    private List<Integer> normalRefQualityScores;
    private List<Integer> normalAltQualityScores;

    private List<String> rejectionReasons = new ArrayList<String>();
    private boolean rejected = false; // summary judgement... keep or reject the site

    public CandidateMutation(GenomeLoc location, char refAllele) {
        this.location = location;
        this.refAllele = refAllele;
    }

    public int getScore() {
        return score;
    }

    public boolean isGermlineAtRisk() {
        return (dbsnpSite && !cosmicSite);
    }

    public int getCountOfNormalsObservedIn() {
        if (panelOfNormalsVC == null) { return 0; }

        // it's either an integer, or a collection of integers
        int count = 0;
        Object o = panelOfNormalsVC.getAttribute("AC");
        if (o == null) { return 0; }
        if (o instanceof String) { return Integer.parseInt((String) o); }
        if (o instanceof Collection) {
            for(String s : ((Collection<String>) o)) {
                count += Integer.parseInt(s);
            }
            return count;
        }
        throw new RuntimeException("Unexpected value processing panel of normals allele count: " + o);
    }

    public void setRejectionReasons(List<String> rejectionReasons) {
        if (rejectionReasons != null && rejectionReasons.size() > 0) {
            setRejected(true);
        }
        this.rejectionReasons = rejectionReasons;
    }

    public void addRejectionReason(String reason) {
        setRejected(true);
        getRejectionReasons().add(reason);
    }




    // -------------------------------------------------------------------------
    // GENERATED CODE BELOW THIS POINT
    // -------------------------------------------------------------------------

    public int[] getStrandContingencyTable() {
        return strandContingencyTable;
    }

    public void setStrandContingencyTable(int[] strandContingencyTable) {
        this.strandContingencyTable = strandContingencyTable;
    }

    public GenomeLoc getLocation() {
        return location;
    }

    public boolean isDbsnpSite() {
        return dbsnpSite;
    }

    public void setDbsnpSite(boolean dbsnpSite) {
        this.dbsnpSite = dbsnpSite;
    }

    public VariantContext getDbsnpVC() {
        return dbsnpVC;
    }

    public void setDbsnpVC(VariantContext dbsnpVC) {
        this.dbsnpVC = dbsnpVC;
    }

    public boolean isCovered() {
        return covered;
    }

    public void setCovered(boolean covered) {
        this.covered = covered;
    }

    public String getSequenceContext() {
        return sequenceContext;
    }

    public void setSequenceContext(String sequenceContext) {
        this.sequenceContext = sequenceContext;
    }

    public char getRefAllele() {
        return refAllele;
    }

    public char getAltAllele() {
        return altAllele;
    }

    public void setAltAllele(char altAllele) {
        this.altAllele = altAllele;
    }

    public boolean isRejected() {
        return rejected;
    }

    public void setRejected(boolean rejected) {
        this.rejected = rejected;
    }

    public double getInitialTumorLod() {
        return initialTumorLod;
    }

    public void setInitialTumorLod(double initialTumorLod) {
        this.initialTumorLod = initialTumorLod;
    }

    public double getInitialNormalLod() {
        return initialNormalLod;
    }

    public void setInitialNormalLod(double initialNormalLod) {
        this.initialNormalLod = initialNormalLod;
    }

    public double getTumorLodFStar() {
        return tumorLodFStar;
    }

    public void setTumorLodFStar(double tumorLodFStar) {
        this.tumorLodFStar = tumorLodFStar;
    }


    public double getTumorLodFStarForward() {
        return tumorLodFStarForward;
    }

    public void setTumorLodFStarForward(double tumorLodFStarForward) {
        this.tumorLodFStarForward = tumorLodFStarForward;
    }

    public double getTumorLodFStarReverse() {
        return tumorLodFStarReverse;
    }

    public void setTumorLodFStarReverse(double tumorLodFStarReverse) {
        this.tumorLodFStarReverse = tumorLodFStarReverse;
    }

    public double getTumorF() {
        return tumorF;
    }

    public void setTumorF(double tumorF) {
        this.tumorF = tumorF;
    }

    public double getNormalF() {
        return normalF;
    }

    public void setNormalF(double normalF) {
        this.normalF = normalF;
    }

    public int getInitialTumorRefQualitySum() {
        return initialTumorRefQualitySum;
    }

    public void setInitialTumorRefQualitySum(int initialTumorRefQualitySum) {
        this.initialTumorRefQualitySum = initialTumorRefQualitySum;
    }

    public int getInitialTumorAltQualitySum() {
        return initialTumorAltQualitySum;
    }

    public void setInitialTumorAltQualitySum(int initialTumorAltQualitySum) {
        this.initialTumorAltQualitySum = initialTumorAltQualitySum;
    }

    public int getInitialTumorNonRefQualitySum() {
        return initialTumorNonRefQualitySum;
    }

    public void setInitialTumorNonRefQualitySum(int initialTumorNonRefQualitySum) {
        this.initialTumorNonRefQualitySum = initialTumorNonRefQualitySum;
    }

    public int getInitialNormalRefQualitySum() {
        return initialNormalRefQualitySum;
    }

    public void setInitialNormalRefQualitySum(int initialNormalRefQualitySum) {
        this.initialNormalRefQualitySum = initialNormalRefQualitySum;
    }

    public int getInitialNormalAltQualitySum() {
        return initialNormalAltQualitySum;
    }

    public void setInitialNormalAltQualitySum(int initialNormalAltQualitySum) {
        this.initialNormalAltQualitySum = initialNormalAltQualitySum;
    }

    public DiploidGenotype getInitialNormalBestGenotype() {
        return initialNormalBestGenotype;
    }

    public void setInitialNormalBestGenotype(DiploidGenotype initialNormalBestGenotype) {
        this.initialNormalBestGenotype = initialNormalBestGenotype;
    }

    public int getInitialTumorReadDepth() {
        return initialTumorReadDepth;
    }

    public void setInitialTumorReadDepth(int initialTumorReadDepth) {
        this.initialTumorReadDepth = initialTumorReadDepth;
    }

    public int getInitialNormalReadDepth() {
        return initialNormalReadDepth;
    }

    public void setInitialNormalReadDepth(int initialNormalReadDepth) {
        this.initialNormalReadDepth = initialNormalReadDepth;
    }

    public String getTumorSampleName() {
        return tumorSampleName;
    }

    public void setTumorSampleName(String tumorSampleName) {
        this.tumorSampleName = tumorSampleName;
    }

    public String getNormalSampleName() {
        return normalSampleName;
    }

    public void setNormalSampleName(String normalSampleName) {
        this.normalSampleName = normalSampleName;
    }

    public List<String> getRejectionReasons() {
        return rejectionReasons;
    }

    public int getInitialTumorRefCounts() {
        return initialTumorRefCounts;
    }

    public void setInitialTumorRefCounts(int initialTumorRefCounts) {
        this.initialTumorRefCounts = initialTumorRefCounts;
    }

    public int getInitialTumorAltCounts() {
        return initialTumorAltCounts;
    }

    public void setInitialTumorAltCounts(int initialTumorAltCounts) {
        this.initialTumorAltCounts = initialTumorAltCounts;
    }

    public int getInitialNormalRefCounts() {
        return initialNormalRefCounts;
    }

    public void setInitialNormalRefCounts(int initialNormalRefCounts) {
        this.initialNormalRefCounts = initialNormalRefCounts;
    }

    public int getInitialNormalAltCounts() {
        return initialNormalAltCounts;
    }

    public void setInitialNormalAltCounts(int initialNormalAltCounts) {
        this.initialNormalAltCounts = initialNormalAltCounts;
    }

    public int getTumorQ20Count() {
        return tumorQ20Count;
    }

    public void setTumorQ20Count(int tumorQ20Count) {
        this.tumorQ20Count = tumorQ20Count;
    }

    public int getNormalQ20Count() {
        return normalQ20Count;
    }

    public void setNormalQ20Count(int normalQ20Count) {
        this.normalQ20Count = normalQ20Count;
    }

    public int getTumorInsertionCount() {
        return tumorInsertionCount;
    }

    public void setTumorInsertionCount(int tumorInsertionCount) {
        this.tumorInsertionCount = tumorInsertionCount;
    }

    public int getTumorDeletionCount() {
        return tumorDeletionCount;
    }

    public void setTumorDeletionCount(int tumorDeletionCount) {
        this.tumorDeletionCount = tumorDeletionCount;
    }

    public int getTotalReads() {
        return totalReads;
    }

    public void setTotalReads(int totalReads) {
        this.totalReads = totalReads;
    }

    public int getMapQ0Reads() {
        return mapQ0Reads;
    }

    public void setMapQ0Reads(int mapQ0Reads) {
        this.mapQ0Reads = mapQ0Reads;
    }

    public double getContaminationFraction() {
        return contaminationFraction;
    }

    public void setContaminationFraction(double contaminationFraction) {
        this.contaminationFraction = contaminationFraction;
    }

    public double getContaminantLod() {
        return contaminantLod;
    }

    public void setContaminantLod(double contaminantLod) {
        this.contaminantLod = contaminantLod;
    }

    public List<Integer> getTumorAltForwardOffsetsInRead() {
        return tumorAltForwardOffsetsInRead;
    }

    public void setTumorAltForwardOffsetsInRead(List<Integer> tumorAltForwardOffsetsInRead) {
        this.tumorAltForwardOffsetsInRead = tumorAltForwardOffsetsInRead;
    }

    public List<Integer> getTumorAltReverseOffsetsInRead() {
        return tumorAltReverseOffsetsInRead;
    }

    public void setTumorAltReverseOffsetsInRead(List<Integer> tumorAltReverseOffsetsInRead) {
        this.tumorAltReverseOffsetsInRead = tumorAltReverseOffsetsInRead;
    }

    public Double getTumorForwardOffsetsInReadMedian() {
        return tumorForwardOffsetsInReadMedian;
    }

    public void setTumorForwardOffsetsInReadMedian(Double tumorForwardOffsetsInReadMedian) {
        this.tumorForwardOffsetsInReadMedian = tumorForwardOffsetsInReadMedian;
    }

    public Double getTumorForwardOffsetsInReadMad() {
        return tumorForwardOffsetsInReadMad;
    }

    public void setTumorForwardOffsetsInReadMad(Double tumorForwardOffsetsInReadMad) {
        this.tumorForwardOffsetsInReadMad = tumorForwardOffsetsInReadMad;
    }

    public Double getTumorReverseOffsetsInReadMedian() {
        return tumorReverseOffsetsInReadMedian;
    }

    public void setTumorReverseOffsetsInReadMedian(Double tumorReverseOffsetsInReadMedian) {
        this.tumorReverseOffsetsInReadMedian = tumorReverseOffsetsInReadMedian;
    }

    public Double getTumorReverseOffsetsInReadMad() {
        return tumorReverseOffsetsInReadMad;
    }

    public void setTumorReverseOffsetsInReadMad(Double tumorReverseOffsetsInReadMad) {
        this.tumorReverseOffsetsInReadMad = tumorReverseOffsetsInReadMad;
    }

    public double getPower() {
        return power;
    }

    public void setPower(double power) {
        this.power = power;
    }

    public double getTumorPower() {
        return tumorPower;
    }

    public void setTumorPower(double tumorPower) {
        this.tumorPower = tumorPower;
    }

    public double getNormalPower() {
        return normalPower;
    }

    public void setNormalPower(double normalPower) {
        this.normalPower = normalPower;
    }

    public boolean isCosmicSite() {
        return cosmicSite;
    }

    public void setCosmicSite(boolean cosmicSite) {
        this.cosmicSite = cosmicSite;
    }

    public boolean isSeenInPanelOfNormals() {
        return (panelOfNormalsVC != null);
    }

    public VariantContext getPanelOfNormalsVC() {
        return panelOfNormalsVC;
    }

    public void setPanelOfNormalsVC(VariantContext panelOfNormalsVC) {
        this.panelOfNormalsVC = panelOfNormalsVC;
    }

    public double getPowerToDetectPositiveStrandArtifact() {
        return powerToDetectPositiveStrandArtifact;
    }

    public void setPowerToDetectPositiveStrandArtifact(double powerToDetectPositiveStrandArtifact) {
        this.powerToDetectPositiveStrandArtifact = powerToDetectPositiveStrandArtifact;
    }

    public double getPowerToDetectNegativeStrandArtifact() {
        return powerToDetectNegativeStrandArtifact;
    }

    public void setPowerToDetectNegativeStrandArtifact(double powerToDetectNegativeStrandArtifact) {
        this.powerToDetectNegativeStrandArtifact = powerToDetectNegativeStrandArtifact;
    }

    public double getNormalPowerWithSNPPrior() {
        return normalPowerWithSNPPrior;
    }

    public void setNormalPowerWithSNPPrior(double normalPowerWithSNPPrior) {
        this.normalPowerWithSNPPrior = normalPowerWithSNPPrior;
    }

    public double getNormalPowerNoSNPPrior() {
        return normalPowerNoSNPPrior;
    }

    public void setNormalPowerNoSNPPrior(double normalPowerNoSNPPrior) {
        this.normalPowerNoSNPPrior = normalPowerNoSNPPrior;
    }

    public int getTumorAltMaxMapQ() {
        return tumorAltMaxMapQ;
    }

    public void setTumorAltMaxMapQ(int tumorAltMaxMapQ) {
        this.tumorAltMaxMapQ = tumorAltMaxMapQ;
    }

    public int getTumorRefMaxMapQ() {
        return tumorRefMaxMapQ;
    }

    public void setTumorRefMaxMapQ(int tumorRefMaxMapQ) {
        this.tumorRefMaxMapQ = tumorRefMaxMapQ;
    }

    public List<Integer> getTumorRefQualityScores() {
        return tumorRefQualityScores;
    }

    public void setTumorRefQualityScores(List<Integer> tumorRefQualityScores) {
        this.tumorRefQualityScores = tumorRefQualityScores;
    }

    public List<Integer> getTumorAltQualityScores() {
        return tumorAltQualityScores;
    }

    public void setTumorAltQualityScores(List<Integer> tumorAltQualityScores) {
        this.tumorAltQualityScores = tumorAltQualityScores;
    }

    public List<Integer> getNormalRefQualityScores() {
        return normalRefQualityScores;
    }

    public void setNormalRefQualityScores(List<Integer> normalRefQualityScores) {
        this.normalRefQualityScores = normalRefQualityScores;
    }

    public List<Integer> getNormalAltQualityScores() {
        return normalAltQualityScores;
    }

    public void setNormalAltQualityScores(List<Integer> normalAltQualityScores) {
        this.normalAltQualityScores = normalAltQualityScores;
    }
}

