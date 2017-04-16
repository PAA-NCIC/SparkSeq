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

import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.StringUtil;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Author: wbc
 */
public class CallStatsGenerator {
    private static final String TAB = "\t";
    private final FormatUtil fmt = new FormatUtil();

    private List<String> headers;

    public CallStatsGenerator(boolean enableQscoreOutput) {
        // set up the output columns (ordered!)
        this.headers = setupOutputColumns(enableQscoreOutput);
    }

    private List<String> setupOutputColumns(boolean enableQscoreOutput) {
        List<String> l = new ArrayList<String>();

        l.add("contig");
        l.add("position");
        l.add("context");
        l.add("ref_allele");
        l.add("alt_allele");
        l.add("tumor_name");
        l.add("normal_name");
        l.add("score");
        l.add("dbsnp_site");
        l.add("covered");
        l.add("power");
        l.add("tumor_power");
        l.add("normal_power");
        l.add("normal_power_nsp");
        l.add("normal_power_wsp");
        l.add("total_reads");
        l.add("map_Q0_reads");
        l.add("init_t_lod");
        l.add("t_lod_fstar");
        l.add("t_lod_fstar_forward");
        l.add("t_lod_fstar_reverse");
        l.add("tumor_f");
        l.add("contaminant_fraction");
        l.add("contaminant_lod");
        l.add("t_q20_count");
        l.add("t_ref_count");
        l.add("t_alt_count");
        l.add("t_ref_sum");
        l.add("t_alt_sum");
        l.add("t_ref_max_mapq");
        l.add("t_alt_max_mapq");
        l.add("t_ins_count");
        l.add("t_del_count");
        l.add("normal_best_gt");
        l.add("init_n_lod");
        l.add("normal_f");
        l.add("n_q20_count");
        l.add("n_ref_count");
        l.add("n_alt_count");
        l.add("n_ref_sum");
        l.add("n_alt_sum");
        l.add("power_to_detect_positive_strand_artifact");
        l.add("power_to_detect_negative_strand_artifact");
        l.add("strand_bias_counts");
        l.add("tumor_alt_fpir_median");
        l.add("tumor_alt_fpir_mad");
        l.add("tumor_alt_rpir_median");
        l.add("tumor_alt_rpir_mad");
        l.add("observed_in_normals_count");

        if (enableQscoreOutput) {
            l.add("tumor_ref_qscores");
            l.add("tumor_alt_qscores");
            l.add("normal_ref_qscores");
            l.add("normal_alt_qscores");
        }

        l.add("failure_reasons");
        l.add("judgement");
        return l;
    }

    public String generateHeader() {
        return StringUtil.join(TAB, this.headers);
    }

    public String generateCallStats(CandidateMutation candidate) {
        HashMap<String, String> d = new HashMap<String, String>();

        String keepString = "REJECT";
        if (!candidate.isRejected()) {
            keepString = "KEEP";
        }

        String siteInfo =
                getSiteInfoString(candidate.isDbsnpSite(), candidate.isCosmicSite());

        String strandInfo =
                getStrandTableString(candidate.getStrandContingencyTable());

        d.put("contig", candidate.getLocation().getContig());
        d.put("position", format(candidate.getLocation().getStart()));
        d.put("context", candidate.getSequenceContext());
        d.put("ref_allele", ""+candidate.getRefAllele());
        d.put("alt_allele", ""+candidate.getAltAllele());
        d.put("tumor_name", candidate.getTumorSampleName());
        d.put("normal_name", candidate.getNormalSampleName());
        d.put("score", format(candidate.getScore()));
        d.put("dbsnp_site", siteInfo);
        d.put("covered", (candidate.isCovered()?"COVERED":"UNCOVERED"));
        d.put("power", format(candidate.getPower()));
        d.put("tumor_power", format(candidate.getTumorPower()));
        d.put("normal_power", format(candidate.getNormalPower()));
        d.put("normal_power_nsp", format(candidate.getNormalPowerNoSNPPrior()));
        d.put("normal_power_wsp", format(candidate.getNormalPowerWithSNPPrior()));
        d.put("total_reads", format(candidate.getTotalReads()));
        d.put("map_Q0_reads", format(candidate.getMapQ0Reads()));
        d.put("init_t_lod", format(candidate.getInitialTumorLod()));
        d.put("t_lod_fstar", format(candidate.getTumorLodFStar()));
        d.put("t_lod_fstar_forward", format(candidate.getTumorLodFStarForward()));
        d.put("t_lod_fstar_reverse", format(candidate.getTumorLodFStarReverse()));
        d.put("tumor_f", format(candidate.getTumorF()));
        d.put("contaminant_fraction", format(candidate.getContaminationFraction()));
        d.put("contaminant_lod", format(candidate.getContaminantLod()));
        d.put("t_q20_count", format(candidate.getTumorQ20Count()));
        d.put("t_ref_count", format(candidate.getInitialTumorRefCounts()));
        d.put("t_alt_count", format(candidate.getInitialTumorAltCounts()));
        d.put("t_ref_sum", format(candidate.getInitialTumorRefQualitySum()));
        d.put("t_alt_sum", format(candidate.getInitialTumorAltQualitySum()));
        d.put("t_ref_max_mapq", format(candidate.getTumorRefMaxMapQ()));
        d.put("t_alt_max_mapq", format(candidate.getTumorAltMaxMapQ()));
        d.put("t_ins_count", format(candidate.getTumorInsertionCount()));
        d.put("t_del_count", format(candidate.getTumorDeletionCount()));
        d.put("normal_best_gt", format(candidate.getInitialNormalBestGenotype().toString()));
        d.put("init_n_lod", format(candidate.getInitialNormalLod()));
        d.put("normal_f", format(candidate.getNormalF()));
        d.put("n_q20_count", format(candidate.getNormalQ20Count()));
        d.put("n_ref_count", format(candidate.getInitialNormalRefCounts()));
        d.put("n_alt_count", format(candidate.getInitialNormalAltCounts()));
        d.put("n_ref_sum", format(candidate.getInitialNormalRefQualitySum()));
        d.put("n_alt_sum", format(candidate.getInitialNormalAltQualitySum()));
        d.put("power_to_detect_positive_strand_artifact", format(candidate.getPowerToDetectPositiveStrandArtifact()));
        d.put("power_to_detect_negative_strand_artifact", format(candidate.getPowerToDetectNegativeStrandArtifact()));
        d.put("strand_bias_counts", format(strandInfo));
        d.put("tumor_alt_fpir_median", candidate.getTumorForwardOffsetsInReadMedian()==null?"n/a":format(candidate.getTumorForwardOffsetsInReadMedian()));
        d.put("tumor_alt_fpir_mad", candidate.getTumorForwardOffsetsInReadMad()==null?"n/a":format(candidate.getTumorForwardOffsetsInReadMad()));
        d.put("tumor_alt_rpir_median", candidate.getTumorReverseOffsetsInReadMedian()==null?"n/a":format(candidate.getTumorReverseOffsetsInReadMedian()));
        d.put("tumor_alt_rpir_mad", candidate.getTumorReverseOffsetsInReadMad()==null?"n/a":format(candidate.getTumorReverseOffsetsInReadMad()));
        d.put("observed_in_normals_count", format(candidate.getCountOfNormalsObservedIn()));
        d.put("tumor_ref_qscores", format(candidate.getTumorRefQualityScores()));
        d.put("tumor_alt_qscores", format(candidate.getTumorAltQualityScores()));
        d.put("normal_ref_qscores", format(candidate.getNormalRefQualityScores()));
        d.put("normal_alt_qscores", format(candidate.getNormalAltQualityScores()));
        d.put("failure_reasons", StringUtil.join(",", candidate.getRejectionReasons().toArray(new String[candidate.getRejectionReasons().size()])));
        d.put("judgement", keepString);

        return generate(d);
    }

    private String generate(HashMap<String, String> d) {
        String[] msg = new String[headers.size()];

        for(int i=0; i<msg.length; i++) {
            String value = d.get(headers.get(i));
            if (value == null) { value = ""; }
            msg[i] = value;
        }
        return StringUtil.join(TAB, msg);
    }

    private String format(String s) { return s; }
    private String format(Integer i) { return fmt.format(i); }
    private String format(Double d) {
        if (d == null) { return "n/a"; }

        String s = fmt.format(d);
        return ("-0".equals(s))?"0":s;
    }

    private String format(List<Integer> ints) {
        if (ints == null || ints.size() == 0) {
            return "n/a";
        }

        return StringUtil.join(",", ints);
    }

    private String getSiteInfoString(boolean isDbsnpSite, boolean isCosmicSite) {
        String siteInfo = "NOVEL";
        if (isDbsnpSite) {
            siteInfo = "DBSNP";
        }
        if (isCosmicSite) {
            siteInfo = "COSMIC";
        }
        if (isDbsnpSite && isCosmicSite) {
            siteInfo = "DBSNP+COSMIC";
        }
        return siteInfo;
    }

    private String getStrandTableString(int[] ci) {
        StringBuilder sb = new StringBuilder();
        sb.append("(");
        sb.append(ci[0]).append(",");
        sb.append(ci[1]).append(",");
        sb.append(ci[2]).append(",");
        sb.append(ci[3]).append(")");
        return sb.toString();
    }

}
