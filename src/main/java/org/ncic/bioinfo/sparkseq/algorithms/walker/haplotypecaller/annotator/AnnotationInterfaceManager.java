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
package org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator;

import org.ncic.bioinfo.sparkseq.algorithms.utils.DeprecatedToolChecks;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.GenotypeAnnotation;
import org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.interfaces.InfoFieldAnnotation;
import org.ncic.bioinfo.sparkseq.exceptions.GATKException;
import org.ncic.bioinfo.sparkseq.exceptions.UserException;
import sun.reflect.annotation.AnnotationType;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;

/**
 * Author: wbc
 */
public class AnnotationInterfaceManager {
    private static final String NULL_ANNOTATION_NAME = "none";
    private static final String NULL_ANNOTATION_GROUP_NAME = "none";
    private static final String packageName = "org.ncic.bioinfo.sparkseq.algorithms.walker.haplotypecaller.annotator.";

    public static List<InfoFieldAnnotation> createInfoFieldAnnotations(List<String> infoFieldAnnotations) {
        List<InfoFieldAnnotation> annotations = new ArrayList<>();
        try {
            for (String name : infoFieldAnnotations) {
                annotations.add((InfoFieldAnnotation) Class.forName(packageName + name).newInstance());
            }
        } catch (Exception e) {
            throw new GATKException("Can't find annotation");
        }
        return annotations;
    }

    public static List<GenotypeAnnotation> createGenotypeAnnotations(List<String> genotypeAnnotations) {
        List<GenotypeAnnotation> annotations = new ArrayList<>();
        try {
            for (String name : genotypeAnnotations) {
                annotations.add((GenotypeAnnotation) Class.forName(packageName + name).newInstance());
            }
        } catch (Exception e) {
            throw new GATKException("Can't find annotation");
        }
        annotations.add(new DepthPerAlleleBySample());
        annotations.add(new DepthPerSampleHC());
        annotations.add(new StrandBiasBySample());
        return annotations;
    }
}
