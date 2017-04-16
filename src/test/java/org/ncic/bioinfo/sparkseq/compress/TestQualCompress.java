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
package org.ncic.bioinfo.sparkseq.compress;

import org.ncic.bioinfo.sparkseq.compress.huffman.BitInputStream;
import org.ncic.bioinfo.sparkseq.compress.huffman.BitOutputStream;
import org.ncic.bioinfo.sparkseq.compress.huffman.CanonicalCode;
import org.ncic.bioinfo.sparkseq.compress.huffman.CodeTree;
import org.ncic.bioinfo.sparkseq.compress.huffman.FrequencyTable;
import org.ncic.bioinfo.sparkseq.compress.huffman.HuffmanCompress;
import org.ncic.bioinfo.sparkseq.compress.huffman.HuffmanDecompress;
import org.ncic.bioinfo.sparkseq.algorithms.walker.AbstractTestCase;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Author: wbc
 */
public class TestQualCompress extends AbstractTestCase {

    public void testBuild() throws IOException {

        // 差值是-94到+94，因此需要189个元素
        int[] statistic = new int[257];
        int[] data = {8, 133, 256, 600, 2121, 3653, 3880, 6160, 5320, 12500, 8653, 1596, 2225, 3380, 4168, 5658, 7257, 8217, 8857, 10968, 12772, 15556, 18880, 22895, 26990, 32692, 41494, 48039, 60265, 74272, 92364, 117290, 168864, 217019, 294182, 507508, 866129, 913998, 7873334, 915057, 866421, 507922, 294456, 216970, 168602, 117199, 92304, 74216, 60236, 47965, 41475, 32743, 27045, 22985, 18967, 15558, 12793, 10995, 8858, 8262, 7311, 5712, 4202, 3403, 2241, 1604, 9607, 14151, 5853, 6961, 4276, 4048, 2304, 617, 247, 117, 7};

        for (int i = 0; i < data.length; i++) {
            statistic[i + 56] = data[i];
        }
        // 使用平滑
        for (int i = 0; i < statistic.length; i++) {
            statistic[i]++;
        }
        FrequencyTable freqs = new FrequencyTable(statistic);
        freqs.increment(256);
        CodeTree code = freqs.buildCodeTree();
        CanonicalCode canonCode = new CanonicalCode(code, 257);
        code = canonCode.toCodeTree();

        byte[] b = {94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94};
        InputStream in = new ByteArrayInputStream(b);
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        BitOutputStream bitOut = new BitOutputStream(out);
        HuffmanCompress.compress(code, in, bitOut);
        bitOut.close();
        byte[] compressed = out.toByteArray();

        InputStream in2 = new ByteArrayInputStream(compressed);
        BitInputStream bitIn = new BitInputStream(in2);
        ByteArrayOutputStream out2 = new ByteArrayOutputStream();
        HuffmanDecompress.decompress(code, bitIn, out2);
        byte[] decompressed = out2.toByteArray();

        assertEquals(new String(b), new String(decompressed));
    }

    public void testQualityCompress() {
        byte[] qual = {63, 63, 63, 63, 64, 64, 64, 64, 62, 35, 36, 71, 73, 73};
        byte[] compressed = QualityCompressTools.compressQual(qual);
        byte[] decompressed = QualityCompressTools.deCompressQual(compressed);
        assertEquals(new String(qual), new String(decompressed));
    }

}
