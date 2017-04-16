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
import org.ncic.bioinfo.sparkseq.exceptions.CompressException;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Author: wbc
 */
public class QualityCompressTools {

    private static CodeTree code = getCode();

    private static final CodeTree getCode() {
        // 差值是-94到+94，因此需要189个元素
        int[] statistic = new int[257];
        //int[] data = {8, 133, 256, 600, 2121, 3653, 3880, 6160, 5320, 12500, 8653, 1596, 2225, 3380, 4168, 5658, 7257, 8217, 8857, 10968, 12772, 15556, 18880, 22895, 26990, 32692, 41494, 48039, 60265, 74272, 92364, 117290, 168864, 217019, 294182, 507508, 866129, 913998, 7873334, 915057, 866421, 507922, 294456, 216970, 168602, 117199, 92304, 74216, 60236, 47965, 41475, 32743, 27045, 22985, 18967, 15558, 12793, 10995, 8858, 8262, 7311, 5712, 4202, 3403, 2241, 1604, 9607, 14151, 5853, 6961, 4276, 4048, 2304, 617, 247, 117, 7};

        int[] data = {175, 3061, 2379, 3744, 2596, 3190, 3130, 2985, 4477, 6139, 4409, 679, 743, 1153, 1543, 2132, 2758, 3516, 6337, 8123, 8639, 9538, 10567, 10981, 12769, 16170, 18882, 25177, 27573, 36628, 50111, 60441, 81285, 93054, 167301, 298000, 413930, 867784, 13896141, 847358, 405231, 280493, 161785, 86856, 72891, 55915, 46484, 33046, 24817, 23020, 17385, 13772, 10836, 10746, 9064, 9063, 8846, 7934, 5475, 2034, 1555, 1159, 892, 670, 454, 364, 182, 141, 84, 59, 38, 45, 3, 16913, 5, 0, 46, 66, 17, 29, 52, 257, 45, 102, 40, 117, 71, 159, 233, 61, 312, 6939, 141, 76, 142, 95, 132, 202, 408, 663, 500, 512, 751, 955, 610, 821, 1433, 1662, 2635, 3692, 143471, 588};
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
        return code;
    }

    public static byte[] compressQual(byte[] qual) {
        if (qual.length < 1) {
            throw new CompressException("Can't compress empty quality string");
        }

        // 如果是一个*，说明是空的base，不用压缩
        if (qual.length == 1 && qual[0] == '*') {
            return qual;
        }

        byte[] deltaArray = new byte[qual.length];
        deltaArray[0] = qual[0];

        for (int i = 1; i < qual.length; i++) {
            deltaArray[i] = (byte) (qual[i] - qual[i - 1] + 94);  // 可能是正负94，修正为非负
        }

        InputStream in = new ByteArrayInputStream(deltaArray);
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        BitOutputStream bitOut = new BitOutputStream(out);
        try {
            HuffmanCompress.compress(code, in, bitOut);
            bitOut.close();
        } catch (IOException e) {
            throw new CompressException();
        }
        return out.toByteArray();
    }

    public static byte[] deCompressQual(byte[] compressedQual) {
        // 如果是一个*，说明是空的base，不用解压缩
        if (compressedQual.length == 1 && compressedQual[0] == '*') {
            return compressedQual;
        }

        InputStream in2 = new ByteArrayInputStream(compressedQual);
        BitInputStream bitIn = new BitInputStream(in2);
        ByteArrayOutputStream out2 = new ByteArrayOutputStream();

        try {
            HuffmanDecompress.decompress(code, bitIn, out2);
        } catch (IOException e) {
            throw new CompressException();
        }

        byte[] deltaArray = out2.toByteArray();
        byte[] qual = new byte[deltaArray.length];
        qual[0] = deltaArray[0];
        for (int i = 1; i < qual.length; i++) {
            qual[i] = (byte) (qual[i - 1] + deltaArray[i] - 94);
        }
        return qual;
    }
}
