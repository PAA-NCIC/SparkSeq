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

import org.ncic.bioinfo.sparkseq.exceptions.PipelineException;

/**
 * Author: wbc
 */
public class BaseCompressTools {
    /**
     * 目前只支持0-255的读长，要支持更长的读长，需要增加前面length的长度
     * qual的取值范围是33-126，因此当base中出现非AGCT时，将qual对应位置替换为127
     * 这样，将base按照二位进行压缩
     * 在最前面记录一个字节长度的base个数
     */
    public static byte[] compressBase(byte[] base, byte[] qual) {
        // 如果是一个*，说明是空的base，不用压缩
        if(base.length == 1 && base[0] == '*') {
            return base;
        }

        int rawLength = base.length;
        int compressedLength = (base.length + 3) / 4;
        if (rawLength > 255) {
            throw new PipelineException("Can't support read longer than 255");
        }
        byte[] result = new byte[compressedLength + 1];   //多一个字节用于存rawLength
        result[0] = (byte) (rawLength - 128);

        int resultIdx = 1;
        byte tmp = 0;
        int flag = 0;
        for (int i = 0; i < rawLength; i++) {
            byte b = base[i];
            switch (b) {
                case 'A':
                case 'a':
                    break;
                case 'G':
                case 'g':
                    tmp += 1;
                    break;
                case 'C':
                case 'c':
                    tmp += 2;
                    break;
                case 'T':
                case 't':
                    tmp += 3;
                    break;
                default:
                    qual[i] = 127;
            }
            flag++;
            if (flag == 4 || i == rawLength - 1) {
                result[resultIdx] = tmp;
                resultIdx++;
                flag = 0;
                tmp = 0;
            } else {
                tmp <<= 2;
            }
        }
        return result;
    }

    public static byte[] decompressBase(byte[] compressedBases, byte[] qual) {
        // 如果是一个*，说明是空的base，不用解压缩
        if(compressedBases.length == 1 && compressedBases[0] == '*') {
            return compressedBases;
        }

        int rawLength = compressedBases[0] + 128;
        byte[] result = new byte[rawLength];
        int compressIdx = 1;
        int baseLeftInTmp = 0;
        byte tmp = 0;
        for (int i = 0; i < rawLength; i++) {
            if (baseLeftInTmp == 0) {
                tmp = compressedBases[compressIdx];
                baseLeftInTmp = Math.min(rawLength - i, 4);
                compressIdx++;
            }
            byte b = (byte) (tmp >>> (baseLeftInTmp * 2 - 2) & 3);
            setBase(b, qual, result, i);
            baseLeftInTmp--;
        }
        return result;
    }

    private static void setBase(byte baseVal, byte[] qual, byte[] result, int idx) {
        switch (baseVal) {
            case 0:
                if(qual[idx] == 127) {
                    result[idx] = 'N';
                    qual[idx] = 33;
                } else {
                    result[idx] = 'A';
                }
                break;
            case 1:
                result[idx] = 'G';
                break;
            case 2:
                result[idx] = 'C';
                break;
            case 3:
                result[idx] = 'T';
                break;
        }
    }
}
