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

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
import com.esotericsoftware.kryo.serializers.JavaSerializer;
import junit.framework.TestCase;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.Serializable;

/**
 * Author: wbc
 */
public class TestKryo extends TestCase {

    public void testKyro() {
        Kryo kryo = new Kryo();
        kryo.setReferences(false);
        //byte[] test = new byte[10];
        //kryo.register(test.getClass(), new JavaSerializer());

        byte[] sequence = "TGGGATGAGAGCATGAGAAGGTGGAGCTAAGGTGGGAGACCGTCTACCCCCGACCCTGTGTGGTGCACTGACCGTGACTCTCTGCACCTTCTCGTGGGGGA".getBytes();
        byte[] quality = "DCDE8EFFFFBEEEEFGGGGEFEGGGFGGGFF6FF/>@A?EEEE=DCAC5C@@@B=8B64A########################################".getBytes();

        byte[] serialSequence = serializationObject(kryo, sequence);
        byte[] serialQual = serializationObject(kryo, quality);
        System.out.println(sequence.length + " " + serialSequence.length);
        System.out.println(quality.length + " " + serialQual.length);
    }

    private byte[] serializationObject(Kryo kryo, byte[] obj) {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        Output output = new Output(baos);
        kryo.writeClassAndObject(output, obj);
        output.flush();
        output.close();

        byte[] b = baos.toByteArray();
        try {
            baos.flush();
            baos.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return b;
    }
}
