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
package org.ncic.bioinfo.sparkseq.compress.huffman;

import java.io.IOException;
import java.util.List;
import java.util.Objects;

/**
 * Encodes symbols and writes to a Huffman-coded bit stream. Not thread-safe.
 * @see HuffmanDecoder
 */
public final class HuffmanEncoder {
	
	// The underlying bit output stream (not null).
	private BitOutputStream output;
	
	/**
	 * The code tree to use in the next {@link#write(int)} operation. Must be given a non-{@code null}
	 * value before calling write(). The tree can be changed after each symbol encoded, as long
	 * as the encoder and decoder have the same tree at the same point in the code stream.
	 */
	public CodeTree codeTree;
	
	
	
	/**
	 * Constructs a Huffman encoder based on the specified bit output stream.
	 * @param out the bit output stream to write to
	 * @throws NullPointerException if the output stream is {@code null}
	 */
	public HuffmanEncoder(BitOutputStream out) {
		Objects.requireNonNull(out);
		output = out;
	}
	
	
	
	/**
	 * Encodes the specified symbol and writes to the Huffman-coded output stream.
	 * @param symbol the symbol to encode, which is non-negative and must be in the range of the code tree
	 * @throws IOException if an I/O exception occurred
	 * @throws NullPointerException if the current code tree is {@code null}
	 * @throws IllegalArgumentException if the symbol value is negative or has no binary code
	 */
	public void write(int symbol) throws IOException {
		if (codeTree == null)
			throw new NullPointerException("Code tree is null");
		List<Integer> bits = codeTree.getCode(symbol);
		for (int b : bits)
			output.write(b);
	}
	
}
