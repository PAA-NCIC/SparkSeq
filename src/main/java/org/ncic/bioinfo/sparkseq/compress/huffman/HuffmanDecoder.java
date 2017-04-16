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
import java.util.Objects;

/**
 * Reads from a Huffman-coded bit stream and decodes symbols. Not thread-safe.
 * @see HuffmanEncoder
 */
public final class HuffmanDecoder {
	
	// The underlying bit input stream (not null).
	private BitInputStream input;
	
	/**
	 * The code tree to use in the next {@link#read()} operation. Must be given a non-{@code null}
	 * value before calling read(). The tree can be changed after each symbol decoded, as long
	 * as the encoder and decoder have the same tree at the same point in the code stream.
	 */
	public CodeTree codeTree;
	
	
	
	/**
	 * Constructs a Huffman decoder based on the specified bit input stream.
	 * @param in the bit input stream to read from
	 * @throws NullPointerException if the input stream is {@code null}
	 */
	public HuffmanDecoder(BitInputStream in) {
		Objects.requireNonNull(in);
		input = in;
	}
	
	
	
	/**
	 * Reads from the input stream to decode the next Huffman-coded symbol.
	 * @return the next symbol in the stream, which is non-negative
	 * @throws IOException if an I/O exception occurred
	 * @throws NullPointerException if the current code tree is {@code null}
	 */
	public int read() throws IOException {
		if (codeTree == null)
			throw new NullPointerException("Code tree is null");
		
		InternalNode currentNode = codeTree.root;
		while (true) {
			int temp = input.readNoEof();
			Node nextNode;
			if      (temp == 0) nextNode = currentNode.leftChild;
			else if (temp == 1) nextNode = currentNode.rightChild;
			else throw new AssertionError("Invalid value from readNoEof()");
			
			if (nextNode instanceof Leaf)
				return ((Leaf)nextNode).symbol;
			else if (nextNode instanceof InternalNode)
				currentNode = (InternalNode)nextNode;
			else
				throw new AssertionError("Illegal node type");
		}
	}
	
}
