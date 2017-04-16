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
package org.ncic.bioinfo.sparkseq.exceptions;

import htsjdk.samtools.SAMRecord;
import org.ncic.bioinfo.sparkseq.algorithms.utils.DocumentedGATKFeature;
import org.ncic.bioinfo.sparkseq.algorithms.utils.GenomeLoc;
import org.ncic.bioinfo.sparkseq.algorithms.utils.HelpConstants;

import java.io.File;

/**
 * Author: wbc
 */
@DocumentedGATKFeature(
        groupName = HelpConstants.DOCS_CAT_USRERR,
        summary = "Errors caused by incorrect user behavior, such as bad files, bad arguments, etc.")
public class UserException extends ReviewedGATKException {

    public UserException(String msg) {
        super(msg);
    }

    public UserException(String msg, Throwable e) {
        super(msg, e);
    }

    private UserException(Throwable e) {
        super("", e);
    } // cannot be called, private access

    protected static String getMessage(Throwable t) {
        String message = t.getMessage();
        return message != null ? message : t.getClass().getName();
    }

    public static class BadInput extends UserException {
        public BadInput(String message) {
            super(String.format("Bad input: %s", message));
        }
    }

    public static class CommandLineException extends UserException {
        public CommandLineException(String message) {
            super(String.format("Invalid command line: %s", message));
        }
    }

    public static class MalformedGenomeLoc extends UserException {
        public MalformedGenomeLoc(String message, GenomeLoc loc) {
            super(String.format("Badly formed genome loc: %s: %s", message, loc));
        }

        public MalformedGenomeLoc(String message) {
            super(String.format("Badly formed genome loc: %s", message));
        }
    }

    public static class ErrorWritingBamFile extends UserException {
        public ErrorWritingBamFile(String message) {
            super(String.format("An error occurred when trying to write the BAM file.  Usually this happens when there is not enough space in the directory to which the data is being written (generally the temp directory) or when your system's open file handle limit is too small.  To tell Java to use a bigger/better file system use -Djava.io.tmpdir=X on the command line.  The exact error was %s", message));
        }
    }

    public static class MalformedBAM extends UserException {
        public MalformedBAM(SAMRecord read, String message) {
            this(read.getFileSource() != null ? read.getFileSource().getReader().toString() : "(none)", message);
        }

        public MalformedBAM(File file, String message) {
            this(file.toString(), message);
        }

        public MalformedBAM(String source, String message) {
            super(String.format("SAM/BAM file %s is malformed: %s", source, message));
        }
    }

    public static class MisencodedBAM extends UserException {
        public MisencodedBAM(SAMRecord read, String message) {
            this(read.getFileSource() != null ? read.getFileSource().getReader().toString() : "(none)", message);
        }

        public MisencodedBAM(String source, String message) {
            super(String.format("SAM/BAM file %s appears to be using the wrong encoding for quality scores: %s; please see the GATK --help documentation for options related to this error", source, message));
        }
    }

    public static class BadArgumentValue extends CommandLineException {
        public BadArgumentValue(String arg, String message) {
            super(String.format("Argument %s has a bad value: %s", arg, message));
        }
    }

    public static class CouldNotReadInputFile extends UserException {
        public CouldNotReadInputFile(String message, Exception e) {
            super(String.format("Couldn't read file because %s caused by %s", message, getMessage(e)));
        }

        public CouldNotReadInputFile(File file) {
            super(String.format("Couldn't read file %s", file.getAbsolutePath()));
        }

        public CouldNotReadInputFile(File file, String message) {
            super(String.format("Couldn't read file %s because %s", file.getAbsolutePath(), message));
        }

        public CouldNotReadInputFile(String file, String message) {
            super(String.format("Couldn't read file %s because %s", file, message));
        }

        public CouldNotReadInputFile(File file, String message, Exception e) {
            super(String.format("Couldn't read file %s because %s with exception %s", file.getAbsolutePath(), message, getMessage(e)));
        }

        public CouldNotReadInputFile(File file, Exception e) {
            this(file, getMessage(e));
        }

        public CouldNotReadInputFile(String message) {
            super(message);
        }
    }

    public static class IncompatibleReadFiltersException extends CommandLineException {
        public IncompatibleReadFiltersException(final String filter1, final String filter2) {
            super(String.format("Two read filters are enabled that are incompatible and cannot be used simultaneously: %s and %s", filter1, filter2));
        }
    }

    public static class ReadMissingReadGroup extends MalformedBAM {
        public ReadMissingReadGroup(final SAMRecord read) {
            super(read, String.format("Read %s is missing the read group (RG) tag, which is required by the GATK.  Please use " + HelpConstants.forumPost("discussion/59/companion-utilities-replacereadgroups to fix this problem"), read.getReadName()));
        }
    }

    public static class MalformedFile extends UserException {
        public MalformedFile(String message) {
            super(String.format("Unknown file is malformed: %s", message));
        }

        public MalformedFile(String message, Exception e) {
            super(String.format("Unknown file is malformed: %s caused by %s", message, getMessage(e)));
        }

        public MalformedFile(File f, String message) {
            super(String.format("File %s is malformed: %s", f.getAbsolutePath(), message));
        }

        public MalformedFile(File f, String message, Exception e) {
            super(String.format("File %s is malformed: %s caused by %s", f.getAbsolutePath(), message, getMessage(e)));
        }

        public MalformedFile(String name, String message) {
            super(String.format("File associated with name %s is malformed: %s", name, message));
        }

        public MalformedFile(String name, String message, Exception e) {
            super(String.format("File associated with name %s is malformed: %s caused by %s", name, message, getMessage(e)));
        }
    }

    public static class DeprecatedAnnotation extends UserException {
        public DeprecatedAnnotation(String annotationName, String version) {
            super(String.format("Annotation %s is no longer available in the GATK; it has been deprecated since version %s", annotationName, version));
        }
    }
}
