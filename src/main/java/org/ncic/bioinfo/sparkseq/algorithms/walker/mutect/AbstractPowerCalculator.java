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

import java.util.HashMap;

/**
 * Author: wbc
 */
public class AbstractPowerCalculator {
    protected HashMap<PowerCacheKey, Double> cache = new HashMap<PowerCacheKey, Double>();
    protected double constantEps;
    protected double constantLodThreshold;

    protected static class PowerCacheKey {
        private int n;
        private double delta;

        public PowerCacheKey(int n, double delta) {
            this.n = n;
            this.delta = delta;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            PowerCacheKey that = (PowerCacheKey) o;

            if (Double.compare(that.delta, delta) != 0) return false;
            if (n != that.n) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result;
            long temp;
            result = n;
            temp = delta != +0.0d ? Double.doubleToLongBits(delta) : 0L;
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            return result;
        }
    }

    protected static double calculateLogLikelihood(int depth, int alts, double eps, double f) {
        double a = (depth-alts) * Math.log10(f*eps + (1d-f)*(1d-eps));
        double b = (alts) * Math.log10(f*(1d-eps) + (1d-f)*eps);
        return (a+b);
    }

}
