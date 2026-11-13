/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.examples.base.core;

import java.util.Random;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;

/**
 *
 * @author Jean Palate
 */
@lombok.experimental.UtilityClass
public class DataBlocks {

    public void main(String[] args) {
        test2();
    }

    private void test1() {
        Random rnd = new Random();
        int n = 10003;
        DataBlock X = DataBlock.make(n);
        DataBlock Y = DataBlock.make(n);
        X.set(i -> rnd.nextDouble());
        Y.set(i -> rnd.nextDouble());
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < 100000; ++i) {
            Y.addAY(1.25, X);
            Y.copy(X);
            double z = X.dot(Y);
        }
        long t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
    }

    private void test2() {
        int m = 50, n = 100;
        DataBlock X = DataBlock.make(m);
        DataBlock Y = DataBlock.make(n);
        FastMatrix A = FastMatrix.make(m, n);
        X.set(i -> 1);
        Y.set(i -> i);
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < 1000000; ++i) {
            A.addXaYt(1.25, X, Y);
        }
        long t1 = System.currentTimeMillis();
        System.out.println(t1 - t0);
        //        System.out.println(A);
    }
}
