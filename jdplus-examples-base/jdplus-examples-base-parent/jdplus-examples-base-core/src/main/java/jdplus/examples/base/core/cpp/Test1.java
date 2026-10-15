/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.examples.base.core.cpp;

import java.util.Random;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.SymmetricMatrix;

/**
 *
 * @author palat
 */
@lombok.experimental.UtilityClass
public class Test1 {

    public void main(String[] args) {
        FastMatrix X = FastMatrix.make(5000, 50);
        Random rnd = new Random();
        X.set((i, j) -> rnd.nextDouble());
        long start = System.currentTimeMillis();
        for (int j = 0; j < 1000; ++j) {
            FastMatrix M = FastMatrix.square(50);
            SymmetricMatrix.XtX(X, M);
            SymmetricMatrix.lcholesky(M.deepClone());
        }

        long end = System.currentTimeMillis();
        System.out.print(end - start);
    }
}
