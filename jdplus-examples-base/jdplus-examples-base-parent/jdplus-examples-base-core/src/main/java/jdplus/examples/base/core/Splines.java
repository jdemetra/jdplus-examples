/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.examples.base.core;

import java.io.IOException;
import java.io.InputStream;
import jdplus.sts.base.core.msts.AtomicModels;
import jdplus.sts.base.core.msts.CompositeModel;
import jdplus.sts.base.core.msts.CompositeModelEstimation;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.math.functions.Optimizer;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.ssf.SsfInitialization;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.MatrixFactory;
import tck.demetra.data.MatrixSerializer;

/**
 *
 * @author Jean Palate
 */
@lombok.experimental.UtilityClass
public class Splines {

    public void main(String[] args) throws IOException {

        InputStream stream = Splines.class.getResourceAsStream("marees.csv");
        Matrix marees = MatrixSerializer.read(stream);
        DoubleSeq data = marees.column(0).extract(1, 10000);
        //        RegularSpline rs = RegularSpline.of(12.46, 6);
        //        RegularSpline rs2 = RegularSpline.of(354.22, 20);
        //        SplineData sd = new SplineData(rs);
        //        SplineData sd2 = new SplineData(rs2);

        CompositeModel model = new CompositeModel();
        model.add(AtomicModels.localLevel("l", 1, true, Double.NaN));
        model.add(AtomicModels.regularSplines("s1", 12.42 * 10, 80, 0, 1, false));
        //        model.add(AtomicModels.regularSplines("s2", 354.22, 5, 0, 1, true));
        model.add(AtomicModels.noise("n", 1, false));
        int len = data.length();
        FastMatrix M = FastMatrix.make(len, 1);
        M.column(0).copy(data);
        CompositeModelEstimation rslt = model.estimate(
                M, false, true, SsfInitialization.Augmented_Robust, Optimizer.LevenbergMarquardt, 1e-9, null);
        System.out.println(MatrixFactory.columnBind(M, rslt.getFastSmoothedComponents(0)));

        //     CompositeSsf ssf = CompositeSsf.builder()
        //                .add(LocalLinearTrend.stateComponent(0.1, 0.1), LocalLinearTrend.defaultLoading())
        //                .add(Noise.of(1), Noise.defaultLoading())
        //                .add(SplineComponent.stateComponent(sd, 1, 0), SplineComponent.loading(sd, 0))
        //                .add(SplineComponent.stateComponent(sd2, 1, 0), SplineComponent.loading(sd2, 0))
        //                .build();
        //
        //
        //
        //        SmoothingOutput rslt = AkfToolkit.robustSmooth(ssf, new SsfData(data), true, true);
        ////        System.out.println(rslt.getSmoothing().getComponent(ssf.componentsPosition()[0]));
        ////        System.out.println(rslt.getSmoothing().getComponent(ssf.componentsPosition()[1]));
        //
        //        ISsfLoading loading = SplineComponent.loading(sd, 0);
        //        ISsfLoading loading2 = SplineComponent.loading(sd2, 0);
        //        DoubleSeq t = rslt.getSmoothing().getComponent(0);
        //        DoubleSeq n = rslt.getSmoothing().getComponent(2);
        //        for (int i = 0; i < rslt.getSmoothing().size(); ++i) {
        //            double z = loading.ZX(i, rslt.getSmoothing().a(i).extract(ssf.componentsPosition()[2],
        // ssf.componentsDimension()[2]));
        //            double z2 = loading2.ZX(i, rslt.getSmoothing().a(i).extract(ssf.componentsPosition()[3],
        // ssf.componentsDimension()[3]));
        //            System.out.print(data.get(i));
        //            System.out.print('\t');
        //            System.out.print(t.get(i));
        //            System.out.print('\t');
        //            System.out.print(n.get(i));
        //            System.out.print('\t');
        //            System.out.print(z);
        //            System.out.print('\t');
        //            System.out.println(z2);
        //        }
        //
    }
}
