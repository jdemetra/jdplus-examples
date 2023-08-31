/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.examples.base.core.ssm;

import java.util.Arrays;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.data.DataBlockStorage;
import jdplus.toolkit.base.core.math.functions.IParametersDomain;
import jdplus.toolkit.base.core.math.functions.ParamValidation;
import jdplus.toolkit.base.core.math.functions.levmar.LevenbergMarquardtMinimizer;
import jdplus.toolkit.base.core.math.functions.ssq.ISsqFunction;
import jdplus.toolkit.base.core.math.functions.ssq.ISsqFunctionPoint;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.ssf.ISsfLoading;
import jdplus.toolkit.base.core.ssf.StateComponent;
import jdplus.toolkit.base.core.ssf.arima.AR1;
import jdplus.toolkit.base.core.ssf.basic.Coefficients;
import jdplus.toolkit.base.core.ssf.basic.Loading;
import jdplus.toolkit.base.core.ssf.benchmarking.SsfCumulator;
import jdplus.toolkit.base.core.ssf.composite.CompositeLoading;
import jdplus.toolkit.base.core.ssf.composite.CompositeState;
import jdplus.toolkit.base.core.ssf.dk.DkToolkit;
import jdplus.toolkit.base.core.ssf.likelihood.DiffuseLikelihood;
import jdplus.toolkit.base.core.ssf.sts.LocalLevel;
import jdplus.toolkit.base.core.ssf.univariate.Ssf;
import jdplus.toolkit.base.core.ssf.univariate.SsfData;
import tck.demetra.data.Data;

/**
 *
 * @author LEMASSO
 */
@lombok.experimental.UtilityClass
public class Chow_Lin {
    
    public void main(String[] args){
        
        // Prepare data
        int freq = 4;
        double yc[] = new double[Data.PCRA.length * freq]; 
        Arrays.fill(yc, Double.NaN);
        for (int i = 0; i < yc.length; i++){
            if ((i+1) % freq == 0) {
                int k = (i+1) / freq - 1;
                yc[i] = Data.PCRA[k];
            }
        }
        FastMatrix Myc = FastMatrix.make(yc.length, 1);
        Myc.column(0).copyFrom(yc, 0);
        SsfData data= new SsfData(yc);
        
        FastMatrix Mx = FastMatrix.make(Data.IND_PCR.length, 1);
        Mx.column(0).copyFrom(Data.IND_PCR, 0);
        
        // Example 1: Chow-Lin with no constant and no indicator - fixed rho
        StateComponent stoch = AR1.of(.9);
        StateComponent c = SsfCumulator.of(stoch, AR1.defaultLoading(), 4, 0);
        Ssf ssf1 = Ssf.of(c, SsfCumulator.defaultLoading(AR1.defaultLoading(), 4, 0));
        DataBlockStorage rslt1 = DkToolkit.fastSmooth(ssf1, data); 
        System.out.println("Example 1:");
        System.out.println(rslt1.item(1));
 
        // Example 2: Chow-Lin with constant and indicator - fixed rho
        StateComponent cst = LocalLevel.stateComponent(0);
        StateComponent stoch2 = AR1.of(.9);
        StateComponent coef = Coefficients.fixedCoefficients(1); //class RegSsf could also be used
        StateComponent cstate = CompositeState.builder()
                .add(cst)
                .add(stoch2)
                .add(coef)
                .build();
        CompositeLoading cloading = new CompositeLoading(new int [] {1,1,1}, new ISsfLoading[] {LocalLevel.defaultLoading(), AR1.defaultLoading(), Loading.regression(Mx)});
        StateComponent c2 = SsfCumulator.of(cstate, cloading, 4, 0);
        Ssf ssf2 = Ssf.of(c2, SsfCumulator.defaultLoading(cloading, 4, 0));
        DataBlockStorage rslt2 = DkToolkit.fastSmooth(ssf2, data); 
        System.out.println("Example 2:");
        System.out.println(rslt2.item(0));
        System.out.println(rslt2.item(1));
        System.out.println(rslt2.item(2));
        System.out.println(rslt2.item(3));
      
        // Example 3: Chow-Lin with constant and indicator - rho unknown
        CL fn = new CL(Data.PCRA, Data.IND_PCR);   
        ISsqFunctionPoint ssq = fn.ssqEvaluate(DoubleSeq.of(.5));
        LevenbergMarquardtMinimizer lm = LevenbergMarquardtMinimizer.builder()
             .functionPrecision(0.000000001)
             .maxIter(30)
             .build();
        boolean ok = lm.minimize(ssq);
        CL_Point result = (CL_Point) lm.getResult();
        Ssf ssf3 = result.ssf();
        DataBlockStorage rslt3 = DkToolkit.fastSmooth(ssf3, data);
        System.out.println("Example 3:");
        System.out.println(rslt3.item(0));
        System.out.println(rslt3.item(1));
        System.out.println(rslt3.item(2));
        System.out.println(rslt3.item(3));
        DoubleSeq p3 = result.getParameters();
        System.out.println("parameter");
        System.out.println(p3);
        
    }  
}

class CL_domain implements IParametersDomain{

    @Override
    public boolean checkBoundaries(DoubleSeq inparams) {
        return Math.abs(inparams.get(0))<=.999;
    }

    @Override
    public double epsilon(DoubleSeq inparams, int idx) {
        return 1e-6;
    }

    @Override
    public int getDim() {
        return 1;
    }

    @Override
    public double lbound(int idx) {
        return -.999;
    }

    @Override
    public double ubound(int idx) {
        return .999;
    }

    @Override
    public ParamValidation validate(DataBlock ioparams) {
        double r = ioparams.get(0);
        if(Math.abs(r) <= .999){
            return ParamValidation.Valid;
        }else{
            r = 1/r;
            if(Math.abs(r) > .999){
                r = .998;
            }
            ioparams.set(0,r);
            return ParamValidation.Changed;
        }
    }
}

class CL implements ISsqFunction{
    
    private final double[] series;
    private final double[] indicator;
    
    public CL(double[] series, double[] quarterlyIndicator){
        this.series = series;
        this.indicator = quarterlyIndicator;
    }
    
    public double[] getY(){
        return series;
    }
    
    public double[] getX(){
        return indicator;
    }
    
    @Override
    public IParametersDomain getDomain() {
        return new CL_domain();
    }

    @Override
    public ISsqFunctionPoint ssqEvaluate(DoubleSeq parameters) {
        return new CL_Point(this, parameters.get(0));
    }
}

class CL_Point implements ISsqFunctionPoint{
    
    private final CL cl;
    private final double rho;
    
    public CL_Point(CL cl, double rho){
        this.cl = cl;
        this.rho = rho;
    }
    
    @Override
    public ISsqFunction getSsqFunction() {
        return cl;
    }

    @Override
    public DoubleSeq getE() {
        SsfData yc= new SsfData(CumulateSeriesByYear(cl.getY(),4));
        DiffuseLikelihood likelihood = DkToolkit.likelihood(ssf(), yc, true, true);
        return likelihood.deviances();
    }

    @Override
    public DoubleSeq getParameters() {
        return DoubleSeq.of(rho);
    }
    
    public SsfData CumulateSeriesByYear(double[] series, int frequency){
        double yc[] = new double[series.length * frequency]; 
        Arrays.fill(yc, Double.NaN);
        for (int i = 0; i < yc.length; i++){
            if ((i+1) % frequency == 0) {
                int k = (i+1) / frequency - 1;
                yc[i] = series[k];
            }
        }
        return new SsfData(yc);
    }
    
    public Ssf ssf(){
        StateComponent cst = LocalLevel.stateComponent(0);
        StateComponent stoch2 = AR1.of(rho);
        StateComponent coef = Coefficients.fixedCoefficients(1); 
        StateComponent cstate = CompositeState.builder()
                .add(cst)
                .add(stoch2)
                .add(coef)
                .build();
        FastMatrix Mx = FastMatrix.make(cl.getX().length, 1);
        Mx.column(0).copyFrom(cl.getX(), 0);
        CompositeLoading cloading = new CompositeLoading(new int [] {1,1,1}, new ISsfLoading[] {LocalLevel.defaultLoading(), AR1.defaultLoading(), Loading.regression(Mx)});
        StateComponent c = SsfCumulator.of(cstate, cloading, 4, 0);
        Ssf ssf = Ssf.of(c, SsfCumulator.defaultLoading(cloading, 4, 0));
        return ssf;
    } 
} 