/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.examples.base.core.ssm;

import jdplus.sts.base.api.SeasonalModel;
import jdplus.sts.base.core.SeasonalComponent;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.data.DataBlockStorage;
import jdplus.toolkit.base.core.math.functions.IParametersDomain;
import jdplus.toolkit.base.core.math.functions.ParamValidation;
import jdplus.toolkit.base.core.math.functions.levmar.LevenbergMarquardtMinimizer;
import jdplus.toolkit.base.core.math.functions.ssq.ISsqFunction;
import jdplus.toolkit.base.core.math.functions.ssq.ISsqFunctionPoint;
import jdplus.toolkit.base.core.ssf.StateComponent;
import jdplus.toolkit.base.core.ssf.composite.CompositeSsf;
import jdplus.toolkit.base.core.ssf.dk.DkToolkit;
import jdplus.toolkit.base.core.ssf.likelihood.DiffuseLikelihood;
import jdplus.toolkit.base.core.ssf.sts.LocalLinearTrend;
import jdplus.toolkit.base.core.ssf.sts.Noise;
import jdplus.toolkit.base.core.ssf.univariate.Ssf;
import jdplus.toolkit.base.core.ssf.univariate.SsfData;
import tck.demetra.data.Data;

/**
 *
 * @author LEMASSO
 */
@lombok.experimental.UtilityClass
public class BasicStructuralModel {
     public void main(String[] args){
     
     SsfData data= new SsfData(Data.US_UNEMPL);
     
     // Example 1: Fixed parameters
     StateComponent llt = LocalLinearTrend.stateComponent(1, 1);
     StateComponent seas = SeasonalComponent.of(SeasonalModel.HarrisonStevens, 4, 1);
     StateComponent n = Noise.of(1);    
     CompositeSsf ssf1 = CompositeSsf.builder()
                .add(llt, LocalLinearTrend.defaultLoading())
                .add(seas, SeasonalComponent.defaultLoading())
                .add(n, Noise.defaultLoading())
                .build();
     DataBlockStorage rslt1 = DkToolkit.fastSmooth(ssf1, data);   
     System.out.println("Example 1:");
     System.out.println(rslt1.item(0));
     System.out.println(rslt1.item(1));
     System.out.println(rslt1.item(2));
     System.out.println(rslt1.item(3));
     System.out.println(rslt1.item(4));
     System.out.println(rslt1.item(5));
     DataBlock Z=DataBlock.make(ssf1.getStateDim());
     ssf1.measurement().loading().Z(0, Z);
     System.out.println("Z");
     System.out.println(Z);
     System.out.println();
     
     
     // Example 2: Estimated parameters
     BSM fn = new BSM(DoubleSeq.of(Data.US_UNEMPL));   
     ISsqFunctionPoint ssq = fn.ssqEvaluate(DoubleSeq.of(.1,.1,.1));
     LevenbergMarquardtMinimizer lm = LevenbergMarquardtMinimizer.builder()
             .functionPrecision(0.001)
             .maxIter(30)
             .build();
     boolean ok = lm.minimize(ssq);
     BSM_Point result = (BSM_Point) lm.getResult();
     Ssf ssf2 = result.ssf();
     DataBlockStorage rslt2 = DkToolkit.fastSmooth(ssf2, data);
     DoubleSeq p2 = result.getParameters();
     System.out.println("Example 2:");
     System.out.println(rslt2.item(0));
     System.out.println(rslt2.item(1));
     System.out.println(rslt2.item(2));
     System.out.println(rslt2.item(3));
     System.out.println(rslt2.item(4));
     System.out.println(rslt2.item(5));
     System.out.println("parameter");
     System.out.println(p2);
     
    }
}    
   
class BSM_Domain implements IParametersDomain{

    @Override
    public boolean checkBoundaries(DoubleSeq inparams) {
        return true;
    }

    @Override
    public double epsilon(DoubleSeq inparams, int idx) {
        return 1e-6;
    }

    @Override
    public int getDim() {
        return 3;
    }

    @Override
    public double lbound(int idx) {
        return Double.NEGATIVE_INFINITY;
    }

    @Override
    public double ubound(int idx) {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public ParamValidation validate(DataBlock ioparams) {
        return ParamValidation.Valid;
    }
}

class BSM implements ISsqFunction{
    
    private final DoubleSeq data; 
    
    public BSM(DoubleSeq data){
        this.data=data;
    }
    
    public DoubleSeq getData(){
        return data;
    }
    
    @Override
    public IParametersDomain getDomain() {
        return new BSM_Domain();
    }

    @Override
    public ISsqFunctionPoint ssqEvaluate(DoubleSeq parameters) {
        return new BSM_Point(this, parameters.get(0), parameters.get(1),parameters.get(2));
    }
}

class BSM_Point implements ISsqFunctionPoint{
    
    private final BSM bsm;
    private final double nstde;
    private final double slstde;
    private final double seasstde;
    // we fix the std error of the trend (level)
    
    public BSM_Point(BSM bsm, double ne, double sle, double sease){
        this.bsm=bsm;
        this.nstde=ne;
        this.slstde=sle;
        this.seasstde=sease;
    }
    
    @Override
    public ISsqFunction getSsqFunction() {
        return bsm;
    }

    @Override
    public DoubleSeq getE() {
       SsfData data= new SsfData(bsm.getData());
       DiffuseLikelihood likelihood = DkToolkit.likelihood(ssf(), data, true, true);
       return likelihood.deviances();
    }

    @Override
    public DoubleSeq getParameters() {
        return DoubleSeq.of(nstde, slstde, seasstde);
    }
    
    public Ssf ssf(){
        StateComponent llt = LocalLinearTrend.stateComponent(.1, slstde*slstde);
        StateComponent seas = SeasonalComponent.of(SeasonalModel.HarrisonStevens, 4, seasstde*seasstde);
        StateComponent n = Noise.of(nstde*nstde);
        CompositeSsf ssf = CompositeSsf.builder()
                .add(llt, LocalLinearTrend.defaultLoading())
                .add(seas, SeasonalComponent.defaultLoading())
                .add(n, Noise.defaultLoading())
                .build();
        return ssf;
    
    }
}
