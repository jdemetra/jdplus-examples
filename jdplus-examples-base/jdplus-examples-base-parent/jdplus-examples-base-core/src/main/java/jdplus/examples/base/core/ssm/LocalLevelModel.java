/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.examples.base.core.ssm;

import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.core.data.DataBlock;
import tck.demetra.data.Data;
import jdplus.toolkit.base.core.ssf.univariate.SsfData;
import jdplus.toolkit.base.core.ssf.StateComponent;
import jdplus.toolkit.base.core.ssf.sts.LocalLevel;
import jdplus.toolkit.base.core.ssf.univariate.Ssf;
import jdplus.toolkit.base.core.data.DataBlockStorage;
import jdplus.toolkit.base.core.math.functions.IParametersDomain;
import jdplus.toolkit.base.core.math.functions.ParamValidation;
import jdplus.toolkit.base.core.math.functions.levmar.LevenbergMarquardtMinimizer;
import jdplus.toolkit.base.core.math.functions.ssq.ISsqFunction;
import jdplus.toolkit.base.core.math.functions.ssq.ISsqFunctionPoint;
import jdplus.toolkit.base.core.ssf.ISsfLoading;
import jdplus.toolkit.base.core.ssf.composite.CompositeSsf;
import jdplus.toolkit.base.core.ssf.dk.DkToolkit;
import jdplus.toolkit.base.core.ssf.likelihood.DiffuseLikelihood;
import jdplus.toolkit.base.core.ssf.sts.Noise;



/**
 *
 * @author LEMASSO
 */
@lombok.experimental.UtilityClass
public class LocalLevelModel {
    public void main(String[] args){
        SsfData data= new SsfData(Data.NILE);
        
        // Example 1: Fixed parameters with both error terms in state equations
        StateComponent ll = LocalLevel.stateComponent(.01);
        StateComponent n = Noise.of(1);
        CompositeSsf ssf1 = CompositeSsf.builder()
                .add(ll, LocalLevel.defaultLoading())
                .add(n, Noise.defaultLoading())
                .build();
        DataBlockStorage rslt1 = DkToolkit.fastSmooth(ssf1, data);
        
        System.out.println("Example 1:");
        System.out.println(rslt1.item(0));
        System.out.println(rslt1.item(1));
        DataBlock Z=DataBlock.make(ssf1.getStateDim());
        ssf1.measurement().loading().Z(0, Z);
        System.out.println("Z");
        System.out.println(Z);
        System.out.println();
        
        // Example 2: Fixed parameters, with error term in the measurement equation
        ISsfLoading loading = LocalLevel.defaultLoading();      
        Ssf ssf2 = Ssf.of(ll, loading, 1);
        DataBlockStorage rslt2 = DkToolkit.fastSmooth(ssf2, data);
        
        System.out.println("Example 2:");
        System.out.println(rslt2.item(0));
        DataBlock Z2=DataBlock.make(ssf2.getStateDim());
        ssf2.measurement().loading().Z(0, Z2);
        System.out.println("Z");
        System.out.println(Z2);
        System.out.println();
        
        // Example 3: Estimated parameters
        LLM fn = new LLM(DoubleSeq.of(Data.NILE)); 
        LevenbergMarquardtMinimizer lm = LevenbergMarquardtMinimizer.builder()
                .functionPrecision(1e-5)
                .maxIter(100)
                .build();
        boolean ok = lm.minimize(fn.ssqEvaluate(DoubleSeq.of(.1)));
        LLM_Point result = (LLM_Point) lm.getResult(); // trick: make a cast here
        Ssf ssf3 = result.ssf();
        DoubleSeq p3 = result.getParameters(); // don't forget that they come from the parameters domain (here square root transformation) 
        DataBlockStorage rslt3 = DkToolkit.fastSmooth(ssf3, data);
        
        System.out.println("Example 3:");
        System.out.println(rslt3.item(0));
        System.out.println(rslt3.item(1));
        System.out.println("parameter");
        System.out.println(p3);
    }
}


class LLM implements ISsqFunction{ 
    
    private final DoubleSeq data;
    
    public LLM(DoubleSeq data){
        this.data=data;
    }
    
    public DoubleSeq getData(){
        return data;
    }
    
    @Override
    public IParametersDomain getDomain() {
        return new LLM_Domain();
    }

    @Override
    public ISsqFunctionPoint ssqEvaluate(DoubleSeq parameters) {
        return new LLM_Point(this, parameters.get(0)); 
    }
}

class LLM_Domain implements IParametersDomain{

    @Override
    public boolean checkBoundaries(DoubleSeq inparams) {
        return true; // always ok in this case
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
        return Double.NEGATIVE_INFINITY;
    }

    @Override
    public double ubound(int idx) {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public ParamValidation validate(DataBlock ioparams) {
        return ParamValidation.Valid; //always valid in this case
    }   
}

class LLM_Point implements ISsqFunctionPoint{
    
    private final LLM llm;
    private final double lstde;
    
    public LLM_Point(LLM llm, double le){
        this.llm=llm;
        this.lstde=le;
    }
    
    @Override
    public ISsqFunction getSsqFunction() {
        return llm;
    }

    @Override
    
    public DoubleSeq getE() {
        SsfData data= new SsfData(llm.getData());
        DiffuseLikelihood likelihood = DkToolkit.likelihood(ssf(), data, true, true);
        return likelihood.deviances();
    }

    @Override
    public DoubleSeq getParameters() {
        return DoubleSeq.of(lstde);
    }
    
    public Ssf ssf(){
        StateComponent ll = LocalLevel.stateComponent(lstde*lstde); // Math.exp(lstde) instead of lstde*lstde for Koopman transformation 
        StateComponent n = Noise.of(1);
        CompositeSsf ssf = CompositeSsf.builder()
                .add(ll, LocalLevel.defaultLoading())
                .add(n, Noise.defaultLoading())
                .build();
        return ssf;
    
    }
}











        