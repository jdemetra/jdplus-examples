/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.examples.base.core.ssm;

import jdplus.toolkit.base.core.data.DataBlockStorage;
import jdplus.toolkit.base.core.ssf.State;
import jdplus.toolkit.base.core.ssf.StateComponent;
import jdplus.toolkit.base.core.ssf.StateInfo;
import jdplus.toolkit.base.core.ssf.UpdateInformation;
import jdplus.toolkit.base.core.ssf.composite.CompositeSsf;
import jdplus.toolkit.base.core.ssf.dk.DefaultDiffuseFilteringResults;
import jdplus.toolkit.base.core.ssf.dk.DkToolkit;
import jdplus.toolkit.base.core.ssf.sts.LocalLevel;
import jdplus.toolkit.base.core.ssf.sts.LocalLinearTrend;
import jdplus.toolkit.base.core.ssf.sts.Noise;
import jdplus.toolkit.base.core.ssf.univariate.DefaultSmoothingResults;
import jdplus.toolkit.base.core.ssf.univariate.IFilteringResults;
import jdplus.toolkit.base.core.ssf.univariate.SsfData;
import jdplus.toolkit.base.core.ssf.univariate.StateFilteringResults;
import tck.demetra.data.Data;

/**
 *
 * @author LEMASSO
 */
@lombok.experimental.UtilityClass
public class HodrickPrescott {
    public void main(String[] args){
        SsfData data= new SsfData(Data.US_UNEMPL);
        
        StateComponent llt = LocalLinearTrend.stateComponent(0, 0.000625);
        StateComponent n = Noise.of(1);
        CompositeSsf ssf = CompositeSsf.builder()
                .add(llt, LocalLinearTrend.defaultLoading())
                .add(n, Noise.defaultLoading())
                .build();      
        
        // Two-Sided HP Filter
        DefaultSmoothingResults hp2 = DkToolkit.smooth(ssf, data, true, false);
        System.out.println("Two-Sided HP Filter");
        System.out.println(hp2.getComponent(0));
        
        // One-Sided HP Filter (!contemporaneous filtering: E(a(t)|y(0)…y(t)))
        StateFilteringResults fr = new StateFilteringResults(StateInfo.Concurrent, true);
        fr.prepare(ssf.getStateDim(), 0, data.length());
        DkToolkit.sqrtFilter(ssf, data, fr, true);
        System.out.println("One-Sided HP Filter");
        System.out.println(fr.getComponent(0));        
    }
}


