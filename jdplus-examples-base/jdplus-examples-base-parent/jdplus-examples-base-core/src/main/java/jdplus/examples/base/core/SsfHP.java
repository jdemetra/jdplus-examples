/*
 * Copyright 2022 National Bank of Belgium
 * 
 * Licensed under the EUPL, Version 1.2 or – as soon they will be approved 
 * by the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 * 
 * https://joinup.ec.europa.eu/software/page/eupl
 * 
 * Unless required by applicable law or agreed to in writing, software 
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and 
 * limitations under the Licence.
 */
package jdplus.examples.base.core;

import tck.demetra.data.Data;
import jdplus.toolkit.base.core.arima.ArimaModel;
import jdplus.toolkit.base.core.math.linearfilters.BackFilter;
import jdplus.toolkit.base.core.ssf.StateComponent;
import jdplus.toolkit.base.core.ssf.arima.SsfArima;
import jdplus.toolkit.base.core.ssf.composite.CompositeSsf;
import jdplus.toolkit.base.core.ssf.dk.DkToolkit;
import jdplus.toolkit.base.core.ssf.sts.Noise;
import jdplus.toolkit.base.core.ssf.univariate.DefaultSmoothingResults;
import jdplus.toolkit.base.core.ssf.univariate.SsfData;

/**
 * Model-based implementation of Hodrick-Prescott, using KF
 *
 * @author palatej
 */
@lombok.experimental.UtilityClass
public class SsfHP {

    public void main(String[] args) {
        // RW of order two

        // AR, D, MA, var
        double lambda = 40;
        ArimaModel rw2 = new ArimaModel(BackFilter.ONE, BackFilter.ofInternal(1, -2, 1), BackFilter.ONE, 1);
        StateComponent signal = SsfArima.stateComponent(rw2);
        StateComponent n = Noise.of(lambda);

        // create a composite state space form
        CompositeSsf ssf = CompositeSsf.builder()
                .add(signal, SsfArima.defaultLoading())
                .add(n, Noise.defaultLoading())
                .build();

        // smoothing using Durbin-Koopman for diffuse initialization
        // and with the specified variances (not estimated)
        SsfData data = new SsfData(Data.NILE);
        DefaultSmoothingResults rslts = DkToolkit.sqrtSmooth(ssf, data, true, true);
        int[] pos = ssf.componentsPosition();
        System.out.println(rslts.getComponent(pos[0]));
        System.out.println(rslts.getComponent(pos[1]));
    }
}
