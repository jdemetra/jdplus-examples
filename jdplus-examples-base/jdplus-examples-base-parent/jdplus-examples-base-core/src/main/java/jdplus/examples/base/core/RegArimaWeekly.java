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

import java.io.IOException;
import java.io.InputStream;
import jdplus.highfreq.base.api.ExtendedAirlineModellingSpec;
import jdplus.highfreq.base.api.ExtendedAirlineSpec;
import jdplus.highfreq.base.core.extendedairline.ExtendedAirlineKernel;
import jdplus.highfreq.base.core.regarima.HighFreqRegArimaModel;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.modelling.TransformationType;
import jdplus.toolkit.base.api.modelling.highfreq.TransformSpec;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.api.timeseries.regression.ModellingContext;
import jdplus.toolkit.base.core.stats.likelihood.LikelihoodStatistics;
import tck.demetra.data.Data;
import tck.demetra.data.MatrixSerializer;

/**
 *
 * @author palatej
 */
@lombok.experimental.UtilityClass
public class RegArimaWeekly {

    public void main(String[] args) throws IOException {

        // reads the data
        Matrix usclaims = usclaims();

        // creates the weekly time series
        TsData usclaims1 = TsData.of(TsPeriod.weekly(1967, 1, 5), usclaims.column(0));
        //TsData usclaims2 = TsData.of(TsPeriod.weekly(1967, 1, 5), usclaims.column(1)); // not used

        // RegArima processing
        // Creates a specification for a weekly series
        ExtendedAirlineModellingSpec spec = ExtendedAirlineModellingSpec.builder()
                .transform(TransformSpec.builder()
                        .function(TransformationType.Log)
                        .build())
                .stochastic(ExtendedAirlineSpec.DEFAULT_W)
                .build();
        ModellingContext context = new ModellingContext();
        ExtendedAirlineKernel kernel = ExtendedAirlineKernel.of(spec, context);
        HighFreqRegArimaModel model = kernel.process(usclaims1, null);
        LikelihoodStatistics statistics = model.getEstimation().getStatistics();
        System.out.println(statistics);

        // same without log transformation
        ExtendedAirlineModellingSpec spec2 = ExtendedAirlineModellingSpec.builder()
                .transform(TransformSpec.builder()
                        .function(TransformationType.None)
                        .build())
                .stochastic(ExtendedAirlineSpec.DEFAULT_W)
                .build();
        kernel = ExtendedAirlineKernel.of(spec2, context);
        model = kernel.process(usclaims1, null);

        statistics = model.getEstimation().getStatistics();
        System.out.println(statistics);
        // logs should be preferred (AIC <)
    }

    public Matrix usclaims() throws IOException {
        InputStream stream = Data.class.getResourceAsStream("usclaims.txt");
        Matrix usclaims = MatrixSerializer.read(stream);
        return usclaims;
    }

}
