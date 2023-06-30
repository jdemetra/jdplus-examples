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
import java.util.ArrayList;
import java.util.List;
import jdplus.highfreq.base.api.DecompositionSpec;
import jdplus.highfreq.base.api.ExtendedAirlineDecompositionSpec;
import jdplus.highfreq.base.api.ExtendedAirlineModellingSpec;
import jdplus.highfreq.base.api.ExtendedAirlineSpec;
import jdplus.highfreq.base.core.extendedairline.ExtendedAirlineKernel;
import jdplus.highfreq.base.core.extendedairline.ExtendedAirlineResults;
import jdplus.highfreq.base.core.extendedairline.decomposition.ExtendedAirlineDecompositionKernel;
import jdplus.highfreq.base.core.regarima.HighFreqRegArimaModel;
import jdplus.sa.base.api.ComponentType;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.modelling.ComponentInformation;
import jdplus.toolkit.base.api.modelling.TransformationType;
import jdplus.toolkit.base.api.modelling.highfreq.OutlierSpec;
import jdplus.toolkit.base.api.modelling.highfreq.TransformSpec;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsDataTable;
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
public class SaWeekly {

    public void main(String[] args) throws IOException {

        // reads the data
        Matrix usclaims = usclaims();

        // creates the weekly time series
        TsData usclaims1 = TsData.of(TsPeriod.weekly(1967, 1, 5), usclaims.column(0));
        //TsData usclaims2 = TsData.of(TsPeriod.weekly(1967, 1, 5), usclaims.column(1)); // not used

        ModellingContext context = new ModellingContext();
        // Not used in this example. See SaDaily

        // RegArima processing
        // Creates a specification for a weekly series
        ExtendedAirlineModellingSpec spec =  ExtendedAirlineModellingSpec.builder()
                .transform(TransformSpec.builder()
                        .function(TransformationType.Log)
                        .build())
                .stochastic(ExtendedAirlineSpec.DEFAULT_W)
                .outlier(OutlierSpec.builder()
                        .criticalValue(6)
                        .ao(true)
                        .build())
                .build();

        DecompositionSpec dspec = DecompositionSpec.builder()
                .periodicities(new double[]{365.25 / 7})
                .build();
        ExtendedAirlineDecompositionSpec allspec = ExtendedAirlineDecompositionSpec.builder()
                .preprocessing(spec)
                .decomposition(dspec)
                .build();
        ExtendedAirlineDecompositionKernel dkernel = new ExtendedAirlineDecompositionKernel(allspec, context);
        ExtendedAirlineResults rslts = dkernel.process(usclaims1, null);
        List<TsData> main = new ArrayList<>();
        main.add(rslts.getFinals().getSeries(ComponentType.Series, ComponentInformation.Value));
        main.add(rslts.getFinals().getSeries(ComponentType.SeasonallyAdjusted, ComponentInformation.Value));
        main.add(rslts.getFinals().getSeries(ComponentType.Trend, ComponentInformation.Value));
        main.add(rslts.getFinals().getSeries(ComponentType.Seasonal, ComponentInformation.Value));
        main.add(rslts.getFinals().getSeries(ComponentType.Irregular, ComponentInformation.Value));
        System.out.println(TsDataTable.of(main));

    }

    public Matrix usclaims() throws IOException {
        InputStream stream = Data.class.getResourceAsStream("usclaims.txt");
        Matrix usclaims = MatrixSerializer.read(stream);
        return usclaims;
    }

}
