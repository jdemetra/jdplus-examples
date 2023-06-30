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

import java.util.Arrays;
import jdplus.toolkit.base.api.data.AggregationType;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsDataTable;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.api.timeseries.TsUnit;
import jdplus.toolkit.base.core.arima.ArimaModel;
import jdplus.toolkit.base.core.math.linearfilters.BackFilter;
import jdplus.toolkit.base.core.ssf.StateComponent;
import jdplus.toolkit.base.core.ssf.arima.SsfArima;
import jdplus.toolkit.base.core.ssf.composite.CompositeSsf;
import jdplus.toolkit.base.core.ssf.dk.DkToolkit;
import jdplus.toolkit.base.core.ssf.sts.Noise;
import jdplus.toolkit.base.core.ssf.univariate.DefaultSmoothingResults;
import jdplus.toolkit.base.core.ssf.univariate.SsfData;
import jdplus.toolkit.base.core.stats.DescriptiveStatistics;
import tck.demetra.data.Data;

/**
 *
 * @author PALATEJ
 */
@lombok.experimental.UtilityClass
public class TimeSeries {

    @lombok.Value
    @lombok.Builder(builderClassName = "Builder")
    class Statistics {

        int n, nmissing;
        double max, min, average, stdev, q25, q50, q75;
    }

    public static void main(String[] arg) {
        TsData retail = TsData.ofInternal(TsPeriod.monthly(1992, 1), Data.RETAIL_BOOKSTORES);
        TsData nile = TsData.ofInternal(TsPeriod.yearly(1900), Data.NILE);

        TsData normalize = normalize(retail);
        TsData pct = pct(retail, 12);
        TsData delta = delta(retail, 1, 2);
        TsData aggregate = aggregate(retail, AggregationType.Average, 4, true);
        TsData aggregateMax = aggregate(retail.drop(0, 1), AggregationType.Max, 4, false);

        TsDataTable table = TsDataTable.of(Arrays.asList(normalize, pct, delta, aggregate, aggregateMax));
        System.out.println(table);

        hodrickPrescott(nile, 40);
        
        Statistics statistics = statistics(nile);
        System.out.println(statistics);
    }

    static TsData normalize(TsData s) {
        return s.normalize();
    }

    static TsData pct(TsData s, int lag) {
        return s.pctVariation(lag);
    }

    static TsData delta(TsData s, int lag, int power) {
        return s.delta(lag, power);
    }

    static TsData aggregate(TsData s, AggregationType type, int freq, boolean complete) {
        return s.aggregate(TsUnit.ofAnnualFrequency(freq), type, complete);
    }

    static void hodrickPrescott(TsData s, int lambda) {
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
        SsfData data = new SsfData(s.getValues());
        DefaultSmoothingResults rslts = DkToolkit.sqrtSmooth(ssf, data, true, true);
        int[] pos = ssf.componentsPosition();

        TsData trend = TsData.of(s.getStart(), rslts.getComponent(pos[0]));
        TsData noise = TsData.of(s.getStart(), rslts.getComponent(pos[1]));

        TsDataTable table = TsDataTable.of(Arrays.asList(s, trend, noise));
        System.out.println(table);
    }
    
    static Statistics statistics(TsData s){
        DescriptiveStatistics ds=DescriptiveStatistics.of(s.getValues());
        double[] quantiles = ds.quantiles(4);
        return Statistics.builder()
                .n(ds.getDataCount())
                .nmissing(ds.getMissingValuesCount())
                .max(ds.getMax())
                .min(ds.getMin())
                .average(ds.getAverage())
                .stdev(ds.getStdev())
                .q25(quantiles[0])
                .q50(quantiles[1])
                .q75(quantiles[2])
                .build();
                
    }

}
