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

import jdplus.toolkit.base.api.data.AggregationType;
import tck.demetra.data.Data;
import jdplus.benchmarking.base.api.univariate.TemporalDisaggregationSpec;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.benchmarking.base.core.univariate.TemporalDisaggregationProcessor;
import jdplus.benchmarking.base.core.univariate.TemporalDisaggregationResults;

/**
 *
 * Chow-Lin
 * @author palatej
 */
@lombok.experimental.UtilityClass
public class TemporalDisaggregation {
    
    public void main(String[] arg){
        TsData y = TsData.ofInternal(TsPeriod.yearly(1978),  Data.PCRA);
        TsData q = TsData.ofInternal(TsPeriod.quarterly(1977, 1),  Data.IND_PCR);
        TemporalDisaggregationSpec spec = TemporalDisaggregationSpec.builder()
                .aggregationType(AggregationType.Sum)
                .residualsModel(TemporalDisaggregationSpec.Model.Ar1)
                .constant(true)
                .build();
        TemporalDisaggregationResults rslt = TemporalDisaggregationProcessor.process(y, new TsData[]{q}, spec);
        System.out.println(rslt.getLikelihood());
        System.out.println(rslt.getDisaggregatedSeries());
    }
    
}
