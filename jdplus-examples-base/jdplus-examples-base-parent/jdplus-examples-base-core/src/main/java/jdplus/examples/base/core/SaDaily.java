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
import java.util.Map;
import jdplus.highfreq.base.api.DecompositionSpec;
import jdplus.highfreq.base.api.ExtendedAirlineDecompositionSpec;
import jdplus.highfreq.base.api.ExtendedAirlineModellingSpec;
import jdplus.highfreq.base.api.ExtendedAirlineSpec;
import jdplus.highfreq.base.core.extendedairline.ExtendedAirlineResults;
import jdplus.highfreq.base.core.extendedairline.decomposition.ExtendedAirlineDecompositionKernel;
import jdplus.sa.base.api.ComponentType;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.modelling.ComponentInformation;
import jdplus.toolkit.base.api.modelling.TransformationType;
import jdplus.toolkit.base.api.modelling.highfreq.HolidaysSpec;
import jdplus.toolkit.base.api.modelling.highfreq.OutlierSpec;
import jdplus.toolkit.base.api.modelling.highfreq.RegressionSpec;
import jdplus.toolkit.base.api.modelling.highfreq.TransformSpec;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsDataTable;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.api.timeseries.calendars.Calendar;
import jdplus.toolkit.base.api.timeseries.calendars.EasterRelatedDay;
import jdplus.toolkit.base.api.timeseries.calendars.FixedDay;
import jdplus.toolkit.base.api.timeseries.calendars.Holiday;
import jdplus.toolkit.base.api.timeseries.calendars.HolidaysOption;
import jdplus.toolkit.base.api.timeseries.regression.ModellingContext;
import jdplus.toolkit.base.core.stats.likelihood.LikelihoodStatistics;
import tck.demetra.data.Data;
import tck.demetra.data.MatrixSerializer;

/**
 *
 * @author palatej
 */
@lombok.experimental.UtilityClass
public class SaDaily {

    public void main(String[] args) throws IOException {

        // reads the data
        Matrix edf = edf();
        // daily time series
        TsData EDF = TsData.of(TsPeriod.daily(1996, 1, 1), edf.column(0));

        // creates the holidays variable
        Holiday[] france = france();
        // creates a modelling context to store the different variables used in the modelling
        ModellingContext context = new ModellingContext();
        context.getCalendars().set("FR", new Calendar(france));

        // build the spec
        ExtendedAirlineModellingSpec spec = ExtendedAirlineModellingSpec.builder()
                .transform(TransformSpec.builder()
                        .function(TransformationType.Log)
                        .build())
                .stochastic(ExtendedAirlineSpec.DEFAULT_WD)
                .outlier(OutlierSpec.builder()
                        .criticalValue(8)
                        .ao(true)
                        .build())
                .regression(RegressionSpec.builder()
                        .holidays(HolidaysSpec.builder()
                                // See modelling context
                                .holidays("FR")
                                // when an holiday falls on a Sunday, we consider it as a Sunday, not as an holiday
                                .holidaysOption(HolidaysOption.Skip)
                                // several regression variables, one for each holiday type (Christmas, New Year...)
                                .single(false)
                                .build())
                        .build())
                .build();

        DecompositionSpec dspec = DecompositionSpec.builder()
                .periodicities(new double[]{7, 365.25})
                .forecastsCount(100)
                .build();
        ExtendedAirlineDecompositionSpec allspec = ExtendedAirlineDecompositionSpec.builder()
                .preprocessing(spec)
                .decomposition(dspec)
                .build();
        ExtendedAirlineDecompositionKernel kernel = new ExtendedAirlineDecompositionKernel(allspec, context);
        ExtendedAirlineResults rslt = kernel.process(EDF, null);
        List<TsData> main = new ArrayList<>();
        main.add(rslt.getFinals().getSeries(ComponentType.Series, ComponentInformation.Value));
        main.add(rslt.getFinals().getSeries(ComponentType.SeasonallyAdjusted, ComponentInformation.Value));
        main.add(rslt.getFinals().getSeries(ComponentType.Trend, ComponentInformation.Value));
        main.add(rslt.getFinals().getSeries(ComponentType.Seasonal, ComponentInformation.Value));
        main.add(rslt.getFinals().getSeries(ComponentType.Irregular, ComponentInformation.Value));
        System.out.println(TsDataTable.of(main));
        Map<String, Class> dictionary = rslt.getDictionary();
        dictionary.keySet().forEach(v -> System.out.println(v));

        
        // for instance (likelihood + seasonal forecast)
        
        // direct function call
        LikelihoodStatistics statistics = rslt.getPreprocessing().getEstimation().getStatistics();
        System.out.println(statistics);
        // or use of the dictionary (generic interface)
        System.out.println(rslt.getData("s_f", TsData.class));

    }

    private void addDefault(List<Holiday> holidays) {
        holidays.add(FixedDay.NEWYEAR);
        holidays.add(FixedDay.MAYDAY);
        holidays.add(FixedDay.ASSUMPTION);
        holidays.add(FixedDay.ALLSAINTSDAY);
        holidays.add(FixedDay.CHRISTMAS);
        holidays.add(EasterRelatedDay.EASTERMONDAY);
        holidays.add(EasterRelatedDay.ASCENSION);
        holidays.add(EasterRelatedDay.WHITMONDAY);
    }

    public Holiday[] france() {
        List<Holiday> holidays = new ArrayList<>();
        addDefault(holidays);
        holidays.add(new FixedDay(5, 8));
        holidays.add(new FixedDay(7, 14));
        holidays.add(FixedDay.ARMISTICE);
        return holidays.stream().toArray(i -> new Holiday[i]);
    }

    public Matrix edf() throws IOException {
        // easily replace by another txt or csv file (series put in columns, no header)
        InputStream stream = Data.class.getResourceAsStream("edf.txt");
        return MatrixSerializer.read(stream);
    }

}
