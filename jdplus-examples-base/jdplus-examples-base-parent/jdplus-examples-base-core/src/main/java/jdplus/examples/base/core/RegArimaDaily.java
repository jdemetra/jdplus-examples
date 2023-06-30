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
import jdplus.highfreq.base.api.ExtendedAirlineModellingSpec;
import jdplus.highfreq.base.api.ExtendedAirlineSpec;
import jdplus.highfreq.base.core.extendedairline.ExtendedAirlineKernel;
import jdplus.highfreq.base.core.regarima.HighFreqRegArimaModel;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.modelling.TransformationType;
import jdplus.toolkit.base.api.modelling.highfreq.HolidaysSpec;
import jdplus.toolkit.base.api.modelling.highfreq.OutlierSpec;
import jdplus.toolkit.base.api.modelling.highfreq.RegressionSpec;
import jdplus.toolkit.base.api.modelling.highfreq.TransformSpec;
import jdplus.toolkit.base.api.processing.ProcessingLog;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.api.timeseries.calendars.Calendar;
import jdplus.toolkit.base.api.timeseries.calendars.EasterRelatedDay;
import jdplus.toolkit.base.api.timeseries.calendars.FixedDay;
import jdplus.toolkit.base.api.timeseries.calendars.Holiday;
import jdplus.toolkit.base.api.timeseries.calendars.HolidaysOption;
import jdplus.toolkit.base.api.timeseries.regression.ModellingContext;
import tck.demetra.data.Data;
import tck.demetra.data.MatrixSerializer;

/**
 *
 * @author palatej
 */
@lombok.experimental.UtilityClass
public class RegArimaDaily {
    public void main(String[] args) throws IOException {

        // reads the data
        Matrix edf = edf();

        // creates the weekly time series
        Holiday[] france = france();
        ModellingContext context=new ModellingContext();
        context.getCalendars().set("FR", new Calendar(france));
        // daily time series
        TsData EDF=TsData.of(TsPeriod.daily(1996, 1, 1) , edf.column(0));
        
        // build the psec
        ExtendedAirlineModellingSpec spec=ExtendedAirlineModellingSpec.builder()
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
                                        .holidays("FR")
                                        .holidaysOption(HolidaysOption.Skip)
                                        .single(false)
                                        .build())
                        .build())
                .build();
        ExtendedAirlineKernel kernel=ExtendedAirlineKernel.of(spec, context);
        HighFreqRegArimaModel rslt = kernel.process(EDF, ProcessingLog.dummy());
        Map<String, Class> dictionary = rslt.getDictionary();
        dictionary.keySet().forEach(v->System.out.println(v));
        
        // for instance, forecasts of the calendar effects for the next year
        System.out.println(rslt.getData("cal_f(365)", TsData.class));
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
