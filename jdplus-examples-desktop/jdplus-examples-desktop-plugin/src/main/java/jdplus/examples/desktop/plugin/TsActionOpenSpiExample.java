package jdplus.examples.desktop.plugin;

import jdplus.toolkit.base.api.timeseries.Ts;
import jdplus.toolkit.desktop.plugin.TsActionOpenSpi;
import jdplus.toolkit.desktop.plugin.notification.MessageUtil;
import lombok.NonNull;
import nbbrd.design.DirectImpl;
import nbbrd.service.ServiceProvider;

@DirectImpl
@ServiceProvider
public final class TsActionOpenSpiExample implements TsActionOpenSpi {

    @Override
    public void open(@NonNull Ts ts) {
        MessageUtil.info("You just opened '" + ts.getName() + "'");
    }

    @Override
    public @NonNull String getName() {
        return "Message box with TS name";
    }
}
