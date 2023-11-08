package jdplus.examples.desktop.plugin;


import jdplus.toolkit.base.api.processing.ProcDocument;
import jdplus.toolkit.base.api.util.Id;
import jdplus.toolkit.base.api.util.LinearId;
import jdplus.toolkit.desktop.plugin.ui.processing.IProcDocumentItemFactory;
import nbbrd.design.DirectImpl;
import nbbrd.service.ServiceProvider;
import org.checkerframework.checker.nullness.qual.NonNull;

import javax.swing.*;

@DirectImpl
@ServiceProvider
public final class ProcDocumentItemFactoryExample implements IProcDocumentItemFactory {

    @Override
    public int getPosition() {
        return 1000;
    }

    @Override
    public @NonNull Class<? extends ProcDocument> getDocumentType() {
        return ProcDocument.class;
    }

    @Override
    public @NonNull Id getItemId() {
        return new LinearId("_ROOT_", "_CHILD_");
    }

    @Override
    public @NonNull JComponent getView(@NonNull ProcDocument procDocument) throws IllegalArgumentException {
        JTextArea result = new JTextArea();
        result.setEditable(false);
        result.setText(getProcDocumentAsString(procDocument));
        return new JScrollPane(result);
    }

    private static String getProcDocumentAsString(ProcDocument<?, ?, ?> procDocument) {
        return "KEY = " + procDocument.getKey() + System.lineSeparator()
                + "SPEC = " + procDocument.getSpecification() + System.lineSeparator()
                + "STATUS = " + procDocument.getStatus();
    }
}
