package jdplus.examples.base.core.claude;

import java.util.Locale;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.function.IntConsumer;
import java.util.stream.IntStream;

/**
 * Décomposition de Cholesky PAR BLOCS et PARALLÈLE, en Java pur, avec les conventions
 * de LAPACK (dpotrf / dpotrs) : matrices stockées par colonnes (column-major),
 * élément (i, j) à l'index i + j * ld.
 *
 * Pour chaque bloc diagonal de nb colonnes :
 *   1) factorisation du bloc diagonal (dpotf2)                    -> séquentiel (petit)
 *   2) résolution triangulaire du panneau (dtrsm)                 -> parallèle
 *   3) mise à jour du bloc restant (dsyrk), par tuiles            -> parallèle
 *
 * Les tâches parallèles écrivent dans des zones disjointes de la matrice et ne lisent que
 * des données non modifiées pendant l'étape en cours : aucune synchronisation n'est
 * nécessaire en dehors de l'attente de fin de chaque étape.
 */
public final class CholeskyParallel {

    /** Taille de bloc par défaut (à ajuster : 32, 64, 128). */
    public static final int DEFAULT_NB = 64;

    /** En dessous de cette taille, tout est exécuté sur le thread appelant. */
    private static final int PARALLEL_MIN_N = 256;

    private static final int ROW_TILE = 128; // tuile en lignes (mise à jour, cas inférieur)
    private static final int COL_TILE = 64; // tuile en colonnes (mise à jour, cas supérieur)
    private static final int PANEL_ROWS = 256; // lignes par tâche pour le panneau (cas inférieur)
    private static final int PANEL_COLS = 32; // colonnes par tâche pour le panneau (cas supérieur)
    private static final int UPDATE_COLS = 32; // colonnes par tâche pour la mise à jour (cas supérieur)

    private CholeskyParallel() {}

    // ------------------------------------------------------------------
    // dpotrf
    // ------------------------------------------------------------------

    /** nb = 64, pool commun de la JVM. */
    public static int dpotrf(boolean lower, int n, double[] a, int lda) {
        return dpotrf(lower, n, a, lda, DEFAULT_NB, ForkJoinPool.commonPool());
    }

    /**
     * Équivalent de LAPACK dpotrf, parallèle.
     * Calcule A = L * L^T (lower = true) ou A = U^T * U (lower = false), sur place.
     * Seul le triangle demandé est lu et écrasé.
     *
     * @param pool pool de threads à utiliser (parallélisme 1 ou null = exécution séquentielle)
     * @return 0 si succès ; k > 0 si le mineur principal d'ordre k n'est pas
     *         strictement positif (matrice non définie positive)
     */
    public static int dpotrf(boolean lower, int n, double[] a, int lda, int nb, ForkJoinPool pool) {
        checkArgs(n, a, lda, "a");
        if (nb < 1) throw new IllegalArgumentException("nb < 1");

        if (nb >= n) {
            return lower ? potf2Lower(a, lda, 0, n) : potf2Upper(a, lda, 0, n);
        }
        ForkJoinPool p = (n < PARALLEL_MIN_N) ? null : pool;
        return lower ? potrfLower(n, a, lda, nb, p) : potrfUpper(n, a, lda, nb, p);
    }

    // ---------------- Cas A = L * L^T ----------------

    private static int potrfLower(int n, double[] a, int lda, int nb, ForkJoinPool pool) {
        for (int j = 0; j < n; j += nb) {
            final int jb = Math.min(nb, n - j);
            final int jj = j + j * lda; // index de A(j, j)

            // 1) Bloc diagonal
            int info = potf2Lower(a, lda, jj, jb);
            if (info != 0) return j + info;

            final int m = n - j - jb; // lignes sous le bloc
            if (m <= 0) continue;

            final int p = (j + jb) + j * lda; // index de A21(0, 0)
            final int t = (j + jb) + (j + jb) * lda; // index de A22(0, 0)

            // 2) Panneau : A21 <- A21 * L11^{-T}. Les lignes sont indépendantes.
            int nChunks = (m + PANEL_ROWS - 1) / PANEL_ROWS;
            parallelFor(pool, nChunks, ch -> {
                int r0 = ch * PANEL_ROWS;
                int r1 = Math.min(r0 + PANEL_ROWS, m);
                for (int c = 0; c < jb; c++) {
                    int pc = p + c * lda;
                    double d = a[jj + c + c * lda];
                    for (int i = r0; i < r1; i++) {
                        a[pc + i] /= d;
                    }
                    for (int c2 = c + 1; c2 < jb; c2++) {
                        double l = a[jj + c2 + c * lda]; // L11(c2, c)
                        int pc2 = p + c2 * lda;
                        for (int i = r0; i < r1; i++) {
                            a[pc2 + i] -= a[pc + i] * l;
                        }
                    }
                }
            });

            // 3) Mise à jour : A22 <- A22 - A21 * A21^T (triangle inférieur).
            //    Chaque tuile de lignes écrit dans ses propres lignes de A22.
            int nTiles = (m + ROW_TILE - 1) / ROW_TILE;
            parallelFor(pool, nTiles, tile -> {
                int i0 = tile * ROW_TILE;
                int i1 = Math.min(i0 + ROW_TILE, m);
                for (int c = 0; c < i1; c++) {
                    int tc = t + c * lda;
                    int is = Math.max(i0, c);
                    for (int k = 0; k < jb; k++) {
                        int pk = p + k * lda;
                        double v = a[pk + c];
                        for (int i = is; i < i1; i++) {
                            a[tc + i] -= a[pk + i] * v;
                        }
                    }
                }
            });
        }
        return 0;
    }

    // ---------------- Cas A = U^T * U ----------------

    private static int potrfUpper(int n, double[] a, int lda, int nb, ForkJoinPool pool) {
        for (int j = 0; j < n; j += nb) {
            final int jb = Math.min(nb, n - j);
            final int jj = j + j * lda;

            // 1) Bloc diagonal
            int info = potf2Upper(a, lda, jj, jb);
            if (info != 0) return j + info;

            final int m = n - j - jb; // colonnes à droite du bloc
            if (m <= 0) continue;

            final int p = j + (j + jb) * lda; // index de A12(0, 0)
            final int t0 = (j + jb) + (j + jb) * lda; // index de A22(0, 0)

            // 2) Panneau : A12 <- U11^{-T} * A12. Les colonnes sont indépendantes.
            int nPanel = (m + PANEL_COLS - 1) / PANEL_COLS;
            parallelFor(pool, nPanel, ch -> {
                int c0 = ch * PANEL_COLS;
                int c1 = Math.min(c0 + PANEL_COLS, m);
                for (int c = c0; c < c1; c++) {
                    int pc = p + c * lda;
                    for (int r = 0; r < jb; r++) {
                        int ur = jj + r * lda; // colonne r de U11
                        double s = a[pc + r];
                        for (int k = 0; k < r; k++) {
                            s -= a[ur + k] * a[pc + k];
                        }
                        a[pc + r] = s / a[ur + r];
                    }
                }
            });

            // 3) Mise à jour : A22 <- A22 - A12^T * A12 (triangle supérieur).
            //    Chaque tâche écrit dans ses propres colonnes de A22.
            int nUpdate = (m + UPDATE_COLS - 1) / UPDATE_COLS;
            parallelFor(pool, nUpdate, ch -> {
                int c0 = ch * UPDATE_COLS;
                int c1 = Math.min(c0 + UPDATE_COLS, m);
                for (int rb = 0; rb < c1; rb += COL_TILE) {
                    int re = Math.min(rb + COL_TILE, m);
                    for (int c = Math.max(rb, c0); c < c1; c++) {
                        int pc = p + c * lda;
                        int tc = t0 + c * lda;
                        int rEnd = Math.min(re - 1, c);
                        for (int r = rb; r <= rEnd; r++) {
                            int pr = p + r * lda;
                            double s = 0.0;
                            for (int k = 0; k < jb; k++) {
                                s += a[pr + k] * a[pc + k];
                            }
                            a[tc + r] -= s;
                        }
                    }
                }
            });
        }
        return 0;
    }

    // ---------------- Factorisation non bloquée d'un bloc (dpotf2) ----------------

    /** Cholesky inférieure non bloquée sur le bloc m x m dont A(0,0) est à l'index base. */
    private static int potf2Lower(double[] a, int lda, int base, int m) {
        for (int j = 0; j < m; j++) {
            int cj = base + j * lda;
            double d = a[cj + j];
            if (!(d > 0.0)) return j + 1; // gère aussi NaN
            d = Math.sqrt(d);
            a[cj + j] = d;
            for (int i = j + 1; i < m; i++) {
                a[cj + i] /= d;
            }
            for (int c = j + 1; c < m; c++) {
                int cc = base + c * lda;
                double ljc = a[cj + c];
                for (int i = c; i < m; i++) {
                    a[cc + i] -= a[cj + i] * ljc;
                }
            }
        }
        return 0;
    }

    /** Cholesky supérieure non bloquée sur le bloc m x m dont A(0,0) est à l'index base. */
    private static int potf2Upper(double[] a, int lda, int base, int m) {
        for (int j = 0; j < m; j++) {
            int cj = base + j * lda;
            double d = a[cj + j];
            if (!(d > 0.0)) return j + 1;
            d = Math.sqrt(d);
            a[cj + j] = d;
            for (int c = j + 1; c < m; c++) {
                a[base + j + c * lda] /= d;
            }
            for (int c = j + 1; c < m; c++) {
                int cc = base + c * lda;
                double ujc = a[cc + j];
                for (int r = j + 1; r <= c; r++) {
                    a[cc + r] -= a[base + j + r * lda] * ujc;
                }
            }
        }
        return 0;
    }

    // ------------------------------------------------------------------
    // dpotrs
    // ------------------------------------------------------------------

    /** Pool commun de la JVM. */
    public static void dpotrs(boolean lower, int n, int nrhs, double[] a, int lda, double[] b, int ldb) {
        dpotrs(lower, n, nrhs, a, lda, b, ldb, ForkJoinPool.commonPool());
    }

    /**
     * Équivalent de LAPACK dpotrs, parallèle sur les seconds membres (colonnes de B).
     * Résout A * X = B, où A = L * L^T (lower = true) ou A = U^T * U (lower = false)
     * a été factorisée par dpotrf. B (n x nrhs) est écrasé par la solution X.
     * Avec un seul second membre, la résolution est séquentielle.
     */
    public static void dpotrs(
            boolean lower, int n, int nrhs, double[] a, int lda, double[] b, int ldb, ForkJoinPool pool) {
        checkArgs(n, a, lda, "a");
        if (nrhs < 0) throw new IllegalArgumentException("nrhs < 0");
        if (ldb < Math.max(1, n)) throw new IllegalArgumentException("ldb < max(1, n)");
        if (b.length < ldb * nrhs) throw new IllegalArgumentException("b trop petit");

        parallelFor(pool, nrhs, k -> solveOne(lower, n, a, lda, b, k * ldb));
    }

    /** Résout pour une colonne de B commençant à l'index o. */
    private static void solveOne(boolean lower, int n, double[] a, int lda, double[] b, int o) {
        if (lower) {
            for (int j = 0; j < n; j++) { // L * y = b
                double bj = b[o + j] / a[j + j * lda];
                b[o + j] = bj;
                for (int i = j + 1; i < n; i++) {
                    b[o + i] -= bj * a[i + j * lda];
                }
            }
            for (int j = n - 1; j >= 0; j--) { // L^T * x = y
                double t = b[o + j];
                for (int i = j + 1; i < n; i++) {
                    t -= a[i + j * lda] * b[o + i];
                }
                b[o + j] = t / a[j + j * lda];
            }
        } else {
            for (int j = 0; j < n; j++) { // U^T * y = b
                double t = b[o + j];
                for (int i = 0; i < j; i++) {
                    t -= a[i + j * lda] * b[o + i];
                }
                b[o + j] = t / a[j + j * lda];
            }
            for (int j = n - 1; j >= 0; j--) { // U * x = y
                double bj = b[o + j] / a[j + j * lda];
                b[o + j] = bj;
                for (int i = 0; i < j; i++) {
                    b[o + i] -= bj * a[i + j * lda];
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Outils
    // ------------------------------------------------------------------

    /**
     * Exécute body(0..count-1) en parallèle dans le pool donné et attend la fin de toutes
     * les tâches. Séquentiel si pool est null, de parallélisme 1, ou si count <= 1.
     */
    private static void parallelFor(ForkJoinPool pool, int count, IntConsumer body) {
        if (count <= 0) return;
        if (pool == null || pool.getParallelism() <= 1 || count == 1) {
            for (int i = 0; i < count; i++) body.accept(i);
            return;
        }
        try {
            pool.submit(() -> IntStream.range(0, count).parallel().forEach(body))
                    .get();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Calcul interrompu", e);
        } catch (ExecutionException e) {
            Throwable c = e.getCause();
            if (c instanceof RuntimeException) throw (RuntimeException) c;
            if (c instanceof Error) throw (Error) c;
            throw new RuntimeException(c);
        }
    }

    /** Met à zéro le triangle non utilisé, pour obtenir une vraie matrice L (ou U). */
    public static void zeroOtherTriangle(boolean lower, int n, double[] a, int lda) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                if ((lower && i < j) || (!lower && i > j)) {
                    a[i + j * lda] = 0.0;
                }
            }
        }
    }

    private static void checkArgs(int n, double[] a, int lda, String name) {
        if (n < 0) throw new IllegalArgumentException("n < 0");
        if (lda < Math.max(1, n)) throw new IllegalArgumentException("lda < max(1, n)");
        if (a.length < lda * n) throw new IllegalArgumentException(name + " trop petit");
    }

    // ------------------------------------------------------------------
    // Test : référence non bloquée, bloc séquentiel, bloc parallèle
    // ------------------------------------------------------------------

    /** Meilleur temps (ms) sur 3 exécutions ; le dernier facteur calculé est copié dans out. */
    private static double bestTimeMs(boolean lower, int n, double[] a, int nb, ForkJoinPool pool, double[] out) {
        double best = Double.MAX_VALUE;
        for (int rep = 0; rep < 3; rep++) {
            double[] f = a.clone();
            long t0 = System.nanoTime();
            int info = dpotrf(lower, n, f, n, nb, pool);
            double ms = (System.nanoTime() - t0) / 1e6;
            if (info != 0) throw new ArithmeticException("Non définie positive, info = " + info);
            best = Math.min(best, ms);
            System.arraycopy(f, 0, out, 0, f.length);
        }
        return best;
    }

    private static double maxDiff(boolean lower, int n, double[] x, double[] y) {
        double d = 0;
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                if ((lower && i >= j) || (!lower && i <= j)) {
                    d = Math.max(d, Math.abs(x[i + j * n] - y[i + j * n]));
                }
            }
        }
        return d;
    }

    public static void main(String[] args) {
        int n = 1500;
        int threads = Runtime.getRuntime().availableProcessors();
        Random rnd = new Random(42);

        // Matrice symétrique à diagonale strictement dominante => définie positive
        double[] a = new double[n * n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i <= j; i++) {
                double v = rnd.nextDouble() * 2 - 1;
                a[i + j * n] = v;
                a[j + i * n] = v;
            }
        }
        for (int i = 0; i < n; i++) a[i + i * n] = n;

        // 4 seconds membres pour tester la résolution parallèle
        int nrhs = 4;
        double[] b = new double[n * nrhs];
        for (int i = 0; i < b.length; i++) b[i] = rnd.nextDouble();

        ForkJoinPool seqPool = new ForkJoinPool(1);
        ForkJoinPool parPool = new ForkJoinPool(threads);
        System.out.println("n = " + n + ", threads = " + threads + ", nb = " + DEFAULT_NB);

        try {
            for (boolean lower : new boolean[] {true, false}) {
                double[] fRef = new double[n * n];
                double[] fSeq = new double[n * n];
                double[] fPar = new double[n * n];

                double tRef = bestTimeMs(lower, n, a, n, null, fRef); // non bloqué
                double tSeq = bestTimeMs(lower, n, a, DEFAULT_NB, seqPool, fSeq);
                double tPar = bestTimeMs(lower, n, a, DEFAULT_NB, parPool, fPar);

                // Résolution parallèle sur nrhs seconds membres, vérification du résidu
                double[] x = b.clone();
                dpotrs(lower, n, nrhs, fPar, n, x, n, parPool);
                double maxRes = 0;
                for (int k = 0; k < nrhs; k++) {
                    for (int i = 0; i < n; i++) {
                        double s = 0;
                        for (int j = 0; j < n; j++) s += a[i + j * n] * x[j + k * n];
                        maxRes = Math.max(maxRes, Math.abs(s - b[i + k * n]));
                    }
                }

                System.out.println(lower ? "Inférieur (L)" : "Supérieur (U)");
                System.out.printf(
                        Locale.FRANCE,
                        "  non bloqué : %.1f ms | bloc séquentiel : %.1f ms | bloc parallèle : %.1f ms%n",
                        tRef,
                        tSeq,
                        tPar);
                System.out.printf(Locale.FRANCE, "  accélération parallèle vs bloc séquentiel : x%.2f%n", tSeq / tPar);
                System.out.println("  écart parallèle vs non bloqué = " + maxDiff(lower, n, fPar, fRef));
                System.out.println("  écart parallèle vs bloc séq.  = " + maxDiff(lower, n, fPar, fSeq));
                System.out.println("  résidu max |A x - b|          = " + maxRes);
            }
        } finally {
            seqPool.shutdown();
            parPool.shutdown();
        }
    }
}
