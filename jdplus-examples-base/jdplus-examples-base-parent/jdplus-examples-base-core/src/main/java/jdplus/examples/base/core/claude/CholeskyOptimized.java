package jdplus.examples.base.core.claude;

import java.util.Locale;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.function.IntConsumer;
import java.util.stream.IntStream;

/**
 * Décomposition de Cholesky par blocs, parallèle, avec micro-noyau en registres, en Java pur.
 * Conventions de LAPACK (dpotrf / dpotrs) : matrices stockées par colonnes (column-major),
 * élément (i, j) à l'index i + j * ld.
 *
 * Pour chaque bloc diagonal de nb colonnes :
 *   1) factorisation du bloc diagonal (dpotf2)                       -> séquentiel (petit)
 *   2) résolution triangulaire du panneau (dtrsm) + PACKING du panneau
 *      dans un tampon contigu, par micro-panneaux de 4 lignes/colonnes -> parallèle
 *   3) mise à jour du bloc restant (dsyrk) avec un micro-noyau 4x4 dont les 16 sommes
 *      partielles restent dans des variables locales (registres)      -> parallèle
 *
 * Le packing rend les accès du noyau contigus (pas de saut de lda entre deux itérations
 * de k), et le noyau 4x4 fait 16 multiplications-additions pour 8 lectures, au lieu de
 * 1 lecture + 1 lecture/écriture de C par multiplication-addition dans la version simple.
 *
 * Aucune dépendance. Math.fma n'est volontairement pas utilisé : sans instruction FMA
 * matérielle, il est très lent en Java.
 */
public final class CholeskyOptimized {

    /** Taille de bloc par défaut. Le micro-noyau profite d'un bloc plus grand : essayer 96 ou 128. */
    public static final int DEFAULT_NB = 64;

    /** En dessous de cette taille, tout est exécuté sur le thread appelant. */
    private static final int PARALLEL_MIN_N = 256;

    private static final int PANEL_ROWS = 256; // lignes par tâche pour le panneau (cas L) - multiple de 4
    private static final int PANEL_COLS = 32; // colonnes par tâche pour le panneau (cas U) - multiple de 4

    private CholeskyOptimized() {}

    // ------------------------------------------------------------------
    // dpotrf
    // ------------------------------------------------------------------

    /** nb = 64, pool commun de la JVM. */
    public static int dpotrf(boolean lower, int n, double[] a, int lda) {
        return dpotrf(lower, n, a, lda, DEFAULT_NB, ForkJoinPool.commonPool());
    }

    /**
     * Équivalent de LAPACK dpotrf.
     * Calcule A = L * L^T (lower = true) ou A = U^T * U (lower = false), sur place.
     * Seul le triangle demandé est lu et écrasé.
     *
     * @param pool pool de threads (null ou parallélisme 1 = exécution séquentielle,
     *             le micro-noyau est utilisé dans tous les cas)
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
        // Tampon de packing : au plus ceil(n/4) micro-panneaux de 4 * nb valeurs
        double[] ap = new double[((n + 3) / 4) * 4 * nb];
        return lower ? potrfLower(n, a, lda, nb, p, ap) : potrfUpper(n, a, lda, nb, p, ap);
    }

    // ---------------- Cas A = L * L^T ----------------

    private static int potrfLower(int n, double[] a, int lda, int nb, ForkJoinPool pool, double[] ap) {
        final int par = (pool == null) ? 1 : Math.max(1, pool.getParallelism());

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

            // 2) Panneau : A21 <- A21 * L11^{-T}, puis packing des lignes traitées.
            //    Micro-panneau q = lignes 4q..4q+3 ; ap[q*4*jb + k*4 + r] = A21(4q + r, k),
            //    complété par des zéros au-delà de la ligne m-1.
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
                for (int q = r0 / 4; q * 4 < r1; q++) {
                    int base = q * 4 * jb;
                    for (int k = 0; k < jb; k++) {
                        int src = p + k * lda;
                        int dst = base + k * 4;
                        for (int s = 0; s < 4; s++) {
                            int row = q * 4 + s;
                            ap[dst + s] = (row < m) ? a[src + row] : 0.0;
                        }
                    }
                }
            });

            // 3) Mise à jour : A22 <- A22 - A21 * A21^T (triangle inférieur).
            //    Chaque tâche traite une tranche de micro-panneaux de lignes [qa, qb).
            final int nq = (m + 3) / 4;
            final int tileMp = Math.max(2, Math.min(16, (nq + 4 * par - 1) / (4 * par)));
            int nTiles = (nq + tileMp - 1) / tileMp;
            parallelFor(pool, nTiles, tile -> {
                double[] tmp = new double[16];
                int qa = tile * tileMp;
                int qb = Math.min(qa + tileMp, nq);
                for (int qc = 0; qc < qb; qc++) {
                    int colStart = qc * 4;
                    int nc = Math.min(4, m - colStart);
                    int ib = qc * 4 * jb;
                    for (int qi = Math.max(qa, qc); qi < qb; qi++) {
                        int rowStart = qi * 4;
                        int mr = Math.min(4, m - rowStart);
                        microKernel(
                                ap,
                                qi * 4 * jb,
                                ib,
                                jb,
                                a,
                                t + rowStart + colStart * lda,
                                lda,
                                mr,
                                nc,
                                (qi == qc) ? DIAG_LOWER : DIAG_NONE,
                                tmp);
                    }
                }
            });
        }
        return 0;
    }

    // ---------------- Cas A = U^T * U ----------------

    private static int potrfUpper(int n, double[] a, int lda, int nb, ForkJoinPool pool, double[] ap) {
        final int par = (pool == null) ? 1 : Math.max(1, pool.getParallelism());

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

            // 2) Panneau : A12 <- U11^{-T} * A12, puis packing des colonnes traitées.
            //    Micro-panneau q = colonnes 4q..4q+3 ; ap[q*4*jb + k*4 + s] = A12(k, 4q + s),
            //    complété par des zéros au-delà de la colonne m-1.
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
                for (int q = c0 / 4; q * 4 < c1; q++) {
                    int base = q * 4 * jb;
                    for (int s = 0; s < 4; s++) {
                        int col = q * 4 + s;
                        if (col < m) {
                            int src = p + col * lda;
                            for (int k = 0; k < jb; k++) {
                                ap[base + k * 4 + s] = a[src + k];
                            }
                        } else {
                            for (int k = 0; k < jb; k++) {
                                ap[base + k * 4 + s] = 0.0;
                            }
                        }
                    }
                }
            });

            // 3) Mise à jour : A22 <- A22 - A12^T * A12 (triangle supérieur).
            //    Chaque tâche traite une tranche de micro-panneaux de colonnes [qa, qb).
            final int nq = (m + 3) / 4;
            final int tileMp = Math.max(2, Math.min(16, (nq + 4 * par - 1) / (4 * par)));
            int nTiles = (nq + tileMp - 1) / tileMp;
            parallelFor(pool, nTiles, tile -> {
                double[] tmp = new double[16];
                int qa = tile * tileMp;
                int qb = Math.min(qa + tileMp, nq);
                for (int qr = 0; qr < qb; qr++) {
                    int rowStart = qr * 4;
                    int mr = Math.min(4, m - rowStart);
                    int ia = qr * 4 * jb;
                    for (int qc = Math.max(qa, qr); qc < qb; qc++) {
                        int colStart = qc * 4;
                        int nc = Math.min(4, m - colStart);
                        microKernel(
                                ap,
                                ia,
                                qc * 4 * jb,
                                jb,
                                a,
                                t0 + rowStart + colStart * lda,
                                lda,
                                mr,
                                nc,
                                (qr == qc) ? DIAG_UPPER : DIAG_NONE,
                                tmp);
                    }
                }
            });
        }
        return 0;
    }

    // ---------------- Micro-noyau 4x4 ----------------

    private static final int DIAG_NONE = 0; // bloc hors diagonale : on écrit tout le bloc utile
    private static final int DIAG_LOWER = 1; // bloc diagonal, triangle inférieur : on écrit r >= s
    private static final int DIAG_UPPER = 2; // bloc diagonal, triangle supérieur : on écrit r <= s

    /**
     * Calcule le produit 4x4  P(r, s) = sum_k ap[ia + 4k + r] * ap[ib + 4k + s]
     * (sommes partielles dans 16 variables locales) puis retranche P du bloc de A situé à
     * l'index cIdx : A(cIdx + r + s * lda) -= P(r, s).
     *
     * @param mr   nombre de lignes valides du bloc (1..4)
     * @param nc   nombre de colonnes valides du bloc (1..4)
     * @param diag DIAG_NONE, DIAG_LOWER ou DIAG_UPPER : seuls les éléments du triangle voulu
     *             sont écrits, l'autre triangle de A n'est jamais modifié
     * @param tmp  tampon de 16 valeurs pour les blocs de bord ou diagonaux
     */
    private static void microKernel(
            double[] ap,
            int ia,
            int ib,
            int jb,
            double[] a,
            int cIdx,
            int lda,
            int mr,
            int nc,
            int diag,
            double[] tmp) {
        double c00 = 0, c10 = 0, c20 = 0, c30 = 0;
        double c01 = 0, c11 = 0, c21 = 0, c31 = 0;
        double c02 = 0, c12 = 0, c22 = 0, c32 = 0;
        double c03 = 0, c13 = 0, c23 = 0, c33 = 0;

        int pa = ia;
        int pb = ib;
        for (int k = 0; k < jb; k++) {
            double a0 = ap[pa], a1 = ap[pa + 1], a2 = ap[pa + 2], a3 = ap[pa + 3];

            double b0 = ap[pb];
            c00 += a0 * b0;
            c10 += a1 * b0;
            c20 += a2 * b0;
            c30 += a3 * b0;
            double b1 = ap[pb + 1];
            c01 += a0 * b1;
            c11 += a1 * b1;
            c21 += a2 * b1;
            c31 += a3 * b1;
            double b2 = ap[pb + 2];
            c02 += a0 * b2;
            c12 += a1 * b2;
            c22 += a2 * b2;
            c32 += a3 * b2;
            double b3 = ap[pb + 3];
            c03 += a0 * b3;
            c13 += a1 * b3;
            c23 += a2 * b3;
            c33 += a3 * b3;

            pa += 4;
            pb += 4;
        }

        if (mr == 4 && nc == 4 && diag == DIAG_NONE) {
            // Cas courant : bloc complet hors diagonale
            int i1 = cIdx + lda;
            int i2 = i1 + lda;
            int i3 = i2 + lda;
            a[cIdx] -= c00;
            a[cIdx + 1] -= c10;
            a[cIdx + 2] -= c20;
            a[cIdx + 3] -= c30;
            a[i1] -= c01;
            a[i1 + 1] -= c11;
            a[i1 + 2] -= c21;
            a[i1 + 3] -= c31;
            a[i2] -= c02;
            a[i2 + 1] -= c12;
            a[i2 + 2] -= c22;
            a[i2 + 3] -= c32;
            a[i3] -= c03;
            a[i3 + 1] -= c13;
            a[i3 + 2] -= c23;
            a[i3 + 3] -= c33;
            return;
        }

        // Blocs de bord ou diagonaux : écriture avec masque
        tmp[0] = c00;
        tmp[1] = c10;
        tmp[2] = c20;
        tmp[3] = c30;
        tmp[4] = c01;
        tmp[5] = c11;
        tmp[6] = c21;
        tmp[7] = c31;
        tmp[8] = c02;
        tmp[9] = c12;
        tmp[10] = c22;
        tmp[11] = c32;
        tmp[12] = c03;
        tmp[13] = c13;
        tmp[14] = c23;
        tmp[15] = c33;

        for (int s = 0; s < nc; s++) {
            int rs;
            int re;
            if (diag == DIAG_LOWER) {
                rs = s;
                re = mr;
            } else if (diag == DIAG_UPPER) {
                rs = 0;
                re = Math.min(mr, s + 1);
            } else {
                rs = 0;
                re = mr;
            }
            int col = cIdx + s * lda;
            for (int r = rs; r < re; r++) {
                a[col + r] -= tmp[r + 4 * s];
            }
        }
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
    // Tests et mesures
    // ------------------------------------------------------------------

    /** Matrice symétrique à diagonale strictement dominante => définie positive. */
    private static double[] randomSpd(int n, long seed) {
        Random rnd = new Random(seed);
        double[] a = new double[n * n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i <= j; i++) {
                double v = rnd.nextDouble() * 2 - 1;
                a[i + j * n] = v;
                a[j + i * n] = v;
            }
        }
        for (int i = 0; i < n; i++) a[i + i * n] = n + 1;
        return a;
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

    /** Meilleur temps (ms) sur reps exécutions ; le dernier facteur calculé est copié dans out. */
    private static double bestTimeMs(
            boolean lower, int n, double[] a, int nb, ForkJoinPool pool, int reps, double[] out) {
        double best = Double.MAX_VALUE;
        for (int rep = 0; rep < reps; rep++) {
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

    /** Vérifie le noyau sur de nombreuses tailles (bords non multiples de 4) et tailles de bloc. */
    private static void selfTest() {
        int[] sizes = {1, 2, 3, 5, 17, 63, 64, 65, 100, 257, 300, 513};
        int[] nbs = {8, 16, 64};
        ForkJoinPool pool = new ForkJoinPool(Math.max(2, Runtime.getRuntime().availableProcessors()));
        double worst = 0;
        try {
            for (int n : sizes) {
                double[] a = randomSpd(n, 1000 + n);
                for (boolean lower : new boolean[] {true, false}) {
                    double[] ref = a.clone();
                    if (dpotrf(lower, n, ref, n, n, null) != 0) { // nb >= n : non bloqué
                        throw new IllegalStateException("référence en échec, n = " + n);
                    }
                    for (int nb : nbs) {
                        double[] f = a.clone();
                        int info = dpotrf(lower, n, f, n, nb, pool);
                        if (info != 0) {
                            throw new IllegalStateException("info = " + info + " (n=" + n + ", nb=" + nb + ")");
                        }
                        double d = maxDiff(lower, n, f, ref);
                        worst = Math.max(worst, d);
                        if (d > 1e-9) {
                            throw new IllegalStateException(
                                    "écart " + d + " (n=" + n + ", nb=" + nb + ", " + (lower ? "L" : "U") + ")");
                        }
                    }
                }
            }

            // Matrice non définie positive : le pivot k doit être signalé (info = k + 1)
            int n = 300;
            int k = 150;
            for (boolean lower : new boolean[] {true, false}) {
                double[] d = new double[n * n];
                for (int i = 0; i < n; i++) d[i + i * n] = 1.0;
                d[k + k * n] = -1.0;
                int info = dpotrf(lower, n, d, n, 64, pool);
                if (info != k + 1) {
                    throw new IllegalStateException("info attendu " + (k + 1) + ", obtenu " + info);
                }
            }
        } finally {
            pool.shutdown();
        }
        System.out.println("Auto-test OK (écart max vs non bloqué = " + worst + ")\n");
    }

    private static void benchmark() {
        int n = 601; // volontairement non multiple de 4
        int threads = Runtime.getRuntime().availableProcessors();
        double gflop = (double) n * n * n / 3.0 / 1e9;

        double[] a = randomSpd(n, 42);
        Random rnd = new Random(7);
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

                double tRef = bestTimeMs(lower, n, a, n, null, 1, fRef); // non bloqué
                double tSeq = bestTimeMs(lower, n, a, DEFAULT_NB, seqPool, 3, fSeq); // noyau, 1 thread
                double tPar = bestTimeMs(lower, n, a, DEFAULT_NB, parPool, 3, fPar); // noyau, parallèle

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
                        Locale.FRANCE, "  non bloqué          : %8.1f ms (%.2f GFLOP/s)%n", tRef, gflop / (tRef / 1e3));
                System.out.printf(
                        Locale.FRANCE, "  noyau, 1 thread     : %8.1f ms (%.2f GFLOP/s)%n", tSeq, gflop / (tSeq / 1e3));
                System.out.printf(
                        Locale.FRANCE,
                        "  noyau, %2d threads   : %8.1f ms (%.2f GFLOP/s), accélération x%.2f%n",
                        threads,
                        tPar,
                        gflop / (tPar / 1e3),
                        tSeq / tPar);
                System.out.println("  écart vs non bloqué = " + maxDiff(lower, n, fPar, fRef));
                System.out.println("  résidu max |A x - b| = " + maxRes);

                // Balayage de la taille de bloc (version parallèle)
                StringBuilder sb = new StringBuilder("  balayage de nb :");
                double[] tmp = new double[n * n];
                for (int nb : new int[] {32, 64, 96, 128, 192}) {
                    double t = bestTimeMs(lower, n, a, nb, parPool, 2, tmp);
                    sb.append(String.format(Locale.ROOT, " nb=%d: %.0f ms |", nb, t));
                }
                System.out.println(sb);
            }
        } finally {
            seqPool.shutdown();
            parPool.shutdown();
        }
    }

    public static void main(String[] args) {
        selfTest();
        benchmark();
    }
}
