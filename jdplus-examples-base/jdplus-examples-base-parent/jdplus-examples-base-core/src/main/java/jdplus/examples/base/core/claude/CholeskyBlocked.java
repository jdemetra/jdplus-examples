package jdplus.examples.base.core.claude;

import java.util.Locale;
import java.util.Random;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.SymmetricMatrix;

/**
 * Décomposition de Cholesky PAR BLOCS et résolution de systèmes, en Java pur,
 * avec les conventions de LAPACK (dpotrf / dpotrs) : matrices stockées par
 * colonnes (column-major), élément (i, j) à l'index i + j * ld.
 *
 * Algorithme de dpotrf : right-looking par blocs de nb colonnes. pour chaque
 * bloc diagonal : 1) factorisation non bloquée du bloc diagonal (dpotf2) 2)
 * résolution triangulaire du panneau (dtrsm) 3) mise à jour du bloc restant,
 * par tuiles (dsyrk)
 */
public final class CholeskyBlocked {

    /**
     * Taille de bloc par défaut (à ajuster selon le cache de la machine : 32,
     * 64, 128).
     */
    public static final int DEFAULT_NB = 32;

    private static final int ROW_TILE = 256; // tuile en lignes pour la mise à jour (cas inférieur)
    private static final int COL_TILE = 64; // tuile en colonnes pour la mise à jour (cas supérieur)

    private CholeskyBlocked() {}

    // ------------------------------------------------------------------
    // dpotrf
    // ------------------------------------------------------------------
    public static int dpotrf(boolean lower, int n, double[] a, int lda) {
        return dpotrf(lower, n, a, lda, DEFAULT_NB);
    }

    /**
     * Équivalent de LAPACK dpotrf avec taille de bloc nb. Calcule A = L * L^T
     * (lower = true) ou A = U^T * U (lower = false), sur place. Seul le
     * triangle demandé est lu et écrasé.
     *
     * @return 0 si succès ; k > 0 si le mineur principal d'ordre k n'est pas
     * strictement positif (matrice non définie positive)
     */
    public static int dpotrf(boolean lower, int n, double[] a, int lda, int nb) {
        checkArgs(n, a, lda, "a");
        if (nb < 1) {
            throw new IllegalArgumentException("nb < 1");
        }

        if (nb >= n) {
            return lower ? potf2Lower(a, lda, 0, n) : potf2Upper(a, lda, 0, n);
        }
        return lower ? potrfLower(n, a, lda, nb) : potrfUpper(n, a, lda, nb);
    }

    // ---------------- Cas A = L * L^T ----------------
    private static int potrfLower(int n, double[] a, int lda, int nb) {
        for (int j = 0; j < n; j += nb) {
            int jb = Math.min(nb, n - j);
            int jj = j + j * lda; // index de A(j, j)

            // 1) Bloc diagonal
            int info = potf2Lower(a, lda, jj, jb);
            if (info != 0) {
                return j + info;
            }

            int m = n - j - jb; // nombre de lignes sous le bloc
            if (m <= 0) {
                continue;
            }

            int p = (j + jb) + j * lda; // index de A21(0, 0)

            // 2) Panneau : A21 <- A21 * L11^{-T}
            for (int c = 0; c < jb; c++) {
                int pc = p + c * lda;
                double d = a[jj + c + c * lda];
                for (int i = 0; i < m; i++) {
                    a[pc + i] /= d;
                }
                for (int c2 = c + 1; c2 < jb; c2++) {
                    double l = a[jj + c2 + c * lda]; // L11(c2, c)
                    int pc2 = p + c2 * lda;
                    for (int i = 0; i < m; i++) {
                        a[pc2 + i] -= a[pc + i] * l;
                    }
                }
            }

            // 3) Mise à jour : A22 <- A22 - A21 * A21^T (triangle inférieur), par tuiles de lignes
            int t = (j + jb) + (j + jb) * lda; // index de A22(0, 0)
            for (int i0 = 0; i0 < m; i0 += ROW_TILE) {
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
            }
        }
        return 0;
    }

    // ---------------- Cas A = U^T * U ----------------
    private static int potrfUpper(int n, double[] a, int lda, int nb) {
        for (int j = 0; j < n; j += nb) {
            int jb = Math.min(nb, n - j);
            int jj = j + j * lda;

            // 1) Bloc diagonal
            int info = potf2Upper(a, lda, jj, jb);
            if (info != 0) {
                return j + info;
            }

            int m = n - j - jb; // nombre de colonnes à droite du bloc
            if (m <= 0) {
                continue;
            }

            int p = j + (j + jb) * lda; // index de A12(0, 0)

            // 2) Panneau : A12 <- U11^{-T} * A12 (substitution avant, colonne par colonne)
            for (int c = 0; c < m; c++) {
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

            // 3) Mise à jour : A22 <- A22 - A12^T * A12 (triangle supérieur), par tuiles de colonnes
            int t0 = (j + jb) + (j + jb) * lda; // index de A22(0, 0)
            for (int rb = 0; rb < m; rb += COL_TILE) {
                int re = Math.min(rb + COL_TILE, m);
                for (int c = rb; c < m; c++) {
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
        }
        return 0;
    }

    // ---------------- Factorisation non bloquée d'un bloc (dpotf2) ----------------
    /**
     * Cholesky inférieure non bloquée sur le bloc m x m dont A(0,0) est à
     * l'index base.
     */
    private static int potf2Lower(double[] a, int lda, int base, int m) {
        //        int dinc = lda + 1;
        //        for (int j = 0, cj = base, jj = 0; j < m; j++, cj += lda, jj += dinc) {
        //            int cjmax = cj + m;
        //            double d = a[jj];
        //            if (!(d > 0.0)) {
        //                return j + 1; // gère aussi NaN
        //            }
        //            for (int c = j + 1; c < m; c++) {
        //                int cc = base + c * lda;
        //                double ljc = a[cj + c] / d;
        //                for (int i = c; i < m; i++) {
        //                    a[cc + i] -= a[cj + i] * ljc;
        //                }
        //            }
        //            d = Math.sqrt(d);
        //            a[jj] = d;
        //            for (int i = jj + 1; i < cjmax; i++) {
        //                a[i] /= d;
        //            }
        //        }
        //        int dinc = lda + 1;
        //        for (int j = 0, jj = base, cj = base, rj = base; j < m; j++, jj += dinc, cj += lda, rj++) {
        //            // column j
        //            // diagonal element
        //            double d = a[jj];
        //            for (int ji = rj; ji < jj; ji += lda) {
        //                double x = a[ji];
        //                d -= x * x;
        //            }
        //            if (!(d > 0.0)) {
        //                return j + 1; // gère aussi NaN
        //            }
        //
        //            d = Math.sqrt(d);
        //            a[jj] = d;
        //            for (int c = 0, cc = base; c < j; c++, cc += lda) {
        //                double temp = a[cc + j];
        //                for (int i = j + 1; i < m; i++) {
        //                    a[cj + i] -= a[cc + i] * temp;
        //                }
        //            }
        //            double scale = 1 / d;
        //            for (int i = j + 1; i < m; i++) {
        //                a[cj + i] *= scale;
        //            }
        //        }

        int dinc = 1 + lda, end = base + m * dinc;
        for (int idiag = base, irow = base, cend = base + m; idiag != end; ++irow, idiag += dinc, cend += lda) {
            // compute aii;
            double aii = a[idiag];
            for (int j = irow; j != idiag; j += lda) {
                double x = a[j];
                aii -= x * x;
            }
            if (aii < 0) { // negative
                return 1;
            } else {
                aii = Math.sqrt(aii);
                a[idiag] = aii;
                // compute elements i+1 : n of column i
                for (int jx = irow; jx != idiag; jx += lda) {
                    double temp = a[jx];
                    if (temp != 0) {
                        for (int ia = jx + 1, iy = idiag + 1; iy < cend; ++ia, ++iy) {
                            a[iy] -= temp * a[ia];
                        }
                    }
                }
                double scale = 1 / aii;
                for (int iy = idiag + 1; iy < cend; ++iy) {
                    a[iy] *= scale;
                }
            }
        }
        return 0;
    }

    /**
     * Cholesky supérieure non bloquée sur le bloc m x m dont A(0,0) est à
     * l'index base.
     */
    private static int potf2Upper(double[] a, int lda, int base, int m) {
        for (int j = 0; j < m; j++) {
            int cj = base + j * lda;
            double d = a[cj + j];
            if (!(d > 0.0)) {
                return j + 1;
            }
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
    /**
     * Équivalent de LAPACK dpotrs. Résout A * X = B, où A = L * L^T (lower =
     * true) ou A = U^T * U (lower = false) a été factorisée par dpotrf. B (n x
     * nrhs) est écrasé par la solution X.
     */
    public static void dpotrs(boolean lower, int n, int nrhs, double[] a, int lda, double[] b, int ldb) {
        checkArgs(n, a, lda, "a");
        if (nrhs < 0) {
            throw new IllegalArgumentException("nrhs < 0");
        }
        if (ldb < Math.max(1, n)) {
            throw new IllegalArgumentException("ldb < max(1, n)");
        }
        if (b.length < ldb * nrhs) {
            throw new IllegalArgumentException("b trop petit");
        }

        for (int k = 0; k < nrhs; k++) {
            int o = k * ldb;

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
    }

    /**
     * Met à zéro le triangle non utilisé, pour obtenir une vraie matrice L (ou
     * U).
     */
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
        if (n < 0) {
            throw new IllegalArgumentException("n < 0");
        }
        if (lda < Math.max(1, n)) {
            throw new IllegalArgumentException("lda < max(1, n)");
        }
        if (a.length < lda * n) {
            throw new IllegalArgumentException(name + " trop petit");
        }
    }

    // ------------------------------------------------------------------
    // Test : version par blocs vs version non bloquée, résidu, temps
    // ------------------------------------------------------------------
    public static void main(String[] args) {
        int N = 50, K = 100000;
        Random rnd = new Random(42);

        // Matrice symétrique à diagonale strictement dominante => définie positive
        double[] a = new double[N * N];
        for (int j = 0; j < N; j++) {
            for (int i = 0; i <= j; i++) {
                double v = rnd.nextDouble() * 2 - 1;
                a[i + j * N] = v;
                a[j + i * N] = v;
            }
        }
        for (int i = 0; i < N; i++) {
            a[i + i * N] = N;
        }

        double[] b = new double[N];
        for (int i = 0; i < N; i++) {
            b[i] = rnd.nextDouble();
        }

        for (int k = 0; k < 3; ++k) {
            for (boolean lower : new boolean[] {false, true}) {
                double[] fBlocked = null;
                double[] fRef = null;
                double[] fJD = null;

                int info1 = 0, info2 = 0;

                long t0 = System.nanoTime();
                for (int l = 0; l < K; ++l) {
                    fBlocked = a.clone();
                    info1 = dpotrf(lower, N, fBlocked, N, 32);
                }
                long t1 = System.nanoTime();
                for (int l = 0; l < K; ++l) {
                    fRef = a.clone();
                    info2 = dpotrf(lower, N, fRef, N, N); // nb >= n : version non bloquée
                }
                long t2 = System.nanoTime();
                for (int l = 0; l < K; ++l) {
                    fJD = a.clone();
                    FastMatrix L = FastMatrix.builder(fJD)
                            .nrows(N)
                            .ncolumns(N)
                            .columnIncrement(N)
                            .build();
                    if (lower) {
                        SymmetricMatrix.lcholesky(L);
                    } else {
                        SymmetricMatrix.ucholesky(L);
                    }
                }
                long t3 = System.nanoTime();

                if (info1 != 0 || info2 != 0) {
                    throw new ArithmeticException("Non définie positive : " + info1 + " / " + info2);
                }
                if (k == 2) {
                    // Écart entre factorisations (triangle utile uniquement)
                    double maxDiff = 0;
                    for (int j = 0; j < N; j++) {
                        for (int i = 0; i < N; i++) {
                            if ((lower && i >= j) || (!lower && i <= j)) {
                                maxDiff = Math.max(maxDiff, Math.abs(fBlocked[i + j * N] - fRef[i + j * N]));
                            }
                        }
                    }
                    double maxDiff2 = 0;
                    for (int j = 0; j < N; j++) {
                        for (int i = 0; i < N; i++) {
                            if ((lower && i >= j) || (!lower && i <= j)) {
                                maxDiff2 = Math.max(maxDiff2, Math.abs(fRef[i + j * N] - fJD[i + j * N]));
                            }
                        }
                    }

                    // Résidu de la résolution A x = b
                    double[] x = b.clone();
                    dpotrs(lower, N, 1, fRef, N, x, N);
                    double maxRes = 0;
                    for (int i = 0; i < N; i++) {
                        double s = 0;
                        for (int j = 0; j < N; j++) {
                            s += a[i + j * N] * x[j];
                        }
                        maxRes = Math.max(maxRes, Math.abs(s - b[i]));
                    }

                    System.out.println((lower ? "Inférieur (L)" : "Supérieur (U)"));
                    System.out.printf(
                            Locale.FRANCE,
                            "  bloc : %.1f ms | non bloqué : %.1f ms | JD+ : %.1f ms",
                            (t1 - t0) / 1e6,
                            (t2 - t1) / 1e6,
                            (t3 - t2) / 1e6);
                    System.out.println("  écart bloc vs non bloqué = " + maxDiff);
                    System.out.println("  écart JD vs non bloqué = " + maxDiff2);
                    System.out.println("  résidu max |A x - b|     = " + maxRes);
                }
            }
        }
    }
}
