import ABC3.Skeleton.GaloisRep.WeilFunctionField
import ABC3.Found.GaloisRep.RootFromIdeal

/-!
# スケルトン —— **`div(f_P ∘ [n])` の因子計算**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★層 3 が 2 節点に還元された

第 139 ブロック(`Found`)で最後の一歩

    J^n = (μ f_P)   かつ   J = (g)     ⟹   ∃ g', g'^n = μ f_P

が取れた。★したがって `g_P` の存在に残るのは**この 2 つだけ**である:

| 節点 | 内容 | 見積もり |
|---|---|---|
| `exists_fractionalIdeal_pow` (D1) | `∃ J, J^n = (μ f_P)`——因子が `n` で割れる | 20-50 |
| `fractionalIdeal_isPrincipal` (D2) | その `J` が単項——Abel–Jacobi | 10-25 |

## ★★★★★★見通しが変わった——**不分岐性は要らない**

当初は「`[n]` は不分岐だから `ord_Q(μ f) = ord_{[n]Q}(f)`」を経由する予定で、
不分岐性の証明(`deg[n] = n²` と `#E[n] = n²` が要る)を勘定に入れていた。
★**それは不要である。**一般の付値論で

    ord_Q(μ f) = e_Q · ord_{[n]Q}(f)

であり、`div(f_P) = n(P) − n(O)` から `ord_{[n]Q}(f_P) ∈ {n, −n, 0}`。
★★したがって **`e_Q` が何であっても `n` で割れる**。これが D1 である。

★★★D2 でも `e_Q` は消える。平行移動 `τ_T`(`T ∈ E[n]`)が `μ` を保つ
(`[n]∘τ_T = [n]`)ので `e_Q` は各ファイバー上で一定であり、共通因子 `e` を括り出すと

    Σ_{Q ∈ [n]⁻¹(P)} Q = n²Q₀ = n·P = 0,      Σ_{T ∈ E[n]} T = 0

から類は 0 になる。★★★★**`e` の値を知る必要はない。**

## ★★★★因子を幾何で読む準備は済んでいる

| 材料 | ブロック |
|---|---|
| `IsDedekindDomain F[W]`(因子機構が開く) | 第 137 |
| `HeightOneSpectrum F[W] ↔ 曲線上の点` | 第 138 |
| `(f_P) = XYIdeal_P^n` | 第 126 |
| `μ`(群法則で固定された `[n]` の引き戻し) | 第 118・125 |
| 平行移動 `τ_T`(体の自己同型) | 第 120・124 |
| 単元は定数 | 第 128 |
-/

namespace ABC3.Skeleton.GaloisRep

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve WeierstrassCurve.Affine nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F]

/-- ★★★★★**D1 —— `(μ f_P)` は `n` 乗である**(因子が `n` で割れる)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`ord_Q(μ f_P) = e_Q · ord_{[n]Q}(f_P)` で右の因子が `{n, −n, 0}` のいずれかだから、
**分岐指数 `e_Q` を知らなくても** `n` で割れる。 -/
theorem exists_fractionalIdeal_pow (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x y : F}
    (h : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n) (hP : n • (Point.some x y h) = 0)
    (μ : W.CoordinateRing →+* W.FunctionField)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    ∃ J : FractionalIdeal W.CoordinateRing⁰ W.FunctionField,
      J ^ n = FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP) := by
  sorry

/-- ★★★★★**D2 —— その `J` は単項である**(Abel–Jacobi)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Point.toClass` の単射性(mathlib)で、類が 0 であることから単項性が出る。
★★類の計算では `Σ_{Q ∈ [n]⁻¹(P)} Q = 0` と `Σ_{T ∈ E[n]} T = 0` を使う。 -/
theorem fractionalIdeal_isPrincipal (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x y : F}
    (h : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n) (hP : n • (Point.some x y h) = 0)
    (μ : W.CoordinateRing →+* W.FunctionField)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    (J : FractionalIdeal W.CoordinateRing⁰ W.FunctionField)
    (hJ : J ^ n = FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) :
    ∃ g : W.FunctionField, J = FractionalIdeal.spanSingleton W.CoordinateRing⁰ g := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_fractionalIdeal_pow.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——div(f_P ∘ [n]) が n で割れること)",
    sectionId := "genell-thm-3-8" }

def exists_fractionalIdeal_pow.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1 の証明(div(f_P∘[n]) の計算)"
      (.absent "mathlib に Weil 対およびその構成要素は 0 件(2026-08-20、`WeilPairing|weil_pairing` で全文検索して 0 件)") 19,
    .implicitStep
      "★★★★★★**不分岐性は要らない**(2026-08-20 の見通し変更)。一般の付値論で `ord_Q(μ f) = e_Q · ord_{[n]Q}(f)` であり、`div(f_P) = n(P) − n(O)` から右の因子は `{n, −n, 0}` のいずれか。★したがって `e_Q` が何であっても `n` で割れる。`deg[n] = n²` も `#E[n] = n²` も**この節点には不要**である" 19,
    .implicitStep
      "★★★★因子機構は **mathlib に完備**である——`IsDedekindDomain.HeightOneSpectrum`・`intValuation`・`Ideal.finprod_heightOneSpectrum_factorization`・adic 付値。`IsDedekindDomain F[W]` は**第 137 ブロックで取れた**ので、そのまま使える(0 ブロック)" 19,
    .implicitStep
      "★★★因子を幾何の言葉(曲線上の点)で読む対応は**第 138 ブロックで取れた**——`HeightOneSpectrum F[W] ↔ { (c,y₀) : W.Equation c y₀ }`(0 ブロック)" 19,
    .implicitStep
      "★★★`ord_Q(μ f) = e_Q · ord_{v'}(f)` の段。`μ` を体の埋め込み `μ̃ : F(W) → F(W)` に延ばし、`v` の下にある `v' := μ̃⁻¹(m_v) ∩ F[W]` を取る。mathlib の `IsDedekindDomain.HeightOneSpectrum` と `Valuation.comap` で書けるはずだが、拡大の分岐指数 `e` を扱う API の在庫は未実測(15-40 ブロック)" 19,
    .implicitStep
      "★★`(f_P) = XYIdeal_P^n` は **第 126 ブロックで済**(`xyIdeal_pow_isPrincipal_integral`)。これが `ord_{v'}(f_P) ∈ {n, 0}` を与える(0 ブロック)" 19,
    .implicitStep
      "★`μ` は第 118・125 ブロックで**群法則によって固定された形で**得られている(`exists_mulByNHom_charZero`)(0 ブロック)" 19,
    .implicitStep
      "★★★因子が `n` で割れることから `J` を作る段。Dedekind 環では分数イデアルが素因子の自由アーベル群なので、指数を `n` で割ったものが `J` になる。mathlib の `FractionalIdeal` の因子分解 API の在庫は未実測(5-15 ブロック)" 19 ]

def fractionalIdeal_isPrincipal.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——因子の類が自明であること)",
    sectionId := "genell-thm-3-8" }

def fractionalIdeal_isPrincipal.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.3.5(Abel–Jacobi)"
      (.absent "mathlib に Weil 対およびその構成要素は 0 件(2026-08-20)") 19,
    .implicitStep
      "★★★★★★**分岐指数 `e_Q` はここでも消える**——平行移動 `τ_T`(`T ∈ E[n]`)は `[n]∘τ_T = [n]` により `μ` を保つので `e_Q` は各ファイバー上で一定。共通因子 `e` を括り出すと `Σ_{Q ∈ [n]⁻¹(P)} Q = n²Q₀ = n·P = 0` と `Σ_{T ∈ E[n]} T = 0` から類は 0 になる。★`e` の値を知る必要はない" 19,
    .implicitStep
      "★★★★`Point.toClass` の**単射性は mathlib にある**(`Point.toClass_injective`、2026-08-20 実測)。類が 0 なら `J` が単項であることはここから出る(2-5 ブロック)" 19,
    .implicitStep
      "★★★平行移動 `τ_T` が体の自己同型であることは **第 120・124 ブロックで済**(`translateAlgEquiv`・`exists_translateAut_all`)(0 ブロック)" 19,
    .implicitStep
      "★★`Σ_{T ∈ E[n]} T = 0` の段。`E[n] ≅ (ℤ/n)²`((G1) で取得済)から `Σ_{(a,b)} (a,b) = (n·Σa, n·Σb) = 0`。有限アーベル群の全元の和に関する mathlib の在庫は未実測(3-8 ブロック)" 19,
    .implicitStep
      "★★★ファイバー `[n]⁻¹(P) = {Q₀ + T : T ∈ E[n]}` の記述と、`Σ_{Q ∈ ファイバー} Q = n²Q₀ + Σ_T T = n·P + 0 = 0` の計算(5-12 ブロック)" 19,
    .implicitStep
      "★★アフィン座標環の類群は `Pic⁰(E)` である(無限遠点 `O` を落とすと次数が消える)。`toClass` が `(R) − (O)` の類を与えることを使って、因子の類を点の和として読む(3-8 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep
