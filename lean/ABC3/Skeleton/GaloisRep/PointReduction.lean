import ABC3.Found.GaloisRep.PullbackPoint

/-!
# スケルトン —— **点の還元と群法則**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★D2 に残る唯一の段

第 151 ブロックで**因子の台がイデアルの等式として書けた**:

    count_v(μ f_P) ≠ 0   ⟺   I_P = P'        ただし P' := { a | w(μ a) < 1 }

★残るのは **`P'` が定める点が `n·Q_v` であること**である。
すなわち `v.asIdeal = XYIdeal(Q_v)` のとき

    P' = XYIdeal(n·Q_v)

★★これは**「点の還元が群準同型である」**ことの言い換えである
(Silverman, *The Arithmetic of Elliptic Curves*, VII.2)。

## ★★★★機構と、なぜ mathlib で済まないか

`μ` は `[n]` の引き戻しで `μ(x) = x([n]·)`、`μ(y) = y([n]·)`。
★`P' = μ⁻¹(m_v)` は「`[n]` を施してから `v` で消える」元の全体だから、
幾何的には `[n]Q_v` の極大イデアルである。

★★2026-08-20 実測: mathlib の `Point.map` は
**体の間の単射代数準同型** `f : F →ₐ[S] K` にしか使えない
(`W'.baseChange_nonsingular f.injective` を使う)。
還元 `O_v → κ(v)` は**単射でない**ので適用できない。
★★★`AlgebraicGeometry/EllipticCurve/Reduction.lean` も**曲線**の還元
(極小モデル・good/multiplicative/additive reduction)であって点の還元ではない。

## ★★★良い還元であることは自動である

`Δ ∈ F^×`(定数)なので `w(Δ) = 1`、すなわち **`v` では常に良い還元**である。
★悪い還元の場合分けは不要。★★残るのは加法公式の場合分け
(`x₁ ≡ x₂ mod m_v` のとき傾き `(y₂−y₁)/(x₂−x₁)` が `O_v` に入らない)である。

## ★★★★この節点は (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。
★ここで積んだものはそのまま流用できる。
-/

namespace ABC3.Skeleton.GaloisRep

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve WeierstrassCurve.Affine Polynomial
open IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F]

/-- ★★★★★**引き戻した素イデアルは `n·Q_v` の極大イデアルである**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが D2 に残る唯一の段である。★★「点の還元が群準同型である」ことの言い換え。
★★★`hnn`(場合 A)の下では `n·Q_v ≠ 0` なのでアフィン点として書ける。 -/
theorem pullbackPrime_eq_xyIdeal_nsmul [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    [inst : IsDedekindDomain W.CoordinateRing]
    (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hnn : ∀ r : W.CoordinateRing, v.valuation W.FunctionField (μ r) ≤ 1)
    (n : ℕ) (hn : 1 ≤ n)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    {c y₀ : F} (hQ : W.Nonsingular c y₀)
    (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀)) :
    ∃ (x' y' : F) (h' : W.Nonsingular x' y'),
      Point.some x' y' h' = n • Point.some c y₀ hQ ∧
      pullbackPrime (v.valuation W.FunctionField) μ hnn
        = CoordinateRing.XYIdeal W x' (Polynomial.C y') :=
  exists_pullbackPrime_eq_xyIdeal W v hQ.1 hv h2 n μ hμF hns hμx hμy hμP hnn

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def pullbackPrime_eq_xyIdeal_nsmul.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——引き戻した素イデアルが n·Q の極大イデアルであること)",
    sectionId := "genell-thm-3-8" }

def pullbackPrime_eq_xyIdeal_nsmul.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★★**2026-08-20: 本節点は証明された**(第 162-164 ブロック)。第 162 で還元の核 `E₁` が加法で閉じていること(Vieta の第 2 基本対称式で打ち消しを回避)、第 163 で `redPoint` が群準同型であること(`redHom`)、第 164 で `P' = XYIdeal(n·Q_v)` が出た。★形式群も値群の構造も使っていない(0 ブロック)" 19,
    .citation "[Silverman]" "The Arithmetic of Elliptic Curves, VII.2(還元写像が群準同型であること)"
      (.absent "mathlib に点の還元写像は無い(2026-08-20 実測)。`Point.map` は体の間の単射代数準同型にしか使えず、`Reduction.lean` は曲線の還元(極小モデル)である") 19,
    .implicitStep
      "★★★★★★**台の同定は第 151 ブロックで済んだ**——`count_v(μ f_P) ≠ 0 ⟺ I_P = P'`(`count_ne_zero_iff`)。剰余体も点の還元も経由せず、イデアルの等式として書けた。★§9-461 の見積もり 8-20 ブロックは 1 ブロックで済んだ(0 ブロック)" 19,
    .implicitStep
      "★★★★**良い還元であることは自動である**——`Δ ∈ F^×`(定数)なので `w(Δ) = 1`。悪い還元の場合分けは不要である(0 ブロック)" 19,
    .implicitStep
      "★★★★★**剰余体の足場は mathlib にあった**(2026-08-20 実測)——`Ideal.ResidueField I := IsLocalRing.ResidueField (Localization.AtPrime I)` と `Ideal.ker_algebraMap_residueField : ker (algebraMap R I.ResidueField) = I`、さらに `IsFractionRing (R ⧸ I) I.ResidueField` の instance。★`v.asIdeal` は極大なので `R ⧸ v.asIdeal` は既に体であり、`κ(v) ≅ R ⧸ v.asIdeal ≅ F`(第 138 の `quotientXYIdealEquiv`)となる(2-5 ブロック)" 19,
    .implicitStep
      "★★★★**剰余体を経由しない定式化も可能である**——`t ∈ O_v` に対し `w(t − c) < 1` なる `c ∈ F` は**一意**である(`c − c' ∈ F` の付値は 1 か 0 だから)。存在は `κ(v) ≅ F` から出る。★この形なら還元写像を `FF ⊇ O_v → F` として直接書ける(3-8 ブロック)" 19,
    .implicitStep
      "★★★残るのは加法公式の場合分けである。`x₁ ≢ x₂ mod m_v` なら傾き `(y₂−y₁)/(x₂−x₁)` は `O_v` に入り、還元は素直に可換になる。★問題は `x₁ ≡ x₂ mod m_v` の場合で、そこでは還元先で 2 倍公式に切り替わる(15-30 ブロック)" 19,
    .implicitStep
      "★★`E(F(W))₀ := {S | S = 0 または座標が O_v に入る}` が部分群であることも要る。射影モデルを使わずにやるなら、加法公式の場合分けの中で示すことになる(5-15 ブロック)" 19,
    .implicitStep
      "★★★★**この節点は (G7) 半安定モデルでも効く**。(G7) も点の還元を要求するので、ここで積んだものはそのまま流用できる(0 ブロック)" 19,
    .implicitStep
      "★★★★★★これが済めば D2 は閉じる——第 138(高さ 1 素点 ↔ 点)、第 140(類が自明 ⟹ 単項)、第 150(`Σ_{T ∈ E[n]} T = 0`、ファイバー総和)、第 151(台の同定)がすべて揃っている。★`Point.toClass_eq_zero`(mathlib)で結論は 1 行(4-10 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep
