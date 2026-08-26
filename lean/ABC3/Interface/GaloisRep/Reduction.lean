import ABC3.Interface.GaloisRep.Representation
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.AlgebraicGeometry.EllipticCurve.Reduction
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.GroupTheory.QuotientGroup.Basic

/-!
# Galois 表現のスケルトン(3/3)—— **還元・Tate 曲線・Faltings 高さ**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.15–p.17。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★本ファイルが受ける 3 件

| # | 受けるもの | mathlib の在庫(2026-08-17 実測) |
|---|---|---|
| G6 | Tate 曲線と局所高さ `v(q_E)` | ★**無い**。原文の典拠は [FC] Ch. III, Cor. 7.3 |
| G7 | 半安定還元と `𝓞_L` 上のモデル | ★`EllipticCurve/Reduction.lean` は**ある**が Néron モデルは無い |
| G8 | Faltings 高さ `ht^Falt = deg(ω_E)` | ★**無い**。Arakelov 側 (D2) に従属 |

★★★**G6 が `Definition 3.3` そのもの**であり、`Lemma 3.2` が要求するものである。
★★G8 は **Arakelov 側と Galois 側の合流点**である——
`deg(ω_E)` は算術直線束の次数だから、(D2) が要る。
-/

namespace ABC3.Interface.GaloisRep

open ABC3.Meta WeierstrassCurve NumberField

/-! ## ★★★G6 —— Tate 曲線と局所高さ -/

/-- **(G6)** 乗法還元をもつ楕円曲線の **Tate 母数 `q_E`** と**局所高さ** `v_K(q_E)`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★★**原文の `Definition 3.3` そのものである。**
★`Lemma 3.2`(局所の階数 1 部分群)がこの上で述べられる。

★★★**2026-08-17: 退化検査で作り直した。**以前は `LocalField : Type` /
`Curve : LocalField → Type` と**世界ごと posit** していたので、
`PUnit` と定数 `1` で埋まってしまった(実測済)。
★**mathlib の `WeierstrassCurve` と `Valuation` に接地し**、
**Tate 一意化そのもの**を要求する形に改めた。 -/
structure TateCurveData where
  /-- ★★**Tate 母数** `q_E ∈ Kˣ`。

  ★★★**述語は posit しない**——mathlib の
  `WeierstrassCurve.HasSplitMultiplicativeReduction` に接地する。
  ★これで「`HasSplitMultRed := False`」という**空虚 witness も不可能**になる。 -/
  tateParam : {R : Type} → [CommRing R] → [IsDomain R] → [IsDiscreteValuationRing R] →
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R] →
    {K : Type} → [Field K] → [Algebra R K] → [IsFractionRing R K] →
    (W : WeierstrassCurve K) → [W.IsElliptic] → [W.IsMinimal R] →
    W.HasSplitMultiplicativeReduction R → Kˣ
  /-- ★★★**Tate 一意化** `E(K) ≅ Kˣ / q^ℤ`。

  ★★★**これが内容の本体である。**`q` が実際に一意化していることを要求する。 -/
  uniformization : ∀ {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [DecidableEq K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R),
    Nonempty (W.toAffine.Point ≃+
      Additive (Kˣ ⧸ Subgroup.zpowers (tateParam W h : Kˣ)))
  /-- ★★★**局所高さ** `v_K(q_E) ∈ ℤ_{>0}`(原文 `Definition 3.3`)。

  ★★★★**2026-08-26 の訂正(逸脱の記録)**——以前は付値
  `v : Kˣ →* Multiplicative ℤ` を**任意に受けていた**。しかしそれだと
  `v := 1`(自明な準同型、`toAdd (1 q) = 0`)を入れたとき
  `localHeight_eq` が `localHeight = 0` を、`localHeight_pos` が `0 < localHeight` を
  要求して**衝突する**——すなわち**構造そのものが充足不能**だった。
  ★離散付値環の**正規化付値**(mathlib の `IsDiscreteValuationRing.maximalIdeal` の
  adic 付値)に固定して直す。★★これは**弱めるのではなく強める**訂正である
  (`localHeight` が `q` から一意に決まる)。

  ★★★★**2026-08-26 の第 2 の訂正**——`Δ ≠ 0`(`IsElliptic`)を要求する。
  mathlib の `HasMultiplicativeReduction` は `v(Δ) < 1` としか言わず、
  **`Δ = 0` を排除しない**(`v(0) = 0 < 1`)。実際 `y² + xy = x³` は
  `Δ = 0`・`c₄ = 1`・接線の 2 次式 `X² + X` が分解、で**すべての仮定を満たす**が、
  これは**結節三次曲線**であり滑らかな点の群は `Kˣ` そのもの——
  Tate 母数は `q = 1` にあたるので**局所高さが `0`** になり `localHeight_pos` と衝突する。 -/
  localHeight : {R : Type} → [CommRing R] → [IsDomain R] → [IsDiscreteValuationRing R] →
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R] →
    {K : Type} → [Field K] → [Algebra R K] → [IsFractionRing R K] →
    (W : WeierstrassCurve K) → [W.IsElliptic] → [W.IsMinimal R] →
    W.HasSplitMultiplicativeReduction R → ℕ
  /-- ★★★**局所高さは Tate 母数の付値そのもの**——定義を型に出す。

  ★正規化付値では `v(q) = ofAdd (−localHeight)`(乗法的な書き方)である。
  ★★これがあるので `localHeight := 1` のような定数 witness は通らない。 -/
  localHeight_eq : ∀ {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R),
    IsDedekindDomain.HeightOneSpectrum.valuation K (IsDiscreteValuationRing.maximalIdeal R)
        ((tateParam W h : Kˣ) : K)
      = (Multiplicative.ofAdd (-(localHeight W h : ℤ)) : Multiplicative ℤ)
  /-- ★局所高さは正(原文が `∈ Z_{>0}` と書いている)。 -/
  localHeight_pos : ∀ {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R), 0 < localHeight W h

def TateCurveData.waiting : WaitingFor :=
  { what := "(G6) Tate 曲線 E(K̄) ≅ K̄^x/q^Z、Tate 母数 q_E、および局所高さ v_K(q_E)(= [GenEll] Definition 3.3)"
    trackB := "Found/GaloisRep — ★mathlib に Tate 曲線は無い(2026-08-16 実測)。★★**FLT にも無い**——`FLT/TateCurve/TateCurve.lean` は 20 行の入口宣言だけ(2026-08-17、clone して計数)。★★★原文の典拠は [FC](Faltings-Chai, *Degenerations of Abelian Varieties*)Chapter III, Corollary 7.3 で、これも mathlib に無い" }

/-! ## ★★G7 —— 半安定還元とモデル -/

/-- **半安定**であること —— 離散付値環 `R` の上で**良還元または乗法還元**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★★**posit ではなく mathlib の還元型による定義である。**
`Reduction.lean` の三分律(良・乗法・加法)のうち加法還元を除いたものが半安定。 -/
def IsSemiStableAt (R : Type) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsMinimal R] : Prop :=
  W.HasGoodReduction R ∨ W.HasMultiplicativeReduction R

/-- ★半安定とは「加法還元でない」ことである(三分律から)。 -/
theorem isSemiStableAt_iff_not_additive (R : Type) [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsMinimal R] :
    IsSemiStableAt R W ↔ ¬ W.HasAdditiveReduction R := by
  constructor
  · rintro (h | h)
    · exact h.not_hasAdditiveReduction
    · exact h.not_hasAdditiveReduction
  · intro h
    rcases WeierstrassCurve.hasGoodReduction_or_hasMultiplicativeReduction_or_hasAdditiveReduction
      R (W := W) with hg | hm | ha
    · exact Or.inl hg
    · exact Or.inr hm
    · exact absurd ha h

/-- **(G7)** 半安定還元と、`𝓞_L` 上への延長が与える **Néron 微分** `ω_E`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★原文は「`E` は `L` のすべての有限素点で半安定還元をもつ」を繰り返し仮定し、
そこから `𝓞_L` 上の 1 次元半アーベルスキームへの延長を得る。

★★★**その延長が使われるのは `ω_E` を作るためだけ**である
(`Proposition 3.4` の `ht^Falt = deg(ω_E)`)。ゆえに `ω_E` を型に出す。

★★★**`ω_E` は階数 1** —— これが退化(`PUnit` は階数 0)を殺す。 -/
structure SemistableModelData where
  /-- 台となる Tate 曲線。 -/
  toTateCurveData : TateCurveData
  /-- ★★**Néron 微分が定める分数イデアル** `ω_E ⊆ L`。

  ★★★★**2026-08-26 の訂正(4 つ目の穴を塞ぐ)**——以前は
  `omega : … → Type` と**階数 1 だけ**を課していた。
  `Check/GaloisRep/OmegaNondegenerate.lean` が示したとおり、それは
  **`ω_E := 𝓞_L`(曲線を無視する定数)で満たせてしまう**。
  ★ここでは `ω_E` を **`L` の分数イデアル**として持ち、
  **変数変換での変わり方**を課す——これが「その欄と入力データの両方を言及する条件」である。

  ★★不変微分 `ω = dx/(2y + a₁x + a₃)` は変数変換で `ω' = u·ω` と変わるので、
  Néron 微分を `ω` で割った係数は `u⁻¹` 倍され、分数イデアルも `u⁻¹` 倍される。 -/
  omegaFrac : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L →
    Submodule (𝓞 L) L
  /-- ★分数イデアルであること(有限生成)——`L` 全体は f.g. でないので、これが効く。

  ★★★★**2026-08-26 の訂正**——`[W.IsElliptic]` を付けた。`Δ = 0` の曲線に
  Néron 微分は無い(第 304 で `TateCurveData` に入れたのと同じ理由)。
  ★欄そのものは全曲線で定義しておき、**性質だけを楕円曲線に限る**。 -/
  omegaFrac_fg : ∀ (L : Type) [Field L] [NumberField L] (W : WeierstrassCurve L) [W.IsElliptic],
    (omegaFrac L W).FG
  /-- ★零でない(階数 1 の代わり)。 -/
  omegaFrac_ne_bot : ∀ (L : Type) [Field L] [NumberField L] (W : WeierstrassCurve L)
    [W.IsElliptic], omegaFrac L W ≠ ⊥
  /-- ★★★★**変数変換で `u⁻¹` 倍される**——これが `ω_E` を曲線に結びつける条件。 -/
  omegaFrac_variableChange : ∀ (L : Type) [Field L] [NumberField L] (W : WeierstrassCurve L)
    [W.IsElliptic] (C : WeierstrassCurve.VariableChange L) (x : L),
    x ∈ omegaFrac L (C • W) ↔ ((C.u : L) * x) ∈ omegaFrac L W
  /-- 大域の半安定性(原文が繰り返し仮定するもの)。 -/
  SemiStable : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L → Prop
  /-- ★★★**その意味を固定する**——各離散付値環の上で半安定であること。

  ★これがあるので `SemiStable := True` のような自由な posit にはならない。 -/
  semiStable_iff : ∀ (L : Type) [Field L] [NumberField L] (W : WeierstrassCurve L),
    SemiStable L W ↔
      ∀ (R : Type) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
        [Algebra R L] [IsFractionRing R L] [W.IsMinimal R], IsSemiStableAt R W

def SemistableModelData.waiting : WaitingFor :=
  { what := "(G7) 半安定還元(★mathlib の還元型で定義済——posit ではない)と、𝓞_L 上への延長が与える Néron 微分の分数イデアル omega_E(★2026-08-26: 変数変換での変わり方を課して曲線に結びつけた)"
    trackB := "Found/GaloisRep — ★mathlib は `AlgebraicGeometry/EllipticCurve/Reduction.lean` を持つ(2026-08-17 実測)が、**Néron モデル・半アーベルスキームは無い**。★★原文は延長の存在を [FC] Chapter I, Proposition 2.7 に帰している" }

/-! ## ★★★G8 —— Faltings 高さ(Arakelov 側との合流点) -/

/-- **(G8)** **Faltings 高さ** `ht^Falt(E) = deg(ω_E)`。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★**ここが Arakelov 側と Galois 側の合流点である。**
`ω_E` は `Spec 𝓞_L` 上の算術直線束なので、その次数を取るには
**(D2)(`APic(Spec 𝓞_F)` と `deg_F`)が要る**。

★`Lemma 3.5` / `Lemma 3.7` / `Proposition 3.4` はすべてこの上で述べられる。 -/
structure FaltingsHeightData where
  /-- 台となる半安定モデル。 -/
  toSemistableModelData : SemistableModelData
  /-- ★★**Faltings 高さ** `ht^Falt(E) = deg(ω_E)`。 -/
  htFalt : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L → ℝ
  /-- ★★★★★**`ht^Falt` は変数変換で変わらない**。

  ## ★★★★2026-08-20 の訂正(§9-404)

  以前は `htFalt_congr`——「omega_E と omega_E' が加群として同型なら
  `ht^Falt(E) = ht^Falt(E')`」——を課していたが、
  ★**これは本物の Neron 微分を排除していた**。

  ★★`L = Q` では整数環が `Z` なので omega_E は**常に階数 1 の自由加群**であり、
  どの 2 つの曲線でも同型になる。★★★したがって旧条件は
  「`Q` 上の Faltings 高さは定数」を強制するが、
  `degInf_ge_localHeight`(局所高さは非有界)と `prop_3_4` に**矛盾**する。

  ★★★★真の原因:**`deg` は計量込みの算術直線束の次数**であり、
  加群としての同型だけでは `ht^Falt` は決まらない。
  ★`ht^Falt = deg(omega_E)` を正しく固めるには **(D3) の計量**が要る
  ——未塗りの穴として §9-404 に記録する。

  ★★代わりに課すのは**真であり、かつ効く**条件である。
  ★★★例えば `htFalt := log|Delta|`(変数変換で `u^12` 倍される)を落とす。

  ## ★★★★★★★★ 2026-08-26 の確認——**この 2 本ではまだ足りない**

  ★第 329 で witness を組んだところ、**`htFalt := deg∞/12` がこの 2 本を満たす**
  ことが分かった(`Check/GaloisRep/HtFaltNotPinned.lean` の `htFalt_not_pinned`)。
  ★★すなわち `prop_3_4` は**恒等的に成り立つ形**で埋まり、
  **`Proposition 3.4` の数学的内容は witness では証明されない**。
  ★★★★★これが界面の欠陥 #6 であり、**初めての「弱すぎる」型**である
  (#1 #2 #4 #5 は充足不能、#3 は弱すぎるで塞いだ)。
  ★★★★★★**塞ぐにはアルキメデス素点での計量が要る**——
  周期束の共体積が入って初めて `ht^Falt` は固定される。
  ★その最下流の葉は `Skeleton/GenEll/Uniformization.lean` の `exists_periodPair`
  (任意の `E/ℂ` が `ℂ/Λ` であること)である。 -/
  htFalt_variableChange : ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
    (C : WeierstrassCurve.VariableChange L), htFalt L (C • E) = htFalt L E
  /-- 無限遠因子の次数 `deg∞`。 -/
  degInf : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L → ℝ
  /-- ★`deg∞` は非負(局所高さの和だから)。 -/
  degInf_nonneg : ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L),
    0 ≤ degInf L E
  /-- ★★★**`deg∞` は局所高さの和である**——(G6) に縛りつける。

  原文 (GenEll p.18):
  > First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

  ★★★**これが「全部 0」の退化を殺す。**乗法還元をもつ素点が 1 つでもあれば
  `deg∞ > 0` が強制される。★`d = [L:ℚ]` の正規化は原文どおり。

  ## ★★★★★2026-08-26 の訂正(4 つ目の充足不能)

  以前は `Lv` を **`L` の任意の拡大**として量化していた。しかしそれだと
  **分岐が効いて偽になる**:`Lv` を `L_v` の分岐次数 `e` の拡大に取ると
  局所高さ `v(q_E)` は `e` 倍になり(`q_E` は同じ元、付値が `e` 倍)、
  右辺はいくらでも大きくなる。★左辺は `L` と `E` だけで決まるから、
  **`e` を大きくすれば必ず破れる**。

  ★★★★訂正:`R` が **`L` の素点 `p` の上にあって不分岐**であること——
  すなわち `L` の元の付値が `p` の付値と一致すること——を仮定に加える。
  ★これは原文の「`v` は `E_L` の局所高さ」(= `L` の素点)という読みに戻す訂正である。 -/
  degInf_ge_localHeight : ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
    (Lv : Type) [Field Lv] [Algebra L Lv]
    (R : Type) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    [Algebra R Lv] [IsFractionRing R Lv] [(E.baseChange Lv).IsElliptic]
    [(E.baseChange Lv).IsMinimal R]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    (p : IsDedekindDomain.HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L,
      (IsDedekindDomain.HeightOneSpectrum.valuation Lv
          (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
        = (IsDedekindDomain.HeightOneSpectrum.valuation L p) x),
    (Module.finrank ℚ L : ℝ) * degInf L E
      ≥ (toSemistableModelData.toTateCurveData.localHeight (E.baseChange Lv) h : ℝ)
        * Real.log 2
  /-- ★★★**`Proposition 3.4`** —— `deg∞` は `ht^Falt` で上から抑えられる。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any -/
  prop_3_4 : ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L]
    (E : WeierstrassCurve L),
    toSemistableModelData.SemiStable L E →
    degInf L E / (12 * (1 + ε)) ≤ htFalt L E + C

def FaltingsHeightData.waiting : WaitingFor :=
  { what := "(G8) Faltings 高さ ht^Falt(E) = deg(omega_E) と、Proposition 3.4(deg_infty を ht^Falt で抑える)"
    trackB := "Found/GaloisRep — ★★★**Arakelov 側との合流点**。`omega_E` は Spec 𝓞_L 上の算術直線束なので、次数を取るには **(D2)** が要る。★`ADiv` / `deg_F` は `Found/GenEll/ArithDiv.lean` に実装済み(sorry 0)なので、(D1)(D2) が入れば `omega_E` の構成だけが残る。★★`Lemma 3.5` / `Lemma 3.7` / `Proposition 3.4` はすべてこの上で述べられる" }

/-! ## ★出典の紐付け(`.src`) -/

def TateCurveData.src : Source :=
  { paper := "GenEll", pdfPage := 15, item := "Definition 3.3(Tate 母数と局所高さの定式化のみ)",
    sectionId := "genell-def-3-3" }

def SemistableModelData.src : Source :=
  { paper := "GenEll", pdfPage := 15, item := "Definition 3.3(半安定還元と 𝓞_L 上のモデルのみ)",
    sectionId := "genell-def-3-3" }

def FaltingsHeightData.src : Source :=
  { paper := "GenEll", pdfPage := 17, item := "Proposition 3.4(Faltings 高さの定式化のみ)",
    sectionId := "genell-prop-3-4" }

end ABC3.Interface.GaloisRep
