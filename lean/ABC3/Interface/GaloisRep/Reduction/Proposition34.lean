/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Interface.GaloisRep.Representation
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.AlgebraicGeometry.EllipticCurve.Reduction
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.GroupTheory.QuotientGroup.Basic
import Mathlib.NumberTheory.NumberField.InfinitePlace.Embeddings
import Mathlib.NumberTheory.ModularForms.Discriminant
import Mathlib.NumberTheory.ModularForms.EisensteinSeries.Basic
import ABC3.Interface.GaloisRep.Reduction.Definition33

/-!
# Reduction —— `[GenEll] Proposition 3.4` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Interface.GaloisRep

open ABC3.Meta WeierstrassCurve NumberField

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
  /-- ★★★★★**アルキメデス素点でのノルム** `‖Δ‖_arch(E^σ)`。

  ## ★★★★★★★★ 2026-08-26 の訂正(欠陥 #6 を塞ぐ)

  以前は `htFalt` に `htFalt_variableChange` と `prop_3_4` の 2 本しか課しておらず、
  ★`htFalt := deg∞/12` で満たせてしまった(§9-659 以降の記録)。
  ★★塞ぐには**アルキメデス素点での計量**が要る——それがこの欄である。

  ★★★`ω_E` のアルキメデス norm は周期束 `Λ_σ` の共体積を通じて
  `‖Δ‖_arch = |Δ(Λ_σ)|·covol(Λ_σ)⁶` で与えられ、これは**相似で不変**な量である。 -/
  archNorm : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L → (L →+* ℂ) → ℝ
  /-- ★楕円曲線では正である。 -/
  archNorm_pos : ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic]
    (σ : L →+* ℂ), 0 < archNorm L E σ
  /-- ★★★★★★**`archNorm` の同定**——これが `archNorm` を完全に固定する。

  `j(τ) = j(E^σ)` なる `τ` について `archNorm = 4096·π¹²·‖Δ(τ)‖·(Im τ)⁶`
  (`Δ` はモジュラー判別式、右辺は `Δ` の **Petersson ノルム**の定数倍)。

  ★★`j : ℍ → ℂ` は全射なので、そのような `τ` は必ず存在する。
  ★★★したがってこの条件は `archNorm` を**一意に決める**——定数で埋める逃げ道は無い。 -/
  archNorm_eq : ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic]
    (σ : L →+* ℂ) (τ : UpperHalfPlane),
    ModularForm.E₄ τ ^ 3 / ModularForm.discriminant τ = (E.map σ).j →
    archNorm L E σ = 4096 * Real.pi ^ 12 * (‖ModularForm.discriminant τ‖ * τ.im ^ 6)
  /-- ★★★★★★★**`ht^Falt` の定義式**——これが欠陥 #6 を塞ぐ。

  真の Faltings 高さは

      12·ht^Falt(E) = deg∞(E) − (1/d)·Σ_{σ} log( (2π)¹²·‖Δ‖_arch(E^σ) )

  である(`d = [L:ℚ]`)。★有限素点側は `deg∞`、アルキメデス側が第 2 項である。

  ★★★★**これで `htFalt` の欄は自由でなくなる**——`degInf` と `archNorm` から
  一意に決まる。★`archNorm` は `archNorm_eq` で固定されているので、
  **`htFalt := deg∞/12` はもはや通らない**(アルキメデス項が残る)。 -/
  htFalt_eq : ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic],
    12 * htFalt L E
      = degInf L E
        - (∑ σ : (L →+* ℂ), Real.log ((2 * Real.pi) ^ 12 * archNorm L E σ))
            / (Module.finrank ℚ L : ℝ)
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

def FaltingsHeightData.src : Source :=
  { paper := "GenEll", pdfPage := 17, item := "Proposition 3.4(Faltings 高さの定式化のみ)",
    sectionId := "genell-prop-3-4" }

end ABC3.Interface.GaloisRep
