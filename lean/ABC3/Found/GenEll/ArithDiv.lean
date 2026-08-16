import Mathlib.NumberTheory.NumberField.InfinitePlace.Basic
import Mathlib.RingTheory.DedekindDomain.Ideal.Lemmas
import Mathlib.RingTheory.DedekindDomain.AdicValuation
import Mathlib.RingTheory.DedekindDomain.FiniteAdeleRing
import Mathlib.RingTheory.Ideal.Norm.AbsNorm
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import ABC3.Meta.Claim

/-!
# [GenEll] §1 の算術因子 —— `ADiv(F)` と `deg_F`(層 A)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.4。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.4):
> An arithmetic divisor on F is defined to be a finite formal sum

## ★★なぜここだけ先に `Found/` に置けるのか —— 層の境界

§1 が最終的に要るのは**高さ**だが、原文の定義は

```
ht_M̄(x) ≝ deg_F(x_F^* M̄) ∈ ℝ
```

である(p.4、目視)。すなわち **一般の `X` 上の算術直線束は
「`APic(Spec 𝒪_F)` の元を作る」ためにだけ要る**。そして原文自身が

```
ADiv(F)/APrc(F) ⥲ APic(Spec(𝒪_F))
```

という同型を述べている。★**右辺はスキーム論だが、左辺は数体の素点だけで書ける。**

| 層 | 中身 | mathlib(2026-08-16 実測) | 置き場所 |
|---|---|---|---|
| **A** | `ADiv(F)` / `deg_F` | ★材料あり(`InfinitePlace` 70 件 / `HeightOneSpectrum` 301 件) | ★**本ファイル** |
| B | スキーム上の直線束・`APic(X)`・引き戻し | **無い**(`LineBundle` 0 件) | `Interface/` |
| C | `X^arc`(解析化)・hermitian 計量・`ι_X` | **無い**(`analytification` / `GAGA` / `complex analytic space` すべて 0 件) | `Interface/` |
| D | 同型 `ADiv/APrc ≅ APic` と `ht` の定義 | — | `Skeleton/` |

★**層 C が「もう一段下」である**——hermitian 計量を載せる先の `X^arc` そのものが無い。
**Arakelov 理論は `L0:基礎` の最下段ではなく、その下にさらに層がある。**
詳細は `ResearchPaper/1_Structured/Arithmetic Elliptic Curves in General Position/definition-1-1.html`。

## ★この形にした理由(原文の指定どおり)

原文 (GenEll p.4):
> — where cv ∈ Z if v ∈ V(F )non and cv ∈ R if v ∈ V(F )arc [cf. [Szp], §1.1]. Here, if

★**係数の型が素点によって違う**(非アルキメデスは `ℤ`、アルキメデスは `ℝ`)。
したがって `ADiv(F)` は自由アーベル群ではなく、**2 つの `Finsupp` の直積**である。
これは飾りではない——`deg_F` が `log(q_v)`(離散)と `1`(連続)を混ぜる理由そのものである。
-/

namespace ABC3.Found.GenEll

open NumberField

section ArithmeticDivisor

variable (F : Type*) [Field F] [NumberField F]

/-- 有限素点(= 整数環の高さ 1 素イデアル)。原文の `𝕍(F)^non`。 -/
abbrev FinitePlace := IsDedekindDomain.HeightOneSpectrum (𝓞 F)

/-- **算術因子** `ADiv(F)`。

原文 (GenEll p.4):
> — where cv ∈ Z if v ∈ V(F )non and cv ∈ R if v ∈ V(F )arc [cf. [Szp], §1.1]. Here, if

★係数の型が素点の種類で違うので、直積で書く。「有限形式和」は `Finsupp` である。 -/
def ADiv : Type _ := (FinitePlace F →₀ ℤ) × (InfinitePlace F →₀ ℝ)

noncomputable instance : AddCommGroup (ADiv F) :=
  inferInstanceAs (AddCommGroup ((FinitePlace F →₀ ℤ) × (InfinitePlace F →₀ ℝ)))

noncomputable instance : Inhabited (ADiv F) := ⟨0⟩

variable {F}

/-- 有限素点の成分。 -/
def ADiv.fin (a : ADiv F) : FinitePlace F →₀ ℤ := a.1

/-- 無限素点の成分。 -/
def ADiv.arc (a : ADiv F) : InfinitePlace F →₀ ℝ := a.2

omit [NumberField F] in
@[simp] theorem ADiv.fin_add (a b : ADiv F) : (a + b).fin = a.fin + b.fin := rfl

omit [NumberField F] in
@[simp] theorem ADiv.arc_add (a b : ADiv F) : (a + b).arc = a.arc + b.arc := rfl

omit [NumberField F] in
@[simp] theorem ADiv.fin_zero : (0 : ADiv F).fin = 0 := rfl

omit [NumberField F] in
@[simp] theorem ADiv.arc_zero : (0 : ADiv F).arc = 0 := rfl

/-- **effective**(原文 "if all of the cv are ≥ 0")。 -/
def ADiv.IsEffective (a : ADiv F) : Prop :=
  (∀ v, 0 ≤ a.fin v) ∧ (∀ v, 0 ≤ a.arc v)

/-- 有限素点 `v` の剰余体の濃度 `q_v`。原文「q_v for the cardinality of the residue field of F_v」。 -/
noncomputable def residueCard (v : FinitePlace F) : ℕ := Ideal.absNorm v.asIdeal

end ArithmeticDivisor

section Degree

variable {F : Type*} [Field F] [NumberField F]

/-- **次数写像** `deg_F`。

原文 (GenEll p.4):
> [cf. [Szp], §1.1] determines a homomorphism APic(Spec(OF )) → R, which we shall

原文は `𝕍(F)^non ∋ v ↦ log(q_v)`、`𝕍(F)^arc ∋ v ↦ 1` で定める。
★**離散側は `log(q_v)`、連続側は `1`** ——この非対称が `ADiv` の型の非対称と対応する。 -/
noncomputable def deg (a : ADiv F) : ℝ :=
  a.fin.sum (fun v n => (n : ℝ) * Real.log (residueCard v)) + a.arc.sum (fun _ r => r)

@[simp] theorem deg_zero : deg (0 : ADiv F) = 0 := by
  simp only [deg, ADiv.fin_zero, ADiv.arc_zero, Finsupp.sum_zero_index, add_zero]

theorem deg_add (a b : ADiv F) : deg (a + b) = deg a + deg b := by
  have hfin : (a.fin + b.fin).sum (fun v n => (n : ℝ) * Real.log (residueCard v))
      = a.fin.sum (fun v n => (n : ℝ) * Real.log (residueCard v))
        + b.fin.sum (fun v n => (n : ℝ) * Real.log (residueCard v)) :=
    Finsupp.sum_add_index' (by intro v; simp) (by intro v m n; push_cast; ring)
  have harc : (a.arc + b.arc).sum (fun _ r => r)
      = a.arc.sum (fun _ r => r) + b.arc.sum (fun _ r => r) :=
    Finsupp.sum_add_index' (by intro v; simp) (by intro v m n; ring)
  simp only [deg, ADiv.fin_add, ADiv.arc_add, hfin, harc]
  ring

/-- `deg` を加法準同型として束ねたもの。 -/
noncomputable def degHom : ADiv F →+ ℝ where
  toFun := deg
  map_zero' := deg_zero
  map_add' := deg_add

@[simp] theorem degHom_apply (a : ADiv F) : degHom a = deg a := rfl

/-- **正規化した次数** `deg_F ≝ (1/[F:ℚ]) · deg_F`(原文の下線つき `deg`)。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★正規化する理由は**有限次拡大で不変にする**ためであり、
その不変性が `ht` を `X(ℚ̄)` の上で well-defined にする。
★**不変性そのものは本ファイルでは示していない**——`Skeleton/GenEll/Heights.lean` の
`deg_normalized_base_change` を参照(そこが層 D)。 -/
noncomputable def degNormalized (a : ADiv F) : ℝ :=
  deg a / (Module.finrank ℚ F : ℝ)

theorem degNormalized_add (a b : ADiv F) :
    degNormalized (a + b) = degNormalized a + degNormalized b := by
  simp only [degNormalized, deg_add]
  ring

@[simp] theorem degNormalized_zero : degNormalized (0 : ADiv F) = 0 := by
  simp [degNormalized]

end Degree

section Principal

variable {F : Type*} [Field F] [NumberField F]

open IsDedekindDomain

/-- **`ord_v(f)`** —— 原文「we shall write ord_v(−) : F_v → Z for the order defined by v」。

`v.valuation` は `ℤᵐ⁰ = WithZero (Multiplicative ℤ)` に値を取るので、
`f ≠ 0` を使って `ℤ` へ落とす。 -/
noncomputable def ordv (v : FinitePlace F) (f : Fˣ) : ℤ :=
  Multiplicative.toAdd
    (WithZero.unzero ((Valuation.ne_zero_iff (v.valuation F)).2 (Units.ne_zero f)))

theorem ordv_eq_zero_iff (v : FinitePlace F) (f : Fˣ) :
    ordv v f = 0 ↔ v.valuation F (f : F) = 1 := by
  constructor
  · intro h
    have := WithZero.coe_unzero ((Valuation.ne_zero_iff (v.valuation F)).2 (Units.ne_zero f))
    rw [← this]
    have h1 : WithZero.unzero ((Valuation.ne_zero_iff (v.valuation F)).2 (Units.ne_zero f)) = 1 := by
      simp only [ordv] at h
      have h2 := congrArg Multiplicative.ofAdd h
      rwa [ofAdd_toAdd, ofAdd_zero] at h2
    rw [h1]
    rfl
  · intro h
    simp only [ordv]
    have : WithZero.unzero ((Valuation.ne_zero_iff (v.valuation F)).2 (Units.ne_zero f)) = 1 := by
      apply WithZero.coe_injective
      rw [WithZero.coe_unzero, h]
      rfl
    rw [this]
    rfl

/-- ★**`ord_v(f) ≠ 0` となる `v` は有限個**。

`f` と `f⁻¹` の台(mathlib の `HeightOneSpectrum.Support`)の合併に含まれる。
★これが「有限形式和」であることの中身であり、`Finsupp` を作る根拠である。 -/
theorem ordv_finite_support (f : Fˣ) :
    {v : FinitePlace F | ordv v f ≠ 0}.Finite := by
  have hsub : {v : FinitePlace F | ordv v f ≠ 0}
      ⊆ HeightOneSpectrum.Support (𝓞 F) (f : F)
        ∪ HeightOneSpectrum.Support (𝓞 F) ((f⁻¹ : Fˣ) : F) := by
    intro v hv
    by_contra hmem
    simp only [Set.mem_union, not_or] at hmem
    apply hv
    rw [ordv_eq_zero_iff]
    have h1 : v.valuation F (f : F) ≤ 1 := not_lt.mp hmem.1
    have h2 : v.valuation F ((f⁻¹ : Fˣ) : F) ≤ 1 := not_lt.mp hmem.2
    have hmul : v.valuation F (f : F) * v.valuation F ((f⁻¹ : Fˣ) : F) = 1 := by
      rw [← Valuation.map_mul]
      simp
    exact le_antisymm h1 (by
      by_contra hlt
      rw [not_le] at hlt
      have : v.valuation F (f : F) * v.valuation F ((f⁻¹ : Fˣ) : F) < 1 :=
        mul_lt_one_of_lt_of_le hlt h2
      exact absurd hmul (ne_of_lt this))
  exact Set.Finite.subset
    ((HeightOneSpectrum.Support.finite (𝓞 F) (f : F)).union
      (HeightOneSpectrum.Support.finite (𝓞 F) ((f⁻¹ : Fˣ) : F))) hsub

/-- **主算術因子** `ADiv(f)`。

原文 (GenEll p.4):
> An arithmetic divisor on F is defined to be a finite formal sum

原文の定義:
`ADiv(f) ≝ Σ_{v∈𝕍(F)^non} ord_v(f)·v − Σ_{v∈𝕍(F)^arc} [F_v : ℝ]·log(|f|_v)·v`。

★アルキメデス側の `[F_v : ℝ]` は `InfinitePlace.mult`(実素点で 1、複素素点で 2)である。
★アルキメデス素点は有限個(`Fintype (InfinitePlace F)`)なので `Finsupp` は自動で作れる。 -/
noncomputable def principalADiv (f : Fˣ) : ADiv F :=
  (Finsupp.onFinset (ordv_finite_support f).toFinset (fun v => ordv v f)
     (by intro v hv; simpa using hv),
   Finsupp.onFinset Finset.univ
     (fun v : InfinitePlace F => -((v.mult : ℝ) * Real.log (v (f : F))))
     (by intro v _; exact Finset.mem_univ v))

/-- **主算術因子のなす部分群** `APrc(F) ⊆ ADiv(F)`。

原文 (GenEll p.4):
> An arithmetic divisor on F is defined to be a finite formal sum

★原文は「the principal arithmetic divisors determine a subgroup APrc(F) ⊆ ADiv(F)」と述べる。
**部分群であること自体は `principalADiv` が準同型であることから出る**が、
ここでは**生成する部分群**として定義する——両者は一致する(原文の主張)。
★`principalADiv` が準同型であることは**まだ示していない**。示せば
`AddMonoidHom.range` と一致することが言える。 -/
noncomputable def APrc (F : Type*) [Field F] [NumberField F] : AddSubgroup (ADiv F) :=
  AddSubgroup.closure (Set.range (principalADiv : Fˣ → ADiv F))

theorem principalADiv_mem_APrc (f : Fˣ) : principalADiv f ∈ APrc F :=
  AddSubgroup.subset_closure ⟨f, rfl⟩

end Principal

/-! ## ★出典の紐付け(`.src`)

★これは番号付き項目ではなく **§1 の地の文**である。`sectionId` は
`1_Structured/…/definition-1-1.html` の該当 `<section>` を指す。 -/

def ADiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(算術因子 ADiv(F))",
    sectionId := "genell-adiv" }

def deg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(次数写像 deg_F)",
    sectionId := "genell-deg" }

end ABC3.Found.GenEll
