/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithPerfFactorial
import ABC3.Found.Divisor.CartierFrobenioid

/-!
# 算術因子のデータから model Frobenioid へ(`Example 6.3` の `Theorem 5.2, (ii)` の段)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

## ★★幾何版との対応

`CartierFrobenioid.lean`(`Example 6.1`)の `CartierDatum` は
`Γ_A ≤ ℤ[D_A]` を並べたものだった。★算術ではアルキメデス成分が `ℝ≥0` なので
係数を `ℝ` にした `ArithDatum` を置く。

| | 幾何(`Example 6.1`) | 算術(`Example 6.3`) |
|---|---|---|
| 係数 | `ℤ` | **`ℝ`** |
| `Φ(A)` | `effSub (grp A)` | `effR (grp A)` |
| perf-factorial の前提 | `Q`-Cartier | 各素点で**離散か全体**＋**座標ごとに閉じる** |
| `realScale` | 空虚 | **空虚でない**(アルキメデス素点) |

★★組み立て(`phi` が `monoid on 𝒟` であること、`ModelData` を作って
`Theorem 5.2, (ii)` を当てること)は**幾何版とまったく同じ形**である。
-/

namespace ABC3.Found.Divisor

open CategoryTheory ABC3.Found.FrdI

universe v u w

variable (D : Type u) [Category.{v} D]

/-! ## ★1. 算術因子のデータ -/

/-- ★★★**算術因子のデータ** —— `Example 6.3` の `V(L)` と `Φ(L)^gp` を
`𝒟` 上に並べたもの(`CartierDatum` の実係数版)。 -/
structure ArithDatum where
  /-- `V(L)` —— 素点の型。 -/
  primes : D → Type w
  /-- 因子の**引き戻し**(`arithExtend`)。 -/
  pull : ∀ {A B : D}, (B ⟶ A) → ((primes A →₀ ℝ) →+ (primes B →₀ ℝ))
  pull_id : ∀ (A : D) (x : primes A →₀ ℝ), pull (𝟙 A) x = x
  pull_comp : ∀ {A B E : D} (α : B ⟶ A) (β : E ⟶ B) (x : primes A →₀ ℝ),
    pull (β ≫ α) x = pull β (pull α x)
  /-- `Φ(L)^gp` —— 算術因子の群。 -/
  grp : ∀ A : D, AddSubgroup (primes A →₀ ℝ)
  /-- 引き戻しは算術因子を算術因子へ移す。 -/
  pull_mem : ∀ {A B : D} (α : B ⟶ A) {x : primes A →₀ ℝ}, x ∈ grp A → pull α x ∈ grp B
  /-- 引き戻しは有効性を保つ。 -/
  pull_nonneg : ∀ {A B : D} (α : B ⟶ A) {x : primes A →₀ ℝ}, 0 ≤ x → 0 ≤ pull α x
  /-- 引き戻しは**単射**。 -/
  pull_inj : ∀ {A B : D} (α : B ⟶ A), Function.Injective (pull α)
  /-- 各素点で正の元がある。 -/
  gen : ∀ A : D, IsGenSubgroupR (grp A)
  /-- 各素点で局所群が離散か全体(`ord(O_v^▷) ≅ ℤ≥0 / ℝ≥0`)。 -/
  locMono : ∀ A : D, IsLocallyMonoprimeR (grp A)
  /-- 座標ごとに閉じている。 -/
  coord : ∀ A : D, IsCoordwiseR (grp A)

variable {D}

namespace ArithDatum

variable (Δ : ArithDatum.{v, u, w} D)

/-- ★`Φ(A) = Γ_A ∩ ℝ≥0[V_A]` への引き戻しの制限。 -/
noncomputable def mapHom {A B : D} (α : B ⟶ A) : effR (Δ.grp A) →+ effR (Δ.grp B) :=
  AddMonoidHom.codRestrict ((Δ.pull α).comp (effR (Δ.grp A)).subtype) _
    (fun x => mem_effR.mpr ⟨Δ.pull_mem α (mem_effR.mp x.2).1,
      fun s => Finsupp.le_def.mp (Δ.pull_nonneg α (Finsupp.le_def.mpr (mem_effR.mp x.2).2)) s⟩)

@[simp] theorem mapHom_coe {A B : D} (α : B ⟶ A) (x : effR (Δ.grp A)) :
    ((Δ.mapHom α x : effR (Δ.grp B)) : Δ.primes B →₀ ℝ)
      = Δ.pull α ((x : effR (Δ.grp A)) : Δ.primes A →₀ ℝ) := rfl

theorem mapHom_injective {A B : D} (α : B ⟶ A) : Function.Injective (Δ.mapHom α) := by
  intro x y h
  exact Subtype.ext (Δ.pull_inj α
    (congrArg (fun t : effR (Δ.grp B) => (t : Δ.primes B →₀ ℝ)) h))

/-- ★`Φ : 𝒟ᵒᵖ ⥤ 𝔐𝔬𝔫`。 -/
noncomputable def phiFunctor : Dᵒᵖ ⥤ AddCommMonCat.{w} where
  obj A := AddCommMonCat.of (effR (Δ.grp A.unop))
  map {A B} f := AddCommMonCat.ofHom (Δ.mapHom f.unop)
  map_id A := by
    refine AddCommMonCat.hom_ext ?_
    ext x : 2
    exact Δ.pull_id _ _
  map_comp {A B E} f g := by
    refine AddCommMonCat.hom_ext ?_
    ext x : 2
    exact Δ.pull_comp f.unop g.unop _

/-! ## ★2. `monoid on 𝒟` -/

/-- ★★★★**`Φ` は `monoid on 𝒟`**(実係数版)。 -/
noncomputable def phi (hD : IsOfFSMType D) : MonoidOn.{v, u, w} D where
  functor := Δ.phiFunctor
  charInj {A B} α :=
    ⟨Δ.mapHom_injective α,
      charMap_injective_of_sharp (isSharp_effR (Δ.grp B)) (Δ.mapHom_injective α)⟩
  fsmIso {A B} α hα := by
    haveI : IsIso α := hD B A α hα
    haveI : IsIso (Δ.phiFunctor.map α.op) := inferInstance
    exact bijective_of_isIso_addCommMon (Δ.phiFunctor.map α.op)

/-- ★★★**`Φ` は divisorial**。 -/
theorem phi_isDivisorialOn (hD : IsOfFSMType D) : (Δ.phi hD).IsDivisorialOn :=
  fun A => isDivisorial_effR (Δ.grp A)

/-- ★★★★各 `Φ(A)` は **perf-factorial**。

★★幾何版と違い、アルキメデス素点で `M_p ≃+ ℝ≥0` なので
`Definition 2.4, (i)` の `realScale` は**空虚でない**
(`ArithPerfFactorial.lean` の `image_iotaR_univ_of_full`)。 -/
theorem phi_isPerfFactorial (hD : IsOfFSMType D) (A : D) :
    IsPerfFactorial ((Δ.phi hD).val A) :=
  isPerfFactorial_effR (Δ.coord A) (Δ.locMono A)

/-! ## ★3. `Theorem 5.2` へ -/

/-- ★★★★**`Example 6.3` の入力** —— `Δ` と有理関数の単系 `B` から `ModelData` を組む。 -/
noncomputable def modelData (hD : IsOfFSMType D) (bmon : MonoidOn.{v, u, w} D)
    (divB : ∀ A : D, bmon.val A →+ Gp ((Δ.phi hD).val A))
    (divB_nat : ∀ {A B : D} (f : A ⟶ B) (x : bmon.val B),
      divB A (bmon.map f x) = (Δ.phi hD).gpMapOn f (divB B x)) :
    ModelData.{v, u, w} D where
  phi := Δ.phi hD
  bmon := bmon
  divB := divB
  divB_nat := divB_nat

/-- ★★★★★**`Example 6.3` の仮定** —— `Theorem 5.2` の `Hyp` がそのまま出る。 -/
theorem modelHyp (hD : IsOfFSMType D) (bmon : MonoidOn.{v, u, w} D)
    (divB : ∀ A : D, bmon.val A →+ Gp ((Δ.phi hD).val A))
    (divB_nat : ∀ {A B : D} (f : A ⟶ B) (x : bmon.val B),
      divB A (bmon.map f x) = (Δ.phi hD).gpMapOn f (divB B x))
    (hbg : ∀ A : D, IsGroupLike (bmon.val A))
    (htot : IsTotallyEpimorphic D) (hcon : IsConnected D) :
    ModelData.Hyp (Δ.modelData hD bmon divB divB_nat) where
  divisorial := Δ.phi_isDivisorialOn hD
  bmonGroupLike := hbg
  totEpiD := htot
  connectedD := hcon

/-- ★★★★★★**[FrdI] `Example 6.3` の結論** —— 算術因子のデータから
**Frobenioid** `C_{F̄/F}` ができる(`Theorem 5.2, (ii)`)。 -/
theorem arithFrobenioid (hD : IsOfFSMType D) (bmon : MonoidOn.{v, u, w} D)
    (divB : ∀ A : D, bmon.val A →+ Gp ((Δ.phi hD).val A))
    (divB_nat : ∀ {A B : D} (f : A ⟶ B) (x : bmon.val B),
      divB A (bmon.map f x) = (Δ.phi hD).gpMapOn f (divB B x))
    (hbg : ∀ A : D, IsGroupLike (bmon.val A))
    (htot : IsTotallyEpimorphic D) (hcon : IsConnected D) :
    Frobenioid (ModelData.modelPre (Δ.modelHyp hD bmon divB divB_nat hbg htot hcon)) :=
  ModelData.model_frobenioid _

/-- ★★★★★**`Example 6.3`** —— できた Frobenioid は **isotropic 型**。 -/
theorem arithFrobenioid_isotropicType (hD : IsOfFSMType D) (bmon : MonoidOn.{v, u, w} D)
    (divB : ∀ A : D, bmon.val A →+ Gp ((Δ.phi hD).val A))
    (divB_nat : ∀ {A B : D} (f : A ⟶ B) (x : bmon.val B),
      divB A (bmon.map f x) = (Δ.phi hD).gpMapOn f (divB B x))
    (hbg : ∀ A : D, IsGroupLike (bmon.val A))
    (htot : IsTotallyEpimorphic D) (hcon : IsConnected D) :
    IsOfIsotropicType
      (ModelData.modelPre (Δ.modelHyp hD bmon divB divB_nat hbg htot hcon)) :=
  ModelData.model_isotropicType _

/-- ★★★★★**`Example 6.3`** —— できた Frobenioid は
**birationally Frobenius-normalized 型**。 -/
theorem arithFrobenioid_biratFrobNormalizedType (hD : IsOfFSMType D)
    (bmon : MonoidOn.{v, u, w} D)
    (divB : ∀ A : D, bmon.val A →+ Gp ((Δ.phi hD).val A))
    (divB_nat : ∀ {A B : D} (f : A ⟶ B) (x : bmon.val B),
      divB A (bmon.map f x) = (Δ.phi hD).gpMapOn f (divB B x))
    (hbg : ∀ A : D, IsGroupLike (bmon.val A))
    (htot : IsTotallyEpimorphic D) (hcon : IsConnected D) :
    IsOfBirationallyFrobeniusNormalizedType _
      (ModelData.modelPre (Δ.modelHyp hD bmon divB divB_nat hbg htot hcon))
      (ModelData.model_frobenioid _) :=
  ModelData.model_isOfBiratFrobNormalizedType _

end ArithDatum

/-! ### ★出典の紐付け -/

/-- ★★locator —— `Example 6.3` の「`Theorem 5.2, (ii)` を当てる」段
(★**条つき**: `ArithDatum` の各フィールドを `𝒟 = B(G)⁰` の上で実現する段
——Galois 対応で `B(G)⁰` の対象を部分体に読み替える——は未実装)。 -/
def ArithDatum.arithFrobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — 算術因子のデータから model Frobenioid",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
