/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.CartierPerfFactorial
import ABC3.Found.FrdI.Thm52NatIso

/-!
# Cartier 因子のデータから Frobenioid へ —— `Example 6.1` の骨組み

★★★★[FrdI] `Example 6.1` は、幾何(`V`・`K̄`・`D_K`)から

* `𝒟 = B(G)⁰`(底の圏。`Sec6GaloisCat.lean` で**閉じた**)
* `Φ(L) ⊆ ℤ≥0[D_L]`(Cartier 有効因子の単系。**perf-factorial**、`Cartier*.lean`)
* `B(L) ⊆ L^×`(零点・極が `D_L` に載る有理関数)
* `Div_B : B → Φ^gp`

を作り、**`Theorem 5.2, (ii)` を当てるだけ**である。

★本ファイルは**その「当てるだけ」の部分を実際に当てる**。幾何が与えるものを
`CartierDatum` として型に出し、そこから

1. `Φ` が **`monoid on 𝒟`** であること(`phi`)
2. `Φ` が **divisorial** であること(`phi_isDivisorialOn`)
3. 各 `Φ(A)` が **perf-factorial** であること(`phi_isPerfFactorial`)
4. `B`・`Div_B` を足すと **`Theorem 5.2` の `ModelData` と `Hyp`** になること
5. したがって **Frobenioid** ができ、**isotropic 型**かつ
   **birationally Frobenius-normalized 型**であること

を出す。

## ★`monoid on 𝒟` の 2 条件がどこから来るか

`Definition 1.1, (ii)` の (a)(b) は、この設定では**ただで出る**:

* **(a) characteristically injective** —— `Φ(A) = Γ_A ∩ ℤ≥0[D_A]` は **sharp** なので、
  引き戻しが単射でありさえすれば `M^char` の側も単射になる
  (`charMap_injective_of_sharp`)。★幾何では「引き戻しの単射性」は
  `V[L] → V[M]` が支配的であることから来る。
* **(b) FSM 射で同型** —— **`𝒟` が FSM 型**(`B(G)⁰` はそう。`finSubOp_isOfFSMType`)
  なら FSM 射はすべて同型で、関手は同型を同型に送る。

★**したがって残っている幾何は `CartierDatum` の各フィールドの実現だけ**である ——
`D_L`(`V[L]` の素因子)、引き戻し、Cartier 因子の群 `Γ_L`、
そして `D_K` が `K`-`Q`-Cartier であること。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w

/-! ## ★0. 補助 -/

/-- ★`AddCommMonCat` の同型は台の全単射を与える。 -/
theorem bijective_of_isIso_addCommMon {X Y : AddCommMonCat.{w}} (f : X ⟶ Y) [IsIso f] :
    Function.Bijective (AddCommMonCat.Hom.hom f) := by
  have hli : ∀ x, (inv f).hom (f.hom x) = x := by
    intro x
    have h := congrArg (fun t : X ⟶ X => AddCommMonCat.Hom.hom t x) (IsIso.hom_inv_id f)
    simpa using h
  have hri : ∀ y, f.hom ((inv f).hom y) = y := by
    intro y
    have h := congrArg (fun t : Y ⟶ Y => AddCommMonCat.Hom.hom t y) (IsIso.inv_hom_id f)
    simpa using h
  exact ⟨Function.LeftInverse.injective hli, Function.RightInverse.surjective hri⟩

variable (D : Type u) [Category.{v} D]

/-! ## ★1. Cartier 因子のデータ -/

/-- ★★★**Cartier 因子のデータ** —— [FrdI] `Example 6.1` の幾何が与えるもの。

`primes A = D_L`、`grp A = `(`V[L]` 上の Cartier 因子で台が `D_L` に入るもの)、
`pull α = `(`V[L] → V[M]` に沿う因子の引き戻し)。 -/
structure CartierDatum where
  /-- `D_L` —— `V[L]` 上の素因子のうち `D_K` の上にあるもの。 -/
  primes : D → Type w
  /-- 因子の**引き戻し**。`α : B ⟶ A` に対し `ℤ[D_A] → ℤ[D_B]`。 -/
  pull : ∀ {A B : D}, (B ⟶ A) → ((primes A →₀ ℤ) →+ (primes B →₀ ℤ))
  pull_id : ∀ (A : D) (x : primes A →₀ ℤ), pull (𝟙 A) x = x
  pull_comp : ∀ {A B E : D} (α : B ⟶ A) (β : E ⟶ B) (x : primes A →₀ ℤ),
    pull (β ≫ α) x = pull β (pull α x)
  /-- `Γ_L` —— `V[L]` 上の Cartier 因子で台が `D_L` に入るもの。 -/
  grp : ∀ A : D, AddSubgroup (primes A →₀ ℤ)
  /-- 引き戻しは Cartier 性を保つ。 -/
  pull_mem : ∀ {A B : D} (α : B ⟶ A) {x : primes A →₀ ℤ}, x ∈ grp A → pull α x ∈ grp B
  /-- 引き戻しは有効性を保つ。 -/
  pull_nonneg : ∀ {A B : D} (α : B ⟶ A) {x : primes A →₀ ℤ}, 0 ≤ x → 0 ≤ pull α x
  /-- 引き戻しは**単射**(幾何では `V[L] → V[M]` が支配的であることから)。 -/
  pull_inj : ∀ {A B : D} (α : B ⟶ A), Function.Injective (pull α)
  /-- `D_K` が `K`-`Q`-Cartier であること。 -/
  qc : ∀ A : D, IsQCartierSubgroup (grp A)

variable {D}

namespace CartierDatum

variable (Δ : CartierDatum.{v, u, w} D)

/-- ★`Φ(A) = Γ_A ∩ ℤ≥0[D_A]` への引き戻しの制限。 -/
noncomputable def mapHom {A B : D} (α : B ⟶ A) : effSub (Δ.grp A) →+ effSub (Δ.grp B) :=
  AddMonoidHom.codRestrict ((Δ.pull α).comp (effSub (Δ.grp A)).subtype) _
    (fun x => mem_effSub.mpr ⟨Δ.pull_mem α (mem_effSub.mp x.2).1,
      fun s => Finsupp.le_def.mp (Δ.pull_nonneg α (Finsupp.le_def.mpr (mem_effSub.mp x.2).2)) s⟩)

@[simp] theorem mapHom_coe {A B : D} (α : B ⟶ A) (x : effSub (Δ.grp A)) :
    ((Δ.mapHom α x : effSub (Δ.grp B)) : Δ.primes B →₀ ℤ)
      = Δ.pull α ((x : effSub (Δ.grp A)) : Δ.primes A →₀ ℤ) := rfl

theorem mapHom_injective {A B : D} (α : B ⟶ A) : Function.Injective (Δ.mapHom α) := by
  intro x y h
  exact Subtype.ext (Δ.pull_inj α
    (congrArg (fun t : effSub (Δ.grp B) => (t : Δ.primes B →₀ ℤ)) h))

/-- ★`Φ : 𝒟ᵒᵖ ⥤ 𝔐𝔬𝔫`。 -/
noncomputable def phiFunctor : Dᵒᵖ ⥤ AddCommMonCat.{w} where
  obj A := AddCommMonCat.of (effSub (Δ.grp A.unop))
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

/-- ★★★★**`Φ` は `monoid on 𝒟`**。

(a) は `Φ(A)` が sharp かつ引き戻しが単射だから、
(b) は `𝒟` が FSM 型(FSM 射がすべて同型)だから。 -/
noncomputable def phi (hD : IsOfFSMType D) : MonoidOn.{v, u, w} D where
  functor := Δ.phiFunctor
  charInj {A B} α :=
    ⟨Δ.mapHom_injective α,
      charMap_injective_of_sharp (isSharp_effSub (Δ.grp B)) (Δ.mapHom_injective α)⟩
  fsmIso {A B} α hα := by
    haveI : IsIso α := hD B A α hα
    haveI : IsIso (Δ.phiFunctor.map α.op) := inferInstance
    exact bijective_of_isIso_addCommMon (Δ.phiFunctor.map α.op)

/-- ★★★**`Φ` は divisorial**。 -/
theorem phi_isDivisorialOn (hD : IsOfFSMType D) : (Δ.phi hD).IsDivisorialOn :=
  fun A => isDivisorial_effSub (Δ.grp A)

/-- ★★★各 `Φ(A)` は **perf-factorial**。 -/
theorem phi_isPerfFactorial (hD : IsOfFSMType D) (A : D) :
    IsPerfFactorial ((Δ.phi hD).val A) :=
  isPerfFactorial_effSub (Δ.qc A)

/-! ## ★3. `Theorem 5.2` へ -/

/-- ★★★★**`Example 6.1` の入力** —— `Δ` と有理関数の単系 `B` から `ModelData` を組む。 -/
noncomputable def modelData (hD : IsOfFSMType D) (bmon : MonoidOn.{v, u, w} D)
    (divB : ∀ A : D, bmon.val A →+ Gp ((Δ.phi hD).val A))
    (divB_nat : ∀ {A B : D} (f : A ⟶ B) (x : bmon.val B),
      divB A (bmon.map f x) = (Δ.phi hD).gpMapOn f (divB B x)) :
    ModelData.{v, u, w} D where
  phi := Δ.phi hD
  bmon := bmon
  divB := divB
  divB_nat := divB_nat

/-- ★★★★★**`Example 6.1` の仮定** —— `Theorem 5.2` の `Hyp` がそのまま出る。 -/
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

/-- ★★★★★★**[FrdI] `Example 6.1` の結論** —— Cartier 因子のデータから
**Frobenioid** ができる(`Theorem 5.2, (ii)`)。 -/
theorem cartierFrobenioid (hD : IsOfFSMType D) (bmon : MonoidOn.{v, u, w} D)
    (divB : ∀ A : D, bmon.val A →+ Gp ((Δ.phi hD).val A))
    (divB_nat : ∀ {A B : D} (f : A ⟶ B) (x : bmon.val B),
      divB A (bmon.map f x) = (Δ.phi hD).gpMapOn f (divB B x))
    (hbg : ∀ A : D, IsGroupLike (bmon.val A))
    (htot : IsTotallyEpimorphic D) (hcon : IsConnected D) :
    Frobenioid (ModelData.modelPre (Δ.modelHyp hD bmon divB divB_nat hbg htot hcon)) :=
  ModelData.model_frobenioid _

/-- ★★★★★**`Example 6.1`** —— できた Frobenioid は **isotropic 型**。 -/
theorem cartierFrobenioid_isotropicType (hD : IsOfFSMType D) (bmon : MonoidOn.{v, u, w} D)
    (divB : ∀ A : D, bmon.val A →+ Gp ((Δ.phi hD).val A))
    (divB_nat : ∀ {A B : D} (f : A ⟶ B) (x : bmon.val B),
      divB A (bmon.map f x) = (Δ.phi hD).gpMapOn f (divB B x))
    (hbg : ∀ A : D, IsGroupLike (bmon.val A))
    (htot : IsTotallyEpimorphic D) (hcon : IsConnected D) :
    IsOfIsotropicType
      (ModelData.modelPre (Δ.modelHyp hD bmon divB divB_nat hbg htot hcon)) :=
  ModelData.model_isotropicType _

/-- ★★★★★**`Example 6.1`** —— できた Frobenioid は
**birationally Frobenius-normalized 型**。 -/
theorem cartierFrobenioid_biratFrobNormalizedType (hD : IsOfFSMType D)
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

end CartierDatum

/-! ### ★出典の紐付け -/

/-- ★locator —— `Example 6.1` の「`Theorem 5.2, (ii)` を当てる」段
(★**条つき**: `CartierDatum` の各フィールドの幾何的な実現、
すなわち `V[L]` の構成と `D_L`・`Γ_L`・引き戻しは未実装)。 -/
def CartierDatum.cartierFrobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — Cartier 因子のデータから model Frobenioid",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.FrdI
