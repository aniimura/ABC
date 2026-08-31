/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.AlgebraicGeometry.ProjectiveSpectrum.Proper
import Mathlib.RingTheory.MvPolynomial.Homogeneous
import ABC3.Meta.Claim

/-!
# 射影空間 `ℙⁿ_R` とその**固有性**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★`ample-and-projective-embedding` の段 C1

台帳は段 C を「`ℙⁿ_S` と `O(1)`（点の関子を含む）」とし、**absent** としていた。
★2026-08-27 の実測で、**`ℙⁿ_R` そのものは mathlib の部品だけで書ける**と判った:

| 要るもの | mathlib |
|---|---|
| `Proj 𝒜`（スキーム） | `AlgebraicGeometry.Proj` **ある** |
| `MvPolynomial` の次数付け | `MvPolynomial.gradedAlgebra` **ある**（instance ではない） |
| 構造射 `Proj 𝒜 ⟶ Spec 𝒜₀` | `Proj.toSpecZero` **ある** |
| **固有性** | `instance [Algebra.FiniteType (𝒜 0) A] : IsProper (Proj.toSpecZero 𝒜)` **ある** |

★★したがって段 C は 2 つに割れる:
* **C1（本ファイル、閉）**: `ℙⁿ_R` と `ℙⁿ_R → Spec R` の固有性
* **C2（開）**: `O(1)`（twisting sheaf）と点の関子——★mathlib に twisting sheaf は
  **1 件も無い**（2026-08-27 実測）

## ★★★★埋めた 2 つの隙間

1. ★`Algebra.FiniteType ↥(𝒜 0) (MvPolynomial (Fin (n+1)) R)` は**instance ではない**。
   変数 `X i` が生成することを書く（`mvPolyGradeZero_finiteType`）。
2. ★★`𝒜 0` は `R` そのものではなく `↥(homogeneousSubmodule _ R 0)` である。
   `homogeneousSubmodule_zero : 𝒜 0 = 1` から環同型 `R ≃+* ↥(𝒜 0)` を作る
   （`gradeZeroEquiv`）。★これで構造射を **`Spec R` へ**向けられる。

## ★消費側との関係（明示）

`Proposition 1.4, (iv)` が実際に消費するのは `NorthcottCoord.lean` の
`northcott_of_projModel` であり、それが要るのは**同次座標**であって
`ℙⁿ` そのものではない。
★★本ファイルは「`L^{⊗n}` が閉埋め込み `X ↪ ℙᴺ` を与える」の**行き先**を作った段である。
★★★残るのは `O(1)`・very ample・Serre の定理・座標と高さの比較である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MvPolynomial

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★★次数 0 の部分は `R` である -/

/-- ★★**`𝒜 0 ≃+* R`** —— 次数 0 の斉次多項式は定数である。

★`MvPolynomial.homogeneousSubmodule_zero : 𝒜 0 = 1` と
`Submodule.mem_one` で全射性が出る。単射性は `MvPolynomial.C_injective`。 -/
noncomputable def gradeZeroEquiv (n : ℕ) (R : Type) [CommRing R] :
    R ≃+* ↥(MvPolynomial.homogeneousSubmodule (Fin (n + 1)) R 0) := by
  refine RingEquiv.ofBijective
    ({ toFun := fun r => ⟨C r, isHomogeneous_C _ _⟩
       map_one' := Subtype.ext (by simp)
       map_mul' := fun _ _ => Subtype.ext (by simp)
       map_zero' := Subtype.ext (by simp)
       map_add' := fun _ _ => Subtype.ext (by simp) } : R →+* _) ⟨?_, ?_⟩
  · intro a b hab
    exact MvPolynomial.C_injective _ _ (congrArg Subtype.val hab)
  · rintro ⟨p, hp⟩
    have hp' : p ∈ (1 : Submodule R (MvPolynomial (Fin (n + 1)) R)) := by
      rwa [← MvPolynomial.homogeneousSubmodule_zero]
    rw [Submodule.mem_one] at hp'
    obtain ⟨r, hr⟩ := hp'
    exact ⟨r, Subtype.ext (by simpa [MvPolynomial.algebraMap_eq] using hr)⟩

instance gradeZeroEquiv_isIso (n : ℕ) (R : Type) [CommRing R] :
    IsIso (CommRingCat.ofHom (gradeZeroEquiv n R).toRingHom) :=
  ⟨CommRingCat.ofHom (gradeZeroEquiv n R).symm.toRingHom, by ext x; simp, by ext x; simp⟩

/-! ## ★多項式環は次数 0 の部分の上で有限型 -/

/-- ★**`MvPolynomial (Fin (n+1)) R` は `𝒜 0` の上で有限型**——変数 `X i` が生成する。

★★これが `IsProper (Proj.toSpecZero 𝒜)` を呼ぶための唯一の欠けた instance である
（mathlib は `[Algebra.FiniteType (𝒜 0) A]` を要求するが、
`MvPolynomial` についてはそれを instance にしていない）。 -/
instance mvPolyGradeZero_finiteType (n : ℕ) (R : Type) [CommRing R] :
    Algebra.FiniteType
      ↥(MvPolynomial.homogeneousSubmodule (Fin (n + 1)) R 0)
      (MvPolynomial (Fin (n + 1)) R) := by
  classical
  refine ⟨⟨Finset.univ.image (fun i : Fin (n + 1) => (X i : MvPolynomial (Fin (n + 1)) R)), ?_⟩⟩
  rw [Algebra.eq_top_iff]
  intro p
  induction p using MvPolynomial.induction_on with
  | C r =>
    exact Subalgebra.algebraMap_mem _
      (⟨C r, isHomogeneous_C _ _⟩ : ↥(homogeneousSubmodule (Fin (n + 1)) R 0))
  | add p q hp hq => exact Subalgebra.add_mem _ hp hq
  | mul_X p i hp => exact Subalgebra.mul_mem _ hp (Algebra.subset_adjoin (by simp))

/-! ## ★★★★★★★★`ℙⁿ_R` -/

/-- ★★★★★★★★**射影空間 `ℙⁿ_R`** —— 斉次多項式環の `Proj`。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite. -/
noncomputable def projSpace (n : ℕ) (R : Type) [CommRing R] : Scheme.{0} :=
  Proj (MvPolynomial.homogeneousSubmodule (Fin (n + 1)) R)

/-- ★★★★★★**構造射 `ℙⁿ_R ⟶ Spec R`**。

★`Proj.toSpecZero` の行き先は `Spec 𝒜₀` なので、`gradeZeroEquiv` で `Spec R` へ移す。 -/
noncomputable def projSpaceOverSpec (n : ℕ) (R : Type) [CommRing R] :
    projSpace n R ⟶ Spec (CommRingCat.of R) :=
  Proj.toSpecZero _ ≫ Spec.map (CommRingCat.ofHom (gradeZeroEquiv n R).toRingHom)

/-- ★★★★★★★★★**`ℙⁿ_R` は `Spec R` 上固有**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★機構は mathlib の `IsProper (Proj.toSpecZero 𝒜)`（[Stacks] 01MF、付値判定法）に
同型を 1 つ合成するだけである。

★★`unfold` ではなく **`show`** で形を合わせる——`unfold` の後は
インスタンスの合流に失敗する（`tools\lean-idioms.md`）。 -/
instance projSpaceOverSpec_isProper (n : ℕ) (R : Type) [CommRing R] :
    IsProper (projSpaceOverSpec n R) := by
  show IsProper (Proj.toSpecZero (MvPolynomial.homogeneousSubmodule (Fin (n + 1)) R)
      ≫ Spec.map (CommRingCat.ofHom (gradeZeroEquiv n R).toRingHom))
  exact IsProper.mk

/-- ★★★`ℙⁿ_R` は分離的である（絶対の意味で）。 -/
instance projSpace_isSeparated (n : ℕ) (R : Type) [CommRing R] :
    Scheme.IsSeparated (projSpace n R) :=
  inferInstanceAs (Scheme.IsSeparated (Proj (MvPolynomial.homogeneousSubmodule (Fin (n + 1)) R)))

/-! ### ★出典の紐付け(`.src`)

★★`Proposition 1.4, (iv)` の証明が使う語彙（射影埋め込みの**行き先**）を作る配管である。
埋め込みそのものではない。 -/

def projSpace.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——射影空間 ℙⁿ_R)",
    sectionId := "genell-prop-1-4" }

def projSpaceOverSpec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——構造射 ℙⁿ_R ⟶ Spec R)",
    sectionId := "genell-prop-1-4" }

def projSpaceOverSpec_isProper.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ℙⁿ_R は Spec R 上固有。埋め込みは含まない)",
    sectionId := "genell-prop-1-4" }

def projSpaceOverSpec_isProper.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "AlgebraicGeometry.Proj(斉次環の Proj がスキームであること)"
      (.inMathlib "AlgebraicGeometry.Proj") 6,
    .citation "[mathlib]" "IsProper (Proj.toSpecZero 𝒜)([Stacks] 01MF、付値判定法)"
      (.inMathlib "AlgebraicGeometry.Proj.toSpecZero") 6,
    .citation "[mathlib]" "MvPolynomial.gradedAlgebra(多項式環の次数付け)"
      (.inMathlib "MvPolynomial.gradedAlgebra") 6,
    .implicitStep
      ("★mathlib が要求する Algebra.FiniteType (𝒜 0) A は MvPolynomial について " ++
       "instance になっていないので、変数 X i が生成することを書いた") 6,
    .implicitStep
      ("★★残っているのは段 C2 —— O(1)(twisting sheaf)と点の関子である。" ++
       "mathlib に twisting sheaf は 1 件も無い(2026-08-27 実測)。" ++
       "★★★さらに段 D(ample の定義)・段 E(Serre の定理)・" ++
       "座標と高さの比較が要る") 6 ]

end ABC3.Found.GenEll
