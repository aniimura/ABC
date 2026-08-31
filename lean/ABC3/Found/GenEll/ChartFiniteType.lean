/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AmpleDef
import ABC3.Meta.Claim
import Mathlib.AlgebraicGeometry.Morphisms.FiniteType

/-!
# ★★★★★★段 E3b —— チャートの座標環は有限生成である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★これは何か —— 段 E3b

`f : X ⟶ Spec A` が**有限型**で `V` がアフィン開なら、
★`Γ(X, V)` は `A` 上**有限生成な代数**である。

★★段 E3c（`A⁰_{x_i} → Γ(X, X_{s_i})` の全射性）は「生成元それぞれに段 E3a を当て、
指数を最大値で揃える」という形で進むので、**まず生成元が有限個であること**が要る。
それが本ファイルである。

## ★★★機構 —— mathlib の 1 行

    `Scheme.Hom.finiteType_appLE`
      : `LocallyOfFiniteType f` なら、アフィン開 `U ⊆ Y`・`V ⊆ f⁻¹U` について
        `Γ(Y,U) ⟶ Γ(X,V)` は有限型

★`U = ⊤`（`Spec A` の全体、アフィン）と取るだけである。
★★あとは `Algebra.FiniteType.out`（`⊤ : Subalgebra` が有限生成）で
**具体の有限生成集合**を取り出す（段 E3c が使う形）。

## ★測定の記録

★底は `A` ではなく **`Γ(Spec A, ⊤)`** である（`ΓSpecIso` で移り合うが**定義的には等しくない**）。
★★本ファイルは `Γ(Spec A, ⊤)` のまま述べる——同型を挟むと式が読めなくなるうえ、
消費側（段 E3c）は `Γ` の側で働くからである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MonoidalCategory ABC3.Found.Arakelov

/-! ## ★有限型の環準同型からは有限生成集合が取れる -/

/-- ★**有限型の環準同型からは具体の有限生成集合が取れる**（`Algebra.FiniteType.out`）。 -/
theorem exists_finset_adjoin_eq_top {R S : Type} [CommRing R] [CommRing S] (φ : R →+* S)
    (h : φ.FiniteType) :
    ∃ T : Finset S, @Algebra.adjoin R S _ _ φ.toAlgebra (T : Set S) = ⊤ := by
  letI := φ.toAlgebra
  obtain ⟨T, hT⟩ := (h.out : (⊤ : Subalgebra R S).FG)
  exact ⟨T, hT⟩

variable {X : Scheme.{0}} {A : CommRingCat.{0}}

/-! ## ★★★★★★段 E3b -/

/-- ★★★★★**アフィン開の座標環は有限型である**（段 E3b）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★機構は mathlib の `Scheme.Hom.finiteType_appLE` に `U = ⊤` を入れるだけである。 -/
theorem finiteType_appLE_top (f : X ⟶ Spec A) [LocallyOfFiniteType f]
    {V : X.Opens} (hV : IsAffineOpen V) :
    (f.appLE ⊤ V (by simp)).hom.FiniteType :=
  Scheme.Hom.finiteType_appLE f (isAffineOpen_top _) hV _

/-- ★★★★★★**チャートの座標環は有限生成である** —— 段 E3b の使う形。

`f : X ⟶ Spec A` が有限型で `V` がアフィン開なら、
★`Γ(X, V)` を `Γ(Spec A, ⊤)` 上生成する**有限集合**がある。

★★段 E3c はこの生成元それぞれに段 E3a-3（大域化、`§9-831`）を当て、
指数を最大値で揃えることで `A⁰_{x_i} → Γ(X, X_{s_i})` の全射性を出す。 -/
theorem exists_finset_generating (f : X ⟶ Spec A) [LocallyOfFiniteType f]
    {V : X.Opens} (hV : IsAffineOpen V) :
    ∃ T : Finset (Γ(X, V) : Type),
      @Algebra.adjoin (Γ(Spec A, (⊤ : (Spec A).Opens)) : Type) (Γ(X, V) : Type) _ _
        (f.appLE ⊤ V (by simp)).hom.toAlgebra (T : Set (Γ(X, V) : Type)) = ⊤ :=
  exists_finset_adjoin_eq_top _ (finiteType_appLE_top f hV)

/-- ★★★★★★★**`ample` が与えるチャートの座標環は有限生成である**。

`IsAmple M` は各点のまわりに「`X_s` がアフィン」となる切断 `s` を与える（`§9-817`）。
★そこへ段 E3b を当てた形である——**段 E3c が直接使う**。 -/
theorem exists_finset_generating_nonVanishing (f : X ⟶ Spec A) [LocallyOfFiniteType f]
    (M : X.PresheafOfModules) (hM : IsAmple M) (x : X) :
    ∃ (n : ℕ), 0 < n ∧ ∃ s : ((presheafTensorPow M n).obj (op ⊤) : Type),
      x ∈ nonVanishing (presheafTensorPow M n) s ∧
      IsAffineOpen (nonVanishing (presheafTensorPow M n) s) ∧
      ∃ T : Finset (Γ(X, nonVanishing (presheafTensorPow M n) s) : Type),
        @Algebra.adjoin (Γ(Spec A, (⊤ : (Spec A).Opens)) : Type)
          (Γ(X, nonVanishing (presheafTensorPow M n) s) : Type) _ _
          (f.appLE ⊤ (nonVanishing (presheafTensorPow M n) s) (by simp)).hom.toAlgebra
          (T : Set (Γ(X, nonVanishing (presheafTensorPow M n) s) : Type)) = ⊤ := by
  obtain ⟨n, hn, s, hxs, haff⟩ := hM.2 x
  exact ⟨n, hn, s, hxs, haff, exists_finset_generating f haff⟩

/-! ## ★出典の紐付け(`.src`) -/

def finiteType_appLE_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3b——アフィン開の座標環は有限型)",
    sectionId := "genell-prop-1-4" }

def exists_finset_generating.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートの座標環は有限生成)",
    sectionId := "genell-prop-1-4" }

def exists_finset_generating_nonVanishing.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(ample が与えるチャートの座標環は有限生成)",
    sectionId := "genell-prop-1-4" }

def exists_finset_generating.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Scheme.Hom.finiteType_appLE(有限型射のアフィン開への制限)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.finiteType_appLE") 2,
    .citation "[mathlib]" "Algebra.FiniteType.out(⊤ : Subalgebra が有限生成)"
      (.inMathlib "Algebra.FiniteType.out") 2,
    .implicitStep
      ("★底は A ではなく **Γ(Spec A, ⊤)** である(ΓSpecIso で移り合うが定義的には等しくない)。" ++
       "★★本ファイルは Γ(Spec A, ⊤) のまま述べる——同型を挟むと式が読めなくなるうえ、" ++
       "消費側(段 E3c)は Γ の側で働くからである") 4,
    .implicitStep
      ("★★段 E3c はこの生成元それぞれに段 E3a-3(大域化、§9-831)を当て、" ++
       "指数を最大値で揃えることで A⁰_{x_i} → Γ(X, X_{s_i}) の全射性を出す。" ++
       "★そこは本ファイルに無い") 7 ]

end ABC3.Found.GenEll
