/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalAwayHom
import ABC3.Found.GenEll.ChartSurjectiveCoords
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★段 E3c（大域版・完成形）—— `⊓ V` なしの全射判定（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★これは何か

`§9-838` は段 E3c を「有限個の座標条件」に落としたが、**`⊓ V` が付いた形**であった。
`§9-841` でその欠陥を測り、`§9-842` で `⊓ V` の無いチャート写像
`globalAwayHom : A⁰_{x_i} →+* Γ(X, X_{s_i})` を作った。

★本ファイルは**その大域版で全射判定を言い直す**:

> `f : X ⟶ Spec A` が有限型、`X_{s_i}` がアフィン、`φ` が底環の像を覆うとき、
> ★**ある有限集合 `T`** があって、`T` の各元が `s_j/s_i` の形なら
> `A⁰_{x_i} → Γ(X, X_{s_i})` は**全射**である。

★★これが段 E3d（`§9-834`）の消費する形である
——`IsClosedImmersion.of_surjective_of_isAffine` はチャート `X_{s_i}` **全体**の `Γ` を見る。

## ★★★機構 —— `§9-833` はそのまま使える

`§9-833` の `exists_finset_surjective_criterion`（有限生成 ⟹ 全射判定）は
**環準同型の形しか見ていない**ので、`⊓ V` の有無に依らない。
★取り替えるのは 2 本だけである:

| 取り替え | 旧（`⊓ V` 付き） | 新（大域） |
|---|---|---|
| 比の判定 | `mem_range_awayToSectionHom_of_mul_eq`（`§9-837`） | `mem_range_globalAwayHom`（`§9-842`） |
| 底環の像 | `range_appLE_subset`（`§9-838`） | `range_appLE_subset_global`（本ファイル） |

## ★残っている段（明示）

★★**座標の選び方**は依然として残る——`§9-839`・`§9-840` が作る `t` を座標族に並べる帳簿である。
★★★`§9-839` は `trivValue` の言葉で条件を出しているので、
**`globalRatio` の言葉へ翻訳する段**が要る（`globalRatio_unique` で `t` を同定する）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}} {A : CommRingCat.{0}}

/-! ## ★★★★係数環の像は自動で入る（大域版） -/

/-- ★★★★**定数の像は `φ r` の制限である**（大域版）。 -/
theorem globalAwayHom_awayConst (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) (r : R) :
    globalAwayHom M hM φ s i (awayConst N R i r)
      = X.presheaf.map (homOfLE (le_top : nonVanishing M (s i) ≤ ⊤)).op (φ r) := by
  show globalHomogHom M hM φ s i (MvPolynomial.C r) = _
  rw [globalHomogHom, MvPolynomial.eval₂Hom_C]
  rfl

/-- ★★★★★**底環の像は `A⁰_{x_i}` の像に含まれる**（大域版）。 -/
theorem range_appLE_subset_global (f : X ⟶ Spec A)
    (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (hφ : ∀ r : (Γ(Spec A, (⊤ : (Spec A).Opens)) : Type),
      f.appLE ⊤ ⊤ (by simp) r ∈ Set.range φ) :
    Set.range (f.appLE ⊤ (nonVanishing M (s i)) (by simp)).hom
      ⊆ Set.range (globalAwayHom M hM φ s i) := by
  rintro _ ⟨r, rfl⟩
  obtain ⟨r', hr'⟩ := hφ r
  refine ⟨awayConst N R i r', ?_⟩
  rw [globalAwayHom_awayConst, hr']
  exact appLE_res f ⊤ ⊤ (nonVanishing M (s i)) (by simp) le_top r

/-! ## ★★★★★★★★★★段 E3c（大域版・完成形） -/

/-- ★★★★★★★★★★**段 E3c（大域版・完成形）** —— `⊓ V` なしの全射判定。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`f : X ⟶ Spec A` が有限型、`X_{s_i}` がアフィン、`φ` が底環の像を覆うとき、
★**ある有限集合 `T`** があって

> `T` の各元 `g` が `g = s_j/s_i` の形になる座標 `j` を持つ

なら `A⁰_{x_i} → Γ(X, X_{s_i})` は**全射**である。

★★これが段 E3d（`§9-834`）に**そのまま渡せる形**である
——`IsClosedImmersion.of_surjective_of_isAffine` が要求するのは
チャート `X_{s_i}` **全体**の `Γ` についての全射性だからである。 -/
theorem exists_finset_surjective_globalAwayHom (f : X ⟶ Spec A) [LocallyOfFiniteType f]
    (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (haff : IsAffineOpen (nonVanishing M (s i)))
    (hφ : ∀ r : (Γ(Spec A, (⊤ : (Spec A).Opens)) : Type),
      f.appLE ⊤ ⊤ (by simp) r ∈ Set.range φ) :
    ∃ T : Finset (Γ(X, nonVanishing M (s i)) : Type),
      (∀ g ∈ T, ∃ j : Fin (N + 1), g = globalRatio M hM (s j) (s i)) →
        Function.Surjective (globalAwayHom M hM φ s i) := by
  obtain ⟨T, hT⟩ := exists_finset_surjective_criterion f haff
  refine ⟨T, fun hgen => hT _ (range_appLE_subset_global f M hM φ s i hφ) ?_⟩
  intro g hg
  obtain ⟨j, hj⟩ := hgen g hg
  exact mem_range_globalAwayHom M hM φ s i j g hj

/-! ## ★出典の紐付け(`.src`) -/

def globalAwayHom_awayConst.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(定数の像は φ r の制限——大域版)",
    sectionId := "genell-prop-1-4" }

def range_appLE_subset_global.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(底環の像は A⁰_{x_i} の像に含まれる——大域版)",
    sectionId := "genell-prop-1-4" }

def exists_finset_surjective_globalAwayHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3c 大域版完成形——⊓ V なしの全射判定)",
    sectionId := "genell-prop-1-4" }

def exists_finset_surjective_globalAwayHom.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_finset_surjective_criterion(有限生成 ⟹ 全射判定、§9-833)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_finset_surjective_criterion") 2,
    .citation "[ABC3]" "mem_range_globalAwayHom(大域の比は像に入る、§9-842)"
      (.inProject "ABC3" "ABC3.Found.GenEll.mem_range_globalAwayHom") 2,
    .citation "[ABC3]" "IsClosedImmersion.of_surjective_of_isAffine の消費側(段 E3d、§9-834)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isClosedImmersion_restrict_of_surjective") 2,
    .implicitStep
      ("★§9-833 の全射判定は**環準同型の形しか見ていない**ので ⊓ V の有無に依らない。" ++
       "取り替えたのは 2 本だけである——比の判定(§9-837 → §9-842)と底環の像(§9-838 → 本ファイル)") 3,
    .implicitStep
      ("★★**座標の選び方**は依然として残る——§9-839・§9-840 が作る t を座標族に並べる帳簿である。" ++
       "★★★§9-839 は trivValue の言葉で条件を出しているので、" ++
       "**globalRatio の言葉へ翻訳する段**が要る(globalRatio_unique で t を同定する)") 6 ]

end ABC3.Found.GenEll
