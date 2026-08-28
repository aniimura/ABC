/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlueBump
import ABC3.Found.GenEll.AmpleDef
import ABC3.Meta.Claim

/-!
# ★★★★★★★★大域の比 `s/t ∈ Γ(X, X_t)` —— 「`⊓ V`」を外す道具（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★これは何か —— 2026-08-28 に見つかった段

`§9-815` の `awayToSectionHom` は

    `A⁰_{x_i} →+* Γ(X, X_{s_i} ⊓ V)`

という形である。★**`⊓ V` が付いている**——`trivValue` を使うので自明化する開 `V` に閉じ込められる。

★★しかし段 E3d（`§9-834`）が消費するのは **`Γ(X, X_{s_i})`（`⊓ V` なし）** である
——`IsClosedImmersion.of_surjective_of_isAffine` はチャート **`X_{s_i}` 全体**の `Γ` を見るからである。

★★★**これは 2026-08-28 に測って初めて見えた段である**。台帳の段 E1c・E3c は
どちらも `⊓ V` の形で書かれていて、そこが噛み合っていなかった。

## ★★★★★★正しい道 —— 貼り合わせた比を使う

`M` を自明化しなくても、**比 `s/t` は `X_t` 全体で定まる**。
★段 D4（`§9-AmpleDef` の `exists_glued_ratio`）が既にそれを持っていた
——`sectionRatio` の族は重なりで一致する（`sectionRatio_agree`、段 D3）ので、
`𝒪_X` が層であることから貼り合う。

★★本ファイルはその貼り合わせを**名前のある対象** `globalRatio` にし、
チャートでの値（`globalRatio| = sectionRatio`）と**一意性**を付けたものである。

## ★★★残っている段（明示）

★★★★`homogValue` / `homogRatio` / `awayToSection` を **`globalRatio` で組み直す**段が要る。
★機構は同じ（多項式を切断で評価し、`s_i` の冪で割る）が、
**`trivValue` の代わりに `globalRatio` を使う**ので `⊓ V` が消える。
★★そうすれば段 E3c の全射性が `Γ(X, X_{s_i})` について言え、段 E3d に直結する。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★自分自身への制限は恒等 -/

/-- ★**同じ開への制限は恒等である**（構造層）。 -/
theorem res_self {A : X.Opens} (h : A ≤ A) (r : (Γ(X, A) : Type)) :
    X.presheaf.map (homOfLE h).op r = r := by
  rw [Subsingleton.elim (homOfLE h) (𝟙 A), op_id, X.presheaf.map_id]
  rfl

/-! ## ★★★★★★★★大域の比 -/

/-- ★★★★★★★★**大域の比 `s/t ∈ Γ(X, X_t)`**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★段 D4（`exists_glued_ratio`）が与える貼り合わせに名前を付け、
`iSup_trivIndex`（自明化つき開は `X_t` を覆う）で `Γ(X, X_t)` へ移したものである。

★★**`⊓ V` が付かない**のが要点である——これが段 E3d の消費する形である。 -/
noncomputable def globalRatio (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t : (M.obj (op ⊤) : Type)) : (Γ(X, nonVanishing M t) : Type) :=
  X.presheaf.map (homOfLE (le_of_eq (iSup_trivIndex M hM t).symm)).op
    (exists_glued_ratio M s t).choose

/-- ★★★★★★**大域の比はチャートで `sectionRatio` に戻る**。

★これが `globalRatio` が「比」であることの中身である。 -/
theorem globalRatio_res (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t : (M.obj (op ⊤) : Type)) (i : TrivIndex M) :
    X.presheaf.map (homOfLE (inf_le_left :
        nonVanishing M t ⊓ i.1 ≤ nonVanishing M t)).op (globalRatio M hM s t)
      = sectionRatio M i.1 i.2 s t := by
  have hg := (exists_glued_ratio M s t).choose_spec.1 i
  rw [globalRatio, res_trans]
  show X.presheaf.map (homOfLE _).op (exists_glued_ratio M s t).choose = _
  rw [show (homOfLE (le_trans (inf_le_left : nonVanishing M t ⊓ i.1 ≤ nonVanishing M t)
      (le_of_eq (iSup_trivIndex M hM t).symm)) :
      (nonVanishing M t ⊓ i.1) ⟶ (⨆ j : TrivIndex M, nonVanishing M t ⊓ j.1))
    = TopologicalSpace.Opens.leSupr (fun j : TrivIndex M => nonVanishing M t ⊓ j.1) i from
      Subsingleton.elim _ _]
  exact hg

/-- ★★★★★★★**チャートでの値が比なら大域の比である**（一意性）。

★これで「別の作り方をした大域切断が `globalRatio` に等しい」ことを言える
——段 E3c を `⊓ V` なしで組み直すときに要る。 -/
theorem globalRatio_unique (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t : (M.obj (op ⊤) : Type)) (r : (Γ(X, nonVanishing M t) : Type))
    (hr : ∀ i : TrivIndex M,
      X.presheaf.map (homOfLE (inf_le_left :
          nonVanishing M t ⊓ i.1 ≤ nonVanishing M t)).op r
        = sectionRatio M i.1 i.2 s t) :
    r = globalRatio M hM s t := by
  have hle : (⨆ j : TrivIndex M, nonVanishing M t ⊓ j.1) ≤ nonVanishing M t :=
    le_of_eq (iSup_trivIndex M hM t)
  have hglue : TopCat.Presheaf.IsGluing X.sheaf.obj
      (fun j : TrivIndex M => nonVanishing M t ⊓ j.1) (ratioFamily M s t)
      (X.presheaf.map (homOfLE hle).op r) := by
    intro i
    show X.presheaf.map (TopologicalSpace.Opens.leSupr _ i).op
      (X.presheaf.map (homOfLE hle).op r) = _
    rw [show TopologicalSpace.Opens.leSupr
        (fun j : TrivIndex M => nonVanishing M t ⊓ j.1) i
      = homOfLE (le_iSup (fun j : TrivIndex M => nonVanishing M t ⊓ j.1) i) from
        Subsingleton.elim _ _, res_trans]
    show X.presheaf.map (homOfLE _).op r = _
    exact hr i
  have heq := (exists_glued_ratio M s t).choose_spec.2 _ hglue
  rw [globalRatio, ← heq, res_trans, res_self]

/-! ## ★出典の紐付け(`.src`) -/

def globalRatio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(大域の比 s/t ∈ Γ(X, X_t)——⊓ V なし)",
    sectionId := "genell-prop-1-4" }

def globalRatio_res.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(大域の比はチャートで sectionRatio に戻る)",
    sectionId := "genell-prop-1-4" }

def globalRatio_unique.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートでの値が比なら大域の比である)",
    sectionId := "genell-prop-1-4" }

def globalRatio.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_glued_ratio(比は X_t の上へ貼り合う、段 D4)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_glued_ratio") 2,
    .citation "[ABC3]" "iSup_trivIndex(自明化つき開は X_t を覆う)"
      (.inProject "ABC3" "ABC3.Found.GenEll.iSup_trivIndex") 2,
    .implicitStep
      ("★★★**2026-08-28 に測って見えた段**: §9-815 の awayToSectionHom は " ++
       "A⁰_{x_i} →+* Γ(X, X_{s_i} ⊓ V) という形で、**⊓ V が付いている**" ++
       "(trivValue を使うので自明化する開に閉じ込められる)。" ++
       "★しかし段 E3d(§9-834)が消費するのは Γ(X, X_{s_i})(⊓ V なし)である" ++
       "——IsClosedImmersion.of_surjective_of_isAffine はチャート全体の Γ を見るから。" ++
       "★★台帳の段 E1c・E3c はどちらも ⊓ V の形で書かれていて、そこが噛み合っていなかった") 7,
    .implicitStep
      ("★★★★homogValue / homogRatio / awayToSection を **globalRatio で組み直す**段が要る。" ++
       "機構は同じ(多項式を切断で評価し s_i の冪で割る)が、trivValue の代わりに " ++
       "globalRatio を使うので ⊓ V が消える。★そうすれば段 E3c の全射性が " ++
       "Γ(X, X_{s_i}) について言え、段 E3d に直結する") 7 ]

end ABC3.Found.GenEll
