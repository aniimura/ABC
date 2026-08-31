/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalRatioUnit
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★比のコサイクル則 `s/u = (s/t)·(t/u)` —— 段 E1d の数学（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★これは何か —— 「チャートが重なりで一致する」の中身

段 E1d（チャート射 `X_{s_i} ⟶ ℙᴺ_R` を貼る）で示すべきことは
**2 つのチャート射が `X_{s_i} ⊓ X_{s_j}` の上で一致する**ことである。

★★★★その数学的な中身は 1 行に書ける:

    `s_k/s_j = (s_k/s_i)·(s_i/s_j)`   （`X_{s_i} ⊓ X_{s_j}` の上で）

——**比のコサイクル則**である。★本ファイルがそれを取る。

## ★★★機構 —— `sectionRatio_mul` を 3 回

局所（自明化 `V` を固定）では `(s/t)·t = s`（`sectionRatio_mul`）を 3 回使い、
`u` の値が単元（`isUnit_trivValue_res`）なので消せる:

    `(s/u)·u = s`、`(s/t)·t = s`、`(t/u)·u = t`
    ⟹ `u·(s/u) = s = (s/t)·t = (s/t)·(t/u)·u = u·((s/t)·(t/u))`

★大域へは `TopCat.Sheaf.eq_of_locally_eq'` で上げる（`§9-907` と同じ被覆）。

## ★被覆の補題を切り出した

`§9-907` では `W = X_t ⊓ X_s` の被覆を証明の中でその都度作っていたが、
★**`W ≤ X_t` でありさえすれば同じ議論が通る**ので `trivIndex_covers` として切り出した。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★被覆の補題 -/

/-- ★★★**`X_t` に含まれる開は自明化つき開で覆われる**。

★`iSup_trivIndex`（`X_t` 用）から**直接は出ない**——`⨆_j V_j` が `⊤` とは限らないからである。
経由するのは `W ≤ X_t ≤ ⨆_j (X_t ⊓ V_j) ≤ ⨆_j V_j` で、
そのあと `Opens` が frame であること（`inf_iSup_eq`）を使う。 -/
theorem trivIndex_covers (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (t : (M.obj (op ⊤) : Type)) (W : X.Opens) (hW : W ≤ nonVanishing M t) :
    W ≤ ⨆ j : TrivIndex M, W ⊓ j.1 := by
  have h1 : (⨆ j : TrivIndex M, nonVanishing M t ⊓ j.1) = nonVanishing M t :=
    iSup_trivIndex M hM t
  have h3 : W ≤ ⨆ j : TrivIndex M, (j.1 : X.Opens) := by
    refine le_trans hW ?_
    rw [← h1]
    exact iSup_mono fun j => inf_le_right
  have h2 : W ⊓ (⨆ j : TrivIndex M, (j.1 : X.Opens)) = ⨆ j : TrivIndex M, W ⊓ j.1 :=
    inf_iSup_eq _ _
  rw [← h2, inf_eq_left.mpr h3]

/-! ## ★★★★★★★★局所のコサイクル則 -/

/-- ★★★★★★★★**`s/u = (s/t)·(t/u)`**（自明化を 1 つ固定した形）。

★`sectionRatio_mul`（`(s/t)·t = s`）を 3 回使い、
`u` の値が単元（`isUnit_trivValue_res`）なので消せる。 -/
theorem sectionRatio_cocycle (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s t u : (M.obj (op ⊤) : Type))
    (h1 : (nonVanishing M t ⊓ nonVanishing M u) ⊓ V ≤ nonVanishing M t ⊓ V)
    (h2 : (nonVanishing M t ⊓ nonVanishing M u) ⊓ V ≤ nonVanishing M u ⊓ V) :
    X.presheaf.map (homOfLE h2).op (sectionRatio M V e s u)
      = X.presheaf.map (homOfLE h1).op (sectionRatio M V e s t)
        * X.presheaf.map (homOfLE h2).op (sectionRatio M V e t u) := by
  set Z : X.Opens := (nonVanishing M t ⊓ nonVanishing M u) ⊓ V with hZ
  have hZV : Z ≤ V := inf_le_right
  set p := X.presheaf.map (homOfLE hZV).op (trivValue M V e s) with hp
  set q := X.presheaf.map (homOfLE hZV).op (trivValue M V e t) with hq
  set r := X.presheaf.map (homOfLE hZV).op (trivValue M V e u) with hr
  set A := X.presheaf.map (homOfLE h2).op (sectionRatio M V e s u) with hA
  set B := X.presheaf.map (homOfLE h1).op (sectionRatio M V e s t) with hB
  set C := X.presheaf.map (homOfLE h2).op (sectionRatio M V e t u) with hC
  have eA : A * r = p := by
    have := congrArg (X.presheaf.map (homOfLE h2).op) (sectionRatio_mul M V e s u)
    rw [map_mul] at this
    rw [hA, hr, hp, ← res_trans h2 (inf_le_right : nonVanishing M u ⊓ V ≤ V),
      ← res_trans h2 (inf_le_right : nonVanishing M u ⊓ V ≤ V)]
    exact this
  have eB : B * q = p := by
    have := congrArg (X.presheaf.map (homOfLE h1).op) (sectionRatio_mul M V e s t)
    rw [map_mul] at this
    rw [hB, hq, hp, ← res_trans h1 (inf_le_right : nonVanishing M t ⊓ V ≤ V),
      ← res_trans h1 (inf_le_right : nonVanishing M t ⊓ V ≤ V)]
    exact this
  have eC : C * r = q := by
    have := congrArg (X.presheaf.map (homOfLE h2).op) (sectionRatio_mul M V e t u)
    rw [map_mul] at this
    rw [hC, hr, hq, ← res_trans h2 (inf_le_right : nonVanishing M u ⊓ V ≤ V),
      ← res_trans h2 (inf_le_right : nonVanishing M u ⊓ V ≤ V)]
    exact this
  have hru : IsUnit r := by
    have hu := isUnit_trivValue_res M V e u
    have h3 := hu.map (X.presheaf.map (homOfLE h2).op).hom
    rw [show (X.presheaf.map (homOfLE h2).op).hom
        (X.presheaf.map (homOfLE (inf_le_right : nonVanishing M u ⊓ V ≤ V)).op
          (trivValue M V e u)) = r from by rw [hr, res_trans]] at h3
    exact h3
  refine hru.mul_left_cancel ?_
  calc r * A = A * r := by ring
    _ = p := eA
    _ = B * q := eB.symm
    _ = B * (C * r) := by rw [eC]
    _ = r * (B * C) := by ring

/-! ## ★★★★★★★★★★★★大域のコサイクル則 -/

/-- ★★★★★★★★★★★★**`s/u = (s/t)·(t/u)`**（大域の比、`X_t ⊓ X_u` の上で）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★これが段 E1d の**数学的な中身**である
——「2 つのチャートの座標が重なりで遷移函数で移り合う」ことそのもの。 -/
theorem globalRatio_cocycle (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t u : (M.obj (op ⊤) : Type)) :
    X.presheaf.map (homOfLE (inf_le_right :
        nonVanishing M t ⊓ nonVanishing M u ≤ nonVanishing M u)).op (globalRatio M hM s u)
      = X.presheaf.map (homOfLE (inf_le_left :
          nonVanishing M t ⊓ nonVanishing M u ≤ nonVanishing M t)).op (globalRatio M hM s t)
        * X.presheaf.map (homOfLE (inf_le_right :
          nonVanishing M t ⊓ nonVanishing M u ≤ nonVanishing M u)).op (globalRatio M hM t u) := by
  set W : X.Opens := nonVanishing M t ⊓ nonVanishing M u with hW
  have hsup : W ≤ ⨆ j : TrivIndex M, W ⊓ j.1 :=
    trivIndex_covers M hM t W inf_le_left
  refine TopCat.Sheaf.eq_of_locally_eq' X.sheaf (fun j : TrivIndex M => W ⊓ j.1) W
    (fun j => homOfLE (inf_le_left : W ⊓ j.1 ≤ W)) hsup _ _ ?_
  intro j
  show X.presheaf.map (homOfLE (inf_le_left : W ⊓ j.1 ≤ W)).op _
    = X.presheaf.map (homOfLE (inf_le_left : W ⊓ j.1 ≤ W)).op _
  rw [map_mul, res_trans, res_trans, res_trans]
  rw [← res_trans (show W ⊓ j.1 ≤ nonVanishing M u ⊓ j.1 from
      inf_le_inf_right _ inf_le_right) (inf_le_left : nonVanishing M u ⊓ j.1 ≤ nonVanishing M u),
    ← res_trans (show W ⊓ j.1 ≤ nonVanishing M t ⊓ j.1 from
      inf_le_inf_right _ inf_le_left) (inf_le_left : nonVanishing M t ⊓ j.1 ≤ nonVanishing M t),
    ← res_trans (show W ⊓ j.1 ≤ nonVanishing M u ⊓ j.1 from
      inf_le_inf_right _ inf_le_right) (inf_le_left : nonVanishing M u ⊓ j.1 ≤ nonVanishing M u),
    globalRatio_res, globalRatio_res, globalRatio_res]
  exact sectionRatio_cocycle M j.1 j.2 s t u _ _

/-! ## ★出典の紐付け(`.src`) -/

def trivIndex_covers.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(X_t に含まれる開は自明化つき開で覆われる)",
    sectionId := "genell-prop-1-4" }

def sectionRatio_cocycle.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(比のコサイクル則 s/u = (s/t)·(t/u)——自明化を固定した形)",
    sectionId := "genell-prop-1-4" }

def globalRatio_cocycle.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(比のコサイクル則 s/u = (s/t)·(t/u)——大域の比)",
    sectionId := "genell-prop-1-4" }

def globalRatio_cocycle.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "sectionRatio_mul((s/t)·t = s、段 D3)"
      (.inProject "ABC3" "ABC3.Found.GenEll.sectionRatio_mul") 1,
    .citation "[ABC3]" "globalRatio_res(大域の比はチャートで sectionRatio に戻る、§9-841)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalRatio_res") 1,
    .citation "[mathlib]" "TopCat.Sheaf.eq_of_locally_eq'(局所で等しければ等しい)"
      (.inMathlib "TopCat.Sheaf.eq_of_locally_eq'") 1,
    .implicitStep
      ("★★これが段 E1d(チャート射が重なりで一致する)の**数学的な中身**である" ++
       "——「2 つのチャートの座標が重なりで遷移函数で移り合う」ことそのもの") 3,
    .implicitStep
      ("★配管の段が残る: Away 𝒜 (x_i·x_j) →+* Γ(X, X_{s_i} ⊓ X_{s_j}) を" ++
       "局所化の lift(Away.isLocalization_mul + isUnit_globalRatio_res、§9-907)で作り、" ++
       "Proj.SpecMap_awayMap_awayι で 2 つのチャート射を同じ awayι へ落とす") 4 ]

end ABC3.Found.GenEll
