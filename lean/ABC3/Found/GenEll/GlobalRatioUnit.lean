/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalAwayHom
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★重なりの上で比は単元である —— 段 E1d の入口（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★これは何か —— 段 E1d の残りを開ける鍵

`§9-848`（`GlobalChartToProj.lean`）でチャートごとの射
`X_{s_i} ⟶ ℙᴺ_R` が入り、`§9-836`（`GlueOpens.lean`）で貼り合わせの機構が入った。
★残っているのは**「チャート射が重なりで一致する」**ことだけである。

★★その証明の骨は「`Away 𝒜 (x_i·x_j)` は `Away 𝒜 x_i` を `x_j/x_i` で局所化したもの」
（mathlib の `HomogeneousLocalization.Away.isLocalization_mul`）であり、
局所化の普遍性を使うには **`x_j/x_i` の像が単元であること**が要る。

★★★`globalAwayHom … (x_j/x_i) = s_j/s_i`（`§9-842`）なので、要るのは

    **`X_{s_i} ⊓ X_{s_j}` の上で `s_j/s_i` は単元**

である。★本ファイルがそれを取る。

## ★★★機構 —— 逆元は反対向きの比

`s_j/s_i` の逆元は `s_i/s_j` である。★どちらも大域の比（`§9-841`）として既にあり、
重なりへ制限して掛けると `1` になる。

★★局所（自明化 `V` を 1 つ固定した上）では

    `(s/t)·t = s`（`sectionRatio_mul`）を 2 回使って `(s/t)·(t/s)·t·s = s·t`

で、`t`・`s` の値は自明化の上で単元（`isUnit_trivValue_res`）だから消せる。
★★★大域へは `TopCat.Sheaf.eq_of_locally_eq'` で上げる
——被覆は「自明化つき開集合」（`TrivIndex`）で、`W ≤ ⨆ j, W ⊓ V_j` は
`iSup_trivIndex` と `Opens` の分配律（`inf_iSup_eq`）から出る。

## ★測定の記録

★★★★**`W = X_t ⊓ X_s` が `TrivIndex` の開で覆われること**は
`iSup_trivIndex`（`X_t` 用）から直接は出ない——`⨆_j V_j` が `⊤` とは限らないからである。
★経由するのは `W ≤ X_t ≤ ⨆_j (X_t ⊓ V_j) ≤ ⨆_j V_j` で、
そのあと `Opens` が frame であること（`inf_iSup_eq`）を使う。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★★★★局所の段 —— 自明化を 1 つ固定した上で -/

/-- ★★★★★★**重なりの上で `(s/t)·(t/s) = 1`**（自明化を 1 つ固定した形）。

★機構は `sectionRatio_mul`（`(s/t)·t = s`）を 2 回使うだけである。
★★`t` の値は `X_t ⊓ V` の上で単元（`isUnit_trivValue_res`）なので消せる。 -/
theorem sectionRatio_mul_inv (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s t : (M.obj (op ⊤) : Type))
    (hle1 : (nonVanishing M t ⊓ nonVanishing M s) ⊓ V ≤ nonVanishing M t ⊓ V)
    (hle2 : (nonVanishing M t ⊓ nonVanishing M s) ⊓ V ≤ nonVanishing M s ⊓ V) :
    X.presheaf.map (homOfLE hle1).op (sectionRatio M V e s t)
      * X.presheaf.map (homOfLE hle2).op (sectionRatio M V e t s) = 1 := by
  set Z : X.Opens := (nonVanishing M t ⊓ nonVanishing M s) ⊓ V with hZ
  have hZV : Z ≤ V := inf_le_right
  set a := X.presheaf.map (homOfLE hle1).op (sectionRatio M V e s t) with ha
  set b := X.presheaf.map (homOfLE hle2).op (sectionRatio M V e t s) with hb
  set p := X.presheaf.map (homOfLE hZV).op (trivValue M V e s) with hp
  set q := X.presheaf.map (homOfLE hZV).op (trivValue M V e t) with hq
  have h1 : a * q = p := by
    have := congrArg (X.presheaf.map (homOfLE hle1).op) (sectionRatio_mul M V e s t)
    rw [map_mul] at this
    rw [ha, hq, hp, ← res_trans hle1 (inf_le_right : nonVanishing M t ⊓ V ≤ V),
      ← res_trans hle1 (inf_le_right : nonVanishing M t ⊓ V ≤ V)]
    exact this
  have h2 : b * p = q := by
    have := congrArg (X.presheaf.map (homOfLE hle2).op) (sectionRatio_mul M V e t s)
    rw [map_mul] at this
    rw [hb, hp, hq, ← res_trans hle2 (inf_le_right : nonVanishing M s ⊓ V ≤ V),
      ← res_trans hle2 (inf_le_right : nonVanishing M s ⊓ V ≤ V)]
    exact this
  have hpu : IsUnit p := by
    have hu := isUnit_trivValue_res M V e s
    have h3 := hu.map (X.presheaf.map (homOfLE hle2).op).hom
    rw [show (X.presheaf.map (homOfLE hle2).op).hom
        (X.presheaf.map (homOfLE (inf_le_right : nonVanishing M s ⊓ V ≤ V)).op
          (trivValue M V e s)) = p from by rw [hp, res_trans]] at h3
    exact h3
  have key : p * (a * b) = p * 1 := by
    calc p * (a * b) = a * (b * p) := by ring
      _ = a * q := by rw [h2]
      _ = p := h1
      _ = p * 1 := by rw [mul_one]
  exact hpu.mul_left_cancel key

/-! ## ★★★★★★★★★★大域の段 -/

/-- ★★★★★★★★★★**重なりの上で `(s/t)·(t/s) = 1`**（大域の比の形）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★被覆は「自明化つき開集合」（`TrivIndex`）で、
`TopCat.Sheaf.eq_of_locally_eq'` で局所から上げる。 -/
theorem globalRatio_mul_globalRatio (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t : (M.obj (op ⊤) : Type)) :
    X.presheaf.map (homOfLE (inf_le_left :
        nonVanishing M t ⊓ nonVanishing M s ≤ nonVanishing M t)).op (globalRatio M hM s t)
      * X.presheaf.map (homOfLE (inf_le_right :
        nonVanishing M t ⊓ nonVanishing M s ≤ nonVanishing M s)).op (globalRatio M hM t s)
      = 1 := by
  set W : X.Opens := nonVanishing M t ⊓ nonVanishing M s with hW
  have hsup : W ≤ ⨆ j : TrivIndex M, W ⊓ j.1 := by
    have h1 : (⨆ j : TrivIndex M, nonVanishing M t ⊓ j.1) = nonVanishing M t :=
      iSup_trivIndex M hM t
    have h3 : W ≤ ⨆ j : TrivIndex M, (j.1 : X.Opens) := by
      refine le_trans (inf_le_left) ?_
      rw [← h1]
      exact iSup_mono fun j => inf_le_right
    have h2 : W ⊓ (⨆ j : TrivIndex M, (j.1 : X.Opens)) = ⨆ j : TrivIndex M, W ⊓ j.1 :=
      inf_iSup_eq _ _
    rw [← h2, inf_eq_left.mpr h3]
  refine TopCat.Sheaf.eq_of_locally_eq' X.sheaf (fun j : TrivIndex M => W ⊓ j.1) W
    (fun j => homOfLE (inf_le_left : W ⊓ j.1 ≤ W)) hsup _ _ ?_
  intro j
  show X.presheaf.map (homOfLE (inf_le_left : W ⊓ j.1 ≤ W)).op _
    = X.presheaf.map (homOfLE (inf_le_left : W ⊓ j.1 ≤ W)).op (1 : (Γ(X, W) : Type))
  rw [map_mul, map_one, res_trans, res_trans]
  rw [← res_trans (show W ⊓ j.1 ≤ nonVanishing M t ⊓ j.1 from
      inf_le_inf_right _ inf_le_left) (inf_le_left : nonVanishing M t ⊓ j.1 ≤ nonVanishing M t),
    ← res_trans (show W ⊓ j.1 ≤ nonVanishing M s ⊓ j.1 from
      inf_le_inf_right _ inf_le_right) (inf_le_left : nonVanishing M s ⊓ j.1 ≤ nonVanishing M s),
    globalRatio_res, globalRatio_res]
  exact sectionRatio_mul_inv M j.1 j.2 s t _ _

/-- ★★★★★★★★★★★**重なりの上で比は単元である** —— 段 E1d の入口。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★これが `HomogeneousLocalization.Away.isLocalization_mul` の普遍性に渡す条件である
——`Away 𝒜 (x_i·x_j)` は `Away 𝒜 x_i` を `x_j/x_i` で局所化したものであり、
その像が単元なら重なりへの環準同型が**一意に**延びる。 -/
theorem isUnit_globalRatio_res (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t : (M.obj (op ⊤) : Type)) :
    IsUnit (X.presheaf.map (homOfLE (inf_le_left :
      nonVanishing M t ⊓ nonVanishing M s ≤ nonVanishing M t)).op (globalRatio M hM s t)) :=
  IsUnit.of_mul_eq_one _ (globalRatio_mul_globalRatio M hM s t)

/-! ## ★出典の紐付け(`.src`) -/

def sectionRatio_mul_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの上で (s/t)·(t/s) = 1——自明化を固定した形)",
    sectionId := "genell-prop-1-4" }

def globalRatio_mul_globalRatio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの上で (s/t)·(t/s) = 1——大域の比)",
    sectionId := "genell-prop-1-4" }

def isUnit_globalRatio_res.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの上で比は単元である——段 E1d の入口)",
    sectionId := "genell-prop-1-4" }

def isUnit_globalRatio_res.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "sectionRatio_mul((s/t)·t = s、段 D3)"
      (.inProject "ABC3" "ABC3.Found.GenEll.sectionRatio_mul") 1,
    .citation "[ABC3]" "globalRatio_res(大域の比はチャートで sectionRatio に戻る、§9-841)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalRatio_res") 1,
    .citation "[mathlib]" "TopCat.Sheaf.eq_of_locally_eq'(局所で等しければ等しい)"
      (.inMathlib "TopCat.Sheaf.eq_of_locally_eq'") 1,
    .implicitStep
      ("★★★★測定: W = X_t ⊓ X_s が TrivIndex の開で覆われることは " ++
       "iSup_trivIndex(X_t 用)から直接は出ない——⨆_j V_j が ⊤ とは限らないからである。" ++
       "経由するのは W ≤ X_t ≤ ⨆_j (X_t ⊓ V_j) ≤ ⨆_j V_j で、" ++
       "そのあと Opens が frame であること(inf_iSup_eq)を使う") 2,
    .implicitStep
      ("★★これで HomogeneousLocalization.Away.isLocalization_mul の普遍性に渡す条件が揃った。" ++
       "★次は Away 𝒜 (x_i·x_j) →+* Γ(X, X_{s_i} ⊓ X_{s_j}) を局所化の lift で作り、" ++
       "Proj.SpecMap_awayMap_awayι で 2 つのチャート射を同じ awayι へ落とす段である") 4 ]

end ABC3.Found.GenEll
