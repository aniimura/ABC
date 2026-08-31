/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalRatioCocycle
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★比が消えない所は `X_{s} ∩ X_{t}` である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★これは何か —— 段 E1d の後の「像の同定」

`§9-911` で `ψ : X ⟶ ℙᴺ_R` が入った。★次に要るのは

    **`ψ⁻¹(D₊(x_i)) = X_{s_i}`**

である——これがあれば `IsImmersion` が**target について局所的**であることから
「チャートごとに埋め込み（`§9-848`）」が「`ψ` が埋め込み」に上がる。

★★その計算の核が本ファイル:

    **`X.basicOpen (s/t) = X_s ⊓ X_t`**

——比 `s/t ∈ Γ(X, X_t)` が消えない所はちょうど `X_s ∩ X_t` である。

## ★★★機構 —— 局所では単元倍で消える

自明化 `V` を固定すると `s/t = trivValue(s)·trivValue(t)⁻¹` なので

    `basicOpen(s/t) = basicOpen(trivValue s) ⊓ (X_t ⊓ V) = (X_s ⊓ V) ⊓ (X_t ⊓ V)`

（`nonVanishing_inf`、段 D2）。★大域へは `basicOpen` が制限と両立すること
（`Scheme.basicOpen_res`）と `Opens` の frame 性（`inf_iSup_eq`）で上げる。

## ★測定の記録

★★★★**束の等式は `ext x; simp; tauto` が速い**。
`⊓` の結合・可換・冪等を `le_antisymm` で手で書くと 10 行になるが、
`TopologicalSpace.Opens.coe_inf` で集合に落として `tauto` に投げると 3 行で済む
（2026-08-28 実測）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★★★★局所の段 -/

/-- ★★★★★★**`basicOpen(s/t) = (X_s ⊓ X_t) ⊓ V`**（自明化を 1 つ固定した形）。

★`s/t = trivValue(s)·trivValue(t)⁻¹` で、単元の `basicOpen` は全体なので消える。 -/
theorem basicOpen_sectionRatio (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s t : (M.obj (op ⊤) : Type)) :
    X.basicOpen (sectionRatio M V e s t)
      = (nonVanishing M s ⊓ nonVanishing M t) ⊓ V := by
  rw [sectionRatio, X.basicOpen_mul,
    X.basicOpen_of_isUnit (isUnit_trivValue_res M V e t).unit⁻¹.isUnit,
    X.basicOpen_res, ← nonVanishing_inf M V e s]
  ext x
  simp only [TopologicalSpace.Opens.coe_inf, Set.mem_inter_iff]
  tauto

/-! ## ★★★★★★★★★★★★大域の段 -/

/-- ★★★★★★★★★★★★**`X.basicOpen (s/t) = X_s ⊓ X_t`**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★比 `s/t ∈ Γ(X, X_t)` が消えない所はちょうど `X_s ∩ X_t` である。
★★これが `ψ⁻¹(D₊(x_i)) = X_{s_i}` の計算の核になる。 -/
theorem basicOpen_globalRatio (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t : (M.obj (op ⊤) : Type)) :
    X.basicOpen (globalRatio M hM s t) = nonVanishing M s ⊓ nonVanishing M t := by
  set A : X.Opens := X.basicOpen (globalRatio M hM s t) with hA
  set B : X.Opens := nonVanishing M s ⊓ nonVanishing M t with hB
  have hAle : A ≤ nonVanishing M t := X.basicOpen_le _
  have hBle : B ≤ nonVanishing M t := inf_le_right
  have key : ∀ j : TrivIndex M, A ⊓ (nonVanishing M t ⊓ j.1) = B ⊓ (nonVanishing M t ⊓ j.1) := by
    intro j
    have h1 : X.basicOpen (X.presheaf.map (homOfLE (inf_le_left :
        nonVanishing M t ⊓ j.1 ≤ nonVanishing M t)).op (globalRatio M hM s t))
        = (nonVanishing M t ⊓ j.1) ⊓ A := X.basicOpen_res _ _
    rw [globalRatio_res M hM s t j, basicOpen_sectionRatio] at h1
    rw [inf_comm A _, ← h1, hB]
    ext x
    simp only [TopologicalSpace.Opens.coe_inf, Set.mem_inter_iff]
    tauto
  have hsup : (⨆ j : TrivIndex M, nonVanishing M t ⊓ j.1) = nonVanishing M t :=
    iSup_trivIndex M hM t
  calc A = A ⊓ nonVanishing M t := (inf_eq_left.mpr hAle).symm
    _ = A ⊓ ⨆ j : TrivIndex M, nonVanishing M t ⊓ j.1 := by rw [hsup]
    _ = ⨆ j : TrivIndex M, A ⊓ (nonVanishing M t ⊓ j.1) := inf_iSup_eq _ _
    _ = ⨆ j : TrivIndex M, B ⊓ (nonVanishing M t ⊓ j.1) := iSup_congr key
    _ = B ⊓ ⨆ j : TrivIndex M, nonVanishing M t ⊓ j.1 := (inf_iSup_eq _ _).symm
    _ = B ⊓ nonVanishing M t := by rw [hsup]
    _ = B := inf_eq_left.mpr hBle

/-- ★★★★**チャートの中で見た形** —— `X_t` の中の基本開集合として。

★`Scheme.Opens.ι_image_basicOpen_topIso_inv`（mathlib）そのものであるが、
右辺が `X_s ⊓ X_t` になるのが本ファイルの内容である。 -/
theorem image_basicOpen_globalRatio (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t : (M.obj (op ⊤) : Type)) :
    (nonVanishing M t).ι ''ᵁ ((nonVanishing M t).toScheme.basicOpen
        ((nonVanishing M t).topIso.inv (globalRatio M hM s t)))
      = nonVanishing M s ⊓ nonVanishing M t := by
  rw [Scheme.Opens.ι_image_basicOpen_topIso_inv, basicOpen_globalRatio]

/-! ## ★出典の紐付け(`.src`) -/

def basicOpen_sectionRatio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(basicOpen(s/t) = (X_s ⊓ X_t) ⊓ V——自明化を固定した形)",
    sectionId := "genell-prop-1-4" }

def basicOpen_globalRatio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(比が消えない所は X_s ∩ X_t である)",
    sectionId := "genell-prop-1-4" }

def image_basicOpen_globalRatio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(比が消えない所——チャートの中で見た形)",
    sectionId := "genell-prop-1-4" }

def basicOpen_globalRatio.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "nonVanishing_inf(X_s ⊓ V = basicOpen(trivValue s)、段 D2)"
      (.inProject "ABC3" "ABC3.Found.GenEll.nonVanishing_inf") 1,
    .citation "[ABC3]" "globalRatio_res(大域の比はチャートで sectionRatio に戻る、§9-841)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalRatio_res") 1,
    .citation "[mathlib]" "Scheme.basicOpen_res・Scheme.basicOpen_mul・basicOpen_of_isUnit"
      (.inMathlib "AlgebraicGeometry.Scheme.basicOpen_res") 1,
    .implicitStep
      ("★★次に要るのは ψ⁻¹(D₊(x_i)) = X_{s_i} である" ++
       "——これがあれば IsImmersion が target について局所的であることから" ++
       "「チャートごとに埋め込み(§9-848)」が「ψ が埋め込み」に上がる。" ++
       "本ファイルはその計算の核である") 4,
    .implicitStep
      ("★★★★測定: 束の等式は ext x; simp; tauto が速い。" ++
       "⊓ の結合・可換・冪等を le_antisymm で手で書くと 10 行になるが、" ++
       "TopologicalSpace.Opens.coe_inf で集合に落として tauto に投げると 3 行で済む") 1 ]

end ABC3.Found.GenEll
