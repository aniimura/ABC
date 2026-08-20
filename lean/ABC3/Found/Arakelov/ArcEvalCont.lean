import ABC3.Found.Arakelov.ArcLiftOpen

/-!
# Arakelov (C3) 第 260 ブロック —— ★★★★★**開集合上の評価は連続**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★連続性が閉じた

第 258 で `genNorm (arcEvalOnTop (c·g)) = ‖evalOn c‖` が出ていた。
★本ブロックで **`evalOn` の連続性**が入り、**`genNorm` の連続性が完成**する。

## ★★機構 —— `appLE` に直すと合成則が使える

    evalOn p V h c = ΓSpecIso ((p.appLE V ⊤ _) c)          (★`rfl`)

★★これで mathlib の `appLE_comp_appLE` が使え、`p = q ≫ V.ι` のとき

    evalOn (q ≫ V.ι) V _ c = evalGlobal q ((V.ι.appLE V ⊤ _) c)

となる——右辺は第 252 でそのまま連続である。

★★★`V.ι.appLE V ⊤ _` は mathlib で **`IsIso`** と登録されている
(`Scheme.Opens.ι_preimage_self` の直後)——`V` 上の関数と `V.toScheme` 上の関数の対応である。

| 定理 | 内容 |
|---|---|
| `evalOn_appLE` | ★`evalOn` は `appLE` で書ける(`rfl`) |
| `appLE_comp` / `appTop_eq` | ★★合成則と `appTop` の言い換え |
| `ringSplit` | ★環の合成の適用は入れ子(`rfl`、**明示束縛子**) |
| `evalOn_comp` | ★★★★`evalOn` は `evalGlobal` に落ちる |
| `continuous_evalOn` | ★★★★★**連続性** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {X : Scheme.{0}} (V : X.Opens) (p : Spec (CommRingCat.of ℂ) ⟶ X)

/-- ★`evalOn` は `appLE` で書ける。 -/
theorem evalOn_appLE (h : p ⁻¹ᵁ V = ⊤) (c : ((X.presheaf.obj (op V)) : Type)) :
    evalOn p V h c
      = (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
          ((p.appLE V ⊤ (le_of_eq h.symm)).hom c) :=
  rfl

/-- ★★合成の `appLE`。 -/
theorem appLE_comp (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme) :
    (V.ι.appLE V ⊤ V.ι_preimage_self.ge) ≫ (q.appLE ⊤ ⊤ (by simp))
      = (q ≫ V.ι).appLE V ⊤ (le_of_eq (comp_preimage_eq_top V q).symm) :=
  Scheme.Hom.appLE_comp_appLE q V.ι V ⊤ ⊤ _ _

/-- ★`appTop` は `appLE ⊤ ⊤`。 -/
theorem appTop_eq (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme) :
    (Scheme.Hom.appTop q) = q.appLE ⊤ ⊤ (by simp) :=
  Scheme.Hom.app_eq_appLE q


/-- ★環の合成の適用は入れ子になる（`rfl`）。 -/
theorem ringSplit {A B C : CommRingCat.{0}} (f : A ⟶ B) (g : B ⟶ C) (x : (A : Type)) :
    ((f ≫ g).hom) x = g.hom (f.hom x) := rfl

/-- ★★★`V` を経由する点では `evalOn` は `V.toScheme` 上の `evalGlobal` である。 -/
theorem evalOn_comp (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme)
    (c : ((X.presheaf.obj (op V)) : Type)) :
    evalOn (q ≫ V.ι) V (comp_preimage_eq_top V q) c
      = evalGlobal q ((V.ι.appLE V ⊤ V.ι_preimage_self.ge).hom c) := by
  have h1 := congrArg (fun (m : _ ⟶ _) => (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
    ((CommRingCat.Hom.hom m) c)) (appLE_comp V q)
  have hsp := congrArg (fun z => (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom z)
    (ringSplit (V.ι.appLE V ⊤ V.ι_preimage_self.ge) (q.appLE ⊤ ⊤ (by simp)) c)
  have h2 := congrArg (fun (m : _ ⟶ _) => (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
    ((CommRingCat.Hom.hom m) ((V.ι.appLE V ⊤ V.ι_preimage_self.ge).hom c)))
    (appTop_eq V q)
  exact ((h2.trans hsp.symm).trans h1).symm


/-- ★★★★★**`V` を経由する点の上で `evalOn` は連続である**。 -/
theorem continuous_evalOn (c : ((X.presheaf.obj (op V)) : Type)) :
    @Continuous _ ℂ (arcTopology V.toScheme) _
      (fun q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme =>
        evalOn (q ≫ V.ι) V (comp_preimage_eq_top V q) c) := by
  letI := arcTopology V.toScheme
  have he : (fun q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme =>
        evalOn (q ≫ V.ι) V (comp_preimage_eq_top V q) c)
      = (fun q => evalGlobal q ((V.ι.appLE V ⊤ V.ι_preimage_self.ge).hom c)) :=
    funext (fun q => evalOn_comp V q c)
  rw [he]
  exact continuous_evalGlobal _


/-! ## ★出典の紐付け(`.src`) -/

def continuous_evalOn.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開集合上の評価は連続)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
