import ABC3.Found.Arakelov.ArcGenNorm

/-!
# Arakelov (C3) 第 258 ブロック —— ★★★★★★**半線形性と正規化ノルムの公式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★§9-300 の三重路を**係数環の選び直し**で解いた

§9-300 で「同じ台に 3 通りの綴りがあり、寄せ先が段ごとに違う」と特定した。
★逃げ道は **`arcFiber` の係数環を選び直す**ことだった:

    arcFiberG p F := ((pullback p).obj F).val.obj (op ⊤)     (係数環は Γ(Spec ℂ, ⊤))

★★これは `arcFiber`(係数環 `ℂ`)と**台が同じ**(`rfl`)だが、
**中間の段がすべて `Γ` 側で書ける**——`ℂ` への変換は最後の 1 回だけになる。

★★★半線形性は `Γ` 側では **3 行**で通った:

| 段 | 使うもの |
|---|---|
| 評価の線形性 | ★第 257 `arcEvalOn_smul` |
| 制限射の半線形性 | ★mathlib `PresheafOfModules.map_smul` |
| 繋ぐ | ★`trans` |

## ★★★★★★★これで連続性が閉じる

    genNorm (arcEvalOnTop (c·g)) = ‖evalOn c‖         (`genNorm_arcEvalOnTop`)

★右辺は**正則関数の絶対値**なので、第 252 の連続性がそのまま効く。

| 定義・定理 | 内容 |
|---|---|
| `arcFiberG` | ★係数環を `Γ` に取ったファイバー(台は `rfl` で同じ) |
| `arcEvalOnTopG` / `arcEvalOnTopG_smul` | ★★★`Γ` 側の評価と半線形性 |
| `evalOn` | ★`U` 上の関数の `p` での値 |
| `arcEvalOnTop_smul` | ★★★★`ℂ` 側の半線形性 |
| `genNorm_arcEvalOnTop` | ★★★★★★**正規化ノルム = 正則関数の絶対値** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {X : Scheme.{0}} (F : X.Modules) (p : Spec (CommRingCat.of ℂ) ⟶ X)

/-- ★係数環を `Γ(Spec ℂ, ⊤)` に取ったファイバー。 -/
noncomputable abbrev arcFiberG :=
  ((Scheme.Modules.pullback p).obj F).val.obj (op ⊤)

/-- ★台は `arcFiber` と同じ。 -/
theorem arcFiberG_carrier : (↥(arcFiberG F p) : Type) = (↥(arcFiber p F) : Type) := rfl

/-- ★★`Γ` 側の評価。 -/
noncomputable def arcEvalOnTopG (U : X.Opens) (h : p ⁻¹ᵁ U = ⊤)
    (s : (F.val.obj (op U) : Type)) : ↥(arcFiberG F p) :=
  ((((Scheme.Modules.pullback p).obj F).val.map
    (homOfLE (le_of_eq h.symm)).op)).hom (arcEvalOn F p U s)

/-- ★★★`Γ` 側での半線形性。 -/
theorem arcEvalOnTopG_smul (U : X.Opens) (h : p ⁻¹ᵁ U = ⊤)
    (c : ((X.presheaf.obj (op U)) : Type)) (s : (F.val.obj (op U) : Type)) :
    arcEvalOnTopG F p U h (c • s)
      = (((Spec (CommRingCat.of ℂ)).presheaf.map (homOfLE (le_of_eq h.symm)).op).hom
          ((p.app U).hom c)) • arcEvalOnTopG F p U h s := by
  have h1 := arcEvalOn_smul F p U c s
  have h2' : (((((Scheme.Modules.pullback p).obj F).val.map
        (homOfLE (le_of_eq h.symm)).op))).hom (arcEvalOn F p U (c • s))
      = (((((Scheme.Modules.pullback p).obj F).val.map
        (homOfLE (le_of_eq h.symm)).op))).hom
          (((p.app U).hom c) • arcEvalOn F p U s) :=
    congrArg (fun y => (((((Scheme.Modules.pullback p).obj F).val.map
      (homOfLE (le_of_eq h.symm)).op))).hom y) h1
  have h3 := PresheafOfModules.map_smul (M := ((Scheme.Modules.pullback p).obj F).val)
    (homOfLE (le_of_eq h.symm)).op ((p.app U).hom c) (arcEvalOn F p U s)
  exact h2'.trans h3


/-- ★`U` 上の関数の `p` での値。 -/
noncomputable def evalOn (U : X.Opens) (h : p ⁻¹ᵁ U = ⊤)
    (c : ((X.presheaf.obj (op U)) : Type)) : (CommRingCat.of ℂ : Type) :=
  (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
    (((Spec (CommRingCat.of ℂ)).presheaf.map (homOfLE (le_of_eq h.symm)).op).hom
      ((p.app U).hom c))

/-- ★★★★ℂ 側での半線形性。 -/
theorem arcEvalOnTop_smul (U : X.Opens) (h : p ⁻¹ᵁ U = ⊤)
    (c : ((X.presheaf.obj (op U)) : Type)) (s : (F.val.obj (op U) : Type)) :
    arcEvalOnTop F p U h (c • s) = (evalOn p U h c) • arcEvalOnTop F p U h s := by
  have hG := arcEvalOnTopG_smul F p U h c s
  have h4 : (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom (evalOn p U h c)
      = ((Spec (CommRingCat.of ℂ)).presheaf.map (homOfLE (le_of_eq h.symm)).op).hom
        ((p.app U).hom c) :=
    congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m)
      (((Spec (CommRingCat.of ℂ)).presheaf.map (homOfLE (le_of_eq h.symm)).op).hom
        ((p.app U).hom c))) (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom_inv_id
  exact hG.trans (congrArg (fun r => r • arcEvalOnTopG F p U h s) h4.symm)


variable (hF : IsLocallyTrivial X F.val)

/-- ★★★★★★正規化ノルムは正則関数の絶対値になる。 -/
theorem genNorm_arcEvalOnTop (V : X.Opens) (g : (F.val.obj (op V) : Type))
    (h : p ⁻¹ᵁ V = ⊤) (hg : arcEvalOnTop F p V h g ≠ 0)
    (c : ((X.presheaf.obj (op V)) : Type)) :
    genNorm F hF V g p h (arcEvalOnTop F p V h (c • g)) = ‖evalOn p V h c‖ := by
  rw [arcEvalOnTop_smul, genNorm_smul, genNorm_self F hF V g p h hg, mul_one]


/-! ## ★出典の紐付け(`.src`) -/

def genNorm_arcEvalOnTop.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——半線形性と正規化ノルムの公式)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
